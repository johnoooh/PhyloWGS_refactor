package main

import (
	"encoding/json"
	"fmt"
	"sort"
	"strings"
)

// treeSignature computes the per-tree topology signature used by
// posterior_trees.py:34-49 to group MCMC samples by topology. Two trees
// with the same signature are considered the same posterior topology.
//
// idx is a map from mutation ID (e.g. "s0") to its canonical index string
// (e.g. "0"). Build it once per archive scan via buildMutationIndexMap.
//
// Algorithmically, this is a faithful port of the Python `descend` +
// `_sort` + `sort_and_merge` pipeline:
//
//   - Initialize per-node genotype to "" for ALL nodes.
//   - DFS from "0". For each node WITH data:
//     gtype = parent_genotype_in_dict + concat(idx[m]+"_" for m in ssms+cnvs)
//     genotype[node] = gtype
//     gtypelist.append(_sort(gtype))
//   - For each node WITHOUT data: skip (genotype[node] stays "").
//   - Final signature = sort_and_merge(gtypelist).
//
// Two correctness subtleties matter for parity:
//
//  1. _sort splits the trailing-"_" gtype on "_", sorts the parts (which
//     puts the empty trailing token first), then re-joins with "_"
//     suffixes — producing a leading underscore, e.g. "0_" -> "_0_".
//
//  2. An empty intermediate node does NOT propagate its parent's
//     genotype down. Its slot in the genotype dict remains the initial
//     "" value, so children read "" from the lookup and start fresh.
func treeSignature(snapshot json.RawMessage, idx map[string]string) (string, error) {
	var s struct {
		Structure      map[string][]int               `json:"structure"`
		MutAssignments map[string]map[string][]string `json:"mut_assignments"`
		Populations    map[string]json.RawMessage     `json:"populations"`
	}
	if err := json.Unmarshal(snapshot, &s); err != nil {
		return "", fmt.Errorf("treeSignature: %w", err)
	}

	// genotype mirrors the Python `genotype` dict: every node observed
	// in the tree starts with "". Use the union of populations keys and
	// mut_assignments keys so we cover every node id even if one map is
	// missing an entry.
	genotype := map[string]string{}
	for k := range s.Populations {
		genotype[k] = ""
	}
	for k := range s.MutAssignments {
		if _, ok := genotype[k]; !ok {
			genotype[k] = ""
		}
	}

	var gtypelist []string

	// Recursive DFS from "0". The Python code passes the actual node
	// reference; the parent lookup happens via the genotype dict, NOT
	// via a recursion argument. Replicating that exactly here.
	var descend func(node string, parent string)
	descend = func(node string, parent string) {
		ma := s.MutAssignments[node]
		ndata := len(ma["ssms"]) + len(ma["cnvs"])

		if ndata > 0 {
			gtype := ""
			if parent != "" {
				// Parent exists -> read parent's genotype from the dict.
				// If the parent was an empty node, this returns "" (the
				// initial value). This is the Python behavior we depend on.
				gtype = genotype[parent]
			}
			for _, m := range ma["ssms"] {
				if v, ok := idx[m]; ok {
					gtype += v + "_"
				}
			}
			for _, m := range ma["cnvs"] {
				if v, ok := idx[m]; ok {
					gtype += v + "_"
				}
			}
			genotype[node] = gtype
			gtypelist = append(gtypelist, sortGenotype(gtype))
		}

		// Recurse into children. Sorted for determinism; the final
		// sort_and_merge makes order irrelevant for the output, but
		// stable iteration helps with debugging and tests.
		children := append([]int(nil), s.Structure[node]...)
		sort.Ints(children)
		for _, c := range children {
			descend(fmt.Sprintf("%d", c), node)
		}
	}

	// Use "" as the sentinel for "no parent" since "0" is a valid node ID.
	descend("0", "")

	return sortAndMerge(gtypelist), nil
}

// sortGenotype mirrors posterior_trees.py:_sort. Splits on "_", sorts the
// resulting tokens lexically (this puts the trailing empty token at the
// front), then re-joins with "_" + trailing "_". E.g. "0_1_" -> tokens
// ["0","1",""] -> sorted ["","0","1"] -> "_0_1_".
func sortGenotype(s string) string {
	parts := strings.Split(s, "_")
	sort.Strings(parts)
	var b strings.Builder
	for _, p := range parts {
		b.WriteString(p)
		b.WriteByte('_')
	}
	return b.String()
}

// sortAndMerge mirrors posterior_trees.py:sort_and_merge. Sorts the list
// of genotype strings and concatenates each followed by ";".
func sortAndMerge(gtypes []string) string {
	sorted := append([]string(nil), gtypes...)
	sort.Strings(sorted)
	var b strings.Builder
	for _, g := range sorted {
		b.WriteString(g)
		b.WriteByte(';')
	}
	return b.String()
}

// buildMutationIndexMap returns mutation-id -> canonical index string for
// use as the `idx` argument of treeSignature. Python's posterior_trees.py
// builds this from `load_data` order: `dict([(datum.name, str(i)) for
// i, datum in enumerate(codes)])`. We derive the same ordering by sorting
// the mutation IDs observed in the snapshot's mut_assignments and
// assigning 0,1,2,... in that order.
//
// Stable across archive reads as long as the same mutation set is present
// in every snapshot — which is the invariant for an MCMC archive over a
// fixed dataset.
func buildMutationIndexMap(snapshot json.RawMessage) (map[string]string, error) {
	var s struct {
		MutAssignments map[string]map[string][]string `json:"mut_assignments"`
	}
	if err := json.Unmarshal(snapshot, &s); err != nil {
		return nil, fmt.Errorf("buildMutationIndexMap: %w", err)
	}
	seen := map[string]struct{}{}
	for _, ma := range s.MutAssignments {
		for _, m := range ma["ssms"] {
			seen[m] = struct{}{}
		}
		for _, m := range ma["cnvs"] {
			seen[m] = struct{}{}
		}
	}
	ids := make([]string, 0, len(seen))
	for id := range seen {
		ids = append(ids, id)
	}
	sort.Strings(ids)
	out := make(map[string]string, len(ids))
	for i, id := range ids {
		out[id] = fmt.Sprintf("%d", i)
	}
	return out, nil
}
