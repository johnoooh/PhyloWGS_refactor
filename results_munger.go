package main

// Go port of pwgsresults/result_munger.py (mskcc/phylowgs): the small-node
// removal path and the shared node-removal machinery. Superclone and
// multiprimary filtering live in results_munger_tree.go.
//
// Iteration order: wherever Python iterates a dict (tree summaries,
// populations, structure, the removed-node set) this port iterates keys
// ascending, which is what CPython 2 does for small non-negative int keys.
// Python asserts and KeyErrors become errors naming the tree.

import (
	"fmt"
	"math"
	"sort"
)

// nodeRemoved marks a removed node in an idx map (Python's None). Node ids
// are never negative, so -1 also satisfies Python 2's `None < sidx`.
const nodeRemoved = -1

type mutDest int

const (
	destBest mutDest = iota
	destClonal
)

func sortedKeys[V any](m map[int]V) []int {
	ks := make([]int, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	sort.Ints(ks)
	return ks
}

// pyMean is np.mean on a 1-D vector: float64 sum in index order over n.
// Empty input gives NaN, as numpy does. numpy sums pairwise for >= 8 elements.
func pyMean[T int | float64](xs []T) float64 {
	sum := 0.0
	for _, x := range xs {
		sum += float64(x)
	}
	return sum / float64(len(xs))
}

// findParent returns the first structure key (ascending) whose child list
// contains idx, like Python's scan of structure.items().
func findParent(structure map[int][]int, idx int) (int, bool) {
	for _, p := range sortedKeys(structure) {
		for _, c := range structure[p] {
			if c == idx {
				return p, true
			}
		}
	}
	return 0, false
}

func containsInt(xs []int, v int) bool {
	for _, x := range xs {
		if x == v {
			return true
		}
	}
	return false
}

// appendOwned appends without writing into a backing array another slice
// may share (Python lists never alias like that).
func appendOwned(dst []string, src ...string) []string {
	return append(dst[:len(dst):len(dst)], src...)
}

// RemoveSmallNodes ports ResultMunger.remove_small_nodes. Note Python calls
// _remove_nodes for every tree even when no node is small, so mut counts of
// every tree are recounted from MutAss.
func (rs *ResultSet) RemoveSmallNodes(minSSMs float64) error {
	for _, tidx := range sortedKeys(rs.Trees) {
		small, err := rs.findSmallNodes(tidx, minSSMs)
		if err != nil {
			return err
		}
		if err := rs.removeNodes(small, tidx, destBest); err != nil {
			return err
		}
	}
	return nil
}

// ssmThreshold: >= 1 is a count (int() truncates); < 1 is a fraction of all
// SSMs, rounded half away from zero (Python 2 round == math.Round).
func ssmThreshold(minSSMs float64, numSSMs int) int {
	if minSSMs >= 1 {
		return int(minSSMs)
	}
	return int(math.Round(minSSMs * float64(numSSMs)))
}

// findSmallNodes ports _find_small_nodes, returning node ids ascending.
// Quirk kept: `parent` persists across iterations, so a node found under no
// parent reuses the previous node's parent (Python NameError if none yet).
// The phi check compares a node only with the previous node in sorted order,
// and only if that previous node is in this node's parent's child list.
func (rs *ResultSet) findSmallNodes(tidx int, minSSMs float64) ([]int, error) {
	ts := rs.Trees[tidx]
	threshold := ssmThreshold(minSSMs, len(rs.Muts.SSMs))
	var small []int
	parent, haveParent := 0, false
	lastPhi, lastIdx, haveLast := 0.0, 0, false
	for _, s := range sortedKeys(ts.Populations) {
		if p, ok := findParent(ts.Structure, s); ok {
			parent, haveParent = p, true
		}
		pop := ts.Populations[s]
		phi := pyMean(pop.CellPrev)
		if s > 0 {
			if !haveParent {
				return nil, fmt.Errorf("tree %d: node %d has no parent", tidx, s)
			}
			if haveLast && containsInt(ts.Structure[parent], lastIdx) && !(phi <= lastPhi) {
				return nil, fmt.Errorf("tree %d: node %d phi %v exceeds preceding sibling %d phi %v", tidx, s, phi, lastIdx, lastPhi)
			}
		}
		lastPhi, lastIdx, haveLast = phi, s, true
		if s == 0 || pop.NumSSMs >= threshold {
			continue
		}
		small = append(small, s)
	}
	return small, nil
}

// removeNodes ports _remove_nodes.
func (rs *ResultSet) removeNodes(nodes []int, tidx int, dest mutDest) error {
	ts := rs.Trees[tidx]
	idxMap := map[int]int{}
	for _, n := range nodes {
		if _, ok := ts.Populations[n]; !ok {
			return fmt.Errorf("tree %d: no population %d to remove", tidx, n)
		}
		delete(ts.Populations, n)
		idxMap[n] = nodeRemoved
	}
	if err := removeFromStructure(tidx, idxMap, ts.Structure); err != nil {
		return err
	}
	if err := renumberNodes(tidx, ts, idxMap); err != nil {
		return err
	}
	if err := rs.reassignMuts(tidx, idxMap, dest); err != nil {
		return err
	}
	return rs.correctMutCounts(tidx)
}

// removeFromStructure ports _remove_nodes_from_tree_structure: removed nodes
// in ascending order; each one's children are spliced into its parent, whose
// list is then re-sorted (and dropped if empty).
func removeFromStructure(tidx int, idxMap map[int]int, structure map[int][]int) error {
	for _, rem := range sortedKeys(idxMap) {
		if idxMap[rem] != nodeRemoved {
			continue
		}
		parent, ok := findParent(structure, rem)
		if !ok {
			return fmt.Errorf("tree %d: could not find parent of %d in %v", tidx, rem, structure)
		}
		kept := []int{}
		for _, c := range structure[parent] {
			if c != rem {
				kept = append(kept, c)
			}
		}
		if ch, ok := structure[rem]; ok {
			kept = append(kept, ch...)
			delete(structure, rem)
		}
		sort.Ints(kept)
		if len(kept) == 0 {
			delete(structure, parent)
		} else {
			structure[parent] = kept
		}
	}
	return nil
}

// renumberNodes ports _renumber_nodes: close gaps in 1..max(remaining),
// moving populations and structure keys down in ascending order, then
// renumber child ids. Child lists are not re-sorted (the map is monotone).
func renumberNodes(tidx int, ts *TreeSummary, idxMap map[int]int) error {
	idxs := sortedKeys(ts.Populations)
	if len(idxs) == 0 {
		return fmt.Errorf("tree %d: no populations left", tidx)
	}
	removed := 0
	for s := 1; s <= idxs[len(idxs)-1]; s++ {
		if _, ok := ts.Populations[s]; !ok {
			removed++
		} else if removed > 0 {
			idxMap[s] = s - removed
		}
	}
	for _, s := range idxs {
		n, ok := idxMap[s]
		if !ok {
			continue
		}
		ts.Populations[n] = ts.Populations[s]
		delete(ts.Populations, s)
		if ch, ok := ts.Structure[s]; ok {
			ts.Structure[n] = ch
			delete(ts.Structure, s)
		}
	}
	for p, children := range ts.Structure {
		out := make([]int, len(children))
		for i, c := range children {
			out[i] = c
			if n, ok := idxMap[c]; ok {
				if n == nodeRemoved {
					return fmt.Errorf("tree %d: removed node %d still a child of %d", tidx, c, p)
				}
				out[i] = n
			}
		}
		ts.Structure[p] = out
	}
	return nil
}

// correctMutCounts ports _correct_mut_counts. A non-root node with no MutAss
// entry is a KeyError in Python, so it is an error here too.
func (rs *ResultSet) correctMutCounts(tidx int) error {
	for _, s := range sortedKeys(rs.Trees[tidx].Populations) {
		if s == 0 {
			continue
		}
		m, ok := rs.MutAss[tidx][s]
		if !ok {
			return fmt.Errorf("tree %d: node %d has no mut_assignments entry", tidx, s)
		}
		pop := rs.Trees[tidx].Populations[s]
		pop.NumSSMs, pop.NumCNVs = len(m.SSMs), len(m.CNVs)
	}
	return nil
}
