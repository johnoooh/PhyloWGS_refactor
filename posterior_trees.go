package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
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

// treeRecord pairs an archive tree index with its computed signature.
type treeRecord struct {
	Idx int
	Sig string
}

// PosteriorGroup is one cluster of trees sharing a topology signature.
type PosteriorGroup struct {
	Signature   string
	TreeIndices []int
	Probability float64
}

func groupTreesBySignature(records []treeRecord) map[string][]int {
	out := map[string][]int{}
	for _, r := range records {
		out[r.Sig] = append(out[r.Sig], r.Idx)
	}
	return out
}

// rankGroups returns groups sorted descending by posterior probability,
// ties broken by lowest first tree index for determinism.
func rankGroups(groups map[string][]int, totalTrees int) []PosteriorGroup {
	out := make([]PosteriorGroup, 0, len(groups))
	for sig, idxs := range groups {
		sort.Ints(idxs)
		out = append(out, PosteriorGroup{
			Signature:   sig,
			TreeIndices: idxs,
			Probability: float64(len(idxs)) / float64(totalTrees),
		})
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Probability != out[j].Probability {
			return out[i].Probability > out[j].Probability
		}
		return out[i].TreeIndices[0] < out[j].TreeIndices[0]
	})
	return out
}

// writePosteriorTreeTeX writes one posterior group's standalone-LaTeX file.
// `representative` is the JSON snapshot of one tree in the group (used for
// structure + node membership); `freqs` provides cellular-prevalence stats
// already aggregated across all trees in the group, keyed by node ID
// (e.g. "0", "1") and shaped [n_aggregated_samples][n_timepoints].
func writePosteriorTreeTeX(outPath string, representative json.RawMessage, g PosteriorGroup, freqs map[string][][]float64) error {
	var s struct {
		Populations map[string]struct {
			CellularPrevalence []float64 `json:"cellular_prevalence"`
			NumSSMs            int       `json:"num_ssms"`
			NumCNVs            int       `json:"num_cnvs"`
		} `json:"populations"`
		Structure      map[string][]int               `json:"structure"`
		MutAssignments map[string]map[string][]string `json:"mut_assignments"`
	}
	if err := json.Unmarshal(representative, &s); err != nil {
		return err
	}

	var b strings.Builder
	b.WriteString(`\documentclass{standalone}` + "\n")
	b.WriteString(`\usepackage{tikz}` + "\n")
	b.WriteString(`\usepackage{multicol}` + "\n")
	b.WriteString(`\usetikzlibrary{fit,positioning}` + "\n")
	b.WriteString(`\begin{document}` + "\n")
	b.WriteString(`\begin{tikzpicture}` + "\n")
	b.WriteString(`\node (a) at (0,0){` + "\n")
	b.WriteString(`\begin{tikzpicture}` + "\n")
	b.WriteString(`[grow=east, ->, level distance=20mm,` +
		`every node/.style={circle, minimum size = 8mm, thick, draw =black,inner sep=2mm},` +
		`every label/.append style={shape=rectangle, yshift=-1mm},` +
		`level 2/.style={sibling distance=50mm},` +
		`level 3/.style={sibling distance=20mm},` +
		`level 4/.style={sibling distance=20mm},` +
		`every edge/.style={-latex, thick}]` + "\n")

	// Tree structure: \node {0} child {\node {1} child {...}}.
	counter := 0
	var emitNode func(node string)
	emitNode = func(node string) {
		counter++
		mine := counter
		b.WriteString(fmt.Sprintf(`\node {%d}`, mine))
		for _, c := range s.Structure[node] {
			b.WriteString("child {")
			emitNode(fmt.Sprintf("%d", c))
			b.WriteString("}")
		}
	}
	b.WriteString("\n")
	emitNode("0")
	b.WriteString(";\n")
	b.WriteString(`\end{tikzpicture}` + "\n")
	b.WriteString(`};` + "\n")

	// Index table.
	b.WriteString(`\node (b) at (a.south)[anchor=north,yshift=-.5cm]{` + "\n")
	b.WriteString(`\begin{tikzpicture}` + "\n")
	b.WriteString(`\node (table){` + "\n")
	// Determine n_timepoints from any population's cellular_prevalence.
	nTP := 0
	for _, p := range s.Populations {
		if len(p.CellularPrevalence) > nTP {
			nTP = len(p.CellularPrevalence)
		}
	}
	b.WriteString(`\begin{tabular}{|c|l|l|`)
	for i := 0; i < nTP; i++ {
		b.WriteString("l|")
	}
	b.WriteString("}\n\\hline\n")
	b.WriteString(fmt.Sprintf(`Node & \multicolumn{1}{|c|}{SSMs} & \multicolumn{1}{|c|}{CNVs} & \multicolumn{%d}{|c|}{Clonal frequencies}\\`, nTP))
	b.WriteString("\n\\hline\n")

	counter = 0
	var emitRow func(node string)
	emitRow = func(node string) {
		counter++
		mine := counter
		ma := s.MutAssignments[node]
		nSSM := len(ma["ssms"])
		nCNV := len(ma["cnvs"])
		f := freqs[node]
		row := fmt.Sprintf("%d & %d & %d", mine, nSSM, nCNV)
		for tp := 0; tp < nTP; tp++ {
			mean, sd := meanSDColumn(f, tp)
			row += fmt.Sprintf(` & %.3f $\pm$ %.3f`, mean, sd)
		}
		row += `\\` + "\n"
		b.WriteString(row)
		for _, c := range s.Structure[node] {
			emitRow(fmt.Sprintf("%d", c))
		}
	}
	emitRow("0")

	b.WriteString("\\hline\n")
	b.WriteString("\\end{tabular}\n};\n")
	b.WriteString(`\end{tikzpicture}` + "\n")
	b.WriteString(`};` + "\n")
	b.WriteString(fmt.Sprintf(`\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: %g};`+"\n", g.Probability))
	b.WriteString(`\end{tikzpicture}` + "\n")
	b.WriteString(`\end{document}` + "\n")

	return os.WriteFile(outPath, []byte(b.String()), 0644)
}

// meanSDColumn computes mean and (population) sd across rows of a 2D matrix.
// Returns (0, 0) for empty input.
func meanSDColumn(rows [][]float64, col int) (float64, float64) {
	if len(rows) == 0 {
		return 0, 0
	}
	var sum, sumSq float64
	n := 0
	for _, r := range rows {
		if col < len(r) {
			sum += r[col]
			sumSq += r[col] * r[col]
			n++
		}
	}
	if n == 0 {
		return 0, 0
	}
	mean := sum / float64(n)
	variance := sumSq/float64(n) - mean*mean
	if variance < 0 {
		variance = 0
	}
	return mean, math.Sqrt(variance)
}

// posteriorTreesConfig parameterizes the runPosteriorTrees driver.
//   - OutDir: directory containing trees.zip; posterior_trees/ output is
//     written under it.
//   - NumTrees: cap on top-K groups to write. 0 means "all groups".
//   - SkipPDF: if true, skip pdflatex invocation entirely (tests/headless).
type posteriorTreesConfig struct {
	OutDir   string
	NumTrees int  // 0 means "all groups"
	SkipPDF  bool // tests / headless servers
}

// runPosteriorTrees is the top-level driver that mirrors Python's
// posterior_trees.py main(): read trees.zip, group by topology signature,
// rank by posterior probability, and write a .tex (+ optional .pdf) per
// group up to NumTrees.
func runPosteriorTrees(cfg posteriorTreesConfig) error {
	zipPath := filepath.Join(cfg.OutDir, "trees.zip")
	r, err := newTreeArchiveReader(zipPath)
	if err != nil {
		return fmt.Errorf("reading %s: %w", zipPath, err)
	}
	defer r.Close()
	if r.NumTrees() == 0 {
		return fmt.Errorf("no trees in archive")
	}

	// Build mutation-id index map from the first tree.
	first, _, err := r.LoadTree(0)
	if err != nil {
		return err
	}
	idx, err := buildMutationIndexMap(first)
	if err != nil {
		return err
	}

	// Compute signatures for every tree.
	records := make([]treeRecord, 0, r.NumTrees())
	snapshots := make([]json.RawMessage, r.NumTrees())
	for i := 0; i < r.NumTrees(); i++ {
		raw, _, err := r.LoadTree(i)
		if err != nil {
			return err
		}
		snapshots[i] = raw
		sig, err := treeSignature(raw, idx)
		if err != nil {
			return err
		}
		records = append(records, treeRecord{Idx: i, Sig: sig})
	}

	groups := groupTreesBySignature(records)
	ranked := rankGroups(groups, len(records))

	limit := cfg.NumTrees
	if limit <= 0 || limit > len(ranked) {
		limit = len(ranked)
	}

	outDir := filepath.Join(cfg.OutDir, "posterior_trees")
	if err := os.MkdirAll(outDir, 0755); err != nil {
		return err
	}

	for rank := 0; rank < limit; rank++ {
		g := ranked[rank]
		// Pull all snapshots for this group, then aggregate phi.
		groupSnaps := make([]json.RawMessage, 0, len(g.TreeIndices))
		for _, tidx := range g.TreeIndices {
			groupSnaps = append(groupSnaps, snapshots[tidx])
		}
		freqs, err := aggregateFreqsByNode(groupSnaps)
		if err != nil {
			return err
		}
		texPath := filepath.Join(outDir, fmt.Sprintf("tree_%d_%g.tex", rank, g.Probability))
		if err := writePosteriorTreeTeX(texPath, snapshots[g.TreeIndices[0]], g, freqs); err != nil {
			return err
		}
		if !cfg.SkipPDF {
			if _, err := exec.LookPath("pdflatex"); err == nil {
				cmd := exec.Command("pdflatex",
					"-interaction=nonstopmode",
					"-output-directory="+outDir,
					texPath,
				)
				_ = cmd.Run() // best-effort, mirrors Python (which also swallows OSError)
			} else {
				fmt.Fprintln(os.Stderr, "pdflatex not available; .tex files written without PDF")
			}
		}
	}
	return nil
}

// aggregateFreqsByNode pools cellular_prevalence vectors across all trees
// in a posterior group, indexed by node id. Caller guarantees all trees
// share the same topology+mutation-set signature, so node IDs are stable.
func aggregateFreqsByNode(trees []json.RawMessage) (map[string][][]float64, error) {
	out := map[string][][]float64{}
	for _, t := range trees {
		var s struct {
			Populations map[string]struct {
				CellularPrevalence []float64 `json:"cellular_prevalence"`
			} `json:"populations"`
		}
		if err := json.Unmarshal(t, &s); err != nil {
			return nil, err
		}
		for node, p := range s.Populations {
			out[node] = append(out[node], p.CellularPrevalence)
		}
	}
	return out, nil
}
