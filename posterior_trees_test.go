package main

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestTreeSignature_TwoNodeChain validates a simple parent->child tree where
// both nodes have data. The expected signature is what the upstream Python
// implementation in posterior_trees.py:34-49 produces:
//
//	node 0 gtype = "0_"        -> _sort -> "_0_"
//	node 1 gtype = "0_" + "1_" = "0_1_" -> _sort -> "_0_1_"
//	sort_and_merge(["_0_", "_0_1_"]) = "_0_;_0_1_;"
//
// Note the leading underscore: Python's _sort splits "0_" into ["0", ""],
// sorts to ["", "0"], then joins with trailing "_" producing "_0_". This is
// faithful to the Python — anything else would produce different
// signatures and break grouping parity.
func TestTreeSignature_TwoNodeChain(t *testing.T) {
	snapshot := json.RawMessage(`{
		"populations":{"0":{},"1":{}},
		"structure":{"0":[1]},
		"mut_assignments":{
			"0":{"ssms":["s0"],"cnvs":[]},
			"1":{"ssms":["s1"],"cnvs":[]}
		}
	}`)
	idx := map[string]string{"s0": "0", "s1": "1"}

	sig, err := treeSignature(snapshot, idx)
	if err != nil {
		t.Fatal(err)
	}
	want := "_0_;_0_1_;"
	if sig != want {
		t.Errorf("sig = %q, want %q", sig, want)
	}
}

// TestTreeSignature_SkipsEmptyNodes encodes a critical correctness property
// of the Python algorithm: when an intermediate node has no data, that node
// does NOT update the genotype dict (initialized to "" for all nodes), so a
// child's genotype lookup of an empty parent returns "" (the initial value),
// not the grandparent's genotype.
//
// Tree: 0 -> 1 -> 2, data: 0=[s0], 1=[], 2=[s2], idx={s0:"0", s2:"1"}.
//
//	node 0: data=[s0], parent=None -> gtype = "0_"   -> _sort -> "_0_"
//	node 1: data=[], skipped, genotype[1] stays ""
//	node 2: data=[s2], parent=node1 -> gtype = "" + "1_" = "1_" -> _sort -> "_1_"
//	sort_and_merge(["_0_", "_1_"]) = "_0_;_1_;"
func TestTreeSignature_SkipsEmptyNodes(t *testing.T) {
	snapshot := json.RawMessage(`{
		"populations":{"0":{},"1":{},"2":{}},
		"structure":{"0":[1],"1":[2]},
		"mut_assignments":{
			"0":{"ssms":["s0"],"cnvs":[]},
			"1":{"ssms":[],"cnvs":[]},
			"2":{"ssms":["s2"],"cnvs":[]}
		}
	}`)
	idx := map[string]string{"s0": "0", "s2": "1"}

	sig, err := treeSignature(snapshot, idx)
	if err != nil {
		t.Fatal(err)
	}
	want := "_0_;_1_;"
	if sig != want {
		t.Errorf("sig = %q, want %q", sig, want)
	}
}

// TestTreeSignature_BranchInheritance validates the parent-genotype lookup
// path with multi-mutation accumulation across a deeper chain.
//
// Tree: 0 -> 1 -> 2 with all nodes carrying data.
//
//	node 0: gtype = "0_"        -> _sort -> "_0_"
//	node 1: parent_gt="0_" + "1_" = "0_1_" -> _sort -> "_0_1_"
//	node 2: parent_gt="0_1_" + "2_" = "0_1_2_" -> _sort -> "_0_1_2_"
func TestTreeSignature_BranchInheritance(t *testing.T) {
	snapshot := json.RawMessage(`{
		"populations":{"0":{},"1":{},"2":{}},
		"structure":{"0":[1],"1":[2]},
		"mut_assignments":{
			"0":{"ssms":["s0"],"cnvs":[]},
			"1":{"ssms":["s1"],"cnvs":[]},
			"2":{"ssms":["s2"],"cnvs":[]}
		}
	}`)
	idx := map[string]string{"s0": "0", "s1": "1", "s2": "2"}

	sig, err := treeSignature(snapshot, idx)
	if err != nil {
		t.Fatal(err)
	}
	want := "_0_;_0_1_;_0_1_2_;"
	if sig != want {
		t.Errorf("sig = %q, want %q", sig, want)
	}
}

func TestGroupAndRank_AssignsCorrectPosteriorProbabilities(t *testing.T) {
	// Three trees: two share signature "A;", one has signature "B;".
	// Expect group A: prob 2/3, group B: prob 1/3.
	groups := groupTreesBySignature([]treeRecord{
		{Idx: 0, Sig: "A;"},
		{Idx: 1, Sig: "A;"},
		{Idx: 2, Sig: "B;"},
	})
	ranked := rankGroups(groups, 3 /* total trees */)
	if len(ranked) != 2 {
		t.Fatalf("expected 2 groups, got %d", len(ranked))
	}
	// Highest prob first.
	if ranked[0].Probability != 2.0/3.0 {
		t.Errorf("rank0 prob = %v, want %v", ranked[0].Probability, 2.0/3.0)
	}
	if ranked[1].Probability != 1.0/3.0 {
		t.Errorf("rank1 prob = %v, want %v", ranked[1].Probability, 1.0/3.0)
	}
	if ranked[0].TreeIndices[0] != 0 || ranked[0].TreeIndices[1] != 1 {
		t.Errorf("rank0 indices = %v, want [0 1]", ranked[0].TreeIndices)
	}
}

func TestWritePosteriorTreeTeX_EmitsValidStandaloneDocument(t *testing.T) {
	tmp := t.TempDir()
	out := filepath.Join(tmp, "tree.tex")

	// Single-tree group, single population, one SSM, two samples (so cell_prev is len-2).
	rep := json.RawMessage(`{
		"populations":{
			"0":{"cellular_prevalence":[0.7,0.6],"num_ssms":1,"num_cnvs":0}
		},
		"structure":{"0":[]},
		"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}
	}`)
	group := PosteriorGroup{
		Signature:   "_0_;",
		TreeIndices: []int{0},
		Probability: 1.0,
	}
	freqs := map[string][][]float64{"0": {{0.7, 0.6}}}

	if err := writePosteriorTreeTeX(out, rep, group, freqs); err != nil {
		t.Fatal(err)
	}

	body, _ := os.ReadFile(out)
	for _, want := range []string{
		`\documentclass{standalone}`,
		`\begin{tikzpicture}`,
		`Posterior probability: 1`,
	} {
		if !strings.Contains(string(body), want) {
			t.Errorf("tex missing %q", want)
		}
	}

	// Guard against a stray \\ (LaTeX line break) leaking into tikz scope.
	if strings.Contains(string(body), `\\`+"\n"+`\node`) {
		t.Errorf("found stray \\\\ before \\node — would break pdflatex")
	}
	if !strings.Contains(string(body), "\n"+`\node {1}`) {
		t.Errorf("expected \\node {1} immediately after a newline; got: %s", body)
	}
}

func TestAggregateFreqsByNode_AveragesAcrossTrees(t *testing.T) {
	// Two trees, same topology. Node 0 has cell_prev=[0.7] in tree A,
	// [0.5] in tree B. Aggregate should record both rows.
	a := json.RawMessage(`{"populations":{"0":{"cellular_prevalence":[0.7],"num_ssms":1,"num_cnvs":0}},"structure":{"0":[]},"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}}`)
	b := json.RawMessage(`{"populations":{"0":{"cellular_prevalence":[0.5],"num_ssms":1,"num_cnvs":0}},"structure":{"0":[]},"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}}`)
	freqs, err := aggregateFreqsByNode([]json.RawMessage{a, b})
	if err != nil {
		t.Fatal(err)
	}
	rows := freqs["0"]
	if len(rows) != 2 {
		t.Fatalf("expected 2 rows for node 0, got %d", len(rows))
	}
}

func TestRunPosteriorTrees_EndToEnd(t *testing.T) {
	tmp := t.TempDir()

	// Build a small trees.zip with 2 identical-topology trees.
	zw, _ := newTreeArchiveWriter(filepath.Join(tmp, "trees.zip"))
	snap := json.RawMessage(`{
		"populations":{"0":{"cellular_prevalence":[0.7],"num_ssms":1,"num_cnvs":0}},
		"structure":{"0":[]},
		"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}
	}`)
	zw.WriteSample(0, -100.0, false, snap)
	zw.WriteSample(1, -101.0, false, snap)
	zw.Close()

	if err := runPosteriorTrees(posteriorTreesConfig{
		OutDir:   tmp,
		NumTrees: 0, // 0 = all
		SkipPDF:  true,
	}); err != nil {
		t.Fatal(err)
	}

	// Expect posterior_trees/tree_0_*.tex to exist.
	entries, _ := os.ReadDir(filepath.Join(tmp, "posterior_trees"))
	hasTex := false
	for _, e := range entries {
		if strings.HasSuffix(e.Name(), ".tex") {
			hasTex = true
		}
	}
	if !hasTex {
		t.Errorf("no .tex file written; got entries: %v", entries)
	}
}
