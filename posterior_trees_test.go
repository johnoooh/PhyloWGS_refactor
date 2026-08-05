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
// This is a RAW-INPUT unit test of treeSignature's algorithm in isolation.
// The real pipeline can no longer produce this exact input at the point
// treeSignature runs on an archived snapshot: since M6, snapshotTree computes
// the topology signature BEFORE postProcessSummary removes empty nodes,
// specifically so this inheritance-break logic fires on real data (an
// archived snapshot has already had its empty nodes spliced out by the time
// it's read back). Kept as a direct algorithm-level test — see
// TestSnapshotTree_StampsSignatureBeforeCleanup for the end-to-end case that
// exercises the real call path.
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

// TestSnapshotTree_StampsSignatureBeforeCleanup is the M6 end-to-end
// regression test. It builds a real TSSB with a genuinely empty
// intermediate node (root -> P[data] -> E[empty] -> C[data]) and runs it
// through the actual snapshotTree call path, proving the stamped
// topology_signature reflects the pre-cleanup tree (where E's emptiness
// breaks genotype inheritance) rather than the post-cleanup one (where E
// has already been spliced out and C would incorrectly inherit P's
// genotype). This is the scenario TestTreeSignature_SkipsEmptyNodes tests
// in isolation; this test proves it actually fires on the real pipeline.
func TestSnapshotTree_StampsSignatureBeforeCleanup(t *testing.T) {
	ssms := []*SSM{{ID: "s0"}, {ID: "s1"}}

	cNode := &Node{ID: 3, Params: []float64{0.6}, Data: []int{1}} // s1
	eNode := &Node{ID: 2, Params: []float64{0.5}}                 // empty
	pNode := &Node{ID: 1, Params: []float64{0.7}, Data: []int{0}} // s0
	rootNode := &Node{ID: 0, Params: []float64{1.0}}

	cTSSB := &TSSBNode{Node: cNode}
	eTSSB := &TSSBNode{Node: eNode, Children: []*TSSBNode{cTSSB}}
	pTSSB := &TSSBNode{Node: pNode, Children: []*TSSBNode{eTSSB}}
	rootTSSB := &TSSBNode{Node: rootNode, Children: []*TSSBNode{pTSSB}}

	tssb := &TSSB{Root: rootTSSB, Data: ssms}
	mutIdx := buildMutationIndexMapFromDataset(tssb.Data, tssb.CNVData)

	raw := snapshotTree(tssb, -10.0, 0, 0, mutIdx)
	if raw == nil {
		t.Fatal("snapshotTree returned nil")
	}

	var parsed struct {
		TopologySignature *string `json:"topology_signature"`
	}
	if err := json.Unmarshal(raw, &parsed); err != nil {
		t.Fatalf("unmarshal snapshot: %v", err)
	}
	if parsed.TopologySignature == nil {
		t.Fatal("topology_signature field missing from snapshot")
	}

	// Pre-cleanup signature: E's emptiness breaks inheritance, so C reads
	// "" from its parent rather than P's genotype.
	want := "_0_;_1_;"
	if *parsed.TopologySignature != want {
		t.Errorf("topology_signature = %q, want %q (pre-cleanup, break-at-empty-node)", *parsed.TopologySignature, want)
	}

	// Sanity: confirm E really was spliced out of the returned (cleaned)
	// structure, so we know this test is exercising the actual pre- vs
	// post-cleanup distinction and not just echoing an untouched tree.
	var cleaned struct {
		Populations map[string]json.RawMessage `json:"populations"`
	}
	if err := json.Unmarshal(raw, &cleaned); err != nil {
		t.Fatalf("unmarshal cleaned snapshot: %v", err)
	}
	if len(cleaned.Populations) != 3 {
		t.Errorf("cleaned snapshot has %d populations, want 3 (empty node E should have been removed)", len(cleaned.Populations))
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
	// Keyed by mutation content ("s0"), not node id — see nodeMutKey.
	freqs := map[string][][]float64{nodeMutKey(map[string][]string{"ssms": {"s0"}}): {{0.7, 0.6}}}

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
	// [0.5] in tree B. Aggregate should record both rows, keyed by the
	// node's mutation set ("s0"), not its node id.
	a := json.RawMessage(`{"populations":{"0":{"cellular_prevalence":[0.7],"num_ssms":1,"num_cnvs":0}},"structure":{"0":[]},"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}}`)
	b := json.RawMessage(`{"populations":{"0":{"cellular_prevalence":[0.5],"num_ssms":1,"num_cnvs":0}},"structure":{"0":[]},"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}}`)
	freqs, err := aggregateFreqsByNode([]json.RawMessage{a, b})
	if err != nil {
		t.Fatal(err)
	}
	rows := freqs[nodeMutKey(map[string][]string{"ssms": {"s0"}})]
	if len(rows) != 2 {
		t.Fatalf("expected 2 rows for mutation set {s0}, got %d", len(rows))
	}
}

// TestAggregateFreqsByNode_KeyedByMutationContentNotNodeID is the M7
// regression test. Two samples share the same topology signature, but the
// two non-root children swap phi rank between samples — exactly the
// scenario that produces different node IDs for the same mutation set,
// since node IDs are assigned per-sample by descending cellular prevalence
// (see summarizePops). Keying by node ID would silently pool the WRONG
// rows together (an order statistic: "whichever sibling ranked higher this
// sample"); keying by mutation content must pool the RIGHT rows.
func TestAggregateFreqsByNode_KeyedByMutationContentNotNodeID(t *testing.T) {
	// Sample 1: node "1" (higher phi) holds s0, node "2" (lower phi) holds s1.
	sample1 := json.RawMessage(`{
		"populations":{
			"0":{"cellular_prevalence":[0.9],"num_ssms":0,"num_cnvs":0},
			"1":{"cellular_prevalence":[0.6],"num_ssms":1,"num_cnvs":0},
			"2":{"cellular_prevalence":[0.3],"num_ssms":1,"num_cnvs":0}
		},
		"structure":{"0":[1,2]},
		"mut_assignments":{
			"1":{"ssms":["s0"],"cnvs":[]},
			"2":{"ssms":["s1"],"cnvs":[]}
		}
	}`)
	// Sample 2: sibling ranks swapped — node "1" (higher phi) now holds s1,
	// node "2" (lower phi) now holds s0.
	sample2 := json.RawMessage(`{
		"populations":{
			"0":{"cellular_prevalence":[0.9],"num_ssms":0,"num_cnvs":0},
			"1":{"cellular_prevalence":[0.55],"num_ssms":1,"num_cnvs":0},
			"2":{"cellular_prevalence":[0.5],"num_ssms":1,"num_cnvs":0}
		},
		"structure":{"0":[1,2]},
		"mut_assignments":{
			"1":{"ssms":["s1"],"cnvs":[]},
			"2":{"ssms":["s0"],"cnvs":[]}
		}
	}`)

	freqs, err := aggregateFreqsByNode([]json.RawMessage{sample1, sample2})
	if err != nil {
		t.Fatal(err)
	}

	s0Key := nodeMutKey(map[string][]string{"ssms": {"s0"}})
	s1Key := nodeMutKey(map[string][]string{"ssms": {"s1"}})

	s0Rows := freqs[s0Key]
	if len(s0Rows) != 2 || s0Rows[0][0] != 0.6 || s0Rows[1][0] != 0.5 {
		t.Errorf("s0 rows = %v, want [[0.6] [0.5]] (its own values across both samples, not the order statistic)", s0Rows)
	}
	s1Rows := freqs[s1Key]
	if len(s1Rows) != 2 || s1Rows[0][0] != 0.3 || s1Rows[1][0] != 0.55 {
		t.Errorf("s1 rows = %v, want [[0.3] [0.55]] (its own values across both samples, not the order statistic)", s1Rows)
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

func TestPosteriorTreesSubcommand_DispatchedFromArgs(t *testing.T) {
	// Build args slice as if invoked from CLI.
	tmp := t.TempDir()
	zw, _ := newTreeArchiveWriter(filepath.Join(tmp, "trees.zip"))
	snap := json.RawMessage(`{"populations":{"0":{"cellular_prevalence":[0.7],"num_ssms":1,"num_cnvs":0}},"structure":{"0":[]},"mut_assignments":{"0":{"ssms":["s0"],"cnvs":[]}}}`)
	zw.WriteSample(0, -100.0, false, snap)
	zw.Close()

	// Direct call to the subcommand handler is enough — TestMain-style argv
	// rewriting is overkill for a unit test.
	if err := posteriorTreesSubcommand([]string{"-no-pdf", tmp}); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(filepath.Join(tmp, "posterior_trees")); err != nil {
		t.Errorf("posterior_trees dir missing: %v", err)
	}
}
