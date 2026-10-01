package main

import (
	"encoding/json"
	"math/rand"
	"testing"
)

// cnvOnlyTree builds root -> A (holds s0) -> [C (holds s1), B (holds only c0)].
// In Python, c0 lives in B's node.data like any SSM, so B is a data-bearing
// node for cull_tree, resample_sticks, resample_stick_orders and
// _summarize_pops.
func cnvOnlyTree(t *testing.T) (*TSSB, *TSSBNode, *TSSBNode, *rand.Rand) {
	t.Helper()
	rng := rand.New(rand.NewSource(1))
	ssms := []*SSM{
		{ID: "s0", A: []int{10}, D: []int{20}, MuR: 0.999, MuV: 0.5, LogBinNormConst: []float64{0}},
		{ID: "s1", A: []int{15}, D: []int{20}, MuR: 0.999, MuV: 0.5, LogBinNormConst: []float64{0}},
	}
	cnvs := []*CNV{{ID: "c0", A: []int{5}, D: []int{15}, LogBinNormConst: []float64{0}}}
	tssb := newTSSB(ssms, cnvs, 1, 1, 0.5, rng)
	a := tssb.Root.Children[0]
	c := tssb.spawnChild(a, 1, rng)
	b := tssb.spawnChild(a, 1, rng)
	a.Children = append(a.Children, c, b)
	a.Sticks = append(a.Sticks, 0.5, 0.5)
	a.Node.Data = []int{0}
	c.Node.Data = []int{1}
	tssb.Assignments[1] = c.Node
	ssms[1].Node = c.Node
	cnvs[0].Node = b.Node
	return tssb, a, b, rng
}

func TestCullTree_KeepsCNVOnlyNode(t *testing.T) {
	tssb, a, b, _ := cnvOnlyTree(t)
	tssb.cullTree()
	if len(a.Children) != 2 || a.Children[1] != b {
		t.Fatalf("CNV-only node was culled; A has %d children", len(a.Children))
	}
}

func TestResampleStickOrders_KeepsCNVOnlyNode(t *testing.T) {
	tssb, _, b, rng := cnvOnlyTree(t)
	tssb.resampleStickOrders(rng)
	for _, n := range tssb.getNodes() {
		if n == b.Node {
			return
		}
	}
	t.Fatal("CNV-only node pruned by resampleStickOrders; c0 now points outside the tree")
}

func TestCountData_IncludesCNVs(t *testing.T) {
	_, a, b, _ := cnvOnlyTree(t)
	cnvCount := map[*Node]int{b.Node: 1}
	if got := countData(a, cnvCount); got != 3 {
		t.Fatalf("countData(A) = %d, want 3 (s0 + s1 + c0)", got)
	}
}

func TestSummarizePops_EmitsCNVs(t *testing.T) {
	tssb, _, _, _ := cnvOnlyTree(t)
	raw, err := json.Marshal(summarizePops(tssb, 0, 0))
	if err != nil {
		t.Fatal(err)
	}
	var out struct {
		MutAss map[string]struct {
			CNVs []string `json:"cnvs"`
		} `json:"mut_assignments"`
	}
	if err := json.Unmarshal(raw, &out); err != nil {
		t.Fatal(err)
	}
	total := 0
	for _, m := range out.MutAss {
		total += len(m.CNVs)
	}
	if total != 1 {
		t.Fatalf("summarizePops assigned %d CNVs, want 1", total)
	}
}
