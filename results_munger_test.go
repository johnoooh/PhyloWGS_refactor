package main

import (
	"fmt"
	"reflect"
	"testing"
)

// mkNode describes one node of a hand-built tree: cellular prevalence and the
// SSM/CNV ids assigned to it (counts are derived from the id lists).
type mkNode struct {
	cp         []float64
	ssms, cnvs []string
}

// ssmIDs returns ids prefix0..prefix(n-1).
func ssmIDs(prefix string, from, n int) []string {
	out := make([]string, n)
	for i := range out {
		out[i] = fmt.Sprintf("%s%d", prefix, from+i)
	}
	return out
}

// mkTree builds a TreeSummary + MutAss from nodes (index = node id; node 0 is
// the root and gets no mutass entry, as Go's summarizePops emits).
func mkTree(structure map[int][]int, nodes []mkNode) (*TreeSummary, map[int]*MutAss) {
	ts := &TreeSummary{LLH: -100, Structure: structure, Populations: map[int]*Population{}}
	ma := map[int]*MutAss{}
	for i, n := range nodes {
		ts.Populations[i] = &Population{CellPrev: n.cp, NumSSMs: len(n.ssms), NumCNVs: len(n.cnvs)}
		if i > 0 {
			ma[i] = &MutAss{SSMs: n.ssms, CNVs: n.cnvs}
		}
	}
	return ts, ma
}

// mkMuts gives every SSM in ids default reads (ref 50 / total 100) and every
// CNV in cnvIDs the same, then applies overrides {id: {ref, total}}.
func mkMuts(ids, cnvIDs []string, over map[string][2]int) *MutList {
	ml := &MutList{SSMs: map[string]*MutListSSM{}, CNVReads: map[string][2][]int{}}
	reads := func(id string) (int, int) {
		if o, ok := over[id]; ok {
			return o[0], o[1]
		}
		return 50, 100
	}
	for _, id := range ids {
		r, t := reads(id)
		ml.SSMs[id] = &MutListSSM{RefReads: []int{r}, TotalReads: []int{t}}
	}
	for _, id := range cnvIDs {
		r, t := reads(id)
		ml.CNVReads[id] = [2][]int{{r}, {t}}
	}
	return ml
}

func oneTreeSet(structure map[int][]int, nodes []mkNode, muts *MutList) *ResultSet {
	ts, ma := mkTree(structure, nodes)
	return &ResultSet{Trees: map[int]*TreeSummary{0: ts}, MutAss: map[int]map[int]*MutAss{0: ma}, Muts: muts}
}

func concat(parts ...[]string) []string {
	out := []string{}
	for _, p := range parts {
		out = append(out, p...)
	}
	return out
}

func assertTree(t *testing.T, rs *ResultSet, tidx int, wantStruct map[int][]int, wantNodes []mkNode) {
	t.Helper()
	wantTS, wantMA := mkTree(wantStruct, wantNodes)
	got := rs.Trees[tidx]
	if !reflect.DeepEqual(got.Structure, wantTS.Structure) {
		t.Errorf("structure = %v, want %v", got.Structure, wantTS.Structure)
	}
	if len(got.Populations) != len(wantTS.Populations) {
		t.Errorf("populations keys = %d, want %d", len(got.Populations), len(wantTS.Populations))
	}
	for i, w := range wantTS.Populations {
		if g := got.Populations[i]; g == nil || !reflect.DeepEqual(*g, *w) {
			t.Errorf("population %d = %+v, want %+v", i, g, w)
		}
	}
	if !reflect.DeepEqual(rs.MutAss[tidx], wantMA) {
		for i := range rs.MutAss[tidx] {
			t.Logf("got mutass %d: %+v", i, rs.MutAss[tidx][i])
		}
		t.Errorf("mutass mismatch")
	}
}

// chain3 is 0 -> 1 -> 2 -> 3 with 20 / 1 / 10 SSMs; node 2's SSM is s20.
func chain3(cp1, cp2, cp3 float64, ref20 int) *ResultSet {
	all := ssmIDs("s", 0, 31)
	return oneTreeSet(map[int][]int{0: {1}, 1: {2}, 2: {3}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{cp1}, ssms: all[:20]},
		{cp: []float64{cp2}, ssms: all[20:21]},
		{cp: []float64{cp3}, ssms: all[21:]},
	}, mkMuts(all, nil, map[string][2]int{"s20": {ref20, 100}}))
}

func TestRemoveSmallNodes_MovesToClosestPhiAndRenumbers(t *testing.T) {
	rs := chain3(0.9, 0.5, 0.3, 60) // s20 implied phi 0.8 -> node 1 (delta .1 vs .5)
	if err := rs.RemoveSmallNodes(5); err != nil {
		t.Fatal(err)
	}
	all := ssmIDs("s", 0, 31)
	assertTree(t, rs, 0, map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.9}, ssms: concat(all[:20], []string{"s20"})},
		{cp: []float64{0.3}, ssms: all[21:]},
	})
}

func TestRemoveSmallNodes_TieGoesToLowerIndex(t *testing.T) {
	rs := chain3(0.75, 0.5, 0.25, 75) // implied phi 0.5: deltas exactly .25 and .25
	if err := rs.RemoveSmallNodes(5); err != nil {
		t.Fatal(err)
	}
	all := ssmIDs("s", 0, 31)
	assertTree(t, rs, 0, map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.75}, ssms: concat(all[:20], []string{"s20"})},
		{cp: []float64{0.25}, ssms: all[21:]},
	})
}

func TestRemoveSmallNodes_FractionRoundsToZero(t *testing.T) {
	all := ssmIDs("s", 0, 40)
	rs := oneTreeSet(map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.9}, ssms: all[:39]},
		{cp: []float64{0.5}, ssms: all[39:]},
	}, mkMuts(all, nil, nil))
	if err := rs.RemoveSmallNodes(0.01); err != nil { // round(0.4) = 0
		t.Fatal(err)
	}
	assertTree(t, rs, 0, map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.9}, ssms: all[:39]},
		{cp: []float64{0.5}, ssms: all[39:]},
	})
}

func TestRemoveSmallNodes_CountOfOne(t *testing.T) {
	rs := chain3(0.9, 0.5, 0.3, 60)
	// 1.0 is a count (keep nodes with >= 1 SSM), not a fraction of 31 SSMs.
	if err := rs.RemoveSmallNodes(1.0); err != nil {
		t.Fatal(err)
	}
	all := ssmIDs("s", 0, 31)
	assertTree(t, rs, 0, map[int][]int{0: {1}, 1: {2}, 2: {3}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.9}, ssms: all[:20]},
		{cp: []float64{0.5}, ssms: all[20:21]},
		{cp: []float64{0.3}, ssms: all[21:]},
	})
}

func TestRemoveSmallNodes_NoBestNode(t *testing.T) {
	all := ssmIDs("s", 0, 11)
	rs := oneTreeSet(map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.0}, ssms: all[:10]},
		{cp: []float64{0.0}, ssms: all[10:]},
	}, mkMuts(all, nil, map[string][2]int{"s10": {0, 100}})) // implied phi min(2,1) = 1
	err := rs.RemoveSmallNodes(5)
	want := "tree 0: no node within phi 1 of mutation s10"
	if err == nil || err.Error() != want {
		t.Fatalf("err = %v, want %q", err, want)
	}
}

func TestRemoveSmallNodes_CNVsMoveToCNVs(t *testing.T) {
	all := ssmIDs("s", 0, 30)
	rs := oneTreeSet(map[int][]int{0: {1}, 1: {2}, 2: {3}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.9}, ssms: all[:20], cnvs: []string{"c1"}},
		{cp: []float64{0.5}, cnvs: []string{"c0"}},
		{cp: []float64{0.3}, ssms: all[20:]},
	}, mkMuts(all, []string{"c0", "c1"}, map[string][2]int{"c0": {60, 100}}))
	if err := rs.RemoveSmallNodes(5); err != nil {
		t.Fatal(err)
	}
	assertTree(t, rs, 0, map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.9}, ssms: all[:20], cnvs: []string{"c1", "c0"}},
		{cp: []float64{0.3}, ssms: all[20:]},
	})
}

// Python's _correct_mut_counts indexes mutass[tree][node] for every non-root
// node; a node with no mutass entry raises KeyError there.
func TestRemoveSmallNodes_MissingMutAssEntryErrors(t *testing.T) {
	rs := chain3(0.9, 0.5, 0.3, 60)
	delete(rs.MutAss[0], 3)
	if err := rs.RemoveSmallNodes(1); err == nil {
		t.Fatal("expected error for node with no mutass entry")
	}
}
