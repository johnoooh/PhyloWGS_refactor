package main

import (
	"reflect"
	"testing"
)

// superTree is 0 -> 1 -> 2 -> 3 with clone/child SSM counts n1/n2 (child's
// SSMs only if childHasSSMs, else one CNV) and node 3 holding 5 SSMs.
func superTree(n1, n2 int, childHasSSMs bool) *ResultSet {
	all := ssmIDs("s", 0, n1+n2+5)
	child := mkNode{cp: []float64{0.9, 0.85}, ssms: all[n1 : n1+n2]}
	if !childHasSSMs {
		child = mkNode{cp: []float64{0.9, 0.85}, cnvs: []string{"c0"}}
	}
	return oneTreeSet(map[int][]int{0: {1}, 1: {2}, 2: {3}}, []mkNode{
		{cp: []float64{1, 1}},
		{cp: []float64{0.95, 0.8}, ssms: all[:n1]},
		child,
		{cp: []float64{0.5, 0.5}, ssms: all[n1+n2:]},
	}, mkMuts(all, []string{"c0"}, nil))
}

func TestRemoveSuperclones_Merges(t *testing.T) {
	rs := superTree(10, 40, true)
	if err := rs.RemoveSuperclones(); err != nil {
		t.Fatal(err)
	}
	all := ssmIDs("s", 0, 55)
	// Python: (cp1*10 + cp2*40)/50 elementwise -> [0.91, 0.84] (checked with numpy).
	assertTree(t, rs, 0, map[int][]int{0: {1}, 1: {2}}, []mkNode{
		{cp: []float64{1, 1}},
		{cp: []float64{0.91, 0.84}, ssms: all[:50]},
		{cp: []float64{0.5, 0.5}, ssms: all[50:]},
	})
}

func assertUnchanged(t *testing.T, build func() *ResultSet) {
	t.Helper()
	rs, orig := build(), build()
	if err := rs.RemoveSuperclones(); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(rs, orig) {
		t.Errorf("tree changed but should have been skipped")
	}
}

func TestRemoveSuperclones_SkipsCNVOnlyChild(t *testing.T) {
	assertUnchanged(t, func() *ResultSet { return superTree(10, 0, false) })
}

func TestRemoveSuperclones_SkipsRatioAbove033(t *testing.T) {
	assertUnchanged(t, func() *ResultSet { return superTree(20, 40, true) })
}

func multiprimaryTree() (*TreeSummary, map[int]*MutAss) {
	return mkTree(map[int][]int{0: {1, 2}}, []mkNode{
		{cp: []float64{1}},
		{cp: []float64{0.5}, ssms: []string{"s0"}},
		{cp: []float64{0.4}, ssms: []string{"s1"}},
	})
}

func TestRemoveSuperclones_SkipsMultiprimary(t *testing.T) {
	assertUnchanged(t, func() *ResultSet {
		ts, ma := multiprimaryTree()
		return &ResultSet{Trees: map[int]*TreeSummary{0: ts}, MutAss: map[int]map[int]*MutAss{0: ma},
			Muts: mkMuts([]string{"s0", "s1"}, nil, nil)}
	})
}

// fiveTrees builds trees 0..4 with LLH -idx; those in multi are multiprimary.
func fiveTrees(multi map[int]bool) *ResultSet {
	rs := &ResultSet{Trees: map[int]*TreeSummary{}, MutAss: map[int]map[int]*MutAss{}}
	for i := 0; i < 5; i++ {
		ts, ma := mkTree(map[int][]int{0: {1}}, []mkNode{{cp: []float64{1}}, {cp: []float64{0.5}, ssms: []string{"s0"}}})
		if multi[i] {
			ts, ma = multiprimaryTree()
		}
		ts.LLH = -float64(i)
		rs.Trees[i], rs.MutAss[i] = ts, ma
	}
	return rs
}

// The brief expected trees 2 and 4 to be renumbered to 1 and 2. The Python
// reference does NOT do that: it deletes the multiprimary trees before its
// renumbering loop, so `tidx in polyidxs` is never true there and the loop is
// a no-op. Verified by running pwgsresults.result_munger on this exact case:
// remaining keys [0, 2, 4].
func TestRemoveMultiprimary_KeepsOriginalIndices(t *testing.T) {
	rs := fiveTrees(map[int]bool{1: true, 3: true})
	if err := rs.RemoveMultiprimaryTrees(0.8); err != nil {
		t.Fatal(err)
	}
	for _, k := range []int{0, 2, 4} {
		if rs.Trees[k] == nil || rs.Trees[k].LLH != -float64(k) || rs.MutAss[k] == nil {
			t.Errorf("tree %d missing or moved", k)
		}
	}
	if len(rs.Trees) != 3 || len(rs.MutAss) != 3 {
		t.Errorf("got %d trees / %d mutass, want 3 / 3", len(rs.Trees), len(rs.MutAss))
	}
}

func TestRemoveMultiprimary_RaisesAtThreshold(t *testing.T) {
	rs := fiveTrees(map[int]bool{0: true, 1: true, 2: true, 3: true})
	err := rs.RemoveMultiprimaryTrees(0.8)
	want := "80% of trees are multiprimary (4 of 5), so not enough to report good posterior."
	if err == nil || err.Error() != want {
		t.Fatalf("err = %v, want %q", err, want)
	}
}
