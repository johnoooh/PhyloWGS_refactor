package main

import (
	"encoding/json"
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
