package main

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestPyFloat_Marshal(t *testing.T) {
	a, b := 0.1, 0.2 // variables: constant 0.1+0.2 folds to exactly 0.3
	cases := []struct {
		in   float64
		want string
	}{
		{1.0, "1.0"},
		{0, "0.0"},
		{-1300, "-1300.0"},
		{0.5, "0.5"},
		{1e-05, "1e-05"},
		{0.0001, "0.0001"},
		{1e16, "1e+16"},
		{123456789012345.0, "123456789012345.0"},
		{a + b, "0.30000000000000004"},
		{math.NaN(), "NaN"},
		{math.Inf(1), "Infinity"},
		{math.Inf(-1), "-Infinity"},
	}
	for _, c := range cases {
		got, err := PyFloat(c.in).MarshalJSON()
		if err != nil {
			t.Fatal(err)
		}
		if string(got) != c.want {
			t.Errorf("PyFloat(%v) = %s, want %s", c.in, got, c.want)
		}
	}
}

func TestSnapshotTree_PythonTypes(t *testing.T) {
	tssb, _, _, _ := cnvOnlyTree(t)
	for _, n := range []*TSSBNode{tssb.Root} {
		n.Node.Params = []float64{1}
	}
	raw := string(snapshotTree(tssb, -1300, 0, 0, map[string]string{}))
	if !strings.Contains(raw, `"cellular_prevalence":[1.0]`) {
		t.Errorf("missing float cellular_prevalence: %s", raw)
	}
	if !strings.Contains(raw, `"llh":-1300.0`) {
		t.Errorf("missing float llh: %s", raw)
	}
	var g struct {
		Pops map[string]struct {
			N json.Number `json:"num_ssms"`
		} `json:"populations"`
	}
	if err := json.Unmarshal([]byte(raw), &g); err != nil {
		t.Fatal(err)
	}
	if len(g.Pops) == 0 {
		t.Fatal("no populations")
	}
	for k, p := range g.Pops {
		if strings.ContainsAny(p.N.String(), ".e") {
			t.Errorf("num_ssms of pop %s not an int: %s", k, p.N)
		}
	}
}

func TestWriteMutList_PythonTypes(t *testing.T) {
	dir := t.TempDir()
	ssms := []*SSM{{ID: "s0", A: []int{1}, D: []int{2}, MuR: 1, MuV: 0.5}}
	cnvs := []*CNV{{ID: "c0", A: []int{1}, D: []int{2},
		PhysicalCNVs: []PhysicalCNV{{Chrom: "1", Start: 1, End: 2, MajorCN: 1, MinorCN: 1, CellPrev: []float64{0, 0.718}}}}}
	if err := writeMutList(dir, ssms, cnvs, "d"); err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile(filepath.Join(dir, "mutlist.json"))
	if err != nil {
		t.Fatal(err)
	}
	s := regexpSpaces(string(b))
	for _, want := range []string{`"expected_ref_in_ref":1.0`, `"name":""`, `"cell_prev":[0.0,0.718]`} {
		if !strings.Contains(s, want) {
			t.Errorf("mutlist.json missing %s:\n%s", want, b)
		}
	}
}

func regexpSpaces(s string) string {
	return strings.NewReplacer(" ", "", "\n", "", "\t", "").Replace(s)
}

func TestPyFloat_NaNRejectedByEncodingJSON(t *testing.T) {
	if _, err := json.Marshal(PyFloat(math.NaN())); err == nil {
		t.Error("json.Marshal(NaN PyFloat) returned nil error")
	}
	if _, err := json.Marshal([]PyFloat{PyFloat(math.Inf(1))}); err == nil {
		t.Error("json.Marshal([+Inf PyFloat]) returned nil error")
	}
}
