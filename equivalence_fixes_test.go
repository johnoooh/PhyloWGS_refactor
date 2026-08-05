package main

import (
	"encoding/json"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

// M5: parsePhysicalCNVs must parse cell_prev into a []float64, mirroring
// Python util2.py:71 ([float(C) for C in cell_prev.split('|')]).
func TestParsePhysicalCNVs_CellPrevFloats(t *testing.T) {
	raw := "chrom=1,start=100,end=200,major_cn=2,minor_cn=0,cell_prev=0.0|0.718"
	got := parsePhysicalCNVs(raw)
	if len(got) != 1 {
		t.Fatalf("expected 1 physical CNV, got %d", len(got))
	}
	want := []float64{0.0, 0.718}
	if !reflect.DeepEqual(got[0].CellPrev, want) {
		t.Errorf("CellPrev = %v, want %v", got[0].CellPrev, want)
	}
}

// M5 + b1: mutlist.json must emit physical_cnvs[].cell_prev as a JSON array of
// numbers (not a string) and a top-level dataset_name, matching the Python
// json_writer schema.
func TestWriteMutList_CellPrevArrayAndDatasetName(t *testing.T) {
	dir := t.TempDir()
	ssms := []*SSM{{ID: "s0", A: []int{10}, D: []int{20}, MuR: 0.999, MuV: 0.5}}
	cnvs := []*CNV{{
		ID: "c0", A: []int{5}, D: []int{15},
		PhysicalCNVs: []PhysicalCNV{{Chrom: "1", Start: 100, End: 200, MajorCN: 2, MinorCN: 0, CellPrev: []float64{0.0, 0.718}}},
	}}

	if err := writeMutList(dir, ssms, cnvs, "DS_X"); err != nil {
		t.Fatalf("writeMutList: %v", err)
	}
	b, err := os.ReadFile(filepath.Join(dir, "mutlist.json"))
	if err != nil {
		t.Fatalf("read mutlist.json: %v", err)
	}

	var parsed map[string]json.RawMessage
	if err := json.Unmarshal(b, &parsed); err != nil {
		t.Fatalf("unmarshal: %v", err)
	}

	var ds string
	if err := json.Unmarshal(parsed["dataset_name"], &ds); err != nil {
		t.Fatalf("dataset_name not a string: %v", err)
	}
	if ds != "DS_X" {
		t.Errorf("dataset_name = %q, want %q", ds, "DS_X")
	}

	// Drill into cnvs.c0.physical_cnvs[0].cell_prev and confirm it decodes as []float64.
	var cnvSection map[string]struct {
		PhysicalCNVs []struct {
			CellPrev []float64 `json:"cell_prev"`
		} `json:"physical_cnvs"`
	}
	if err := json.Unmarshal(parsed["cnvs"], &cnvSection); err != nil {
		t.Fatalf("cell_prev is not a JSON array of numbers: %v", err)
	}
	got := cnvSection["c0"].PhysicalCNVs[0].CellPrev
	want := []float64{0.0, 0.718}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("cell_prev = %v, want %v", got, want)
	}
}
