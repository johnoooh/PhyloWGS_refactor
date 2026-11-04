package main

import (
	"archive/zip"
	"path/filepath"
	"strings"
	"testing"
)

func TestLoadGoResults_Sim(t *testing.T) {
	dir := filepath.Join("testdata", "native_results", "sim")
	rs, err := loadGoResults(filepath.Join(dir, "trees.zip"), filepath.Join(dir, "mutlist.json"))
	if err != nil {
		t.Fatalf("loadGoResults: %v", err)
	}
	zr, err := zip.OpenReader(filepath.Join(dir, "trees.zip"))
	if err != nil {
		t.Fatal(err)
	}
	defer zr.Close()
	want := 0
	for _, f := range zr.File {
		if strings.HasPrefix(f.Name, "tree_") {
			want++
		}
	}
	if want == 0 || len(rs.Trees) != want {
		t.Fatalf("trees = %d, want %d (>0)", len(rs.Trees), want)
	}
	for idx, tr := range rs.Trees {
		if tr.Populations[0] == nil || len(tr.Populations[0].CellPrev) == 0 {
			t.Errorf("tree %d: population 0 missing/empty", idx)
		}
		if _, ok := rs.MutAss[idx]; !ok {
			t.Errorf("tree %d: no MutAss entry", idx)
		}
	}
	if len(rs.MutAss) != len(rs.Trees) {
		t.Errorf("MutAss %d != Trees %d", len(rs.MutAss), len(rs.Trees))
	}
	if len(rs.Muts.SSMs) == 0 || rs.Muts.SSMs["s0"] == nil || rs.Muts.SSMs["s0"].Name != "SIM_0" {
		t.Errorf("ssm s0 not loaded: %+v", rs.Muts.SSMs["s0"])
	}
	if len(rs.Muts.CNVs) == 0 || len(rs.Muts.CNVReads) != len(rs.Muts.CNVs) {
		t.Errorf("cnvs=%d cnvreads=%d", len(rs.Muts.CNVs), len(rs.Muts.CNVReads))
	}
	for id, r := range rs.Muts.CNVReads {
		if len(r[0]) == 0 || len(r[1]) == 0 {
			t.Errorf("cnv %s: empty reads", id)
		}
	}
}

func TestLoadGoResults_NoTrees(t *testing.T) {
	dir := t.TempDir()
	zp := filepath.Join(dir, "trees.zip")
	w, err := newTreeArchiveWriter(zp)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSample(0, -1.5, true, []byte(`{}`)); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	_, err = loadGoResults(zp, filepath.Join("testdata", "native_results", "sim", "mutlist.json"))
	if err == nil || !strings.Contains(err.Error(), "no trees") {
		t.Fatalf("err = %v, want 'no trees'", err)
	}
}

func TestLoadGoResults_MissingMutlist(t *testing.T) {
	zp := filepath.Join("testdata", "native_results", "sim", "trees.zip")
	if _, err := loadGoResults(zp, filepath.Join(t.TempDir(), "nope.json")); err == nil {
		t.Fatal("expected error for missing mutlist")
	}
}

func TestIntKeys_NonInteger(t *testing.T) {
	if _, err := intKeys(map[string][]int{"a": {1}}); err == nil {
		t.Fatal("expected error on non-integer key")
	}
}

func TestLoadGoResults_AllFixtures(t *testing.T) {
	for _, c := range []string{"sim", "linear", "identical", "single", "superclone", "cnvonly", "multiprimary"} {
		dir := filepath.Join("testdata", "native_results", c)
		rs, err := loadGoResults(filepath.Join(dir, "trees.zip"), filepath.Join(dir, "mutlist.json"))
		if err != nil || len(rs.Trees) == 0 {
			t.Errorf("%s: err=%v", c, err)
		}
	}
}
