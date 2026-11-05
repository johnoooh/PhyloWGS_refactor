package main

import (
	"archive/zip"
	"encoding/json"
	"os"
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

func writeTestZip(t *testing.T, entries [][2]string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "trees.zip")
	f, err := os.Create(p)
	if err != nil {
		t.Fatal(err)
	}
	zw := zip.NewWriter(f)
	for _, e := range entries {
		w, err := zw.Create(e[0])
		if err != nil {
			t.Fatal(err)
		}
		w.Write([]byte(e[1]))
	}
	zw.Close()
	f.Close()
	return p
}

const okTree = `{"llh":-1.5,"structure":{"0":[]},"populations":{"0":{"cellular_prevalence":[1],"num_ssms":0,"num_cnvs":0}}}`

func loadWithSimMutlist(t *testing.T, zp string) error {
	_, err := loadGoResults(zp, filepath.Join("testdata", "native_results", "sim", "mutlist.json"))
	return err
}

func TestLoadGoResults_DuplicateIdx(t *testing.T) {
	zp := writeTestZip(t, [][2]string{{"tree_5_-1.5", okTree}, {"tree_5_-2.5", okTree}})
	if err := loadWithSimMutlist(t, zp); err == nil || !strings.Contains(err.Error(), "duplicate") {
		t.Fatalf("err = %v, want duplicate", err)
	}
}

func TestLoadGoResults_MalformedName(t *testing.T) {
	zp := writeTestZip(t, [][2]string{{"tree_0_-1.5", okTree}, {"tree_x_-2.5", okTree}, {"tree_3", okTree}})
	if err := loadWithSimMutlist(t, zp); err == nil || !strings.Contains(err.Error(), "2 unparseable") {
		t.Fatalf("err = %v, want 2 unparseable", err)
	}
}

func TestLoadGoResults_MissingKey(t *testing.T) {
	cases := map[string]string{
		"llh":         `{"structure":{},"populations":{}}`,
		"structure":   `{"llh":-1,"populations":{}}`,
		"populations": `{"llh":-1,"structure":{}}`,
	}
	for key, body := range cases {
		zp := writeTestZip(t, [][2]string{{"tree_7_-1", body}})
		err := loadWithSimMutlist(t, zp)
		if err == nil || !strings.Contains(err.Error(), key) || !strings.Contains(err.Error(), "tree 7") {
			t.Errorf("%s: err = %v", key, err)
		}
	}
	// present-but-empty structure and absent mut_assignments are fine
	zp := writeTestZip(t, [][2]string{{"tree_7_-1", `{"llh":-1,"structure":{},"populations":{}}`}})
	if err := loadWithSimMutlist(t, zp); err != nil {
		t.Errorf("empty structure should load: %v", err)
	}
}

func writeMutlist(t *testing.T, body string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "mutlist.json")
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestLoadGoResults_MutlistMissingSSMs(t *testing.T) {
	zp := filepath.Join("testdata", "native_results", "sim", "trees.zip")
	_, err := loadGoResults(zp, writeMutlist(t, `{}`))
	if err == nil || !strings.Contains(err.Error(), `missing key "ssms"`) {
		t.Fatalf("err = %v, want missing key \"ssms\"", err)
	}
	_, err = loadGoResults(zp, writeMutlist(t, `{"ssms":{}}`))
	if err == nil || !strings.Contains(err.Error(), `missing key "cnvs"`) {
		t.Fatalf("err = %v, want missing key \"cnvs\"", err)
	}
}

func TestLoadGoResults_MutlistWrongSample(t *testing.T) {
	dir := filepath.Join("testdata", "native_results", "sim")
	b, err := os.ReadFile(filepath.Join(dir, "mutlist.json"))
	if err != nil {
		t.Fatal(err)
	}
	var m map[string]json.RawMessage
	if err := json.Unmarshal(b, &m); err != nil {
		t.Fatal(err)
	}
	var ssms map[string]json.RawMessage
	if err := json.Unmarshal(m["ssms"], &ssms); err != nil {
		t.Fatal(err)
	}
	delete(ssms, "s0")
	m["ssms"], _ = json.Marshal(ssms)
	nb, _ := json.Marshal(m)
	_, err = loadGoResults(filepath.Join(dir, "trees.zip"), writeMutlist(t, string(nb)))
	if err == nil || !strings.Contains(err.Error(), "not in mutlist") {
		t.Fatalf("err = %v, want 'not in mutlist'", err)
	}
}
