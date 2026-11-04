package main

import (
	"archive/zip"
	"bytes"
	"compress/gzip"
	"encoding/json"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// goldenFlagSets maps a fixture flag directory to the flags regen.sh passed
// to the Python reference (which always got --include-ssm-names except for
// nonames, whose expected output is default minus ssms/*/name).
var goldenFlagSets = map[string][]string{
	"default":  {"--include-ssm-names"},
	"prod":     {"--max-multiprimary", "1.0", "--include-multiprimary", "--include-ssm-names"},
	"keepsc":   {"--keep-superclones", "--include-ssm-names"},
	"mincount": {"--min-ssms", "5", "--include-ssm-names"},
	"nonames":  {},
}

// minGoldenDirs is the number of case x flag fixture directories regen.sh
// writes; fewer found means the walk is broken and the test would be vacuous.
const minGoldenDirs = 17

func gunzipFile(t *testing.T, path string) []byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	zr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("%s: %v", path, err)
	}
	b, err := io.ReadAll(zr)
	if err != nil {
		t.Fatalf("%s: %v", path, err)
	}
	return b
}

// unzipMutass decodes a mutass.zip into {"<tidx>": entry}, as dump_expected.py does.
func unzipMutass(t *testing.T, path string) []byte {
	t.Helper()
	zr, err := zip.OpenReader(path)
	if err != nil {
		t.Fatalf("%s: %v", path, err)
	}
	defer zr.Close()
	out := map[string]json.RawMessage{}
	for _, f := range zr.File {
		if f.Method != zip.Deflate {
			t.Errorf("%s/%s: method %d, want deflate", path, f.Name, f.Method)
		}
		rc, err := f.Open()
		if err != nil {
			t.Fatal(err)
		}
		b, err := io.ReadAll(rc)
		rc.Close()
		if err != nil {
			t.Fatal(err)
		}
		out[strings.TrimSuffix(f.Name, ".json")] = b
	}
	b, err := json.Marshal(out)
	if err != nil {
		t.Fatal(err)
	}
	return b
}

func findGoldenDirs(t *testing.T, root string) []string {
	t.Helper()
	var dirs []string
	err := filepath.WalkDir(root, func(p string, d os.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() && d.Name() == "expected" {
			dirs = append(dirs, filepath.Dir(p))
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	sort.Strings(dirs)
	if len(dirs) < minGoldenDirs {
		t.Fatalf("found %d golden dirs under %s, want >= %d", len(dirs), root, minGoldenDirs)
	}
	return dirs
}

func checkGoldenDir(t *testing.T, flagDir string) {
	caseDir, fs := filepath.Dir(flagDir), filepath.Base(flagDir)
	flags, ok := goldenFlagSets[fs]
	if !ok {
		t.Fatalf("unknown flag set %q", fs)
	}
	out := t.TempDir()
	summ, muts, mutass := filepath.Join(out, "summ.json.gz"), filepath.Join(out, "muts.json.gz"), filepath.Join(out, "mutass.zip")
	args := append(append([]string(nil), flags...), filepath.Base(caseDir),
		filepath.Join(caseDir, "trees.zip"), filepath.Join(caseDir, "mutlist.json"), summ, muts, mutass)
	err := writeResultsSubcommand(args)
	exp := filepath.Join(flagDir, "expected")
	if want, rerr := os.ReadFile(filepath.Join(exp, "error.txt")); rerr == nil {
		w := strings.TrimSpace(string(want))
		if w == "" {
			t.Fatal("empty error.txt")
		}
		if err == nil || !strings.HasSuffix(err.Error(), w) {
			t.Fatalf("want error ending %q, got %v", w, err)
		}
		return
	}
	if err != nil {
		t.Fatalf("writeResultsSubcommand: %v", err)
	}
	for _, c := range []struct {
		file string
		got  []byte
	}{
		{"summ.json", gunzipFile(t, summ)},
		{"muts.json", gunzipFile(t, muts)},
		{"mutass.json", unzipMutass(t, mutass)},
	} {
		want, rerr := os.ReadFile(filepath.Join(exp, c.file))
		if rerr != nil || len(want) == 0 {
			t.Fatalf("%s: missing or empty expected (%v)", c.file, rerr)
		}
		comparePythonJSON(t, c.file, c.got, want)
	}
}

func TestWriteResults_MatchesPythonReference(t *testing.T) {
	root := filepath.Join("testdata", "native_results")
	dirs := findGoldenDirs(t, root)
	for _, d := range dirs {
		rel, _ := filepath.Rel(root, d)
		t.Run(rel, func(t *testing.T) { checkGoldenDir(t, d) })
	}
	t.Logf("compared %d golden dirs", len(dirs))
}

func TestWriteResultsSubcommand_BadArgs(t *testing.T) {
	c := filepath.Join("testdata", "native_results", "sim")
	out := t.TempDir()
	pos := []string{"sim", filepath.Join(c, "trees.zip"), filepath.Join(c, "mutlist.json"),
		filepath.Join(out, "s.gz"), filepath.Join(out, "m.gz"), filepath.Join(out, "a.zip")}
	if err := writeResultsSubcommand(pos[:5]); err == nil || !strings.Contains(err.Error(), "usage") {
		t.Errorf("5 positionals: want usage error, got %v", err)
	}
	for _, v := range []string{"1.5", "-0.1"} {
		if err := writeResultsSubcommand(append([]string{"--max-multiprimary", v}, pos...)); err == nil {
			t.Errorf("--max-multiprimary %s: want error", v)
		}
	}
	if _, err := os.Stat(filepath.Join(out, "s.gz")); err == nil {
		t.Error("output written despite bad args")
	}
}

// TestWritePythonResults_NoNamesOmitsName guards the --include-ssm-names
// switch directly: name present (even "") only when requested.
func TestWritePythonResults_NoNamesOmitsName(t *testing.T) {
	for _, include := range []bool{true, false} {
		rs := &ResultSet{
			Trees: map[int]*TreeSummary{}, MutAss: map[int]map[int]*MutAss{},
			Muts: &MutList{SSMs: map[string]*MutListSSM{"s0": {RefReads: []int{1}, TotalReads: []int{2}, ExpRefInRef: 0.999, ExpRefInVariant: 0.5}}, CNVs: map[string]json.RawMessage{}},
		}
		var buf bytes.Buffer
		if err := writeMutsJSON(&buf, rs.Muts, "d", include); err != nil {
			t.Fatal(err)
		}
		has := strings.Contains(buf.String(), `"name":""`)
		if has != include {
			t.Errorf("include=%v: name present=%v in %s", include, has, buf.String())
		}
	}
}
