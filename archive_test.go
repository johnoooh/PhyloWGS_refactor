package main

import (
	"archive/zip"
	"encoding/json"
	"io"
	"path/filepath"
	"strings"
	"testing"
)

func TestTreeArchiveWriter_WritesNamedTreeAndBurninEntries(t *testing.T) {
	tmp := t.TempDir()
	zipPath := filepath.Join(tmp, "trees.zip")

	w, err := newTreeArchiveWriter(zipPath)
	if err != nil {
		t.Fatalf("newTreeArchiveWriter: %v", err)
	}

	if err := w.WriteSample(0, -1234.5, false, json.RawMessage(`{"iter":0}`)); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSample(1, -1230.1, false, json.RawMessage(`{"iter":1}`)); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSample(-5, -1300.0, true, json.RawMessage(`{"iter":-5}`)); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	r, err := zip.OpenReader(zipPath)
	if err != nil {
		t.Fatalf("zip.OpenReader: %v", err)
	}
	defer r.Close()

	gotNames := make(map[string]bool)
	for _, f := range r.File {
		gotNames[f.Name] = true
	}

	// Exact-name match — Go's 'g'/-1 omits trailing .0 on integral floats.
	for _, want := range []string{"tree_0_-1234.5", "tree_1_-1230.1", "burnin_-5_-1300"} {
		if !gotNames[want] {
			t.Errorf("missing entry %q in zip; got names: %v", want, gotNames)
		}
	}
}

func TestTreeArchiveWriter_ContentRoundTrips(t *testing.T) {
	tmp := t.TempDir()
	zipPath := filepath.Join(tmp, "trees.zip")
	w, err := newTreeArchiveWriter(zipPath)
	if err != nil {
		t.Fatal(err)
	}
	payload := json.RawMessage(`{"populations":{"0":{"num_ssms":3}},"llh":-99.5}`)
	if err := w.WriteSample(7, -99.5, false, payload); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	r, err := zip.OpenReader(zipPath)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if len(r.File) != 1 {
		t.Fatalf("expected 1 entry, got %d", len(r.File))
	}
	rc, err := r.File[0].Open()
	if err != nil {
		t.Fatal(err)
	}
	defer rc.Close()
	got, _ := io.ReadAll(rc)
	if string(got) != string(payload) {
		t.Errorf("content mismatch:\n got: %s\nwant: %s", got, payload)
	}
}

func TestMutassArchiveWriter_WritesPerTreeJSONEntries(t *testing.T) {
	tmp := t.TempDir()
	zipPath := filepath.Join(tmp, "mutass.zip")
	w, err := newMutassArchiveWriter(zipPath, "TestDataset")
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteTree(0, map[string]map[string][]string{
		"0": {"ssms": {"s1", "s2"}, "cnvs": {}},
		"1": {"ssms": {"s3"}, "cnvs": {"c1"}},
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteTree(1, map[string]map[string][]string{
		"0": {"ssms": {"s1"}, "cnvs": {}},
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	r, err := zip.OpenReader(zipPath)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	got := map[string]string{}
	for _, f := range r.File {
		rc, _ := f.Open()
		b, _ := io.ReadAll(rc)
		rc.Close()
		got[f.Name] = string(b)
	}
	if _, ok := got["0.json"]; !ok {
		t.Errorf("missing 0.json; entries: %v", got)
	}
	if _, ok := got["1.json"]; !ok {
		t.Errorf("missing 1.json")
	}
	// minimal sanity: dataset_name embedded
	if !strings.Contains(got["0.json"], `"dataset_name":"TestDataset"`) {
		t.Errorf("0.json missing dataset_name field: %s", got["0.json"])
	}
	if !strings.Contains(got["0.json"], `"mut_assignments"`) {
		t.Errorf("0.json missing mut_assignments wrapper: %s", got["0.json"])
	}
}
