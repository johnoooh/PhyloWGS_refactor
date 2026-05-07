package main

import (
	"archive/zip"
	"encoding/json"
	"io"
	"path/filepath"
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
