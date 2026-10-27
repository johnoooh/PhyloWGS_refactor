package main

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// writeChainFiles must produce exactly what writeResults writes for the same
// chain, leave no .tmp files, and writeResults must not rewrite a chain whose
// files were already written as it finished.
func TestWriteChainFiles_MatchesWriteResultsAndIsSkipped(t *testing.T) {
	ssms, cnvs := loadTinyDataset(t)
	r := runChain(0, ssms, cnvs, 5, 10, 50, 7, 0)
	if len(r.Trees) == 0 {
		t.Fatal("chain produced no trees")
	}
	names := []string{"chain_0_samples.txt", "chain_0_trees.ndjson"}

	early := t.TempDir()
	if err := writeChainFiles(early, r); err != nil {
		t.Fatal(err)
	}
	late := t.TempDir()
	if err := writeResults(late, []ChainResult{r}, 1.1, "T"); err != nil {
		t.Fatal(err)
	}
	for _, n := range names {
		a, err := os.ReadFile(filepath.Join(early, n))
		if err != nil {
			t.Fatal(err)
		}
		b, err := os.ReadFile(filepath.Join(late, n))
		if err != nil {
			t.Fatal(err)
		}
		if len(a) == 0 || !bytes.Equal(a, b) {
			t.Fatalf("%s: early write (%d bytes) differs from writeResults (%d bytes)", n, len(a), len(b))
		}
	}
	if tmp, _ := filepath.Glob(filepath.Join(early, "*.tmp")); len(tmp) != 0 {
		t.Fatalf("leftover temp files: %v", tmp)
	}

	// Already written: writeResults must leave the files alone.
	sentinel := []byte("sentinel\n")
	for _, n := range names {
		if err := os.WriteFile(filepath.Join(early, n), sentinel, 0644); err != nil {
			t.Fatal(err)
		}
	}
	r.ChainFilesWritten = true
	if err := writeResults(early, []ChainResult{r}, 1.1, "T"); err != nil {
		t.Fatal(err)
	}
	for _, n := range names {
		if got, _ := os.ReadFile(filepath.Join(early, n)); !bytes.Equal(got, sentinel) {
			t.Fatalf("%s was rewritten although ChainFilesWritten was set", n)
		}
	}
}
