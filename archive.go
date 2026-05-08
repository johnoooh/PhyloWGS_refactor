package main

import (
	"archive/zip"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// TreeArchiveWriter writes one zip entry per MCMC sample. Entries are
// named "tree_<idx>_<llh>" (or "burnin_<idx>_<llh>" during burn-in),
// matching the layout used by Python util2.TreeWriter so that
// util2.TreeReader can parse the index and LLH from the filename via
// float(). Filenames are NOT guaranteed byte-identical to Python's
// (Go's strconv.FormatFloat 'g'/-1 omits the trailing .0 on integral
// values) but they round-trip cleanly through float().
//
// Entry content is JSON (caller-provided RawMessage), not pickle, so
// downstream Go consumers don't need a Python runtime.
type TreeArchiveWriter struct {
	f  *os.File
	zw *zip.Writer
}

func newTreeArchiveWriter(path string) (*TreeArchiveWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	return &TreeArchiveWriter{f: f, zw: zip.NewWriter(f)}, nil
}

// WriteSample appends one entry. For non-burnin samples idx must be >= 0;
// burnin samples conventionally use negative idx in Python, but we expose
// the explicit isBurnin flag so the caller doesn't have to encode that
// in the index sign.
func (w *TreeArchiveWriter) WriteSample(idx int, llh float64, isBurnin bool, payload json.RawMessage) error {
	prefix := "tree"
	if isBurnin {
		prefix = "burnin"
	}
	// LLH is encoded with strconv.FormatFloat 'g'/-1, which round-trips
	// through Python's float() but may differ from Python's str(float)
	// for integral values (Go: "-1300", Python: "-1300.0").
	llhStr := strconv.FormatFloat(llh, 'g', -1, 64)
	name := fmt.Sprintf("%s_%d_%s", prefix, idx, llhStr)

	entry, err := w.zw.CreateHeader(&zip.FileHeader{
		Name:   name,
		Method: zip.Deflate,
	})
	if err != nil {
		return err
	}
	_, err = entry.Write(payload)
	return err
}

func (w *TreeArchiveWriter) Close() error {
	if err := w.zw.Close(); err != nil {
		w.f.Close()
		return err
	}
	return w.f.Close()
}

// MutassArchiveWriter mirrors pwgsresults.json_writer.JsonWriter.write_mutass:
// one zip entry per included tree named "<idx>.json" with content
// {"mut_assignments": {populationIdx: {"ssms": [...], "cnvs": [...]}}, "dataset_name": "<name>"}.
type MutassArchiveWriter struct {
	f       *os.File
	zw      *zip.Writer
	dataset string
}

func newMutassArchiveWriter(path, datasetName string) (*MutassArchiveWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	return &MutassArchiveWriter{
		f: f, zw: zip.NewWriter(f), dataset: datasetName,
	}, nil
}

// WriteTree records one tree's mutation assignments. The mutass map shape
// matches Python: {populationIdx: {"ssms": [...], "cnvs": [...]}}.
func (w *MutassArchiveWriter) WriteTree(idx int, mutass map[string]map[string][]string) error {
	entry, err := w.zw.CreateHeader(&zip.FileHeader{
		Name:   fmt.Sprintf("%d.json", idx),
		Method: zip.Deflate,
	})
	if err != nil {
		return err
	}
	payload := struct {
		MutAssignments map[string]map[string][]string `json:"mut_assignments"`
		DatasetName    string                         `json:"dataset_name"`
	}{mutass, w.dataset}
	return json.NewEncoder(entry).Encode(payload)
}

func (w *MutassArchiveWriter) Close() error {
	if err := w.zw.Close(); err != nil {
		w.f.Close()
		return err
	}
	return w.f.Close()
}

// TreeArchiveReader reads back zip entries written by TreeArchiveWriter.
// It exposes only non-burnin samples in idx-ascending order, matching
// Python TreeReader semantics.
type TreeArchiveReader struct {
	r       *zip.ReadCloser
	entries []treeEntry
}

type treeEntry struct {
	idx int
	llh float64
	zf  *zip.File
}

func newTreeArchiveReader(path string) (*TreeArchiveReader, error) {
	r, err := zip.OpenReader(path)
	if err != nil {
		return nil, err
	}
	var entries []treeEntry
	for _, f := range r.File {
		if !strings.HasPrefix(f.Name, "tree_") {
			continue
		}
		// "tree_<idx>_<llh>"
		parts := strings.SplitN(f.Name, "_", 3)
		if len(parts) != 3 {
			continue
		}
		idx, err1 := strconv.Atoi(parts[1])
		llh, err2 := strconv.ParseFloat(parts[2], 64)
		if err1 != nil || err2 != nil {
			continue
		}
		entries = append(entries, treeEntry{idx, llh, f})
	}
	sort.Slice(entries, func(i, j int) bool { return entries[i].idx < entries[j].idx })
	return &TreeArchiveReader{r: r, entries: entries}, nil
}

func (r *TreeArchiveReader) NumTrees() int { return len(r.entries) }

func (r *TreeArchiveReader) LoadTree(i int) (json.RawMessage, float64, error) {
	e := r.entries[i]
	rc, err := e.zf.Open()
	if err != nil {
		return nil, 0, err
	}
	defer rc.Close()
	b, err := io.ReadAll(rc)
	if err != nil {
		return nil, 0, err
	}
	return json.RawMessage(b), e.llh, nil
}

func (r *TreeArchiveReader) Close() error { return r.r.Close() }
