package main

import (
	"archive/zip"
	"encoding/json"
	"fmt"
	"os"
	"strconv"
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
