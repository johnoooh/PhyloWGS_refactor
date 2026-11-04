package main

// Go port of write_results.py main() + pwgsresults/json_writer.py JsonWriter
// for phylowgs-go output (reference: testdata/native_results/write_results_go.py).
// Every float goes through PyFloat so decoded output matches Python 2's json.

import (
	"archive/zip"
	"bytes"
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
)

type popOut struct {
	CellPrev []PyFloat `json:"cellular_prevalence"`
	NumSSMs  int       `json:"num_ssms"`
	NumCNVs  int       `json:"num_cnvs"`
}

type treeOut struct {
	LLH         PyFloat        `json:"llh"`
	Structure   map[int][]int  `json:"structure"`
	Populations map[int]popOut `json:"populations"`
	Linearity   PyFloat        `json:"linearity_index"`
	Branching   PyFloat        `json:"branching_index"`
	Clustering  PyFloat        `json:"clustering_index"`
}

type summOut struct {
	DatasetName   string          `json:"dataset_name"`
	Params        struct{}        `json:"params"`
	Trees         map[int]treeOut `json:"trees"`
	TreeDensities map[int]PyFloat `json:"tree_densities"`
}

type ssmOut struct {
	RefReads        []int   `json:"ref_reads"`
	TotalReads      []int   `json:"total_reads"`
	ExpRefInRef     PyFloat `json:"expected_ref_in_ref"`
	ExpRefInVariant PyFloat `json:"expected_ref_in_variant"`
	Name            *string `json:"name,omitempty"`
}

type mutAssOut struct {
	SSMs []string `json:"ssms"`
	CNVs []string `json:"cnvs"`
}

type mutAssFileOut struct {
	MutAssignments map[int]mutAssOut `json:"mut_assignments"`
	DatasetName    string            `json:"dataset_name"`
}

func nonNil[T any](xs []T) []T {
	if xs == nil {
		return []T{}
	}
	return xs
}

// buildSummaries is JsonWriter.write_summaries' payload: per-tree indices
// plus calc_tree_densities over the remaining trees.
func buildSummaries(rs *ResultSet, dataset string) (*summOut, error) {
	idx, err := allTreeIndices(rs.Trees)
	if err != nil {
		return nil, err
	}
	dens, mode, err := treeDensities(rs.Trees, idx)
	if err != nil {
		return nil, err
	}
	log.Printf("write-results: tree densities computed in %s mode over %d trees", mode, len(dens))
	out := &summOut{DatasetName: dataset, Trees: map[int]treeOut{}, TreeDensities: map[int]PyFloat{}}
	for tidx, ts := range rs.Trees {
		st := make(map[int][]int, len(ts.Structure))
		for n, ch := range ts.Structure {
			st[n] = nonNil(ch)
		}
		pops := make(map[int]popOut, len(ts.Populations))
		for n, p := range ts.Populations {
			pops[n] = popOut{pyFloats(p.CellPrev), p.NumSSMs, p.NumCNVs}
		}
		v := idx[tidx]
		out.Trees[tidx] = treeOut{PyFloat(ts.LLH), st, pops, PyFloat(v[0]), PyFloat(v[1]), PyFloat(v[2])}
		out.TreeDensities[tidx] = PyFloat(dens[tidx])
	}
	return out, nil
}

// writeMutsJSON writes JsonWriter.write_mutlist's payload. name is emitted
// (even if empty) only when includeNames is set, as write_results.py does.
func writeMutsJSON(w io.Writer, ml *MutList, dataset string, includeNames bool) error {
	ssms := make(map[string]ssmOut, len(ml.SSMs))
	for id, s := range ml.SSMs {
		o := ssmOut{nonNil(s.RefReads), nonNil(s.TotalReads), PyFloat(s.ExpRefInRef), PyFloat(s.ExpRefInVariant), nil}
		if includeNames {
			name := s.Name
			o.Name = &name
		}
		ssms[id] = o
	}
	cnvs := ml.CNVs
	if cnvs == nil {
		cnvs = map[string]json.RawMessage{}
	}
	return json.NewEncoder(w).Encode(struct {
		SSMs        map[string]ssmOut          `json:"ssms"`
		CNVs        map[string]json.RawMessage `json:"cnvs"`
		DatasetName string                     `json:"dataset_name"`
	}{ssms, cnvs, dataset})
}

func gzipJSON(write func(io.Writer) error) ([]byte, error) {
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if err := write(zw); err != nil {
		return nil, err
	}
	if err := zw.Close(); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// buildMutassZip is JsonWriter.write_mutass: one deflated <tidx>.json per tree.
func buildMutassZip(rs *ResultSet, dataset string) ([]byte, error) {
	var buf bytes.Buffer
	zw := zip.NewWriter(&buf)
	for _, tidx := range sortedKeys(rs.MutAss) {
		ma := make(map[int]mutAssOut, len(rs.MutAss[tidx]))
		for n, m := range rs.MutAss[tidx] {
			ma[n] = mutAssOut{nonNil(m.SSMs), nonNil(m.CNVs)}
		}
		b, err := json.Marshal(mutAssFileOut{ma, dataset})
		if err != nil {
			return nil, fmt.Errorf("mutass tree %d: %w", tidx, err)
		}
		f, err := zw.CreateHeader(&zip.FileHeader{Name: strconv.Itoa(tidx) + ".json", Method: zip.Deflate})
		if err != nil {
			return nil, err
		}
		if _, err := f.Write(b); err != nil {
			return nil, err
		}
	}
	if err := zw.Close(); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// writePythonResults writes the three write_results.py outputs for an
// already-munged result set. All payloads are built before any file is
// written, so an error (e.g. a one-tree density failure) leaves no output.
func writePythonResults(rs *ResultSet, dataset string, includeSSMNames bool, summPath, mutsPath, mutassPath string) error {
	summ, err := buildSummaries(rs, dataset)
	if err != nil {
		return err
	}
	summGz, err := gzipJSON(func(w io.Writer) error { return json.NewEncoder(w).Encode(summ) })
	if err != nil {
		return fmt.Errorf("summaries: %w", err)
	}
	mutsGz, err := gzipJSON(func(w io.Writer) error { return writeMutsJSON(w, rs.Muts, dataset, includeSSMNames) })
	if err != nil {
		return fmt.Errorf("mutlist: %w", err)
	}
	mutass, err := buildMutassZip(rs, dataset)
	if err != nil {
		return err
	}
	for _, o := range []struct {
		path string
		data []byte
	}{{summPath, summGz}, {mutsPath, mutsGz}, {mutassPath, mutass}} {
		if err := os.WriteFile(o.path, o.data, 0o644); err != nil {
			return err
		}
	}
	return nil
}

const writeResultsUsage = "usage: phylowgs-go write-results [--include-ssm-names] [--min-ssms F] [--include-multiprimary] [--keep-superclones] [--max-multiprimary F] DATASET TREES_ZIP MUTLIST_JSON SUMM_OUT MUTS_OUT MUTASS_OUT"

// writeResultsSubcommand is write_results.py main() for phylowgs-go output.
func writeResultsSubcommand(args []string) error {
	fs := flag.NewFlagSet("write-results", flag.ContinueOnError)
	names := fs.Bool("include-ssm-names", false, "Include SSM names in the muts output")
	minSSMs := fs.Float64("min-ssms", 0.01, "Minimum number or percent of SSMs to retain a subclone")
	inclMP := fs.Bool("include-multiprimary", false, "Include multiprimary trees in result")
	keepSC := fs.Bool("keep-superclones", false, "Keep superclones, often artifacts of tree sampling")
	maxMP := fs.Float64("max-multiprimary", 0.8, "Maximum proportion of multiprimary trees if --include-multiprimary is not set")
	if err := fs.Parse(args); err != nil {
		return err
	}
	pos := fs.Args()
	if len(pos) != 6 {
		return fmt.Errorf("%s", writeResultsUsage)
	}
	if !(*maxMP >= 0 && *maxMP <= 1) {
		return fmt.Errorf("--max-multiprimary %v outside range [0, 1]", *maxMP)
	}
	rs, err := loadGoResults(pos[1], pos[2])
	if err != nil {
		return err
	}
	if err := rs.RemoveSmallNodes(*minSSMs); err != nil {
		return err
	}
	if !*keepSC {
		if err := rs.RemoveSuperclones(); err != nil {
			return err
		}
	}
	if !*inclMP {
		if err := rs.RemoveMultiprimaryTrees(*maxMP); err != nil {
			return err
		}
	}
	return writePythonResults(rs, pos[0], *names, pos[3], pos[4], pos[5])
}
