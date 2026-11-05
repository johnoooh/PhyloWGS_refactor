package main

import (
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Population is one node's summary in a sampled tree.
type Population struct {
	CellPrev         []float64
	NumSSMs, NumCNVs int
}

// MutAss lists the SSM and CNV ids assigned to a node.
type MutAss struct{ SSMs, CNVs []string }

// TreeSummary is the per-tree data the Python ResultMunger consumes.
type TreeSummary struct {
	LLH         float64
	Structure   map[int][]int
	Populations map[int]*Population
}

// MutListSSM is one SSM entry of mutlist.json.
type MutListSSM struct {
	RefReads, TotalReads         []int
	ExpRefInRef, ExpRefInVariant float64
	Name                         string
}

// MutList is mutlist.json. CNVs are passed through untouched; CNVReads holds
// {ref_reads, total_reads} per CNV id for the munger's best-node search.
type MutList struct {
	SSMs     map[string]*MutListSSM
	CNVs     map[string]json.RawMessage
	CNVReads map[string][2][]int
}

// ResultSet is everything the result writer needs.
type ResultSet struct {
	Trees  map[int]*TreeSummary
	MutAss map[int]map[int]*MutAss
	Muts   *MutList
}

type treeJSON struct {
	LLH         *float64              `json:"llh"`
	Structure   map[string][]int      `json:"structure"`
	Populations map[string]popJSON    `json:"populations"`
	MutAss      map[string]mutAssJSON `json:"mut_assignments"`
}

type popJSON struct {
	CellPrev []float64 `json:"cellular_prevalence"`
	NumSSMs  int       `json:"num_ssms"`
	NumCNVs  int       `json:"num_cnvs"`
}

type mutAssJSON struct {
	SSMs []string `json:"ssms"`
	CNVs []string `json:"cnvs"`
}

// intKeys converts JSON string node keys to ints, erroring on non-integers.
func intKeys[V any](m map[string]V) (map[int]V, error) {
	out := make(map[int]V, len(m))
	for k, v := range m {
		i, err := strconv.Atoi(k)
		if err != nil {
			return nil, fmt.Errorf("non-integer node key %q", k)
		}
		out[i] = v
	}
	return out, nil
}

func (t *treeJSON) toSummary() (*TreeSummary, map[int]*MutAss, error) {
	for key, absent := range map[string]bool{"llh": t.LLH == nil, "structure": t.Structure == nil, "populations": t.Populations == nil} {
		if absent {
			return nil, nil, fmt.Errorf("missing key %q", key)
		}
	}
	structure, err := intKeys(t.Structure)
	if err != nil {
		return nil, nil, fmt.Errorf("structure: %w", err)
	}
	pops, err := intKeys(t.Populations)
	if err != nil {
		return nil, nil, fmt.Errorf("populations: %w", err)
	}
	ma, err := intKeys(t.MutAss)
	if err != nil {
		return nil, nil, fmt.Errorf("mut_assignments: %w", err)
	}
	ts := &TreeSummary{LLH: *t.LLH, Structure: structure, Populations: map[int]*Population{}}
	for n, p := range pops {
		ts.Populations[n] = &Population{CellPrev: p.CellPrev, NumSSMs: p.NumSSMs, NumCNVs: p.NumCNVs}
	}
	mas := make(map[int]*MutAss, len(ma))
	for n, m := range ma {
		mas[n] = &MutAss{SSMs: m.SSMs, CNVs: m.CNVs}
	}
	return ts, mas, nil
}

// loadGoResults reads Go's trees.zip and mutlist.json, mirroring the Python
// reference loader (write_results_go.py).
func loadGoResults(treesZip, mutlistPath string) (*ResultSet, error) {
	muts, err := loadMutList(mutlistPath)
	if err != nil {
		return nil, err
	}
	r, err := newTreeArchiveReader(treesZip)
	if err != nil {
		return nil, fmt.Errorf("open %s: %w", treesZip, err)
	}
	defer r.Close()
	if r.NumTrees() == 0 {
		return nil, fmt.Errorf("no trees in %s", treesZip)
	}
	if bad := countUnparseableTrees(r); bad > 0 {
		return nil, fmt.Errorf("%d unparseable tree_ entries in %s", bad, treesZip)
	}
	rs := &ResultSet{Trees: map[int]*TreeSummary{}, MutAss: map[int]map[int]*MutAss{}, Muts: muts}
	for i := 0; i < r.NumTrees(); i++ {
		idx := r.entries[i].idx
		if _, dup := rs.Trees[idx]; dup {
			return nil, fmt.Errorf("duplicate tree index %d in %s", idx, treesZip)
		}
		raw, _, err := r.LoadTree(i)
		if err != nil {
			return nil, fmt.Errorf("tree %d: %w", idx, err)
		}
		var tj treeJSON
		if err := json.Unmarshal(raw, &tj); err != nil {
			return nil, fmt.Errorf("tree %d: %w", idx, err)
		}
		ts, mas, err := tj.toSummary()
		if err != nil {
			return nil, fmt.Errorf("tree %d: %w", idx, err)
		}
		if err := checkMutAssIDs(idx, mas, muts); err != nil {
			return nil, err
		}
		rs.Trees[idx], rs.MutAss[idx] = ts, mas
	}
	return rs, nil
}

// checkMutAssIDs fails on the first SSM/CNV id a tree assigns that the
// mutlist lacks (a mutlist from the wrong sample). Ids are checked in node,
// then list order so the reported id is deterministic.
func checkMutAssIDs(idx int, mas map[int]*MutAss, muts *MutList) error {
	for _, n := range sortedKeys(mas) {
		for _, id := range mas[n].SSMs {
			if _, ok := muts.SSMs[id]; !ok {
				return fmt.Errorf("tree %d: ssm %s not in mutlist", idx, id)
			}
		}
		for _, id := range mas[n].CNVs {
			if _, ok := muts.CNVs[id]; !ok {
				return fmt.Errorf("tree %d: cnv %s not in mutlist", idx, id)
			}
		}
	}
	return nil
}

// countUnparseableTrees counts tree_-prefixed zip entries the archive reader
// silently skipped (bad name split, index or llh).
func countUnparseableTrees(r *TreeArchiveReader) int {
	n := 0
	for _, f := range r.r.File {
		if strings.HasPrefix(f.Name, "tree_") {
			n++
		}
	}
	return n - r.NumTrees()
}

func loadMutList(path string) (*MutList, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var raw struct {
		SSMs *map[string]struct {
			RefReads   []int   `json:"ref_reads"`
			TotalReads []int   `json:"total_reads"`
			RefInRef   float64 `json:"expected_ref_in_ref"`
			RefInVar   float64 `json:"expected_ref_in_variant"`
			Name       string  `json:"name"`
		} `json:"ssms"`
		CNVs *map[string]json.RawMessage `json:"cnvs"`
	}
	if err := json.Unmarshal(b, &raw); err != nil {
		return nil, fmt.Errorf("parse %s: %w", path, err)
	}
	// Python indexes mutlist['ssms'] / ['cnvs'] directly: absent is a KeyError.
	if raw.SSMs == nil {
		return nil, fmt.Errorf("mutlist %s: missing key \"ssms\"", path)
	}
	if raw.CNVs == nil {
		return nil, fmt.Errorf("mutlist %s: missing key \"cnvs\"", path)
	}
	ml := &MutList{
		SSMs:     make(map[string]*MutListSSM, len(*raw.SSMs)),
		CNVs:     *raw.CNVs,
		CNVReads: make(map[string][2][]int, len(*raw.CNVs)),
	}
	if ml.CNVs == nil {
		ml.CNVs = map[string]json.RawMessage{}
	}
	for id, s := range *raw.SSMs {
		ml.SSMs[id] = &MutListSSM{s.RefReads, s.TotalReads, s.RefInRef, s.RefInVar, s.Name}
	}
	for id, c := range *raw.CNVs {
		var rd struct {
			Ref   []int `json:"ref_reads"`
			Total []int `json:"total_reads"`
		}
		if err := json.Unmarshal(c, &rd); err != nil {
			return nil, fmt.Errorf("cnv %s: %w", id, err)
		}
		ml.CNVReads[id] = [2][]int{rd.Ref, rd.Total}
	}
	return ml, nil
}
