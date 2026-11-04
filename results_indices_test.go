package main

import (
	"encoding/json"
	"errors"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

func TestTreeIndices_HandComputed(t *testing.T) {
	ts := &TreeSummary{
		Structure: map[int][]int{0: {1}, 1: {2, 3}},
		Populations: map[int]*Population{
			0: {}, 1: {NumSSMs: 4}, 2: {NumSSMs: 2}, 3: {NumSSMs: 2},
		},
	}
	lin, br, cl, err := treeIndices(ts)
	if err != nil {
		t.Fatal(err)
	}
	if lin != 2*(16.0/56) || br != 8.0/56 || cl != 16.0/56 {
		t.Fatalf("got lin=%v br=%v cl=%v", lin, br, cl)
	}
}

func TestTreeIndices_OneSSM(t *testing.T) {
	ts := &TreeSummary{
		Structure:   map[int][]int{0: {1}},
		Populations: map[int]*Population{0: {}, 1: {NumSSMs: 1}},
	}
	if _, _, _, err := treeIndices(ts); err == nil {
		t.Fatal("want error for a tree with fewer than 2 SSMs")
	}
	_, err := allTreeIndices(map[int]*TreeSummary{7: ts})
	if err == nil || err.Error() != "tree 7: fewer than 2 SSMs" {
		t.Fatalf("got %v", err)
	}
}

// Kill: a population named in the structure but absent from populations
// is a KeyError in Python, not a silent zero.
func TestTreeIndices_MissingPopulation(t *testing.T) {
	ts := &TreeSummary{
		Structure:   map[int][]int{0: {1, 2}},
		Populations: map[int]*Population{0: {}, 1: {NumSSMs: 3}},
	}
	if _, _, _, err := treeIndices(ts); err == nil {
		t.Fatal("want error for missing population")
	}
}

type expectedSumm struct {
	Trees map[string]struct {
		Lin       float64          `json:"linearity_index"`
		Br        float64          `json:"branching_index"`
		Cl        float64          `json:"clustering_index"`
		Structure map[string][]int `json:"structure"`
		Pops      map[string]struct {
			NumSSMs int `json:"num_ssms"`
		} `json:"populations"`
	} `json:"trees"`
	Densities map[string]float64 `json:"tree_densities"`
}

func atoiT(t *testing.T, s string) int {
	t.Helper()
	n, err := strconv.Atoi(s)
	if err != nil {
		t.Fatal(err)
	}
	return n
}

// loadExpectedTrees rebuilds the post-munger trees from a Python summ.json.
func loadExpectedTrees(t *testing.T, dir string) (map[int]*TreeSummary, *expectedSumm) {
	t.Helper()
	b, err := os.ReadFile(filepath.Join(dir, "expected", "summ.json"))
	if err != nil {
		t.Fatal(err)
	}
	var es expectedSumm
	if err := json.Unmarshal(b, &es); err != nil {
		t.Fatal(err)
	}
	trees := map[int]*TreeSummary{}
	for k, tr := range es.Trees {
		ts := &TreeSummary{Structure: map[int][]int{}, Populations: map[int]*Population{}}
		for p, kids := range tr.Structure {
			ts.Structure[atoiT(t, p)] = kids
		}
		for p, pop := range tr.Pops {
			ts.Populations[atoiT(t, p)] = &Population{NumSSMs: pop.NumSSMs}
		}
		trees[atoiT(t, k)] = ts
	}
	if len(trees) == 0 {
		t.Fatalf("%s: no trees", dir)
	}
	return trees, &es
}

func TestTreeDensities_Paths(t *testing.T) {
	cases := []struct{ dir, mode string }{
		{"linear/default", "1d"},
		{"identical/default", "1d"},
		{"sim/default", "2d"},
		{"sim/mincount", "2d"},
		{"superclone/default", "1d"},
		{"multiprimary/prod", "1d"},
		{"cnvonly/default", "zeros"},
	}
	for _, c := range cases {
		trees, es := loadExpectedTrees(t, filepath.Join("testdata", "native_results", c.dir))
		idx, err := allTreeIndices(trees)
		if err != nil {
			t.Fatalf("%s: %v", c.dir, err)
		}
		for k, v := range idx {
			want := es.Trees[strconv.Itoa(k)]
			if v != [3]float64{want.Lin, want.Br, want.Cl} {
				t.Errorf("%s tree %d: indices %v, want %v", c.dir, k, v, want)
			}
		}
		dens, mode, err := treeDensities(trees, idx)
		if err != nil {
			t.Fatalf("%s: %v", c.dir, err)
		}
		if mode != c.mode {
			t.Errorf("%s: mode %q, want %q", c.dir, mode, c.mode)
		}
		if len(dens) != len(es.Densities) {
			t.Fatalf("%s: %d densities, want %d", c.dir, len(dens), len(es.Densities))
		}
		maxRel := 0.0
		for k, want := range es.Densities {
			got, ok := dens[atoiT(t, k)]
			rel := math.Abs(got - want)
			if want != 0 {
				rel /= math.Abs(want)
			}
			if !ok || rel > 1e-9 {
				t.Errorf("%s tree %s: density %v, want %v", c.dir, k, got, want)
			}
			maxRel = math.Max(maxRel, rel)
		}
		t.Logf("%s: mode=%s max rel err=%.3g", c.dir, mode, maxRel)
	}
}

func TestTreeDensities_SingleTreeError(t *testing.T) {
	want, err := os.ReadFile("testdata/native_results/single/default/expected/error.txt")
	if err != nil {
		t.Fatal(err)
	}
	trees := map[int]*TreeSummary{0: {}}
	idx := map[int][3]float64{0: {0.5, 0, 0.5}}
	_, _, err = treeDensities(trees, idx)
	if err == nil || !strings.HasSuffix(err.Error(), strings.TrimSpace(string(want))) {
		t.Fatalf("got %v, want suffix %q", err, want)
	}
}

func TestTreeDensities_MissingIndex(t *testing.T) {
	trees := map[int]*TreeSummary{0: {}, 1: {}}
	if _, _, err := treeDensities(trees, map[int][3]float64{0: {}}); err == nil {
		t.Fatal("want error for a tree without indices")
	}
}

func TestKDELinearAlgebra2x2(t *testing.T) {
	a := []float64{1, 2, 3, 4} // |a21| > |a11|: rows swap, det sign flips
	inv, err := lapackInv(a, 2)
	if err != nil {
		t.Fatal(err)
	}
	for i, want := range []float64{-2, 1, 1.5, -0.5} {
		if math.Abs(inv[i]-want) > 1e-15 {
			t.Errorf("inv[%d] = %v, want %v", i, inv[i], want)
		}
	}
	if det, _ := lapackDet(a, 2); math.Abs(det+2) > 1e-15 {
		t.Errorf("det = %v, want -2", det)
	}
	if det, _ := lapackDet([]float64{1, 2, 2, 4}, 2); det != 0 {
		t.Errorf("singular det = %v, want 0", det)
	}
	if _, err := lapackInv([]float64{1, 2, 2, 4}, 2); !errors.Is(err, errKDEFallback) {
		t.Errorf("singular inv: %v", err)
	}
	for _, m := range [][]float64{{4, 2, 2, 1}, {1, 0, 0, -1}, {-1, 0, 0, 1}} {
		if _, err := choleskyUpper(m, 2); !errors.Is(err, errKDEFallback) {
			t.Errorf("cholesky %v: want fallback, got %v", m, err)
		}
	}
	// Without partial pivoting inv[0] cancels to 0 instead of -1.
	inv, _ = lapackInv([]float64{1e-20, 1, 1, 1}, 2)
	if inv[0] != -1 || inv[1] != 1 || inv[2] != 1 || inv[3] != -1e-20 {
		t.Errorf("tiny-pivot inv = %v", inv)
	}
	if _, err := lapackInv([]float64{math.Inf(1), 0, 0, 1}, 2); err != errNotFinite {
		t.Errorf("non-finite inv: %v", err)
	}
}

// One sample in 2-D: np.cov has fact 0, 0*inf is 'invalid' -> fallback.
func TestKDEOneSampleCovIsFallback(t *testing.T) {
	_, err := kdeSelfDensity([][]float64{{0.5}, {0.2}})
	if !errors.Is(err, errKDEFallback) {
		t.Fatalf("got %v", err)
	}
}

// Ruling 8: when the trees sit on exactly two distinct (X, Y) points the
// 2-D covariance is mathematically rank 1, and Python's own result depends
// on the OpenBLAS kernel (2d with Inf densities under Prescott, 1d under
// Sandybridge/Haswell). It is not a parity target; this pins Go's
// deterministic sequential port (reviewer fuzz case 25).
func TestTreeDensities_RankDeficientIsDeterministic(t *testing.T) {
	trees := map[int]*TreeSummary{}
	idx := map[int][3]float64{}
	want := map[int]float64{}
	for _, k := range []int{1, 4, 16, 22, 31} {
		trees[k] = &TreeSummary{}
		idx[k] = [3]float64{0x1.5ad6b5ad6b5adp-1, 0x1.3d1c13d1c13d2p-5, 0x1.22af122af122bp-2}
		want[k] = 0x1.a64751efb286ep+02
	}
	for _, k := range []int{7, 10, 13, 19, 25, 28, 34, 37} {
		trees[k] = &TreeSummary{}
		idx[k] = [3]float64{0x1.aaaaaaaaaaaabp-2, 0x1.c71c71c71c71cp-3, 0x1.71c71c71c71c7p-2}
		want[k] = 0x1.5067548edee61p+03
	}
	dens, mode, err := treeDensities(trees, idx)
	if err != nil || mode != "1d" {
		t.Fatalf("mode %q err %v, want 1d", mode, err)
	}
	if len(dens) != len(want) {
		t.Fatalf("%d densities, want %d", len(dens), len(want))
	}
	for k, w := range want {
		if dens[k] != w {
			t.Errorf("tree %d: %x, want %x", k, dens[k], w)
		}
	}
}
