package main

import (
	"errors"
	"fmt"
	"math"
)

// Port of pwgsresults/index_calculator.py and calc_tree_densities in
// pwgsresults/json_writer.py.

const (
	relAncDesc = "anc_desc"
	relDescAnc = "desc_anc"
	relCousin  = "cousin"
)

// popRelations is IndexCalculator._determine_pop_relations: ancestor and
// descendant pairs from a preorder walk from root 0, every other pair of
// reached vertices is a cousin pair (both orders).
func popRelations(structure map[int][]int) (map[[2]int]string, error) {
	rel := map[[2]int]string{}
	var verts []int
	seen := map[int]bool{}
	var walk func(v int, anc []int) error
	walk = func(v int, anc []int) error {
		if containsInt(anc, v) {
			return fmt.Errorf("cycle through population %d", v)
		}
		if !seen[v] {
			seen[v] = true
			verts = append(verts, v)
		}
		for _, a := range anc {
			rel[[2]int{a, v}] = relAncDesc
			rel[[2]int{v, a}] = relDescAnc
		}
		next := append(append([]int(nil), anc...), v)
		for _, c := range structure[v] {
			if err := walk(c, next); err != nil {
				return err
			}
		}
		return nil
	}
	if err := walk(0, nil); err != nil {
		return nil, err
	}
	for i, v1 := range verts {
		for _, v2 := range verts[i+1:] {
			if _, ok := rel[[2]int{v1, v2}]; ok {
				continue
			}
			rel[[2]int{v1, v2}] = relCousin
			rel[[2]int{v2, v1}] = relCousin
		}
	}
	return rel, nil
}

// calcIndex is IndexCalculator._calc_index: Σ n_i·n_j over ordered pairs of
// the relation type, normalised by N(N−1). The sum is an integer in Python.
func calcIndex(rel map[[2]int]string, pops map[int]*Population, total int, kind string) (float64, error) {
	index := 0
	for pair, r := range rel {
		if r != kind {
			continue
		}
		p1, ok1 := pops[pair[0]]
		p2, ok2 := pops[pair[1]]
		if !ok1 || !ok2 {
			return 0, fmt.Errorf("population missing for pair %v", pair)
		}
		index += p1.NumSSMs * p2.NumSSMs
	}
	norm := float64(index) / float64(total*(total-1))
	if !(norm >= 0 && norm <= 1) {
		return 0, fmt.Errorf("%s index %v outside [0, 1]", kind, norm)
	}
	return norm, nil
}

// treeIndices returns the linearity, branching and clustering indices of one
// tree, with the Python assertions as errors.
func treeIndices(t *TreeSummary) (lin, br, cl float64, err error) {
	total, cc := 0, 0
	for _, p := range t.Populations {
		total += p.NumSSMs
		cc += p.NumSSMs * (p.NumSSMs - 1)
	}
	if total <= 1 {
		return 0, 0, 0, errors.New("fewer than 2 SSMs")
	}
	rel, err := popRelations(t.Structure)
	if err != nil {
		return 0, 0, 0, err
	}
	low, err := calcIndex(rel, t.Populations, total, relAncDesc)
	if err != nil {
		return 0, 0, 0, err
	}
	lin = 2 * low
	if !(lin >= 0 && lin <= 1) {
		return 0, 0, 0, fmt.Errorf("linearity index %v outside [0, 1]", lin)
	}
	if br, err = calcIndex(rel, t.Populations, total, relCousin); err != nil {
		return 0, 0, 0, err
	}
	cl = float64(cc) / float64(total*(total-1))
	if !(cl >= 0 && cl <= 1) {
		return 0, 0, 0, fmt.Errorf("clustering index %v outside [0, 1]", cl)
	}
	// np.isclose(cl, 1 - lin - br): |a-b| <= atol + rtol*|b|.
	other := 1 - lin - br
	if !(math.Abs(cl-other) <= 1e-8+1e-5*math.Abs(other)) {
		return 0, 0, 0, fmt.Errorf("indices do not sum to 1: cl=%v, 1-lin-br=%v", cl, other)
	}
	return lin, br, cl, nil
}

// allTreeIndices runs treeIndices over every tree; values are {lin, br, cl}.
func allTreeIndices(trees map[int]*TreeSummary) (map[int][3]float64, error) {
	out := make(map[int][3]float64, len(trees))
	for _, k := range sortedKeys(trees) {
		lin, br, cl, err := treeIndices(trees[k])
		if err != nil {
			return nil, fmt.Errorf("tree %d: %w", k, err)
		}
		out[k] = [3]float64{lin, br, cl}
	}
	return out, nil
}

// treeDensities is calc_tree_densities: a 2-D gaussian_kde on
// (clustering, branching/(branching+linearity+1e-4)), falling back to 1-D on
// clustering, then to zeros, on LinAlgError or FloatingPointError. The mode
// is "2d", "1d" or "zeros". Other gaussian_kde errors (ValueError in Python)
// are returned.
func treeDensities(trees map[int]*TreeSummary, idx map[int][3]float64) (map[int]float64, string, error) {
	keys := sortedKeys(trees)
	x := make([]float64, len(keys))
	y := make([]float64, len(keys))
	for i, k := range keys {
		v, ok := idx[k]
		if !ok {
			return nil, "", fmt.Errorf("tree %d: no indices", k)
		}
		lin, br, cl := v[0], v[1], v[2]
		x[i] = cl
		y[i] = br / ((br + lin) + 0.0001)
	}
	dens, err := kdeSelfDensity([][]float64{x, y})
	mode := "2d"
	if errors.Is(err, errKDEFallback) {
		dens, err = kdeSelfDensity([][]float64{x})
		mode = "1d"
		if errors.Is(err, errKDEFallback) {
			dens, err = make([]float64, len(x)), nil
			mode = "zeros"
		}
	}
	if err != nil {
		return nil, "", err
	}
	out := make(map[int]float64, len(keys))
	for i, k := range keys {
		out[k] = dens[i]
	}
	return out, mode, nil
}
