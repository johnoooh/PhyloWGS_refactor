package main

import (
	"errors"
	"fmt"
	"math"
)

// Port of scipy 1.2.3 scipy/stats/kde.py gaussian_kde(data)(data) on numpy
// 1.16.6 linked against OpenBLAS, for d ∈ {1, 2}. Matrices are d×d,
// row-major. float64(...) around products stops Go fusing them into FMAs.

// errKDEFallback marks LinAlgError and FloatingPointError (np.seterr
// invalid='raise'), the exceptions calc_tree_densities catches.
var errKDEFallback = errors.New("gaussian_kde failed")

func kdeFallback(format string, args ...interface{}) error {
	return fmt.Errorf("%w: %s", errKDEFallback, fmt.Sprintf(format, args...))
}

// errNotFinite is scipy.linalg's check_finite ValueError (not caught).
var errNotFinite = errors.New("array must not contain infs or NaNs")

func allFinite(xs []float64) bool {
	for _, x := range xs {
		if math.IsNaN(x) || math.IsInf(x, 0) {
			return false
		}
	}
	return true
}

// pairwiseSum is numpy's pairwise_sum (PW_BLOCKSIZE 128, 8 accumulators).
func pairwiseSum(a []float64) float64 {
	n := len(a)
	switch {
	case n < 8:
		res := 0.0
		for _, x := range a {
			res += x
		}
		return res
	case n <= 128:
		var r [8]float64
		copy(r[:], a[:8])
		i := 8
		for ; i < n-n%8; i += 8 {
			for j := 0; j < 8; j++ {
				r[j] += a[i+j]
			}
		}
		res := ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))
		for ; i < n; i++ {
			res += a[i]
		}
		return res
	}
	n2 := n / 2
	n2 -= n2 % 8
	return pairwiseSum(a[:n2]) + pairwiseSum(a[n2:])
}

// npSum is ndarray.sum over a contiguous axis: identity 0 + pairwise sum.
func npSum(a []float64) float64 { return 0 + pairwiseSum(a) }

// npCov is numpy.cov(data, rowvar=1, bias=False, aweights=w).
func npCov(data [][]float64, w []float64) ([]float64, error) {
	d, n := len(data), len(w)
	wSum := npSum(w)
	ww := make([]float64, n)
	for i := range w {
		ww[i] = w[i] * w[i]
	}
	fact := wSum - npSum(ww)/wSum
	if fact <= 0 {
		fact = 0
	}
	x := make([][]float64, d)
	xw := make([][]float64, d)
	for r := range data {
		prod := make([]float64, n)
		for i, v := range data[r] {
			prod[i] = v * w[i]
		}
		avg := npSum(prod) / wSum
		x[r], xw[r] = make([]float64, n), make([]float64, n)
		for i, v := range data[r] {
			x[r][i] = v - avg
			xw[r][i] = x[r][i] * w[i]
		}
	}
	inv := 1 / fact // np.true_divide(1, 0) is 'divide', not raised
	c := make([]float64, d*d)
	for a := 0; a < d; a++ {
		for b := 0; b < d; b++ {
			s := 0.0 // BLAS ddot / dgemm, sequential over samples
			for i := 0; i < n; i++ {
				s += float64(x[a][i] * xw[b][i])
			}
			c[a*d+b] = s * inv
			// numpy raises 'invalid' only when an op creates a NaN from
			// non-NaN operands; a quiet NaN merely propagating does not
			// raise. For indices in [0, 1] the only NaN here is 0*inf when
			// n == 1 (fact 0). Not a fallback for NaN inputs: do not rely
			// on it as one.
			if math.IsNaN(c[a*d+b]) {
				return nil, kdeFallback("invalid value in multiply")
			}
		}
	}
	return c, nil
}

// getf2 is OpenBLAS getf2 (left-looking LU, idamax pivoting, reciprocal
// scaling) on a d×d matrix. It returns L\U, whether rows 0 and 1 were
// swapped, and whether a pivot was exactly zero (info > 0).
func getf2(a []float64, d int) (lu []float64, swapped, singular bool) {
	lu = append([]float64(nil), a...)
	if d == 1 {
		return lu, false, lu[0] == 0
	}
	if math.Abs(lu[2]) > math.Abs(lu[0]) {
		lu[0], lu[2] = lu[2], lu[0]
		lu[1], lu[3] = lu[3], lu[1]
		swapped = true
	}
	if p := lu[0]; p == 0 {
		singular = true
	} else if math.Abs(p) >= 0x1p-1022 {
		lu[2] *= 1 / p
	} else {
		lu[2] /= p
	}
	lu[3] -= float64(lu[2] * lu[1])
	return lu, swapped, singular || lu[3] == 0
}

// lapackInv is scipy.linalg.inv: getrf, then LAPACK dgetri (unblocked:
// trti2 on U, then solve with L, then undo the column interchange).
func lapackInv(a []float64, d int) ([]float64, error) {
	if !allFinite(a) {
		return nil, errNotFinite
	}
	lu, swapped, singular := getf2(a, d)
	if singular {
		return nil, kdeFallback("singular matrix")
	}
	if d == 1 {
		return []float64{1 / lu[0]}, nil
	}
	v11 := 1 / lu[0]
	v22 := 1 / lu[3]
	v12 := float64(v11*lu[1]) * -v22
	t := -lu[2]
	inv := []float64{v11 + float64(t*v12), v12, float64(t * v22), v22}
	if swapped {
		inv[0], inv[1] = inv[1], inv[0]
		inv[2], inv[3] = inv[3], inv[2]
	}
	return inv, nil
}

// lapackDet is scipy.linalg.det (flinalg det.f): 0 when getrf reports a zero
// pivot, else the U diagonal product with a sign flip per row interchange.
func lapackDet(a []float64, d int) (float64, error) {
	if !allFinite(a) {
		return 0, errNotFinite
	}
	lu, swapped, singular := getf2(a, d)
	if singular {
		return 0, nil
	}
	det := lu[0]
	if swapped {
		det = -1 * lu[0]
	}
	if d == 2 {
		det *= lu[3]
	}
	return det, nil
}

// choleskyUpper is scipy.linalg.cholesky(lower=False) via OpenBLAS potf2_U:
// only the upper triangle is read; U[0][1] is scaled by the reciprocal.
func choleskyUpper(a []float64, d int) ([]float64, error) {
	if !allFinite(a) {
		return nil, errNotFinite
	}
	u := make([]float64, d*d)
	ajj := a[0]
	if ajj <= 0 { // potf2: NaN passes
		return nil, kdeFallback("1-th leading minor not positive definite")
	}
	u[0] = math.Sqrt(ajj)
	if d == 1 {
		return u, nil
	}
	u[1] = a[1] * (1 / u[0])
	ajj = a[3] - float64(u[1]*u[1])
	if ajj <= 0 { // potf2: NaN passes
		return nil, kdeFallback("2-th leading minor not positive definite")
	}
	u[3] = math.Sqrt(ajj)
	return u, nil
}

// kdeSelfDensity is scipy.stats.gaussian_kde(data)(data) with default Scott
// bandwidth and uniform weights; data is d rows of n samples.
func kdeSelfDensity(data [][]float64) ([]float64, error) {
	d, n := len(data), len(data[0])
	if d*n <= 1 {
		return nil, errors.New("`dataset` input should have multiple elements.")
	}
	w := make([]float64, n)
	ww := make([]float64, n)
	for i := range w {
		w[i] = 1 / float64(n)
		ww[i] = w[i] * w[i]
	}
	factor := math.Pow(1/npSum(ww), -1/float64(d+4))
	dcov, err := npCov(data, w)
	if err != nil {
		return nil, err
	}
	dinv, err := lapackInv(dcov, d)
	if err != nil {
		return nil, err
	}
	f2 := factor * factor
	cov2pi := make([]float64, d*d)
	invCov := make([]float64, d*d)
	for i := range dcov {
		cov2pi[i] = (2 * math.Pi) * (dcov[i] * f2)
		invCov[i] = dinv[i] / f2
	}
	det, err := lapackDet(cov2pi, d)
	if err != nil {
		return nil, err
	}
	if det < 0 {
		return nil, kdeFallback("invalid value in sqrt")
	}
	norm := math.Sqrt(det) * float64(n)
	whiten, err := choleskyUpper(invCov, d)
	if err != nil {
		return nil, err
	}
	return kdeEvaluate(data, w, whiten, norm)
}

// kdeEvaluate is gaussian_kde.evaluate with points == dataset (m >= n, so
// it loops over data points and accumulates in index order).
func kdeEvaluate(data [][]float64, w, u []float64, norm float64) ([]float64, error) {
	d, n := len(data), len(w)
	sc := make([][]float64, d) // dot(whitening, dataset)
	for r := range sc {
		sc[r] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		if d == 1 {
			sc[0][i] = u[0] * data[0][i]
			continue
		}
		sc[0][i] = float64(u[0]*data[0][i]) + float64(u[1]*data[1][i])
		sc[1][i] = float64(u[2]*data[0][i]) + float64(u[3]*data[1][i])
	}
	res := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			energy := 0.0
			for r := 0; r < d; r++ {
				diff := sc[r][i] - sc[r][j]
				energy += float64(diff * diff)
			}
			energy /= 2.0
			// Mirrors 'invalid' on inf-inf in diff. numpy does not raise
			// when a quiet NaN merely propagates, and for indices in [0, 1]
			// this is unreachable; it must not be relied on as a fallback
			// for NaN inputs.
			if math.IsNaN(energy) {
				return nil, kdeFallback("invalid value in subtract")
			}
			res[j] += float64(w[i] * math.Exp(-energy))
		}
	}
	for j := range res {
		res[j] = res[j] * float64(n) / norm
	}
	return res, nil
}
