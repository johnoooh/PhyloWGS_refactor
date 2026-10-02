package main

import (
	"math/rand"
	"testing"
)

func denseFromNodes(nodes []*Node, ntps int) [][]float64 {
	d := make([][]float64, ntps)
	for tp := range d {
		d[tp] = make([]float64, len(nodes))
		for k, n := range nodes {
			d[tp][k] = n.Pi1[tp]
		}
	}
	return d
}

// The dense MH path must be bit-identical to the reference precomputed path
// for both co-located (4-state) and non-co-located (2-state) CNV SSMs, and
// paramPost must give the same total with and without densePi.
func TestLogLikelihoodMHDense_BitIdentical(t *testing.T) {
	for seed := int64(1); seed <= 20; seed++ {
		rng := rand.New(rand.NewSource(seed))
		tssb, ssms, _ := buildTestTSSB(t, rng)
		tssb.precomputeMHStates()
		_, nodes := tssb.getMixture()
		// Random proposal in Pi1 so old/new states differ.
		for tp := 0; tp < tssb.NTPS; tp++ {
			alpha := make([]float64, len(nodes))
			for i := range alpha {
				alpha[i] = 1
			}
			p := dirichletSample(alpha, rng)
			for i, n := range nodes {
				n.Pi1[tp] = p[i]
			}
		}
		if !mhStatesAligned(tssb.Data, nodes) {
			t.Fatalf("seed %d: precomputed states not aligned with getMixture order", seed)
		}
		dense := denseFromNodes(nodes, tssb.NTPS)
		saw4, saw2 := false, false
		for _, ssm := range ssms {
			if len(ssm.CNVs) == 0 {
				continue
			}
			if ssm.UseFourStates {
				saw4 = true
			} else {
				saw2 = true
			}
			want := logLikelihoodWithCNVTreeMHPrecomputed(ssm, true)
			got := logLikelihoodMHDense(ssm, dense)
			if got != want {
				t.Fatalf("seed %d SSM %s (four=%v): dense %.17g != reference %.17g", seed, ssm.ID, ssm.UseFourStates, got, want)
			}
		}
		if !saw4 || !saw2 {
			t.Fatalf("seed %d: fixture must exercise both 4-state (%v) and 2-state (%v) SSMs", seed, saw4, saw2)
		}
		if a, b := tssb.paramPost(nodes, true, nil), tssb.paramPost(nodes, true, dense); a != b {
			t.Fatalf("seed %d: paramPost reference %.17g != dense %.17g", seed, a, b)
		}
		// The fixture's 2-state SSMs can have state 1 == state 2, which would
		// hide a 1/2 mix-up in the 3/4 reuse. Make them asymmetric (keeping
		// the computeSSMStates invariant 3==1, 4==2) and re-check.
		for _, ssm := range ssms {
			if len(ssm.CNVs) == 0 || ssm.UseFourStates {
				continue
			}
			for k := range ssm.MHStates {
				s := &ssm.MHStates[k]
				s.Nv1 += 0.5 + float64(k%3)
				s.Nr3, s.Nv3 = s.Nr1, s.Nv1
				s.Nr4, s.Nv4 = s.Nr2, s.Nv2
			}
			want := logLikelihoodWithCNVTreeMHPrecomputed(ssm, true)
			if got := logLikelihoodMHDense(ssm, dense); got != want {
				t.Fatalf("seed %d SSM %s asymmetric: dense %.17g != reference %.17g", seed, ssm.ID, got, want)
			}
		}
	}
}

func TestMHStatesAligned_DetectsMisalignment(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	tssb, _, _ := buildTestTSSB(t, rng)
	tssb.precomputeMHStates()
	_, nodes := tssb.getMixture()
	if len(nodes) < 2 {
		t.Skip("need >=2 nodes")
	}
	swapped := append([]*Node{}, nodes...)
	swapped[0], swapped[1] = swapped[1], swapped[0]
	if mhStatesAligned(tssb.Data, swapped) {
		t.Fatal("misaligned node order not detected")
	}
	if mhStatesAligned(tssb.Data, nodes[:len(nodes)-1]) {
		t.Fatal("length mismatch not detected")
	}
}
