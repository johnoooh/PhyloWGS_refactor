package main

import (
	"math"
	"math/rand"
	"testing"
)

// makeTestSSM creates a minimal SSM for unit tests.
func makeTestSSM(id string, nTP int, rng *rand.Rand) *SSM {
	ssm := &SSM{
		ID:              id,
		Name:            id,
		A:               make([]int, nTP),
		D:               make([]int, nTP),
		MuR:             0.999,
		MuV:             0.5,
		LogBinNormConst: make([]float64, nTP),
	}
	for tp := 0; tp < nTP; tp++ {
		ssm.D[tp] = rng.Intn(400) + 100
		ssm.A[tp] = rng.Intn(ssm.D[tp])
		ssm.LogBinNormConst[tp] = logBinCoeff(ssm.D[tp], ssm.A[tp])
	}
	return ssm
}

// makeTestCNV creates a minimal CNV for unit tests.
func makeTestCNV(id string, nTP int, rng *rand.Rand) *CNV {
	cnv := &CNV{
		ID:              id,
		A:               make([]int, nTP),
		D:               make([]int, nTP),
		LogBinNormConst: make([]float64, nTP),
	}
	for tp := 0; tp < nTP; tp++ {
		cnv.D[tp] = rng.Intn(400) + 100
		cnv.A[tp] = rng.Intn(cnv.D[tp])
		cnv.LogBinNormConst[tp] = logBinCoeff(cnv.D[tp], cnv.A[tp])
	}
	return cnv
}

// buildTestTSSB constructs a small tree with a few forked nodes so that
// computeNGenomes traverses a nontrivial set of nodes. It then manually
// assigns ssms and cnv to specific nodes to exercise the 4 CNV cases.
func buildTestTSSB(t *testing.T, rng *rand.Rand) (*TSSB, []*SSM, *CNV) {
	t.Helper()

	nTP := 3
	nSSM := 6

	ssms := make([]*SSM, nSSM)
	for i := 0; i < nSSM; i++ {
		ssms[i] = makeTestSSM("s"+string(rune('0'+i)), nTP, rng)
	}
	cnv := makeTestCNV("c0", nTP, rng)

	// Attach CNV to the first 3 SSMs (major=2, minor=1)
	for i := 0; i < 3; i++ {
		ssms[i].CNVs = []*CNVRef{{CNV: cnv, MaternalCN: 1, PaternalCN: 2}}
	}

	cnvs := []*CNV{cnv}
	tssb := newTSSB(ssms, cnvs, 25.0, 1.0, 0.25, rng)

	// Manually spawn two more child nodes to grow a tree of depth 2 with branching.
	// Root -> childA (existing from newTSSB) -> grand1, grand2
	rootTSSB := tssb.Root
	childA := rootTSSB.Children[0]
	grand1 := tssb.spawnChild(childA, 2, rng)
	grand2 := tssb.spawnChild(childA, 2, rng)

	// Move ssms[3] and ssms[4] to grand1 and grand2, ssms[5] stays at childA
	oldNode := childA.Node
	// Remove 3 and 4 from oldNode.Data
	remove := func(slice []int, vals ...int) []int {
		out := slice[:0]
	outer:
		for _, v := range slice {
			for _, vv := range vals {
				if v == vv {
					continue outer
				}
			}
			out = append(out, v)
		}
		return out
	}
	oldNode.Data = remove(oldNode.Data, 3, 4)

	grand1.Node.Data = append(grand1.Node.Data, 3)
	ssms[3].Node = grand1.Node
	grand2.Node.Data = append(grand2.Node.Data, 4)
	ssms[4].Node = grand2.Node

	// Put SSM 0 (with CNV) on grand1 so SSM-CNV timing differs from co-location
	oldNode.Data = remove(oldNode.Data, 0)
	grand1.Node.Data = append(grand1.Node.Data, 0)
	ssms[0].Node = grand1.Node

	// Move CNV to grand2 so its node differs from SSM 0's
	cnv.Node = grand2.Node

	// Randomize pi/params so likelihood is non-trivial
	for _, n := range tssb.getNodes() {
		for tp := 0; tp < nTP; tp++ {
			n.Pi[tp] = 0.1 + 0.8*rng.Float64()
			n.Params[tp] = n.Pi[tp]
			n.Pi1[tp] = 0.1 + 0.8*rng.Float64()
			n.Params1[tp] = n.Pi1[tp]
		}
	}

	tssb.setNodeHeights()
	tssb.setNodePaths()
	tssb.rebuildAncestorSets()

	return tssb, ssms, cnv
}

// TestPrecomputedMHStatesEquivalence asserts that after calling
// precomputeMHStates(), logLikelihoodWithCNVTreeMH produces identical results
// to the slow tree-walking path, for both old and new state.
//
// This is the key invariant: the optimization must be a pure refactor, not a
// semantic change. Matches Python's mh.cpp precomputed-state inner loop.
func TestPrecomputedMHStatesEquivalence(t *testing.T) {
	rng := rand.New(rand.NewSource(12345))
	tssb, ssms, _ := buildTestTSSB(t, rng)

	// Baseline: compute LLH via the slow tree-walking path for every SSM,
	// for both newState=false and newState=true.
	baselineOld := make([]float64, len(ssms))
	baselineNew := make([]float64, len(ssms))
	for i, ssm := range ssms {
		baselineOld[i] = logLikelihoodWithCNVTreeMH(ssm, tssb, false)
		baselineNew[i] = logLikelihoodWithCNVTreeMH(ssm, tssb, true)
	}

	// Now precompute MH states. After this call, logLikelihoodWithCNVTreeMH
	// should take a fast O(K) dot-product path for CNV-affected SSMs.
	tssb.precomputeMHStates()

	for i, ssm := range ssms {
		gotOld := logLikelihoodWithCNVTreeMH(ssm, tssb, false)
		gotNew := logLikelihoodWithCNVTreeMH(ssm, tssb, true)

		if math.Abs(gotOld-baselineOld[i]) > 1e-10 {
			t.Errorf("SSM %s: old-state LLH mismatch after precompute\n  baseline: %.15e\n  precomp:  %.15e\n  diff:     %.2e",
				ssm.ID, baselineOld[i], gotOld, math.Abs(gotOld-baselineOld[i]))
		}
		if math.Abs(gotNew-baselineNew[i]) > 1e-10 {
			t.Errorf("SSM %s: new-state LLH mismatch after precompute\n  baseline: %.15e\n  precomp:  %.15e\n  diff:     %.2e",
				ssm.ID, baselineNew[i], gotNew, math.Abs(gotNew-baselineNew[i]))
		}
	}
}
