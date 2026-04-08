package main

import (
	"math"
	"math/rand"
	"testing"
)

// ============================================================================
// BLOCKER 1d: dpAlphaLLH must include the root's main-stick beta term
// ============================================================================
//
// Python's TSSB.dp_alpha_llh (tssb.py:247-255) evaluates
//     betapdfln(root['main'], 1.0, dp_alpha)
// at depth 0 whenever self.min_depth <= depth (default min_depth=0 makes this
// always true). The depth-0 term depends on dp_alpha, so it cannot be treated
// as a constant that cancels in the slice sampler's log-likelihood comparison.
//
// Go's closure dpAlphaLLH in resampleHypers (main.go:2357-2370) previously
// gated on `if depth >= 1`, silently excluding the root term. This test pins
// the fix by asserting the computed LLH equals the full sum including depth 0.

func TestDPAlphaLLHIncludesRootMain(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	ssms := []*SSM{makeTestSSM("s0", 3, rng)}
	cnvs := []*CNV{}
	tssb := newTSSB(ssms, cnvs, 25.0, 1.0, 0.25, rng)

	// Pin root.Main to a known, non-trivial value so the depth-0 term has a
	// well-defined contribution (post-fix we unconditionally compute it; this
	// must not blow up on the 1e-30 sentinel that the MCMC loop uses).
	tssb.Root.Main = 0.3

	// Pin the child's main as well so the depth-1 term is well-defined.
	if len(tssb.Root.Children) != 1 {
		t.Fatalf("expected 1 initial child, got %d", len(tssb.Root.Children))
	}
	child := tssb.Root.Children[0]
	child.Main = 0.6

	alpha := 7.0
	decay := 0.25
	tssb.AlphaDecay = decay

	// Expected: sum over all nodes at depth d >= min_depth (0) of
	//     betaPDFLn(node.Main, 1.0, decay^d * alpha)
	// With min_depth=0, this INCLUDES depth 0 (root).
	want := betaPDFLn(tssb.Root.Main, 1.0, math.Pow(decay, 0)*alpha) +
		betaPDFLn(child.Main, 1.0, math.Pow(decay, 1)*alpha)

	got := tssb.dpAlphaLLH(alpha)
	if math.Abs(got-want) > 1e-12 {
		t.Errorf("dpAlphaLLH(%.3f): want %.15g, got %.15g, diff %.3e",
			alpha, want, got, math.Abs(got-want))
	}

	// Also verify the sampler is sensitive to root.Main — if the fix is wrong
	// and root.Main is dropped, the LLH at alpha1 will equal the LLH at alpha2
	// modulo the child term only. Assert that changing root.Main changes the
	// computed LLH by exactly the expected delta.
	oldMain := tssb.Root.Main
	tssb.Root.Main = 0.7
	got2 := tssb.dpAlphaLLH(alpha)
	expectedDelta := betaPDFLn(0.7, 1.0, alpha) - betaPDFLn(oldMain, 1.0, alpha)
	actualDelta := got2 - got
	if math.Abs(actualDelta-expectedDelta) > 1e-12 {
		t.Errorf("dpAlphaLLH is not sensitive to root.Main:\n"+
			"  expected delta from changing root.Main %.2f → %.2f: %.15g\n"+
			"  actual delta:                                       %.15g\n"+
			"  (buggy behavior would produce delta == 0, since the root term is dropped)",
			oldMain, 0.7, expectedDelta, actualDelta)
	}
}

// ============================================================================
// BLOCKER 2a: resampleStickOrders leftover mass bucket must match Python
// ============================================================================
//
// Python (tssb.py:208-210):
//     sub_weights = hstack([all_weights[sub_indices], 1.0 - sum(all_weights)])
//     sub_weights = sub_weights / sum(sub_weights)
// The "unallocated region" bucket is 1 - sum(ALL stick weights), i.e. the
// true unallocated stick mass.
//
// Go previously computed the bucket as 1 - sum(NOT-yet-ordered weights),
// which conflates the already-ordered children's mass with the genuinely
// unallocated mass and inflates P(spawn) by up to 3x.
//
// This test pins the correct leftover computation via the extracted helper
// stickOrderSubWeights(allWeights, newOrder).

func TestStickOrderSubWeightsMatchesPython(t *testing.T) {
	// Scenario: parent has 3 children with stick weights [0.4, 0.3, 0.2].
	// True unallocated mass = 1 - (0.4+0.3+0.2) = 0.1.
	// Child 0 is already in new_order.
	//
	// Python's leftover bucket = 0.1 (true unallocated), giving normalized
	// sub_weights = [0.3, 0.2, 0.1] / 0.6 = [0.5, 0.333..., 0.1666...].
	// P(spawn) = 0.1666...
	//
	// Go's buggy leftover was 1 - (0.3+0.2) = 0.5, giving normalized
	// sub_weights = [0.3, 0.2, 0.5] / 1.0 = [0.3, 0.2, 0.5]. P(spawn) = 0.5.
	allWeights := []float64{0.4, 0.3, 0.2}
	newOrder := []int{0}

	got := stickOrderSubWeights(allWeights, newOrder)
	if len(got) != 3 {
		t.Fatalf("expected 3 sub_weights (2 remaining children + leftover), got %d", len(got))
	}

	// Expected Python values (normalized):
	//   [0.3, 0.2, 0.1] / sum(0.6) = [0.5, 0.333..., 0.1666...]
	want := []float64{0.5, 1.0 / 3.0, 1.0 / 6.0}
	for i, w := range want {
		if math.Abs(got[i]-w) > 1e-12 {
			t.Errorf("sub_weights[%d]: want %.12f, got %.12f, diff %.3e",
				i, w, got[i], math.Abs(got[i]-w))
		}
	}

	// Specifically verify the leftover bucket is NOT 0.5 (the buggy value).
	if math.Abs(got[2]-0.5) < 1e-9 {
		t.Errorf("leftover bucket is 0.5 (buggy Go value). Python would give 0.1666...")
	}

	// Also verify normalization: sub_weights should sum to 1.
	sum := 0.0
	for _, w := range got {
		sum += w
	}
	if math.Abs(sum-1.0) > 1e-12 {
		t.Errorf("sub_weights should sum to 1, got %.15g", sum)
	}
}

func TestStickOrderSubWeightsNoneOrdered(t *testing.T) {
	// Edge case: nothing in new_order yet (first iteration of outer loop).
	// Both buggy and correct implementations should agree here because
	// totalUsed == sum(allWeights) when subIndices == all indices.
	allWeights := []float64{0.4, 0.3, 0.2}
	newOrder := []int{}

	got := stickOrderSubWeights(allWeights, newOrder)
	if len(got) != 4 {
		t.Fatalf("expected 4 sub_weights (3 children + leftover), got %d", len(got))
	}

	// Python: [0.4, 0.3, 0.2, 0.1] normalized by 1.0 = same values.
	want := []float64{0.4, 0.3, 0.2, 0.1}
	for i, w := range want {
		if math.Abs(got[i]-w) > 1e-12 {
			t.Errorf("sub_weights[%d]: want %.12f, got %.12f", i, w, got[i])
		}
	}
}

func TestStickOrderSubWeightsAllOrdered(t *testing.T) {
	// Edge case: every child is already in new_order. Only the leftover bucket
	// remains. Should normalize to [1.0].
	allWeights := []float64{0.4, 0.3, 0.2}
	newOrder := []int{0, 1, 2}

	got := stickOrderSubWeights(allWeights, newOrder)
	if len(got) != 1 {
		t.Fatalf("expected 1 sub_weight (leftover only), got %d", len(got))
	}
	if math.Abs(got[0]-1.0) > 1e-12 {
		t.Errorf("single leftover bucket should be 1.0, got %.15g", got[0])
	}
}

// ============================================================================
// SUSPICIOUS 5b: spawnChild must use one broadcast scalar across pi dims
// ============================================================================
//
// Python alleles.__init__ (alleles.py:35-37):
//     self.pi = rand(1) * parent.pi
//     parent.pi = parent.pi - self.pi
// Here rand(1) returns a length-1 numpy array. Broadcasting (1,) * (ntps,)
// applies ONE random scalar uniformly across all ntps dimensions, making
// child.pi proportional to parent.pi.
//
// Go previously drew rng.Float64() inside a per-dim for loop, breaking the
// proportionality. Fix: hoist the draw outside the loop.

func TestSpawnChildUsesSingleBroadcastScalar(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	ssms := []*SSM{makeTestSSM("s0", 5, rng)}
	cnvs := []*CNV{}
	tssb := newTSSB(ssms, cnvs, 25.0, 1.0, 0.25, rng)

	// Grab the existing child (created by newTSSB), replace its pi with a
	// distinctive, non-uniform vector so the proportionality test is meaningful.
	parent := tssb.Root.Children[0]
	parent.Node.Pi = []float64{0.10, 0.25, 0.50, 0.15, 0.05}
	origParentPi := make([]float64, len(parent.Node.Pi))
	copy(origParentPi, parent.Node.Pi)

	// Spawn a new child under this parent.
	child := tssb.spawnChild(parent, 2, rng)

	// The child's pi should be proportional to the original parent pi:
	//     child.Pi[i] = frac * origParentPi[i]
	// for some SINGLE frac (the broadcast scalar).
	// Compute the implied frac at each dim and assert they all match.
	if len(child.Node.Pi) != len(origParentPi) {
		t.Fatalf("child.Pi length %d != parent.Pi length %d",
			len(child.Node.Pi), len(origParentPi))
	}

	fracs := make([]float64, len(child.Node.Pi))
	for i := range child.Node.Pi {
		if origParentPi[i] == 0 {
			t.Fatalf("origParentPi[%d] is zero; test setup wrong", i)
		}
		fracs[i] = child.Node.Pi[i] / origParentPi[i]
	}

	// All fracs must be equal (within float tolerance). Buggy Go drew
	// independent randoms, so fracs would differ across dims.
	for i := 1; i < len(fracs); i++ {
		if math.Abs(fracs[i]-fracs[0]) > 1e-12 {
			t.Errorf("child.Pi is NOT proportional to parent.Pi:\n"+
				"  frac[0] = %.15g\n"+
				"  frac[%d] = %.15g\n"+
				"  diff = %.3e\n"+
				"  (Python uses one broadcast scalar; Go previously drew\n"+
				"   rng.Float64() per dim, giving independent fracs.)",
				fracs[0], i, fracs[i], math.Abs(fracs[i]-fracs[0]))
		}
	}

	// Sanity: frac should be in (0, 1).
	if fracs[0] <= 0 || fracs[0] >= 1 {
		t.Errorf("frac out of (0,1): %.15g", fracs[0])
	}

	// Sanity: parent.Pi should be origParentPi - child.Pi = (1-frac) * origParentPi.
	// Assert this proportionality too.
	complementFracs := make([]float64, len(origParentPi))
	for i := range origParentPi {
		complementFracs[i] = parent.Node.Pi[i] / origParentPi[i]
	}
	for i := 1; i < len(complementFracs); i++ {
		if math.Abs(complementFracs[i]-complementFracs[0]) > 1e-12 {
			t.Errorf("post-spawn parent.Pi is not proportional to original: frac[0]=%.15g frac[%d]=%.15g",
				complementFracs[0], i, complementFracs[i])
		}
	}
}
