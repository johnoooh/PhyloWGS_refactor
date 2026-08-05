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

	// Also co-locate SSM 1 with the CNV itself (both at grand2), to
	// exercise the UseFourStates=true (M1) code path — previously no SSM
	// in this fixture had Node == its own CNV's Node.
	oldNode.Data = remove(oldNode.Data, 1)
	grand2.Node.Data = append(grand2.Node.Data, 1)
	ssms[1].Node = grand2.Node

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

// TestPrecomputedMHStatesEquivalence is a broad integration sanity check
// over buildTestTSSB's mixed fixture (co-located and non-co-located
// CNV-affected SSMs, plus plain SSMs). It used to compare the precomputed
// fast path against a slow tree-walking fallback path; that fallback was
// deleted as part of the M1 fix (logLikelihoodWithCNVTreeMHPrecomputed and
// the old fallback computed two DIFFERENT, inequivalent reductions —
// data.py's dynamic nv>0 filter vs mh.hpp's static-collapse + nr+nv>0
// guard — so "equivalence" between them was never actually the right
// invariant to assert; MHStateValid is now a hard precondition instead,
// see logLikelihoodWithCNVTreeMH). The deep numeric checks for the new
// reduction now live in TestComputeSSMStates_StaticCollapse_LOHColocated,
// TestLogLikelihoodWithCNVTreeMHPrecomputed_LOHFixture, and
// TestLogLikelihoodWithCNVTreeMHPrecomputed_TwoStateInvariantUnderCollapse.
// This test just confirms the whole pipeline (a mixed tree with both
// UseFourStates=true and =false SSMs) runs to completion with finite
// values and is deterministic across repeated precompute calls.
func TestPrecomputedMHStatesEquivalence(t *testing.T) {
	rng := rand.New(rand.NewSource(12345))
	tssb, ssms, _ := buildTestTSSB(t, rng)

	tssb.precomputeMHStates()

	sawUseFourStates := false
	firstOld := make([]float64, len(ssms))
	firstNew := make([]float64, len(ssms))
	for i, ssm := range ssms {
		if len(ssm.CNVs) > 0 && ssm.UseFourStates {
			sawUseFourStates = true
		}
		old := logLikelihoodWithCNVTreeMH(ssm, tssb, false)
		nw := logLikelihoodWithCNVTreeMH(ssm, tssb, true)
		if math.IsNaN(old) || math.IsInf(old, 0) {
			t.Errorf("SSM %s: old-state LLH is %v, want finite", ssm.ID, old)
		}
		if math.IsNaN(nw) || math.IsInf(nw, 0) {
			t.Errorf("SSM %s: new-state LLH is %v, want finite", ssm.ID, nw)
		}
		firstOld[i], firstNew[i] = old, nw
	}
	if !sawUseFourStates {
		t.Fatal("fixture has no UseFourStates=true (co-located) SSM — buildTestTSSB should exercise both the M1 4-state and 2-state paths")
	}

	// Determinism: re-running precomputeMHStates without changing Pi/Pi1
	// must reproduce identical values.
	tssb.precomputeMHStates()
	for i, ssm := range ssms {
		old := logLikelihoodWithCNVTreeMH(ssm, tssb, false)
		nw := logLikelihoodWithCNVTreeMH(ssm, tssb, true)
		if old != firstOld[i] {
			t.Errorf("SSM %s: old-state LLH not deterministic across precompute calls: %.15e vs %.15e", ssm.ID, firstOld[i], old)
		}
		if nw != firstNew[i] {
			t.Errorf("SSM %s: new-state LLH not deterministic across precompute calls: %.15e vs %.15e", ssm.ID, firstNew[i], nw)
		}
	}
}

// TestComputeSSMStates_StaticCollapse_LOHColocated hand-derives the
// expected post-collapse MHStates for the canonical M1 scenario: an SSM
// co-located with its own zero-copy-allele (LOH, minor_cn=0) CNV. This is
// the independent oracle — the expected values below are derived directly
// from computeSSMStates' Case 4 branch and the params.py:213-219 static
// collapse rule, on paper, not by running any Go code path and reading off
// its output.
//
// Fixture: two-node tree, root -> nodeX. The SSM and its own CNV both sit
// at nodeX (cp=2 paternal, cm=0 maternal — LOH). Pi[0]=0.3 at root, 0.7 at
// nodeX.
//
// Per-node states before collapse (Case 4, "SSM before CNV: variant
// amplified" branch fires at nodeX since isAncestorOfCached(nodeX, nodeX)
// is true by self-inclusion):
//
//	root (Case 1, diploid reference): Nr1..4=2, Nv1..4=0
//	nodeX (Case 4): Nr1=cm=0,Nv1=cp=2 ; Nr2=cp=2,Nv2=cm=0
//	               Nr3=Nr4=max(0,totalCN-1)=1, Nv3=Nv4=min(1,totalCN)=1
//
// Static collapse (params.py:213-219), aggregated with Pi[0]:
//
//	nv1Total = 0.3*0 + 0.7*2 = 1.4   (nonzero)
//	nv2Total = 0.3*0 + 0.7*0 = 0     (ZERO -> every node's state2 := state1)
//
// So nodeX's (Nr2,Nv2) becomes (Nr1,Nv1) = (0,2), overwriting the
// case-4-computed (2,0). Root's (Nr2,Nv2) was already equal to (Nr1,Nv1)
// so the copy there is a no-op. UseFourStates=true (co-located), so
// states 3/4 are NOT touched by the (!UseFourStates) duplication branch —
// they keep their Case-4 timing-derived values.
func TestComputeSSMStates_StaticCollapse_LOHColocated(t *testing.T) {
	root := &Node{ID: 0, Pi: []float64{0.3}}
	nodeX := &Node{ID: 1, Parent: root, Pi: []float64{0.7}}
	root.Children = []*Node{nodeX}

	cnv := &CNV{ID: "c0", Node: nodeX}
	ssm := &SSM{
		ID: "s0", A: []int{1}, D: []int{1}, MuR: 0.999, MuV: 0.5,
		Node: nodeX,
		CNVs: []*CNVRef{{CNV: cnv, MaternalCN: 0, PaternalCN: 2}},
	}

	computeSSMStates(ssm, []*Node{root, nodeX})

	if !ssm.UseFourStates {
		t.Fatal("expected UseFourStates=true for a co-located SSM+CNV")
	}

	var rootState, xState *SSMNodeState
	for i := range ssm.MHStates {
		switch ssm.MHStates[i].Node {
		case root:
			rootState = &ssm.MHStates[i]
		case nodeX:
			xState = &ssm.MHStates[i]
		}
	}
	if rootState == nil || xState == nil {
		t.Fatal("missing expected node states")
	}

	check := func(label string, got, want float64) {
		t.Helper()
		if got != want {
			t.Errorf("%s = %v, want %v", label, got, want)
		}
	}
	// nodeX state 1 must be untouched by the collapse.
	check("nodeX.Nr1", xState.Nr1, 0)
	check("nodeX.Nv1", xState.Nv1, 2)
	// nodeX state 2 must have collapsed to equal state 1 — this pins the
	// collapse DIRECTION: nv1Total != 0 so state1 survives, and it's
	// state2 (whose aggregate variant mass was zero) that gets
	// overwritten, not the other way around.
	check("nodeX.Nr2 (post-collapse)", xState.Nr2, 0)
	check("nodeX.Nv2 (post-collapse)", xState.Nv2, 2)
	// States 3/4 keep their Case-4 timing-derived values (UseFourStates
	// means the duplication branch is skipped).
	check("nodeX.Nr3", xState.Nr3, 1)
	check("nodeX.Nv3", xState.Nv3, 1)
	check("nodeX.Nr4", xState.Nr4, 1)
	check("nodeX.Nv4", xState.Nv4, 1)

	// Root (Case 1, diploid reference) is untouched: its state1==state2
	// already, so the collapse copy is a no-op there.
	check("root.Nr1", rootState.Nr1, 2)
	check("root.Nv1", rootState.Nv1, 0)
	check("root.Nr2", rootState.Nr2, 2)
	check("root.Nv2", rootState.Nv2, 0)
}

// TestLogLikelihoodWithCNVTreeMHPrecomputed_LOHFixture reuses the fixture
// and hand-derivation from TestComputeSSMStates_StaticCollapse_LOHColocated
// and checks the FULL reduction (aggregation + nr+nv>0 guard + fixed
// log(0.25) prior + logsumexp), not just the collapsed per-node states.
// The expected (nr,nv) pairs below are the Pi-weighted aggregates of the
// hand-derived per-node states from that test:
//
//	(nr1,nv1) = 0.3*(2,0) + 0.7*(0,2) = (0.6, 1.4)
//	(nr2,nv2) = 0.3*(2,0) + 0.7*(0,2) = (0.6, 1.4)   [state2, now == state1]
//	(nr3,nv3) = 0.3*(2,0) + 0.7*(1,1) = (1.3, 0.7)
//	(nr4,nv4) = 0.3*(2,0) + 0.7*(1,1) = (1.3, 0.7)   [== state3]
//
// The expected LLH is computed by an independently-written reduction
// (transcribed from mh.hpp: guard nr+nv>0, fixed log(0.25) prior, logsumexp
// over exactly these four pairs) rather than by calling
// logLikelihoodWithCNVTreeMHPrecomputed itself.
func TestLogLikelihoodWithCNVTreeMHPrecomputed_LOHFixture(t *testing.T) {
	root := &Node{ID: 0, Pi: []float64{0.3}}
	nodeX := &Node{ID: 1, Parent: root, Pi: []float64{0.7}}
	root.Children = []*Node{nodeX}

	cnv := &CNV{ID: "c0", Node: nodeX}
	ssm := &SSM{
		ID:              "s0",
		A:               []int{65},
		D:               []int{100},
		MuR:             0.999,
		MuV:             0.5,
		LogBinNormConst: []float64{logBinCoeff(100, 65)},
		Node:            nodeX,
		CNVs:            []*CNVRef{{CNV: cnv, MaternalCN: 0, PaternalCN: 2}},
	}

	computeSSMStates(ssm, []*Node{root, nodeX})
	if !ssm.MHStateValid {
		t.Fatal("computeSSMStates did not set MHStateValid")
	}

	pairs := [4][2]float64{{0.6, 1.4}, {0.6, 1.4}, {1.3, 0.7}, {1.3, 0.7}}
	const prior = -1.3862943611198906 // math.Log(0.25)
	var lls [4]float64
	for i, p := range pairs {
		nr, nv := p[0], p[1]
		mu := (nr*ssm.MuR + nv*(1-ssm.MuR)) / (nr + nv)
		lls[i] = logBinomialLikelihood(ssm.A[0], ssm.D[0], mu) + prior + ssm.LogBinNormConst[0]
	}
	want := logsumexp(lls[:])

	got := logLikelihoodWithCNVTreeMHPrecomputed(ssm, false)
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("got %.15e, want %.15e (diff %.2e)", got, want, math.Abs(got-want))
	}
}

// TestLogLikelihoodWithCNVTreeMHPrecomputed_TwoStateInvariantUnderCollapse
// checks the algebraic invariant the M1 fix depends on for the common
// (non-co-located) case: when UseFourStates=false, computeSSMStates
// duplicates states 1/2 into 3/4 rather than collapsing them (neither
// nv1Total nor nv2Total is zero in this fixture, so only the duplication
// fires, not the zero-mass collapse). That means the new unconditional
// four-term reduction, logsumexp([b1,b2,b1,b2]) + log(0.25), must equal
// the old two-term form, logsumexp([b1,b2]) + log(0.5), since
// log(2*e^b1+2*e^b2) + log(0.25) = logsumexp([b1,b2]) + log(2) + log(0.25)
// = logsumexp([b1,b2]) + log(0.5). This is the regression guard that the
// M1 fix doesn't perturb the majority of non-CNV-colocated data.
//
// Fixture: root -> P (SSM assigned here) -> G (CNV assigned here, a
// descendant of P, so cnv.Node != ssm.Node -> UseFourStates=false
// globally). At G, ssmInAncestors is true (P is an ancestor of G) and
// cnvNode(=G) is a descendant of ssmNode(=P), so the "SSM before CNV:
// variant amplified" branch fires there with cp=2, cm=1 (both nonzero, so
// neither Nv1Total nor Nv2Total collapses to zero).
func TestLogLikelihoodWithCNVTreeMHPrecomputed_TwoStateInvariantUnderCollapse(t *testing.T) {
	root := &Node{ID: 0, Pi: []float64{0.2}}
	p := &Node{ID: 1, Parent: root, Pi: []float64{0.5}}
	g := &Node{ID: 2, Parent: p, Pi: []float64{0.3}}
	root.Children = []*Node{p}
	p.Children = []*Node{g}

	cnv := &CNV{ID: "c0", Node: g}
	ssm := &SSM{
		ID:              "s0",
		A:               []int{40},
		D:               []int{100},
		MuR:             0.999,
		MuV:             0.5,
		LogBinNormConst: []float64{logBinCoeff(100, 40)},
		Node:            p, // NOT g -> UseFourStates=false
		CNVs:            []*CNVRef{{CNV: cnv, MaternalCN: 1, PaternalCN: 2}},
	}

	computeSSMStates(ssm, []*Node{root, p, g})
	if ssm.UseFourStates {
		t.Fatal("expected UseFourStates=false (SSM and CNV are on different nodes)")
	}

	// Aggregate (nr1,nv1) and (nr2,nv2) directly from the collapsed
	// MHStates, computed by the SAME per-node loop the production
	// function uses (this part is not independent of computeSSMStates —
	// see TestComputeSSMStates_StaticCollapse_LOHColocated for the
	// from-scratch hand derivation of a related fixture). What IS
	// independent here is the reduction formula itself: the two-term
	// reference below is written separately from
	// logLikelihoodWithCNVTreeMHPrecomputed's four-term implementation.
	var nr1, nv1, nr2, nv2 float64
	for i := range ssm.MHStates {
		s := &ssm.MHStates[i]
		pi := s.Node.Pi[0]
		nr1 += pi * s.Nr1
		nv1 += pi * s.Nv1
		nr2 += pi * s.Nr2
		nv2 += pi * s.Nv2
	}
	if nv1 == 0 || nv2 == 0 {
		t.Fatalf("fixture invariant violated: expected both nv1Total and nv2Total nonzero (no collapse), got nv1=%v nv2=%v", nv1, nv2)
	}

	mu1 := (nr1*ssm.MuR + nv1*(1-ssm.MuR)) / (nr1 + nv1)
	mu2 := (nr2*ssm.MuR + nv2*(1-ssm.MuR)) / (nr2 + nv2)
	b1 := logBinomialLikelihood(ssm.A[0], ssm.D[0], mu1) + ssm.LogBinNormConst[0]
	b2 := logBinomialLikelihood(ssm.A[0], ssm.D[0], mu2) + ssm.LogBinNormConst[0]
	want := logsumexp([]float64{b1, b2}) + math.Log(0.5)

	got := logLikelihoodWithCNVTreeMHPrecomputed(ssm, false)
	if math.Abs(got-want) > 1e-12 {
		t.Errorf("4-term reduction diverges from 2-term reference: got %.15e, want %.15e (diff %.2e)", got, want, math.Abs(got-want))
	}
}
