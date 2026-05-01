package main

import (
	"encoding/json"
	"math"
	"reflect"
	"sort"
	"testing"
)

// TestPostProcessSummary_RemovesEmptyPops verifies that empty populations
// (num_ssms=0, num_cnvs=0) are removed and remaining pops are renumbered
// contiguously, matching Python's remove_empty_nodes + _summarize_pops.
func TestPostProcessSummary_RemovesEmptyPops(t *testing.T) {
	// Build a summary with 5 populations:
	//   0 (root, empty — always kept)
	//   1 (has 3 SSMs)
	//   2 (empty internal node, children=[3,4])
	//   3 (has 2 SSMs)
	//   4 (empty leaf)
	//
	// Structure: 0→[1,2], 2→[3,4]
	//
	// After post-processing:
	//   0 (root) kept
	//   1 (was 1, has SSMs) → renumbered 1
	//   2 (was empty) removed, children reparented to 0
	//   3 (was 3, has SSMs) → renumbered 2
	//   4 (was empty leaf) removed
	//
	// Expected structure: 0→[1,2]

	type popSummary struct {
		CellPrev []float64 `json:"cellular_prevalence"`
		NumSSMs  int       `json:"num_ssms"`
		NumCNVs  int       `json:"num_cnvs"`
	}

	pops := map[string]popSummary{
		"0": {CellPrev: []float64{1.0}, NumSSMs: 0, NumCNVs: 0},
		"1": {CellPrev: []float64{0.8}, NumSSMs: 3, NumCNVs: 0},
		"2": {CellPrev: []float64{0.5}, NumSSMs: 0, NumCNVs: 0},
		"3": {CellPrev: []float64{0.3}, NumSSMs: 2, NumCNVs: 0},
		"4": {CellPrev: []float64{0.1}, NumSSMs: 0, NumCNVs: 0},
	}
	structure := map[string][]int{
		"0": {1, 2},
		"2": {3, 4},
	}
	mutAss := map[string]interface{}{
		"1": map[string]interface{}{"ssms": []string{"s0", "s1", "s2"}, "cnvs": []string{}},
		"3": map[string]interface{}{"ssms": []string{"s3", "s4"}, "cnvs": []string{}},
	}

	// Marshal to match summarizePops output format
	popsJSON, _ := json.Marshal(pops)
	structJSON, _ := json.Marshal(structure)
	mutAssJSON, _ := json.Marshal(mutAss)

	var popsIface map[string]interface{}
	var structIface map[string]interface{}
	var mutAssIface map[string]interface{}
	json.Unmarshal(popsJSON, &popsIface)
	json.Unmarshal(structJSON, &structIface)
	json.Unmarshal(mutAssJSON, &mutAssIface)

	summary := map[string]interface{}{
		"chain_id":        0,
		"llh":             -1234.5,
		"num_populations": 2,
		"populations":     popsIface,
		"structure":       structIface,
		"mut_assignments": mutAssIface,
	}

	result := postProcessSummary(summary)

	// Verify num_populations updated
	if got := result["num_populations"]; got != 2 {
		t.Errorf("num_populations = %v, want 2", got)
	}

	// Verify only 3 populations remain (0, 1, 2)
	resultPops := result["populations"].(map[string]interface{})
	if len(resultPops) != 3 {
		t.Errorf("got %d populations, want 3", len(resultPops))
	}
	// Population "0" must exist (root always kept)
	if _, ok := resultPops["0"]; !ok {
		t.Error("root population '0' missing")
	}
	// Old pop "3" should now be "2"
	pop2 := resultPops["2"].(map[string]interface{})
	if numSSMs := pop2["num_ssms"]; numSSMs != float64(2) {
		t.Errorf("pop '2' num_ssms = %v, want 2", numSSMs)
	}

	// Verify structure: 0→[1,2] only
	resultStruct := result["structure"].(map[string]interface{})
	root0 := resultStruct["0"]
	expected0 := []interface{}{float64(1), float64(2)}
	if !reflect.DeepEqual(root0, expected0) {
		t.Errorf("structure[0] = %v, want %v", root0, expected0)
	}
	// No other structure entries (1 and 2 are leaves)
	if len(resultStruct) != 1 {
		t.Errorf("got %d structure entries, want 1", len(resultStruct))
	}

	// Verify mut_assignments renumbered
	resultMutAss := result["mut_assignments"].(map[string]interface{})
	if _, ok := resultMutAss["1"]; !ok {
		t.Error("mut_assignments['1'] missing (was pop 1)")
	}
	if _, ok := resultMutAss["2"]; !ok {
		t.Error("mut_assignments['2'] missing (was pop 3)")
	}
	if _, ok := resultMutAss["3"]; ok {
		t.Error("mut_assignments['3'] should not exist after renumbering")
	}
}

// TestPostProcessSummary_AllNonEmpty verifies no-op when all pops have data.
func TestPostProcessSummary_AllNonEmpty(t *testing.T) {
	pops := map[string]interface{}{
		"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": float64(0), "num_cnvs": float64(0)},
		"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.8}, "num_ssms": float64(5), "num_cnvs": float64(0)},
		"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.3}, "num_ssms": float64(3), "num_cnvs": float64(0)},
	}
	structure := map[string]interface{}{
		"0": []interface{}{float64(1), float64(2)},
	}
	mutAss := map[string]interface{}{
		"1": map[string]interface{}{"ssms": []interface{}{"s0"}, "cnvs": []interface{}{}},
		"2": map[string]interface{}{"ssms": []interface{}{"s1"}, "cnvs": []interface{}{}},
	}

	summary := map[string]interface{}{
		"chain_id":        0,
		"llh":             -100.0,
		"num_populations": 2,
		"populations":     pops,
		"structure":       structure,
		"mut_assignments": mutAss,
	}

	result := postProcessSummary(summary)

	resultPops := result["populations"].(map[string]interface{})
	if len(resultPops) != 3 {
		t.Errorf("got %d populations, want 3 (no-op)", len(resultPops))
	}
}

// ── removeSuperclones tests ─────────────────────────────────────────────────

// buildSummary is a helper that constructs a summary map from typed inputs,
// doing a JSON round-trip to get map[string]interface{} throughout.
func buildSummary(t *testing.T, pops map[string]interface{}, structure map[string]interface{}, mutAss map[string]interface{}) map[string]interface{} {
	t.Helper()
	raw, err := json.Marshal(map[string]interface{}{
		"chain_id":        0,
		"llh":             -100.0,
		"populations":     pops,
		"structure":       structure,
		"mut_assignments": mutAss,
	})
	if err != nil {
		t.Fatal(err)
	}
	var out map[string]interface{}
	if err := json.Unmarshal(raw, &out); err != nil {
		t.Fatal(err)
	}
	return out
}

func getPopSSMs(t *testing.T, summary map[string]interface{}, popIdx string) int {
	t.Helper()
	pops := summary["populations"].(map[string]interface{})
	pop := pops[popIdx].(map[string]interface{})
	return int(pop["num_ssms"].(float64))
}

func getStructChildren(t *testing.T, summary map[string]interface{}, popIdx string) []int {
	t.Helper()
	structure := summary["structure"].(map[string]interface{})
	kids, ok := structure[popIdx]
	if !ok {
		return nil
	}
	var result []int
	for _, k := range kids.([]interface{}) {
		result = append(result, int(k.(float64)))
	}
	return result
}

func getPopCP(t *testing.T, summary map[string]interface{}, popIdx string) []float64 {
	t.Helper()
	pops := summary["populations"].(map[string]interface{})
	pop := pops[popIdx].(map[string]interface{})
	cpRaw := pop["cellular_prevalence"].([]interface{})
	var cp []float64
	for _, v := range cpRaw {
		cp = append(cp, v.(float64))
	}
	return cp
}

func countMutAssSSMs(t *testing.T, summary map[string]interface{}, popIdx string) int {
	t.Helper()
	ma := summary["mut_assignments"].(map[string]interface{})
	entry, ok := ma[popIdx]
	if !ok {
		return 0
	}
	ssms := entry.(map[string]interface{})["ssms"].([]interface{})
	return len(ssms)
}

// TestRemoveSuperclones_Triggers tests the superclone pattern:
// 0→[1], 1→[2], pop1 has ≤33% SSMs of pop2, CPs within 10%.
// Expected: pop2 is removed, its SSMs merged into pop1, pop1 CP becomes
// weighted average, children of pop2 become children of pop1.
func TestRemoveSuperclones_Triggers(t *testing.T) {
	// Structure: 0→[1], 1→[2], 2→[3]
	// Pop 1: 2 SSMs, CP=[0.9, 0.85]
	// Pop 2: 10 SSMs, CP=[0.88, 0.83] (within 10% of pop1, 2/10=0.2 ≤ 0.33)
	// Pop 3: 5 SSMs, CP=[0.5, 0.4]
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0, 1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.9, 0.85}, "num_ssms": 2, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.88, 0.83}, "num_ssms": 10, "num_cnvs": 0},
			"3": map[string]interface{}{"cellular_prevalence": []interface{}{0.5, 0.4}, "num_ssms": 5, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1},
			"1": []interface{}{2},
			"2": []interface{}{3},
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s0", "s1"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11"}, "cnvs": []interface{}{}},
			"3": map[string]interface{}{"ssms": []interface{}{"s12", "s13", "s14", "s15", "s16"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSuperclones(summary)

	// Pop 1 should now have all 12 SSMs (2+10)
	if got := getPopSSMs(t, result, "1"); got != 12 {
		t.Errorf("pop 1 num_ssms = %d, want 12", got)
	}
	if got := countMutAssSSMs(t, result, "1"); got != 12 {
		t.Errorf("pop 1 mut_assignments ssms = %d, want 12", got)
	}

	// Pop 1's CP should be weighted average: (2*[0.9,0.85] + 10*[0.88,0.83]) / 12
	expectedCP := []float64{
		(2*0.9 + 10*0.88) / 12.0,
		(2*0.85 + 10*0.83) / 12.0,
	}
	gotCP := getPopCP(t, result, "1")
	for i, exp := range expectedCP {
		if math.Abs(gotCP[i]-exp) > 1e-6 {
			t.Errorf("pop 1 CP[%d] = %f, want %f", i, gotCP[i], exp)
		}
	}

	// Old pop 2 removed → old pop 3 renumbered to 2
	pops := result["populations"].(map[string]interface{})
	if len(pops) != 3 { // 0, 1, 2
		t.Errorf("got %d populations, want 3", len(pops))
	}
	if got := getPopSSMs(t, result, "2"); got != 5 {
		t.Errorf("pop 2 (was 3) num_ssms = %d, want 5", got)
	}

	// Structure: 0→[1], 1→[2]
	if got := getStructChildren(t, result, "0"); !reflect.DeepEqual(got, []int{1}) {
		t.Errorf("structure[0] = %v, want [1]", got)
	}
	if got := getStructChildren(t, result, "1"); !reflect.DeepEqual(got, []int{2}) {
		t.Errorf("structure[1] = %v, want [2]", got)
	}
}

// TestRemoveSuperclones_NoTrigger_MultiChild tests that superclone removal
// does NOT trigger when the clonal node has multiple children.
func TestRemoveSuperclones_NoTrigger_MultiChild(t *testing.T) {
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.9}, "num_ssms": 2, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.5}, "num_ssms": 10, "num_cnvs": 0},
			"3": map[string]interface{}{"cellular_prevalence": []interface{}{0.3}, "num_ssms": 5, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1},
			"1": []interface{}{2, 3}, // multiple children — no superclone
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s0", "s1"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s2"}, "cnvs": []interface{}{}},
			"3": map[string]interface{}{"ssms": []interface{}{"s3"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSuperclones(summary)
	pops := result["populations"].(map[string]interface{})
	if len(pops) != 4 {
		t.Errorf("got %d populations, want 4 (no-op)", len(pops))
	}
}

// TestRemoveSuperclones_NoTrigger_HighRatio tests that superclone removal
// does NOT trigger when the SSM ratio is > 0.33.
func TestRemoveSuperclones_NoTrigger_HighRatio(t *testing.T) {
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.9}, "num_ssms": 5, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.88}, "num_ssms": 10, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1},
			"1": []interface{}{2},
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s0"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s1"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSuperclones(summary)
	pops := result["populations"].(map[string]interface{})
	if len(pops) != 3 {
		t.Errorf("got %d populations, want 3 (no-op, ratio=0.5 > 0.33)", len(pops))
	}
}

// TestRemoveSuperclones_NoTrigger_CPTooFar tests that superclone removal
// does NOT trigger when CPs differ by more than 10%.
func TestRemoveSuperclones_NoTrigger_CPTooFar(t *testing.T) {
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.9}, "num_ssms": 1, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.5}, "num_ssms": 10, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1},
			"1": []interface{}{2},
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s0"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s1"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSuperclones(summary)
	pops := result["populations"].(map[string]interface{})
	if len(pops) != 3 {
		t.Errorf("got %d populations, want 3 (no-op, CP diff=0.4 > 0.1)", len(pops))
	}
}

// ── removeSmallNodes tests ──────────────────────────────────────────────────

// TestRemoveSmallNodes_RemovesBelowThreshold tests that populations with
// fewer SSMs than ceil(minFrac * totalSSMs) are removed and their mutations
// reassigned to the closest-CP population.
func TestRemoveSmallNodes_RemovesBelowThreshold(t *testing.T) {
	// 4 pops: root(0), pop1(50 SSMs), pop2(1 SSM), pop3(49 SSMs)
	// Total SSMs = 100, threshold at 0.05 = 5
	// Pop2 has 1 SSM < 5 → removed, its SSM reassigned to closest CP
	// Pop2 CP=0.78, Pop1 CP=0.8, Pop3 CP=0.3 → reassign to pop1 (closest)
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.8}, "num_ssms": 50, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.78}, "num_ssms": 1, "num_cnvs": 0},
			"3": map[string]interface{}{"cellular_prevalence": []interface{}{0.3}, "num_ssms": 49, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1, 2, 3},
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s0"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s_small"}, "cnvs": []interface{}{}},
			"3": map[string]interface{}{"ssms": []interface{}{"s1"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSmallNodes(summary, 0.05)

	pops := result["populations"].(map[string]interface{})
	if len(pops) != 3 { // 0, 1, 2 (was 3)
		t.Errorf("got %d populations, want 3", len(pops))
	}

	// Pop1 should now have the reassigned SSM from pop2
	if got := countMutAssSSMs(t, result, "1"); got != 2 { // s0 + s_small
		t.Errorf("pop 1 mut_assignments ssms = %d, want 2", got)
	}
	if got := getPopSSMs(t, result, "1"); got != 2 {
		t.Errorf("pop 1 num_ssms = %d, want 2", got)
	}
}

// TestRemoveSmallNodes_NoOp tests that nothing changes when all pops are above threshold.
func TestRemoveSmallNodes_NoOp(t *testing.T) {
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.8}, "num_ssms": 50, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.3}, "num_ssms": 50, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1, 2},
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s0"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s1"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSmallNodes(summary, 0.01)
	pops := result["populations"].(map[string]interface{})
	if len(pops) != 3 {
		t.Errorf("got %d populations, want 3 (no-op)", len(pops))
	}
}

// TestRemoveSmallNodes_ReparentsChildren tests that when a small internal node
// is removed, its children are reparented to the grandparent.
func TestRemoveSmallNodes_ReparentsChildren(t *testing.T) {
	// Structure: 0→[1], 1→[2,3]
	// Pop1 has 1 SSM (below 5% of 100), pop2 has 60, pop3 has 39
	// Pop1 removed → pop2 and pop3 become children of root
	summary := buildSummary(t,
		map[string]interface{}{
			"0": map[string]interface{}{"cellular_prevalence": []interface{}{1.0}, "num_ssms": 0, "num_cnvs": 0},
			"1": map[string]interface{}{"cellular_prevalence": []interface{}{0.8}, "num_ssms": 1, "num_cnvs": 0},
			"2": map[string]interface{}{"cellular_prevalence": []interface{}{0.5}, "num_ssms": 60, "num_cnvs": 0},
			"3": map[string]interface{}{"cellular_prevalence": []interface{}{0.3}, "num_ssms": 39, "num_cnvs": 0},
		},
		map[string]interface{}{
			"0": []interface{}{1},
			"1": []interface{}{2, 3},
		},
		map[string]interface{}{
			"1": map[string]interface{}{"ssms": []interface{}{"s_small"}, "cnvs": []interface{}{}},
			"2": map[string]interface{}{"ssms": []interface{}{"s0"}, "cnvs": []interface{}{}},
			"3": map[string]interface{}{"ssms": []interface{}{"s1"}, "cnvs": []interface{}{}},
		},
	)

	result := removeSmallNodes(summary, 0.05) // threshold = 5

	pops := result["populations"].(map[string]interface{})
	if len(pops) != 3 { // 0, 1(was 2), 2(was 3)
		t.Errorf("got %d populations, want 3", len(pops))
	}

	// Root should now have 2 children (renumbered 1 and 2)
	rootKids := getStructChildren(t, result, "0")
	sort.Ints(rootKids)
	if !reflect.DeepEqual(rootKids, []int{1, 2}) {
		t.Errorf("structure[0] = %v, want [1, 2]", rootKids)
	}
}

// ============================================================================
// filterChainsByInclusion tests
// ============================================================================

func TestFilterChains_AllIncluded(t *testing.T) {
	// Three chains with similar LLH — all should be included with factor=1.1
	results := []ChainResult{
		{SampleLLH: []float64{-100, -99, -101}},
		{SampleLLH: []float64{-100, -98, -102}},
		{SampleLLH: []float64{-101, -100, -99}},
	}
	included, excluded := filterChainsByInclusion(results, 1.1)
	if len(included) != 3 {
		t.Errorf("expected 3 included, got %d: %v", len(included), included)
	}
	if len(excluded) != 0 {
		t.Errorf("expected 0 excluded, got %d: %v", len(excluded), excluded)
	}
}

func TestFilterChains_OneBadChainExcluded(t *testing.T) {
	// Chain 2 has much worse LLH — should be excluded
	results := []ChainResult{
		{SampleLLH: []float64{-100, -99, -101}},
		{SampleLLH: []float64{-100, -98, -102}},
		{SampleLLH: []float64{-10000, -9999, -10001}}, // terrible chain
	}
	included, excluded := filterChainsByInclusion(results, 1.1)
	if len(included) != 2 {
		t.Errorf("expected 2 included, got %d: %v", len(included), included)
	}
	if len(excluded) != 1 || excluded[0] != 2 {
		t.Errorf("expected excluded=[2], got %v", excluded)
	}
}

func TestFilterChains_FactorInf_IncludesAll(t *testing.T) {
	// factor=Inf should include all chains regardless of quality
	results := []ChainResult{
		{SampleLLH: []float64{-100}},
		{SampleLLH: []float64{-10000}},
	}
	included, excluded := filterChainsByInclusion(results, math.Inf(1))
	if len(included) != 2 {
		t.Errorf("expected 2 included with Inf factor, got %d", len(included))
	}
	if len(excluded) != 0 {
		t.Errorf("expected 0 excluded, got %d", len(excluded))
	}
}

func TestFilterChains_Factor1_OnlyBest(t *testing.T) {
	// factor=1.0 should include only the chain(s) with the best logSumLLH
	results := []ChainResult{
		{SampleLLH: []float64{-100, -99}},   // logSumExp ≈ -98.69
		{SampleLLH: []float64{-200, -199}},   // logSumExp ≈ -198.69
		{SampleLLH: []float64{-100, -99}},    // same as chain 0
	}
	included, excluded := filterChainsByInclusion(results, 1.0)
	// Chains 0 and 2 tie, both should be included
	if len(included) != 2 {
		t.Errorf("expected 2 included (tied best), got %d: %v", len(included), included)
	}
	if len(excluded) != 1 || excluded[0] != 1 {
		t.Errorf("expected excluded=[1], got %v", excluded)
	}
}

func TestFilterChains_SingleChain(t *testing.T) {
	results := []ChainResult{
		{SampleLLH: []float64{-500, -499}},
	}
	included, excluded := filterChainsByInclusion(results, 1.1)
	if len(included) != 1 || included[0] != 0 {
		t.Errorf("single chain should always be included, got included=%v", included)
	}
	if len(excluded) != 0 {
		t.Errorf("expected 0 excluded, got %v", excluded)
	}
}

func TestFilterChains_EmptySampleLLH(t *testing.T) {
	// Chain with no samples should be excluded
	results := []ChainResult{
		{SampleLLH: []float64{-100, -99}},
		{SampleLLH: []float64{}},
	}
	included, excluded := filterChainsByInclusion(results, 1.1)
	if len(included) != 1 || included[0] != 0 {
		t.Errorf("expected only chain 0 included, got %v", included)
	}
	if len(excluded) != 1 || excluded[0] != 1 {
		t.Errorf("expected chain 1 excluded, got %v", excluded)
	}
}
