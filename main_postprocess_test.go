package main

import (
	"encoding/json"
	"reflect"
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
