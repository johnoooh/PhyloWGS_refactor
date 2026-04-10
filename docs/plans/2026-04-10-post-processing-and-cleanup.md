# Post-Processing Integration & Debug Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Wire `removeEmptyNodes` into Go's NDJSON output so post-processed K matches Python's `tree_summaries.json.gz`, optionally add `removeSuperclones`, clean up all 9 investigation debug flags, and write an offline comparison script for the 216 HPC fixtures.

**Architecture:** The Go port already has `removeEmptyNodes()` (main.go:3039-3107) and `summarizePops()` (main.go:3114-3215). We add a thin `postProcessSummary()` that operates on the JSON-level dict returned by `summarizePops` — filtering empty populations and reparenting structure entries — then call it from `snapshotTree()`. Debug flag removal is purely mechanical deletion. The comparison script is standalone Python that reads existing Go NDJSON + Python `trees.zip` and compares K distributions.

**Tech Stack:** Go 1.21 (`/tmp/go/bin/go`), module `phylowgs-go`, Python 3 (for comparison script)

---

### Task 1: Add `postProcessSummary()` function

This function takes the `map[string]interface{}` returned by `summarizePops` and removes empty populations (num_ssms == 0 && num_cnvs == 0, except population "0" which is always kept as the normal-cell root). It also reparents children of removed nodes in the `structure` dict and renumbers all population IDs to be contiguous.

This matches the combined effect of Python's `remove_empty_vertices=True` in `result_generator.py:39` → `util2.py:128-155`.

**Files:**
- Modify: `main.go` (add function after `summarizePops`, before `snapshotTree`)
- Test: `main_postprocess_test.go` (new file)

**Step 1: Write the failing test**

Create `main_postprocess_test.go`:

```go
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
```

**Step 2: Run test to verify it fails**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go test -run TestPostProcessSummary -v -count=1`
Expected: FAIL — `postProcessSummary` undefined

**Step 3: Write minimal implementation**

Add to `main.go` between `summarizePops` (after line 3215) and `snapshotTree` (line 3222):

```go
// postProcessSummary filters empty populations from the JSON-level summary
// returned by summarizePops, reparents their children, and renumbers all
// population IDs to be contiguous starting from 0. This matches the combined
// effect of Python's remove_empty_vertices=True (util2.py:128-155) applied
// at the result_generator level (result_generator.py:39).
//
// Population "0" (the normal-cell root) is always kept even if empty.
func postProcessSummary(summary map[string]interface{}) map[string]interface{} {
	popsRaw := summary["populations"].(map[string]interface{})
	structRaw := summary["structure"].(map[string]interface{})
	mutAssRaw := summary["mut_assignments"].(map[string]interface{})

	// Find max pop index
	maxIdx := 0
	for k := range popsRaw {
		idx, _ := strconv.Atoi(k)
		if idx > maxIdx {
			maxIdx = idx
		}
	}

	// Identify empty populations (except root "0")
	emptySet := make(map[int]bool)
	for k, v := range popsRaw {
		idx, _ := strconv.Atoi(k)
		if idx == 0 {
			continue // Always keep root
		}
		pop := v.(map[string]interface{})
		numSSMs := pop["num_ssms"].(float64)
		numCNVs := pop["num_cnvs"].(float64)
		if numSSMs == 0 && numCNVs == 0 {
			emptySet[idx] = true
		}
	}

	if len(emptySet) == 0 {
		// Nothing to remove — return as-is
		return summary
	}

	// Build parent map from structure
	parentOf := make(map[int]int)
	childrenOf := make(map[int][]int)
	for k, v := range structRaw {
		parentIdx, _ := strconv.Atoi(k)
		kids := v.([]interface{})
		for _, kid := range kids {
			childIdx := int(kid.(float64))
			parentOf[childIdx] = parentIdx
			childrenOf[parentIdx] = append(childrenOf[parentIdx], childIdx)
		}
	}

	// Remove empty nodes from structure: reparent children to grandparent
	// Process in reverse order (leaves first) to handle chains of empty nodes
	for idx := maxIdx; idx >= 1; idx-- {
		if !emptySet[idx] {
			continue
		}
		parent, hasParent := parentOf[idx]
		if !hasParent {
			continue
		}
		// Remove idx from parent's children
		newKids := []int{}
		for _, k := range childrenOf[parent] {
			if k != idx {
				newKids = append(newKids, k)
			}
		}
		// Reparent idx's children to parent
		for _, grandchild := range childrenOf[idx] {
			newKids = append(newKids, grandchild)
			parentOf[grandchild] = parent
		}
		childrenOf[parent] = newKids
		delete(childrenOf, idx)
	}

	// Build renumbering map (contiguous IDs, skipping removed nodes)
	renumber := make(map[int]int) // old → new
	newIdx := 0
	for idx := 0; idx <= maxIdx; idx++ {
		if emptySet[idx] {
			continue
		}
		if _, exists := popsRaw[strconv.Itoa(idx)]; !exists {
			continue
		}
		renumber[idx] = newIdx
		newIdx++
	}

	// Build new populations
	newPops := make(map[string]interface{})
	for oldIdx, newIdx := range renumber {
		newPops[strconv.Itoa(newIdx)] = popsRaw[strconv.Itoa(oldIdx)]
	}

	// Build new structure
	newStruct := make(map[string]interface{})
	for parentIdx, kids := range childrenOf {
		if emptySet[parentIdx] {
			continue
		}
		newParent, ok := renumber[parentIdx]
		if !ok {
			continue
		}
		// Renumber and filter children
		newKids := []interface{}{}
		for _, kid := range kids {
			if emptySet[kid] {
				continue // shouldn't happen after reparenting, but safety
			}
			if newKidIdx, ok := renumber[kid]; ok {
				newKids = append(newKids, float64(newKidIdx))
			}
		}
		if len(newKids) > 0 {
			newStruct[strconv.Itoa(newParent)] = newKids
		}
	}

	// Build new mut_assignments
	newMutAss := make(map[string]interface{})
	for k, v := range mutAssRaw {
		oldIdx, _ := strconv.Atoi(k)
		if newIdx, ok := renumber[oldIdx]; ok {
			newMutAss[strconv.Itoa(newIdx)] = v
		}
	}

	// Count non-empty (non-root) populations
	nonEmpty := 0
	for k, v := range newPops {
		if k == "0" {
			continue
		}
		pop := v.(map[string]interface{})
		if pop["num_ssms"].(float64) > 0 || pop["num_cnvs"].(float64) > 0 {
			nonEmpty++
		}
	}

	// Return new summary
	result := make(map[string]interface{})
	for k, v := range summary {
		result[k] = v
	}
	result["populations"] = newPops
	result["structure"] = newStruct
	result["mut_assignments"] = newMutAss
	result["num_populations"] = nonEmpty

	return result
}
```

**Step 4: Run test to verify it passes**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go test -run TestPostProcessSummary -v -count=1`
Expected: PASS

**Step 5: Commit**

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor
git add main.go main_postprocess_test.go
git commit -m "feat: add postProcessSummary to filter empty pops from NDJSON output

Matches Python's remove_empty_vertices=True (util2.py:128-155).
Removes empty populations (except root), reparents children,
renumbers contiguously."
```

---

### Task 2: Wire `postProcessSummary` into `snapshotTree`

**Files:**
- Modify: `main.go:3222-3233` (`snapshotTree` function)

**Step 1: Write the failing test**

Add to `main_postprocess_test.go`:

```go
// TestSnapshotTreeAppliesPostProcessing is an integration-level check
// that snapshotTree's JSON output has no empty non-root populations.
// This requires a real TSSB, so we build a minimal one.
func TestSnapshotTreeAppliesPostProcessing(t *testing.T) {
	// Create a minimal TSSB with one empty node and one non-empty node.
	rng := rand.New(rand.NewSource(42))
	ssms := []*SSM{
		{ID: "s0", Name: "s0", A: []int{10}, D: []int{20}, MuR: 0.999, MuV: 0.5, LogBinNormConst: []float64{0}},
	}
	tssb := newTSSB(ssms, nil, 25.0, 1.0, 0.25, rng)

	snap := snapshotTree(tssb, -100.0, 0, 0)
	if snap == nil {
		t.Fatal("snapshotTree returned nil")
	}

	var result map[string]interface{}
	if err := json.Unmarshal(snap, &result); err != nil {
		t.Fatalf("failed to unmarshal snapshot: %v", err)
	}

	pops := result["populations"].(map[string]interface{})
	for k, v := range pops {
		if k == "0" {
			continue // root is always kept
		}
		pop := v.(map[string]interface{})
		numSSMs := pop["num_ssms"].(float64)
		numCNVs := pop["num_cnvs"].(float64)
		if numSSMs == 0 && numCNVs == 0 {
			t.Errorf("population %s is empty but was not removed by post-processing", k)
		}
	}
}
```

**Step 2: Run test to verify it fails**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go test -run TestSnapshotTreeAppliesPostProcessing -v -count=1`
Expected: FAIL — empty populations present in output

**Step 3: Modify `snapshotTree` to call `postProcessSummary`**

In `main.go`, change `snapshotTree` (lines 3222-3233) from:

```go
func snapshotTree(tssb *TSSB, llh float64, iteration, chainID int) json.RawMessage {
	summary := summarizePops(tssb, llh, chainID)
	// Stamp iteration into the snapshot so each line is self-contained
	summary["iteration"] = iteration
	data, err := json.Marshal(summary)
```

to:

```go
func snapshotTree(tssb *TSSB, llh float64, iteration, chainID int) json.RawMessage {
	summary := summarizePops(tssb, llh, chainID)
	summary = postProcessSummary(summary)
	// Stamp iteration into the snapshot so each line is self-contained
	summary["iteration"] = iteration
	data, err := json.Marshal(summary)
```

Also update the comment block above `snapshotTree` to remove the note about skipping `removeEmptyNodes`.

**Step 4: Run test to verify it passes**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go test -run TestSnapshotTree -v -count=1`
Expected: PASS

**Step 5: Commit**

```bash
git add main.go main_postprocess_test.go
git commit -m "feat: wire postProcessSummary into snapshotTree

NDJSON output now removes empty populations and renumbers
contiguously, matching Python's tree_summaries.json.gz format."
```

---

### Task 3: Remove 8 debug flags and --trace infrastructure

These 9 items were added during the K-inflation investigation and are no longer needed:

| # | Global variable | Flag name | Lines (declaration) |
|---|----------------|-----------|---------------------|
| 1 | `debugClampAlphaDecay` | `--clamp-alpha-decay` | 29 |
| 2 | `debugMHReseed` | `--mh-reseed` | 30 |
| 3 | `debugMHFloat32` | `--mh-float32` | 31 |
| 4 | `debugMHDumpFile` | `--mh-dump` | 32 |
| 5 | `debugMHAcceptAll` | `--mh-accept-all` | 33 |
| 6 | `debugMHRejectAll` | `--mh-reject-all` | 34 |
| 7 | `debugMHTruncPrec` | `--mh-trunc-prec` | 35 |
| 8 | `debugMHFixedStd` | `--mh-fixed-std` | 36 |
| 9 | `Trace`/`TraceCounters`/`--trace` | `--trace` | 178-197, 202-211 |

**Files:**
- Modify: `main.go`
- Delete: `trace_test.go` (if it only tests trace infrastructure)

**Step 1: Remove global variable declarations (lines 29-36)**

Delete these 8 lines from the `var` block:
```go
debugClampAlphaDecay float64 = -1
debugMHReseed        bool
debugMHFloat32       bool
debugMHDumpFile      string
debugMHAcceptAll     bool
debugMHRejectAll     bool
debugMHTruncPrec     bool
debugMHFixedStd      float64 = -1
```

**Step 2: Remove `Trace` field from TSSB struct (line 181)**

Delete:
```go
Trace TraceCounters
```

**Step 3: Remove `TraceCounters` struct (lines 184-197) and `countAllNodes` (lines 199-211)**

Delete entire struct definition and function.

**Step 4: Remove debug blocks from `metropolis()` (lines 1542-1548, 1572-1721, 1726-1731)**

In `metropolis()`:
- Delete lines 1542-1548 (debugMHReseed + debugMHFloat32 blocks)
- Delete lines 1572-1721 (entire dumpDone/doDump/dumpData infrastructure — the `dumpDone` declaration, all `doDump` conditionals, and the dump JSON write block)
- Delete lines 1726-1731 (debugMHAcceptAll/debugMHRejectAll overrides)

After cleanup, `metropolis()` should have a clean accept/reject that is simply:
```go
doAccept := math.Log(rng.Float64()) < logA
if doAccept {
```

**Step 5: Remove trace counter increments (lines 2386, 2470)**

Delete:
```go
t.Trace.SpawnsInStickOrders++
```
and:
```go
t.Trace.SpawnsInFindNode++
```

**Step 6: Remove `debugClampAlphaDecay` usage (lines 2571-2572)**

Delete:
```go
if debugClampAlphaDecay >= 0 {
    t.AlphaDecay = debugClampAlphaDecay
}
```

**Step 7: Remove trace infrastructure from `runChain()` (lines 2642-2660, 2728-2731, 2780, 2782, 2788-2816)**

Delete:
- Lines 2642-2660: traceW/traceFH declaration and file open block
- Lines 2728-2731: kBeforeCull/kAfterCull variables (keep `tssb.cullTree()` call!)
- Lines 2780, 2782: kBeforeStickOrders/kAfterStickOrders variables
- Lines 2788-2816: entire trace record emission block

**Step 8: Remove `debugMHTruncPrec` block (lines 2748-2756) and `debugMHFixedStd` block (lines 2767-2768)**

Delete the `debugMHTruncPrec` block. For `debugMHFixedStd`, change:
```go
if debugMHFixedStd >= 0 {
    mhStd = debugMHFixedStd
} else {
    if mhAcc < 0.08 && mhStd < 10000 {
```
to just:
```go
if mhAcc < 0.08 && mhStd < 10000 {
```
(remove the outer if/else wrapper)

**Step 9: Remove flag definitions and copying (lines 3249-3268)**

Delete:
```go
traceFile := flag.String("trace", ...)
clampAlphaDecay := flag.Float64("clamp-alpha-decay", ...)
mhReseed := flag.Bool("mh-reseed", ...)
mhFloat32 := flag.Bool("mh-float32", ...)
mhDump := flag.String("mh-dump", ...)
mhAcceptAll := flag.Bool("mh-accept-all", ...)
mhRejectAll := flag.Bool("mh-reject-all", ...)
mhTruncPrec := flag.Bool("mh-trunc-prec", ...)
mhFixedStd := flag.Float64("mh-fixed-std", ...)
```
and the copy block:
```go
debugClampAlphaDecay = *clampAlphaDecay
debugMHReseed = *mhReseed
...
debugMHFixedStd = *mhFixedStd
```

**Step 10: Remove trace file routing in chain launch (lines 3350-3357)**

Change the goroutine body from:
```go
chainTrace := ""
if *traceFile != "" {
    if *chains == 1 {
        chainTrace = *traceFile
    } else {
        chainTrace = fmt.Sprintf("%s.chain%d", *traceFile, chainID)
    }
}
results[chainID] = runChain(chainID, ssms, cnvs, *burnin, *samples, *mhIters, chainSeed, chainTrace)
```
to:
```go
results[chainID] = runChain(chainID, ssms, cnvs, *burnin, *samples, *mhIters, chainSeed)
```

And update `runChain` signature to remove `traceFile string` parameter.

**Step 11: Remove `trace_test.go` if it only tests trace infra**

Check the file — if it tests TraceCounters or trace NDJSON, delete it. If it tests other things, keep the non-trace parts.

**Step 12: Remove unused imports**

After all deletions, check for unused imports (`bufio`, `os` if only used by trace/dump, `strconv` if only used by debugMHTruncPrec). The compiler will flag these.

**Step 13: Run full test suite**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go build ./... && /tmp/go/bin/go test -v -count=1 ./...`
Expected: BUILD SUCCESS, all tests PASS

**Step 14: Commit**

```bash
git add -u
git commit -m "chore: remove 9 investigation debug flags and trace infrastructure

Removes debugClampAlphaDecay, debugMHReseed, debugMHFloat32,
debugMHDumpFile, debugMHAcceptAll, debugMHRejectAll, debugMHTruncPrec,
debugMHFixedStd, and the --trace NDJSON diagnostic system.
These were added during the K-inflation investigation (round 2)
and are no longer needed."
```

---

### Task 4: Write offline comparison script

This Python script reads existing Go NDJSON files and Python `trees.zip` pickles from the HPC results, applies `remove_empty_nodes` to the Python pickles (same as `result_generator.py`), and compares the K distributions across all 216 common fixtures.

**No rerun needed** — the existing Go NDJSON already contains full population data (including empty nodes). The script applies `postProcessSummary`-equivalent logic offline.

**Files:**
- Create: `simulation_validation4/compare_k_postprocessed.py`

**Step 1: Write the comparison script**

```python
#!/usr/bin/env python3
"""Compare post-processed K between Go NDJSON and Python trees.zip.

Reads existing HPC results (no rerun needed):
  - Go: simulation_validation4/results/go-cpu/*/chain_*_trees.ndjson
  - Python: simulation_validation4/results/original-python/*/trees.zip

For Go: counts non-empty populations (num_ssms > 0 or num_cnvs > 0) per tree.
  Since the Go NDJSON includes ALL nodes (empty + non-empty), we filter here.

For Python: unpickles TSSB from trees.zip, applies remove_empty_nodes, then
  counts populations with data. This matches result_generator.py:39.

Outputs a TSV comparing median K_pop for each fixture.
"""
import json
import gzip
import os
import sys
import glob
import statistics
import pickle
import zipfile
import numpy as np
from pathlib import Path

RESULTS_DIR = Path(__file__).parent / "results"
GO_DIR = RESULTS_DIR / "go-cpu"
PY_DIR = RESULTS_DIR / "original-python"


def go_k_from_ndjson(ndjson_path):
    """Count non-empty populations per tree line in Go NDJSON."""
    k_values = []
    with open(ndjson_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            pops = obj.get("populations", {})
            non_empty = 0
            for idx_str, pop in pops.items():
                if idx_str == "0":
                    continue  # Skip root
                if pop.get("num_ssms", 0) > 0 or pop.get("num_cnvs", 0) > 0:
                    non_empty += 1
            k_values.append(non_empty)
    return k_values


def python_k_from_trees_zip(trees_zip_path):
    """Extract K_pop from Python trees.zip pickles after remove_empty_nodes."""
    k_values = []
    try:
        zf = zipfile.ZipFile(trees_zip_path, 'r')
    except (zipfile.BadZipFile, FileNotFoundError):
        return k_values

    # Find tree pickle files
    tree_files = sorted([n for n in zf.namelist() if n.startswith('tree_') and n.endswith('.pkl')])
    if not tree_files:
        # Try alternative naming
        tree_files = sorted([n for n in zf.namelist() if 'tree' in n.lower() and '.pkl' in n.lower()])

    for tf in tree_files:
        try:
            data = zf.read(tf)
            tree = pickle.loads(data, encoding='latin1')
            # Count non-empty nodes via traversal
            non_empty = count_non_empty_after_removal(tree)
            k_values.append(non_empty)
        except Exception:
            continue

    zf.close()
    return k_values


def count_non_empty_after_removal(tssb_root):
    """Count non-empty populations after conceptually applying remove_empty_nodes.

    Instead of modifying the tree, we just count nodes that have data and
    would survive removal (i.e., non-empty nodes, which are always kept).
    """
    count = 0
    def traverse(node):
        nonlocal count
        if hasattr(node, 'get_data') and len(node.get_data()) > 0:
            count += 1
        for child in node.children():
            traverse(child)
    try:
        traverse(tssb_root.root['node'])
    except (AttributeError, TypeError):
        pass
    return count


def main():
    # Find common fixtures
    go_fixtures = set(os.listdir(GO_DIR)) if GO_DIR.exists() else set()
    py_fixtures = set(os.listdir(PY_DIR)) if PY_DIR.exists() else set()
    common = sorted(go_fixtures & py_fixtures)

    if not common:
        print("No common fixtures found!", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(common)} common fixtures", file=sys.stderr)
    print("fixture\tgo_median_k\tpy_median_k\tratio\tgo_n\tpy_n")

    for fixture in common:
        # Go: aggregate across all chains
        go_k_all = []
        ndjson_files = glob.glob(str(GO_DIR / fixture / "chain_*_trees.ndjson"))
        for nf in ndjson_files:
            go_k_all.extend(go_k_from_ndjson(nf))

        # Python: from trees.zip
        py_k_all = []
        trees_zip = PY_DIR / fixture / "trees.zip"
        if trees_zip.exists():
            py_k_all = python_k_from_trees_zip(str(trees_zip))

        if not go_k_all or not py_k_all:
            go_med = statistics.median(go_k_all) if go_k_all else "NA"
            py_med = statistics.median(py_k_all) if py_k_all else "NA"
            print(f"{fixture}\t{go_med}\t{py_med}\tNA\t{len(go_k_all)}\t{len(py_k_all)}")
            continue

        go_med = statistics.median(go_k_all)
        py_med = statistics.median(py_k_all)
        ratio = go_med / py_med if py_med > 0 else float('inf')
        print(f"{fixture}\t{go_med:.1f}\t{py_med:.1f}\t{ratio:.3f}\t{len(go_k_all)}\t{len(py_k_all)}")


if __name__ == "__main__":
    main()
```

**Step 2: Run the comparison**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork && python3 simulation_validation4/compare_k_postprocessed.py > simulation_validation4/k_comparison_postprocessed.tsv`

**Step 3: Analyze results**

```bash
# Quick summary statistics
python3 -c "
import pandas as pd
df = pd.read_csv('simulation_validation4/k_comparison_postprocessed.tsv', sep='\t')
df_valid = df[df['ratio'] != 'NA']
df_valid['ratio'] = df_valid['ratio'].astype(float)
print(f'Fixtures compared: {len(df_valid)}')
print(f'Median ratio (Go/Py): {df_valid[\"ratio\"].median():.3f}')
print(f'Mean ratio: {df_valid[\"ratio\"].mean():.3f}')
print(f'Min/Max ratio: {df_valid[\"ratio\"].min():.3f} / {df_valid[\"ratio\"].max():.3f}')
print(f'Fixtures within 20%: {(df_valid[\"ratio\"].between(0.8, 1.2)).sum()} / {len(df_valid)}')
"
```

Expected: ratio near 1.0 (Go ≈ Python after post-processing)

**Step 4: Commit**

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork
git add simulation_validation4/compare_k_postprocessed.py simulation_validation4/k_comparison_postprocessed.tsv
git commit -m "validation: add offline K comparison script and results

Compares post-processed K between Go NDJSON and Python trees.zip
across 216 HPC fixtures. No rerun needed — reads existing data."
```

---

### Task 5: Verify build and all tests pass

**Step 1: Full build**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go build ./...`
Expected: clean build, no errors

**Step 2: Full test suite**

Run: `cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go test -v -count=1 ./...`
Expected: all tests PASS

**Step 3: Verify NDJSON output format manually**

Run a quick smoke test with a small fixture to confirm NDJSON lines no longer contain empty populations:

```bash
# Pick a small fixture
FIXTURE=$(ls simulation_validation4/results/go-cpu/ | head -1)
NDJSON=$(ls simulation_validation4/results/go-cpu/$FIXTURE/chain_*_trees.ndjson | head -1)
# Check first line for empty pops
python3 -c "
import json
with open('$NDJSON') as f:
    obj = json.loads(f.readline())
pops = obj['populations']
empty = [k for k,v in pops.items() if k != '0' and v['num_ssms'] == 0 and v['num_cnvs'] == 0]
print(f'Total pops: {len(pops)}, Empty non-root: {len(empty)}')
if empty:
    print('WARNING: Empty populations still present — postProcessSummary not applied')
else:
    print('OK: No empty non-root populations')
"
```

Note: Existing NDJSON files were written BEFORE this change. Only new runs will have post-processed output. The comparison script handles this by filtering on-the-fly.
