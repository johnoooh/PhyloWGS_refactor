# mutlist.json Output Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Parse the `physical_cnvs` column from `cnv_data.txt` and write a `mutlist.json` file to the output directory matching the Python `write_results.py` / `result_generator.py` format.

**Architecture:** Two-part change. (1) Extend `loadCNVData` to parse the optional 5th column into a new `PhysicalCNV` slice stored on the `CNV` struct. (2) Add a `writeMutList` function called at the end of `writeResults`, using the best chain's `FinalTree` to traverse SSM/CNV data and emit `mutlist.json`.

**Tech Stack:** Go 1.21, standard library only (`encoding/json`, `os`, `strings`, `strconv`). Build with `/tmp/go/bin/go build ./...` from `PhyloWGS_refactor/`.

---

## Background: Python `mutlist.json` Format

From `phylowgs/pwgsresults/result_generator.py` and `phylowgs/util2.py`, the output JSON has this shape:

```json
{
  "ssms": {
    "s0": {
      "ref_reads": [60000, 40383],
      "total_reads": [60000, 63000],
      "expected_ref_in_ref": 0.999,
      "expected_ref_in_variant": 0.5,
      "name": "optional_gene_name"
    }
  },
  "cnvs": {
    "c0": {
      "ref_reads": [60000, 40383],
      "total_reads": [60000, 63000],
      "physical_cnvs": [
        {
          "chrom": "1",
          "start": 141500000,
          "end": 148899999,
          "major_cn": 2,
          "minor_cn": 1,
          "cell_prev": "0.0|0.718"
        }
      ],
      "ssms": [
        {"ssm_id": "s2", "maternal_cn": 1, "paternal_cn": 2}
      ]
    }
  }
}
```

Key points:
- `expected_ref_in_ref` = SSM's `MuR` field
- `expected_ref_in_variant` = SSM's `MuV` field
- `physical_cnvs` is an array of segments parsed from the 5th TSV column
- `ssms` within each CNV = list of `{ssm_id, maternal_cn, paternal_cn}` (the `AffectedSSMs` data with CN values)
- `name` field on SSMs: include it (SSM's `Name` field; may be empty string)

## Background: `physical_cnvs` Column Format

From `phylowgs/parser/test/outputs/multisamp_cnvs/cnv_data.txt`:
```
chrom=1,start=141500000,end=148899999,major_cn=2,minor_cn=1,cell_prev=0.0|0.718
```
Multiple segments are separated by `;`. Each segment is a comma-separated list of `key=value` pairs. The `cell_prev` value contains `|` separators (one per timepoint) — keep it as a raw string.

## Background: SSM-to-CNV CN Values

In `loadCNVData` (line 437-448), when an SSM ref is parsed from the `ssms` column (`s2,1,2`), `parts[1]` → `MaternalCN` and `parts[2]` → `PaternalCN`. These same values need to appear in the `mutlist.json` CNV's `ssms` list. The `CNV` struct currently stores only `AffectedSSMs []string` (SSM IDs). We need to also store the per-SSM CN values there.

---

## Task 1: Add `PhysicalCNV` Struct and Extend `CNV` Struct

**Files:**
- Modify: `PhyloWGS_refactor/main.go:79–87` (CNV struct block)

**Step 1: Add the new struct and extend CNV**

In `main.go`, after line 77 (closing brace of `CNVRef` struct), add a `PhysicalCNV` struct. Then extend the `CNV` struct.

Replace this block (lines 79–87):
```go
// CNV represents a Copy Number Variation
type CNV struct {
	ID              string
	A               []int
	D               []int
	LogBinNormConst []float64
	Node            *Node
	AffectedSSMs    []string // SSM IDs affected by this CNV
}
```

With:
```go
// PhysicalCNV holds genomic annotation for one physical CNV segment.
// Parsed from the physical_cnvs column of cnv_data.txt.
// Fields map directly to the Python cnv_logical_physical_mapping output.
type PhysicalCNV struct {
	Chrom    string `json:"chrom"`
	Start    int    `json:"start"`
	End      int    `json:"end"`
	MajorCN  int    `json:"major_cn"`
	MinorCN  int    `json:"minor_cn"`
	CellPrev string `json:"cell_prev"` // raw "0.0|0.718" string, one value per timepoint
}

// CNVSSMLink holds the SSM ID and copy numbers for one SSM affected by a CNV.
// Used when writing mutlist.json.
type CNVSSMLink struct {
	SSMID      string
	MaternalCN int
	PaternalCN int
}

// CNV represents a Copy Number Variation
type CNV struct {
	ID              string
	A               []int
	D               []int
	LogBinNormConst []float64
	Node            *Node
	AffectedSSMs    []string     // SSM IDs affected by this CNV
	SSMLinks        []CNVSSMLink // SSM IDs with their CN values (for mutlist.json)
	PhysicalCNVs    []PhysicalCNV
}
```

**Step 2: Build to confirm no syntax errors**

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go build ./... 2>&1
```
Expected: no output (clean build).

**Step 3: Commit**

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor
git add main.go
git commit -m "feat: add PhysicalCNV struct and extend CNV with SSMLinks + PhysicalCNVs fields"
```

---

## Task 2: Parse `physical_cnvs` Column in `loadCNVData`

**Files:**
- Modify: `PhyloWGS_refactor/main.go:432–451` (SSM ref parsing block inside `loadCNVData`)

**Step 1: Extend SSM ref parsing to also populate `SSMLinks`**

In the SSM ref parsing loop (lines 432–451), after `cnv.AffectedSSMs = append(...)`, add:
```go
cnv.SSMLinks = append(cnv.SSMLinks, CNVSSMLink{
    SSMID:      ssmID,
    MaternalCN: maternalCN,
    PaternalCN: paternalCN,
})
```

The full updated block (lines 432–451) becomes:
```go
// Parse SSM references
if len(fields) > 3 && fields[3] != "" {
    ssmRefs := strings.Split(fields[3], ";")
    for _, ref := range ssmRefs {
        parts := strings.Split(ref, ",")
        if len(parts) >= 3 {
            ssmID := parts[0]
            maternalCN, _ := strconv.Atoi(parts[1])
            paternalCN, _ := strconv.Atoi(parts[2])
            cnv.AffectedSSMs = append(cnv.AffectedSSMs, ssmID)
            cnv.SSMLinks = append(cnv.SSMLinks, CNVSSMLink{
                SSMID:      ssmID,
                MaternalCN: maternalCN,
                PaternalCN: paternalCN,
            })
            if ssm, ok := ssmMap[ssmID]; ok {
                ssm.CNVs = append(ssm.CNVs, &CNVRef{
                    CNV:        cnv,
                    MaternalCN: maternalCN,
                    PaternalCN: paternalCN,
                })
            }
        }
    }
}
```

**Step 2: Add `physical_cnvs` column parsing**

Immediately after the SSM ref parsing block (after the closing `}` of the `if len(fields) > 3` block, around line 451), add:

```go
// Parse physical_cnvs column (5th column, optional)
if len(fields) > 4 && fields[4] != "" {
    cnv.PhysicalCNVs = parsePhysicalCNVs(fields[4])
}
```

**Step 3: Add `parsePhysicalCNVs` helper function**

Add this new function anywhere before `loadCNVData` (e.g., right before it at line 383):

```go
// parsePhysicalCNVs parses the physical_cnvs column from cnv_data.txt.
// Format: "chrom=1,start=141500000,end=148899999,major_cn=2,minor_cn=1,cell_prev=0.0|0.718"
// Multiple segments are separated by ";".
func parsePhysicalCNVs(raw string) []PhysicalCNV {
    var result []PhysicalCNV
    segments := strings.Split(raw, ";")
    for _, seg := range segments {
        seg = strings.TrimSpace(seg)
        if seg == "" {
            continue
        }
        var p PhysicalCNV
        for _, kv := range strings.Split(seg, ",") {
            idx := strings.Index(kv, "=")
            if idx < 0 {
                continue
            }
            key := kv[:idx]
            val := kv[idx+1:]
            switch key {
            case "chrom":
                p.Chrom = val
            case "start":
                p.Start, _ = strconv.Atoi(val)
            case "end":
                p.End, _ = strconv.Atoi(val)
            case "major_cn":
                p.MajorCN, _ = strconv.Atoi(val)
            case "minor_cn":
                p.MinorCN, _ = strconv.Atoi(val)
            case "cell_prev":
                p.CellPrev = val
            }
        }
        result = append(result, p)
    }
    return result
}
```

**Step 4: Build**

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go build ./... 2>&1
```
Expected: clean build.

**Step 5: Commit**

```bash
git add main.go
git commit -m "feat: parse physical_cnvs column and SSMLinks in loadCNVData"
```

---

## Task 3: Write `writeMutList` Function

**Files:**
- Modify: `PhyloWGS_refactor/main.go` — add new function before `writeResults` (before line 2225)

**Step 1: Add `writeMutList` function**

Insert the following function at line 2224 (immediately before `func writeResults`):

```go
// writeMutList writes mutlist.json to outDir, matching the Python
// write_results.py / result_generator.py output format.
// ssms and cnvs come from the best chain's FinalTree.
func writeMutList(outDir string, ssms []*SSM, cnvs []*CNV) error {
    // Build SSM map for looking up per-CNV SSM annotations
    _ = ssms // used via direct field access below

    // --- SSM section ---
    type ssmEntry struct {
        RefReads              []int   `json:"ref_reads"`
        TotalReads            []int   `json:"total_reads"`
        ExpectedRefInRef      float64 `json:"expected_ref_in_ref"`
        ExpectedRefInVariant  float64 `json:"expected_ref_in_variant"`
        Name                  string  `json:"name,omitempty"`
    }
    ssmsOut := make(map[string]ssmEntry, len(ssms))
    for _, s := range ssms {
        ssmsOut[s.ID] = ssmEntry{
            RefReads:             s.A,
            TotalReads:           s.D,
            ExpectedRefInRef:     s.MuR,
            ExpectedRefInVariant: s.MuV,
            Name:                 s.Name,
        }
    }

    // --- CNV section ---
    type cnvSSMEntry struct {
        SSMID      string `json:"ssm_id"`
        MaternalCN int    `json:"maternal_cn"`
        PaternalCN int    `json:"paternal_cn"`
    }
    type cnvEntry struct {
        RefReads     []int          `json:"ref_reads"`
        TotalReads   []int          `json:"total_reads"`
        PhysicalCNVs []PhysicalCNV  `json:"physical_cnvs"`
        SSMs         []cnvSSMEntry  `json:"ssms"`
    }
    cnvsOut := make(map[string]cnvEntry, len(cnvs))
    for _, c := range cnvs {
        ssmLinks := make([]cnvSSMEntry, 0, len(c.SSMLinks))
        for _, link := range c.SSMLinks {
            ssmLinks = append(ssmLinks, cnvSSMEntry{
                SSMID:      link.SSMID,
                MaternalCN: link.MaternalCN,
                PaternalCN: link.PaternalCN,
            })
        }
        physCNVs := c.PhysicalCNVs
        if physCNVs == nil {
            physCNVs = []PhysicalCNV{} // emit [] not null
        }
        cnvsOut[c.ID] = cnvEntry{
            RefReads:     c.A,
            TotalReads:   c.D,
            PhysicalCNVs: physCNVs,
            SSMs:         ssmLinks,
        }
    }

    out := struct {
        SSMs map[string]ssmEntry `json:"ssms"`
        CNVs map[string]cnvEntry `json:"cnvs"`
    }{
        SSMs: ssmsOut,
        CNVs: cnvsOut,
    }

    data, err := json.MarshalIndent(out, "", "  ")
    if err != nil {
        return err
    }
    return os.WriteFile(filepath.Join(outDir, "mutlist.json"), data, 0644)
}
```

**Step 2: Call `writeMutList` from `writeResults`**

In `writeResults`, after the `best_tree.json` block (after line 2324, inside the `if tssb := results[bestChainIdx].FinalTree; tssb != nil {` block), add a call to `writeMutList`. The updated block becomes:

```go
if tssb := results[bestChainIdx].FinalTree; tssb != nil {
    removeEmptyNodes(tssb.Root, nil)
    treeSummary := summarizePops(tssb, bestOverallLLH, results[bestChainIdx].ChainID)
    treeData, _ := json.MarshalIndent(treeSummary, "", "  ")
    treePath := filepath.Join(outDir, "best_tree.json")
    if err := os.WriteFile(treePath, treeData, 0644); err != nil {
        return err
    }
    if err := writeMutList(outDir, tssb.Data, tssb.CNVData); err != nil {
        return err
    }
}
```

**Step 3: Build**

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor && /tmp/go/bin/go build ./... 2>&1
```
Expected: clean build.

**Step 4: Commit**

```bash
git add main.go
git commit -m "feat: add writeMutList; write mutlist.json matching Python result_generator output"
```

---

## Task 4: Verify on Existing Fixture

**Files:**
- Read: `/tmp/sim_fixtures/` — any existing fixture with CNVs (look for `cnv_data.txt` with content)

**Step 1: Find a fixture with CNVs**

```bash
grep -rl "^c[0-9]" /tmp/sim_fixtures/ --include="cnv_data.txt" | head -5
```

If `/tmp/sim_fixtures/` is empty or has no CNVs, use the reference fixture from the Python parser test suite:

```bash
ls /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/phylowgs/parser/test/outputs/multisamp_cnvs/
```

**Step 2: Run the Go port on the fixture**

```bash
FIXTURE=/tmp/sim_fixtures/<fixture_dir>
/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor/phylowgs_go \
    -ssm $FIXTURE/ssm_data.txt \
    -cnv $FIXTURE/cnv_data.txt \
    -o /tmp/mutlist_verify \
    -b 50 -s 100 -j 1
```

Adjust binary name if needed (check `ls PhyloWGS_refactor/`).

**Step 3: Inspect `mutlist.json`**

```bash
cat /tmp/mutlist_verify/mutlist.json | python3 -m json.tool | head -60
```

Verify:
- Top-level keys are `ssms` and `cnvs`
- Each SSM entry has `ref_reads`, `total_reads`, `expected_ref_in_ref`, `expected_ref_in_variant`
- Each CNV entry has `ref_reads`, `total_reads`, `physical_cnvs` (array, not null), `ssms` (array of `{ssm_id, maternal_cn, paternal_cn}`)
- If using the `multisamp_cnvs` fixture: `c0` should have `physical_cnvs[0].chrom = "1"`, `start = 141500000`, and `ssms = [{ssm_id: "s2", maternal_cn: 1, paternal_cn: 2}]`

**Step 4: Cross-check against Python output (if Python env available)**

If the Singularity container is available:
```bash
singularity exec phylowgs.sif python phylowgs/write_results.py \
    -i $FIXTURE/cnv_data.txt $FIXTURE/ssm_data.txt \
    <trees.zip> /tmp/python_results
diff <(python3 -m json.tool /tmp/mutlist_verify/mutlist.json) \
     <(python3 -m json.tool /tmp/python_results/mutlist.json)
```
Expected: no diff (modulo key ordering).

**Step 5: Commit (if any fixes needed from verification)**

```bash
git add main.go
git commit -m "fix: correct mutlist.json field from verification"
```

---

## Summary of Changes

| Location | Change |
|----------|--------|
| `main.go:79` (before CNV struct) | Add `PhysicalCNV` struct, `CNVSSMLink` struct |
| `main.go:87` (CNV struct) | Add `SSMLinks []CNVSSMLink` and `PhysicalCNVs []PhysicalCNV` fields |
| `main.go:383` (before `loadCNVData`) | Add `parsePhysicalCNVs(raw string) []PhysicalCNV` helper |
| `main.go:441` (inside SSM ref loop) | Populate `cnv.SSMLinks` alongside `cnv.AffectedSSMs` |
| `main.go:451` (after SSM ref block) | Parse `fields[4]` into `cnv.PhysicalCNVs` |
| `main.go:2224` (before `writeResults`) | Add `writeMutList(outDir, ssms, cnvs)` function |
| `main.go:2324` (inside `writeResults`) | Call `writeMutList` after writing `best_tree.json` |
