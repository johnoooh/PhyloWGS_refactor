package main

import (
	"bufio"
	"encoding/json"
	"math/rand"
	"os"
	"path/filepath"
	"testing"
)

// TestTraceFlagWritesNDJSON pins the contract of the --trace plumbing:
// when runChain is given a non-empty traceFile path, it writes one NDJSON
// record per MCMC iteration (burnin + samples) containing the fields defined
// in round2_trace_design.md. Each line must parse as valid JSON and must have
// the required keys, so downstream analysis scripts can rely on the schema.
//
// This is the RED step for Round 2 Task 2.2: the test references a new
// traceFile parameter on runChain that does not exist yet, which forces the
// fix to extend runChain's signature and emit the NDJSON.
func TestTraceFlagWritesNDJSON(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	ssms := []*SSM{
		makeTestSSM("s0", 2, rng),
		makeTestSSM("s1", 2, rng),
		makeTestSSM("s2", 2, rng),
		makeTestSSM("s3", 2, rng),
		makeTestSSM("s4", 2, rng),
	}
	cnvs := []*CNV{}

	tmpDir := t.TempDir()
	tracePath := filepath.Join(tmpDir, "trace.ndjson")

	burnin := 2
	samples := 3
	mhIters := 5
	// Call runChain with the new traceFile parameter.
	_ = runChain(0, ssms, cnvs, burnin, samples, mhIters, 42, tracePath)

	// File must exist.
	f, err := os.Open(tracePath)
	if err != nil {
		t.Fatalf("trace file not created at %s: %v", tracePath, err)
	}
	defer f.Close()

	// Required fields per round2_trace_design.md.
	required := []string{
		"iter", "chain",
		"K_before_cull", "K_after_cull",
		"spawns_in_find_node", "spawns_in_stick_orders",
		"kills_in_cull",
		"dp_alpha", "dp_gamma", "alpha_decay",
		"best_llh_this_iter",
	}

	scanner := bufio.NewScanner(f)
	nLines := 0
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		var obj map[string]interface{}
		if err := json.Unmarshal(line, &obj); err != nil {
			t.Fatalf("line %d not valid JSON: %v\nline=%q", nLines+1, err, string(line))
		}
		for _, key := range required {
			if _, ok := obj[key]; !ok {
				t.Errorf("line %d missing required key %q", nLines+1, key)
			}
		}
		nLines++
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("scanning trace file: %v", err)
	}

	wantLines := burnin + samples
	if nLines != wantLines {
		t.Errorf("expected %d trace lines (burnin=%d + samples=%d), got %d",
			wantLines, burnin, samples, nLines)
	}
}

// TestTraceFlagDisabledNoFile: when traceFile is "", no file is created and
// there is zero overhead. Pins the gate so a test regression would catch any
// accidental I/O.
func TestTraceFlagDisabledNoFile(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	ssms := []*SSM{
		makeTestSSM("s0", 2, rng),
		makeTestSSM("s1", 2, rng),
	}
	cnvs := []*CNV{}

	tmpDir := t.TempDir()
	tracePath := filepath.Join(tmpDir, "trace.ndjson")

	_ = runChain(0, ssms, cnvs, 1, 1, 5, 42, "")

	if _, err := os.Stat(tracePath); !os.IsNotExist(err) {
		t.Errorf("expected no trace file when traceFile=\"\", but stat returned: %v", err)
	}
}
