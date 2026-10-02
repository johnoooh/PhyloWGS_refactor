package main

import (
	"archive/zip"
	"encoding/json"
	"os"
	"path/filepath"
	"reflect"
	"testing"
	"time"
)

func loadTinyDataset(t *testing.T) ([]*SSM, []*CNV) {
	t.Helper()
	tmp := t.TempDir()
	ssmPath := filepath.Join(tmp, "ssm.txt")
	cnvPath := filepath.Join(tmp, "cnv.txt")
	if err := os.WriteFile(ssmPath, []byte(
		"id\tgene\ta\td\tmu_r\tmu_v\n"+
			"s0\tg\t90\t100\t0.999\t0.5\n"+
			"s1\tg\t88\t100\t0.999\t0.5\n"+
			"s2\tg\t70\t100\t0.999\t0.5\n"+
			"s3\tg\t72\t100\t0.999\t0.5\n",
	), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(cnvPath, []byte(
		"cnv\ta\td\tssms\tphysical_cnvs\n"+
			"c0\t900\t1000\ts2,1,2;s3,1,2\tchrom=1,start=1,end=100,major_cn=2,minor_cn=1,cell_prev=0.2\n",
	), 0644); err != nil {
		t.Fatal(err)
	}
	ssms, err := loadSSMData(ssmPath)
	if err != nil {
		t.Fatal(err)
	}
	cnvs, err := loadCNVData(cnvPath, ssms)
	if err != nil {
		t.Fatal(err)
	}
	return ssms, cnvs
}

// A generous timeout that is never hit must not perturb the chain: the budget
// check draws no random numbers and touches no state.
func TestRunChain_UnhitTimeoutIsBitIdentical(t *testing.T) {
	ssms, cnvs := loadTinyDataset(t)
	a := runChain(0, ssms, cnvs, 5, 10, 50, 7, 0)
	b := runChain(0, ssms, cnvs, 5, 10, 50, 7, time.Hour)
	if a.TimedOut || b.TimedOut {
		t.Fatalf("unexpected timeout: a=%v b=%v", a.TimedOut, b.TimedOut)
	}
	if !reflect.DeepEqual(a.SampleLLH, b.SampleLLH) || !reflect.DeepEqual(a.BurninLLH, b.BurninLLH) {
		t.Fatalf("timeout=1h changed the trajectory:\n%v\n%v", a.SampleLLH, b.SampleLLH)
	}
}

func TestRunChain_TimeoutStopsChain(t *testing.T) {
	ssms, cnvs := loadTinyDataset(t)
	r := runChain(3, ssms, cnvs, 5, 1000000, 50, 7, time.Nanosecond)
	if !r.TimedOut {
		t.Fatal("expected TimedOut")
	}
	if n := len(r.BurninLLH) + len(r.SampleLLH); n > 1 {
		t.Fatalf("chain ran %d iterations past a 1ns budget", n)
	}
}

func fakeChain(id int, llh []float64, timedOut bool) ChainResult {
	return ChainResult{ChainID: id, SampleLLH: llh, TimedOut: timedOut}
}

func TestSelectChains_NoTimeoutMatchesFilter(t *testing.T) {
	rs := []ChainResult{
		fakeChain(0, []float64{-100}, false),
		fakeChain(1, []float64{-105}, false),
		fakeChain(2, []float64{-500}, false),
	}
	i1, e1 := filterChainsByInclusion(rs, 1.1)
	i2, e2 := selectChains(rs, 1.1)
	if !reflect.DeepEqual(i1, i2) || !reflect.DeepEqual(e1, e2) {
		t.Fatalf("selectChains diverged from filterChainsByInclusion: %v/%v vs %v/%v", i1, e1, i2, e2)
	}
}

// A timed-out chain is excluded even when its LLH is the best, and it does not
// set the inclusion threshold for the completed chains.
func TestSelectChains_ExcludesTimedOutAndRecomputesThreshold(t *testing.T) {
	rs := []ChainResult{
		fakeChain(0, []float64{-1000}, false),
		fakeChain(1, []float64{-10}, true), // best LLH but timed out
		fakeChain(2, []float64{-1050}, false),
	}
	inc, exc := selectChains(rs, 1.1)
	if !reflect.DeepEqual(inc, []int{0, 2}) || !reflect.DeepEqual(exc, []int{1}) {
		t.Fatalf("got included=%v excluded=%v", inc, exc)
	}
}

func TestSelectChains_AllTimedOutFallsBack(t *testing.T) {
	rs := []ChainResult{
		fakeChain(0, []float64{-100}, true),
		fakeChain(1, []float64{-500}, true),
	}
	inc, _ := selectChains(rs, 1.1)
	if !reflect.DeepEqual(inc, []int{0}) {
		t.Fatalf("got included=%v", inc)
	}
}

// Regression for the straggler failure mode: one chain that cannot finish must
// not prevent the completed chains' results from being written.
func TestWriteResults_WithTimedOutChain(t *testing.T) {
	ssms, cnvs := loadTinyDataset(t)
	done := runChain(0, ssms, cnvs, 2, 5, 20, 11, 0)
	// Simulate a straggler cut off after producing some samples (a real
	// timeout mid-sampling keeps its partial Trees/SampleLLH).
	stuck := runChain(1, ssms, cnvs, 2, 4, 20, 12, 0)
	stuck.TimedOut = true
	if len(stuck.Trees) == 0 || done.TimedOut {
		t.Fatalf("setup: stuck has %d trees, done.TimedOut=%v", len(stuck.Trees), done.TimedOut)
	}
	out := t.TempDir()
	if err := writeResults(out, []ChainResult{done, stuck}, 1.1, "T"); err != nil {
		t.Fatal(err)
	}
	raw, err := os.ReadFile(filepath.Join(out, "summary.json"))
	if err != nil {
		t.Fatal(err)
	}
	var summ struct {
		TimedOut []int `json:"timed_out_chains"`
	}
	if err := json.Unmarshal(raw, &summ); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(summ.TimedOut, []int{1}) {
		t.Fatalf("summary timed_out_chains=%v", summ.TimedOut)
	}
	zr, err := zip.OpenReader(filepath.Join(out, "trees.zip"))
	if err != nil {
		t.Fatal(err)
	}
	defer zr.Close()
	if len(zr.File) != len(done.Trees) {
		t.Fatalf("trees.zip has %d entries, want %d (completed chain only)", len(zr.File), len(done.Trees))
	}
	if _, err := os.Stat(filepath.Join(out, "best_tree.json")); err != nil {
		t.Fatal(err)
	}
}

// Default output is unchanged: summary.json has no timed_out_chains key when
// no chain timed out.
func TestWriteResults_NoTimeoutKeyByDefault(t *testing.T) {
	ssms, cnvs := loadTinyDataset(t)
	done := runChain(0, ssms, cnvs, 2, 3, 20, 11, 0)
	out := t.TempDir()
	if err := writeResults(out, []ChainResult{done}, 1.1, "T"); err != nil {
		t.Fatal(err)
	}
	raw, _ := os.ReadFile(filepath.Join(out, "summary.json"))
	var m map[string]interface{}
	if err := json.Unmarshal(raw, &m); err != nil {
		t.Fatal(err)
	}
	if _, ok := m["timed_out_chains"]; ok {
		t.Fatal("timed_out_chains present without any timeout")
	}
}
