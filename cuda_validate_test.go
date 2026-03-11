//go:build cuda
// +build cuda

package main

import (
	"math"
	"math/rand"
	"phylowgs-go/cuda"
	"testing"
)

// TestGPUvsCPUAccuracy validates that GPU and CPU likelihood computations
// produce numerically identical results (within floating point tolerance).
func TestGPUvsCPUAccuracy(t *testing.T) {
	// Skip if CUDA not available
	if !cuda.Available() {
		t.Skip("CUDA not available")
	}

	// Create test SSMs with random data
	rng := rand.New(rand.NewSource(42))
	nSSM := 100
	nTP := 5

	ssms := make([]*SSM, nSSM)
	for i := 0; i < nSSM; i++ {
		ssm := &SSM{
			ID:              string(rune('A' + i)),
			Name:            string(rune('A' + i)),
			A:               make([]int, nTP),
			D:               make([]int, nTP),
			MuR:             0.999,
			MuV:             0.5,
			LogBinNormConst: make([]float64, nTP),
		}
		for tp := 0; tp < nTP; tp++ {
			ssm.D[tp] = rng.Intn(500) + 50 // 50-550 depth
			ssm.A[tp] = rng.Intn(ssm.D[tp])
			ssm.LogBinNormConst[tp] = logBinCoeff(ssm.D[tp], ssm.A[tp])
		}
		ssms[i] = ssm
	}

	// Initialize GPU engine
	engine, err := cuda.NewGPUEngine(nSSM, nTP)
	if err != nil {
		t.Fatalf("Failed to create GPU engine: %v", err)
	}
	defer engine.Close()

	// Convert SSMs to interface
	ssmInterfaces := make([]cuda.SSMDataInterface, len(ssms))
	for i, s := range ssms {
		ssmInterfaces[i] = s
	}

	// Upload static data
	if err := engine.UploadStaticData(ssmInterfaces); err != nil {
		t.Fatalf("Failed to upload static data: %v", err)
	}

	// Test with various phi values
	testCases := []string{
		"uniform_0.5",
		"random",
		"low",
		"high",
		"edge_cases",
	}

	for _, tc := range testCases {
		t.Run(tc, func(t *testing.T) {
			// Generate phi values
			phiBuf := make([]float64, nSSM*nTP)
			phiPerSSM := make([][]float64, nSSM)

			for i := 0; i < nSSM; i++ {
				phiPerSSM[i] = make([]float64, nTP)
				for tp := 0; tp < nTP; tp++ {
					var phi float64
					switch tc {
					case "uniform_0.5":
						phi = 0.5
					case "random":
						phi = rng.Float64()
					case "low":
						phi = rng.Float64() * 0.1
					case "high":
						phi = 0.9 + rng.Float64()*0.1
					case "edge_cases":
						// Mix of edge cases
						switch tp % 5 {
						case 0:
							phi = 0
						case 1:
							phi = 1
						case 2:
							phi = 1e-15
						case 3:
							phi = 1 - 1e-15
						case 4:
							phi = 0.5
						}
					}
					phiBuf[i*nTP+tp] = phi
					phiPerSSM[i][tp] = phi
				}
			}

			// Compute CPU likelihoods
			cpuLLH := make([]float64, nSSM)
			for i, ssm := range ssms {
				cpuLLH[i] = ssm.logLikelihoodNoCNV(phiPerSSM[i])
			}

			// Compute GPU likelihoods
			gpuLLH, err := engine.ComputeBatchFast(phiBuf, nSSM, nTP)
			if err != nil {
				t.Fatalf("GPU computation failed: %v", err)
			}

			// Compare results
			maxAbsErr := 0.0
			maxRelErr := 0.0
			for i := 0; i < nSSM; i++ {
				absErr := math.Abs(gpuLLH[i] - cpuLLH[i])
				var relErr float64
				if cpuLLH[i] != 0 {
					relErr = absErr / math.Abs(cpuLLH[i])
				}

				if absErr > maxAbsErr {
					maxAbsErr = absErr
				}
				if relErr > maxRelErr {
					maxRelErr = relErr
				}

				// Check relative error (1e-10 tolerance)
				if relErr > 1e-10 {
					t.Errorf("SSM %d: relative error %.2e exceeds tolerance\n  CPU: %.15e\n  GPU: %.15e",
						i, relErr, cpuLLH[i], gpuLLH[i])
				}
			}

			t.Logf("Max absolute error: %.2e", maxAbsErr)
			t.Logf("Max relative error: %.2e", maxRelErr)
		})
	}
}

// TestGPUvsCPUMCMC runs MCMC iterations with both CPU and GPU paths
// and verifies they produce consistent likelihood trajectories.
func TestGPUvsCPUMCMC(t *testing.T) {
	if !cuda.Available() {
		t.Skip("CUDA not available")
	}

	// Load test data
	ssms, err := loadSSMData("/home/john/.openclaw/workspace/phylowgs-optimized/ssm_data.txt")
	if err != nil {
		t.Skipf("Test data not available: %v", err)
	}
	cnvs, err := loadCNVData("/home/john/.openclaw/workspace/phylowgs-optimized/cnv_data.txt", ssms)
	if err != nil {
		t.Skipf("CNV data not available: %v", err)
	}

	// Run short MCMC with CPU
	seed := int64(12345)
	rng := rand.New(rand.NewSource(seed))

	// Deep copy SSMs for CPU run
	ssmsCPU := deepCopySSMs(ssms)
	tssbCPU := newTSSB(ssmsCPU, cnvs, 25.0, 1.0, 0.25, rng)

	// Disable GPU for CPU run
	origUseGPU := useGPU
	useGPU = false

	cpuLLHs := make([]float64, 20)
	for i := 0; i < 20; i++ {
		tssbCPU.resampleAssignments(rng)
		tssbCPU.cullTree()
		tssbCPU.setNodeHeights()
		tssbCPU.setNodePaths()
		tssbCPU.mapDatumToNode()
		tssbCPU.metropolis(100, 100.0, rng)
		tssbCPU.resampleSticks(rng)
		cpuLLHs[i] = tssbCPU.completeDataLogLikelihood()
	}

	// Now run GPU version with same seed
	useGPU = true

	// Initialize GPU engine
	nTP := len(ssms[0].A)
	engine, err := cuda.NewGPUEngine(len(ssms), nTP)
	if err != nil {
		t.Fatalf("Failed to create GPU engine: %v", err)
	}
	gpuEngine = engine
	defer func() {
		gpuEngine.Close()
		gpuEngine = nil
		useGPU = origUseGPU
	}()

	// Reset RNG to same seed
	rng = rand.New(rand.NewSource(seed))

	// Deep copy SSMs for GPU run
	ssmsGPU := deepCopySSMs(ssms)
	tssbGPU := newTSSB(ssmsGPU, cnvs, 25.0, 1.0, 0.25, rng)

	// Initialize GPU for TSSB
	if err := tssbGPU.initGPUForTSSB(); err != nil {
		t.Fatalf("Failed to init GPU for TSSB: %v", err)
	}

	gpuLLHs := make([]float64, 20)
	for i := 0; i < 20; i++ {
		tssbGPU.resampleAssignments(rng)
		tssbGPU.cullTree()
		tssbGPU.setNodeHeights()
		tssbGPU.setNodePaths()
		tssbGPU.mapDatumToNode()
		tssbGPU.metropolis(100, 100.0, rng)
		tssbGPU.resampleSticks(rng)
		gpuLLHs[i] = tssbGPU.completeDataLogLikelihood()
	}

	// Compare likelihood trajectories
	maxRelErr := 0.0
	for i := 0; i < 20; i++ {
		relErr := math.Abs(gpuLLHs[i]-cpuLLHs[i]) / math.Abs(cpuLLHs[i])
		if relErr > maxRelErr {
			maxRelErr = relErr
		}
		if relErr > 1e-10 {
			t.Errorf("Iteration %d: LLH mismatch\n  CPU: %.15e\n  GPU: %.15e\n  RelErr: %.2e",
				i, cpuLLHs[i], gpuLLHs[i], relErr)
		}
	}
	t.Logf("Max relative error across iterations: %.2e", maxRelErr)
}

func deepCopySSMs(ssms []*SSM) []*SSM {
	result := make([]*SSM, len(ssms))
	for i, ssm := range ssms {
		result[i] = &SSM{
			ID:              ssm.ID,
			Name:            ssm.Name,
			A:               append([]int{}, ssm.A...),
			D:               append([]int{}, ssm.D...),
			MuR:             ssm.MuR,
			MuV:             ssm.MuV,
			LogBinNormConst: append([]float64{}, ssm.LogBinNormConst...),
		}
	}
	return result
}

// BenchmarkGPUvsCPU benchmarks GPU vs CPU likelihood computation
func BenchmarkCPULikelihood(b *testing.B) {
	rng := rand.New(rand.NewSource(42))
	nSSM := 1000
	nTP := 5

	ssms := make([]*SSM, nSSM)
	phiPerSSM := make([][]float64, nSSM)
	for i := 0; i < nSSM; i++ {
		ssm := &SSM{
			A:               make([]int, nTP),
			D:               make([]int, nTP),
			MuR:             0.999,
			MuV:             0.5,
			LogBinNormConst: make([]float64, nTP),
		}
		phiPerSSM[i] = make([]float64, nTP)
		for tp := 0; tp < nTP; tp++ {
			ssm.D[tp] = rng.Intn(500) + 50
			ssm.A[tp] = rng.Intn(ssm.D[tp])
			ssm.LogBinNormConst[tp] = logBinCoeff(ssm.D[tp], ssm.A[tp])
			phiPerSSM[i][tp] = rng.Float64()
		}
		ssms[i] = ssm
	}

	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		for i, ssm := range ssms {
			_ = ssm.logLikelihoodNoCNV(phiPerSSM[i])
		}
	}
}

func BenchmarkGPULikelihood(b *testing.B) {
	if !cuda.Available() {
		b.Skip("CUDA not available")
	}

	rng := rand.New(rand.NewSource(42))
	nSSM := 1000
	nTP := 5

	ssms := make([]*SSM, nSSM)
	phiBuf := make([]float64, nSSM*nTP)
	for i := 0; i < nSSM; i++ {
		ssm := &SSM{
			A:               make([]int, nTP),
			D:               make([]int, nTP),
			MuR:             0.999,
			MuV:             0.5,
			LogBinNormConst: make([]float64, nTP),
		}
		for tp := 0; tp < nTP; tp++ {
			ssm.D[tp] = rng.Intn(500) + 50
			ssm.A[tp] = rng.Intn(ssm.D[tp])
			ssm.LogBinNormConst[tp] = logBinCoeff(ssm.D[tp], ssm.A[tp])
			phiBuf[i*nTP+tp] = rng.Float64()
		}
		ssms[i] = ssm
	}

	engine, err := cuda.NewGPUEngine(nSSM, nTP)
	if err != nil {
		b.Fatalf("Failed to create GPU engine: %v", err)
	}
	defer engine.Close()

	ssmInterfaces := make([]cuda.SSMDataInterface, len(ssms))
	for i, s := range ssms {
		ssmInterfaces[i] = s
	}

	if err := engine.UploadStaticData(ssmInterfaces); err != nil {
		b.Fatalf("Failed to upload static data: %v", err)
	}

	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		_, _ = engine.ComputeBatchFast(phiBuf, nSSM, nTP)
	}
}
