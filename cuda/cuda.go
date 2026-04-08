//go:build cuda
// +build cuda

// Package cuda provides GPU-accelerated likelihood computation for PhyloWGS
package cuda

/*
#cgo CFLAGS: -I${SRCDIR}
#cgo LDFLAGS: -L${SRCDIR} -L/usr/local/cuda/lib64 -lcuda_bridge -lcudart -lstdc++ -Wl,-rpath,${SRCDIR} -Wl,-rpath,/usr/local/cuda/lib64
#include "cuda_bridge.h"
#include <stdlib.h>
#include <string.h>
*/
import "C"
import (
	"fmt"
	"runtime"
	"sync"
	"unsafe"
)

// GPULikelihoodEngine manages GPU resources for batched likelihood computation
type GPULikelihoodEngine struct {
	nSSM        int
	nTimepoints int

	// C-allocated buffers (safe to pass to CGo)
	cPhi         *C.double
	cAlt         *C.int
	cDepth       *C.int
	cMuR         *C.double
	cMuV         *C.double
	cLogBinConst *C.double
	cOutLLH      *C.double

	// CNV data buffers
	cHasCNV  *C.int
	cMajorCN *C.int
	cMinorCN *C.int

	staticUploaded bool
	mu             sync.Mutex
}

// SSMDataInterface defines the interface for SSM data access
// Matches the SSM struct in main.go
type SSMDataInterface interface {
	GetA() []int
	GetD() []int
	GetMuR() float64
	GetMuV() float64
	GetLogBinNormConst() []float64
	// CNV data (optional - return false for HasCNV if no CNV)
	HasCNV() bool
	GetMajorCN() int // paternal CN
	GetMinorCN() int // maternal CN
}

// Available returns true if CUDA is available on this system
func Available() bool {
	return C.cuda_available() == 1
}

// NewGPUEngine creates a new GPU likelihood engine
// nSSM: maximum number of SSMs to process in a batch
// nTimepoints: number of timepoints per SSM (must be constant)
func NewGPUEngine(nSSM, nTimepoints int) (*GPULikelihoodEngine, error) {
	if !Available() {
		return nil, fmt.Errorf("CUDA not available")
	}

	ret := C.cuda_init(C.int(nSSM), C.int(nTimepoints))
	if ret != 0 {
		return nil, fmt.Errorf("failed to initialize CUDA")
	}

	sizeSSMxTP := C.size_t(nSSM * nTimepoints)
	sizeSSM := C.size_t(nSSM)

	e := &GPULikelihoodEngine{
		nSSM:        nSSM,
		nTimepoints: nTimepoints,
	}

	// Allocate C memory for buffers
	e.cPhi = (*C.double)(C.malloc(sizeSSMxTP * C.sizeof_double))
	e.cAlt = (*C.int)(C.malloc(sizeSSMxTP * C.sizeof_int))
	e.cDepth = (*C.int)(C.malloc(sizeSSMxTP * C.sizeof_int))
	e.cMuR = (*C.double)(C.malloc(sizeSSM * C.sizeof_double))
	e.cMuV = (*C.double)(C.malloc(sizeSSM * C.sizeof_double))
	e.cLogBinConst = (*C.double)(C.malloc(sizeSSMxTP * C.sizeof_double))
	e.cOutLLH = (*C.double)(C.malloc(sizeSSM * C.sizeof_double))

	// CNV buffers
	e.cHasCNV = (*C.int)(C.malloc(sizeSSM * C.sizeof_int))
	e.cMajorCN = (*C.int)(C.malloc(sizeSSM * C.sizeof_int))
	e.cMinorCN = (*C.int)(C.malloc(sizeSSM * C.sizeof_int))

	if e.cPhi == nil || e.cAlt == nil || e.cDepth == nil ||
		e.cMuR == nil || e.cMuV == nil || e.cLogBinConst == nil || e.cOutLLH == nil ||
		e.cHasCNV == nil || e.cMajorCN == nil || e.cMinorCN == nil {
		e.Close()
		return nil, fmt.Errorf("failed to allocate host buffers")
	}

	// Set a finalizer to clean up if the engine is garbage collected
	runtime.SetFinalizer(e, func(e *GPULikelihoodEngine) {
		e.Close()
	})

	return e, nil
}

// UploadStaticData uploads the static SSM data (A, D, muR, muV, logBinConst) to the GPU.
// Call this once before running MCMC iterations.
// Only phi values change during MCMC, so static data can be uploaded once.
func (e *GPULikelihoodEngine) UploadStaticData(ssms []SSMDataInterface) error {
	e.mu.Lock()
	defer e.mu.Unlock()

	nSSM := len(ssms)
	if nSSM > e.nSSM {
		return fmt.Errorf("too many SSMs: %d > %d", nSSM, e.nSSM)
	}
	if nSSM == 0 {
		return fmt.Errorf("no SSMs provided")
	}

	nTP := len(ssms[0].GetA())
	if nTP > e.nTimepoints {
		return fmt.Errorf("too many timepoints: %d > %d", nTP, e.nTimepoints)
	}

	// Convert C pointers to slices for easier access
	// These are backed by C memory, not Go memory
	altSlice := unsafe.Slice(e.cAlt, e.nSSM*e.nTimepoints)
	depthSlice := unsafe.Slice(e.cDepth, e.nSSM*e.nTimepoints)
	muRSlice := unsafe.Slice(e.cMuR, e.nSSM)
	muVSlice := unsafe.Slice(e.cMuV, e.nSSM)
	logBinConstSlice := unsafe.Slice(e.cLogBinConst, e.nSSM*e.nTimepoints)
	hasCNVSlice := unsafe.Slice(e.cHasCNV, e.nSSM)
	majorCNSlice := unsafe.Slice(e.cMajorCN, e.nSSM)
	minorCNSlice := unsafe.Slice(e.cMinorCN, e.nSSM)

	// Pack static data
	hasCNVData := false
	for i, ssm := range ssms {
		A := ssm.GetA()
		D := ssm.GetD()
		logBinConst := ssm.GetLogBinNormConst()

		muRSlice[i] = C.double(ssm.GetMuR())
		muVSlice[i] = C.double(ssm.GetMuV())

		// CNV data
		if ssm.HasCNV() {
			hasCNVSlice[i] = C.int(1)
			majorCNSlice[i] = C.int(ssm.GetMajorCN())
			minorCNSlice[i] = C.int(ssm.GetMinorCN())
			hasCNVData = true
		} else {
			hasCNVSlice[i] = C.int(0)
			majorCNSlice[i] = C.int(1)
			minorCNSlice[i] = C.int(1)
		}

		base := i * e.nTimepoints
		for tp := 0; tp < nTP; tp++ {
			altSlice[base+tp] = C.int(A[tp])
			depthSlice[base+tp] = C.int(D[tp])
			logBinConstSlice[base+tp] = C.double(logBinConst[tp])
		}
	}

	// Create batch struct - include CNV data if any SSMs have CNVs
	batch := C.LikelihoodBatch{
		n_ssm:         C.int(nSSM),
		n_timepoints:  C.int(nTP),
		phi:           nil, // Not needed for static upload
		alt_counts:    e.cAlt,
		total_depth:   e.cDepth,
		mu_r:          e.cMuR,
		mu_v:          e.cMuV,
		log_bin_const: e.cLogBinConst,
		has_cnv:       nil,
		major_cn:      nil,
		minor_cn:      nil,
	}

	if hasCNVData {
		batch.has_cnv = e.cHasCNV
		batch.major_cn = e.cMajorCN
		batch.minor_cn = e.cMinorCN
	}

	ret := C.cuda_upload_static(&batch)
	if ret != 0 {
		return fmt.Errorf("failed to upload static data to GPU")
	}

	e.staticUploaded = true
	return nil
}

// ComputeBatchFast computes log-likelihoods using only phi values.
// Requires UploadStaticData to have been called first.
// phis is a flat array: phis[i*nTimepoints + tp] = phi for SSM i at timepoint tp
// Returns per-SSM log-likelihoods.
func (e *GPULikelihoodEngine) ComputeBatchFast(phis []float64, nSSM, nTimepoints int) ([]float64, error) {
	e.mu.Lock()
	defer e.mu.Unlock()

	if !e.staticUploaded {
		return nil, fmt.Errorf("static data not uploaded - call UploadStaticData first")
	}

	if nSSM > e.nSSM || nTimepoints > e.nTimepoints {
		return nil, fmt.Errorf("batch size exceeds allocated buffers")
	}

	expected := nSSM * nTimepoints
	if len(phis) < expected {
		return nil, fmt.Errorf("phi array too small: %d < %d", len(phis), expected)
	}

	// Copy phi values to C buffer
	phiSlice := unsafe.Slice(e.cPhi, e.nSSM*e.nTimepoints)
	for i := 0; i < expected; i++ {
		phiSlice[i] = C.double(phis[i])
	}

	ret := C.cuda_compute_phi_only(
		e.cPhi,
		C.int(nSSM),
		C.int(nTimepoints),
		e.cOutLLH,
	)

	if ret != 0 {
		return nil, fmt.Errorf("CUDA computation failed")
	}

	// Copy results to Go slice
	result := make([]float64, nSSM)
	outSlice := unsafe.Slice(e.cOutLLH, e.nSSM)
	for i := 0; i < nSSM; i++ {
		result[i] = float64(outSlice[i])
	}

	return result, nil
}

// ComputeBatch computes log-likelihoods for all SSMs given their phi values.
// This is the full batch computation that uploads all data each time.
// For MCMC use cases, prefer UploadStaticData + ComputeBatchFast.
func (e *GPULikelihoodEngine) ComputeBatch(phis [][]float64, ssms []SSMDataInterface) ([]float64, error) {
	e.mu.Lock()
	defer e.mu.Unlock()

	nSSM := len(ssms)
	if nSSM > e.nSSM {
		return nil, fmt.Errorf("too many SSMs: %d > %d", nSSM, e.nSSM)
	}
	if nSSM == 0 {
		return nil, fmt.Errorf("no SSMs provided")
	}

	nTP := len(ssms[0].GetA())
	if nTP > e.nTimepoints {
		return nil, fmt.Errorf("too many timepoints: %d > %d", nTP, e.nTimepoints)
	}

	if len(phis) != nSSM {
		return nil, fmt.Errorf("phi count mismatch: %d != %d", len(phis), nSSM)
	}

	// Get slices backed by C memory
	phiSlice := unsafe.Slice(e.cPhi, e.nSSM*e.nTimepoints)
	altSlice := unsafe.Slice(e.cAlt, e.nSSM*e.nTimepoints)
	depthSlice := unsafe.Slice(e.cDepth, e.nSSM*e.nTimepoints)
	muRSlice := unsafe.Slice(e.cMuR, e.nSSM)
	muVSlice := unsafe.Slice(e.cMuV, e.nSSM)
	logBinConstSlice := unsafe.Slice(e.cLogBinConst, e.nSSM*e.nTimepoints)

	// Pack data into C buffers
	for i, ssm := range ssms {
		A := ssm.GetA()
		D := ssm.GetD()
		logBinConst := ssm.GetLogBinNormConst()
		phi := phis[i]

		muRSlice[i] = C.double(ssm.GetMuR())
		muVSlice[i] = C.double(ssm.GetMuV())

		base := i * e.nTimepoints
		for tp := 0; tp < nTP; tp++ {
			phiSlice[base+tp] = C.double(phi[tp])
			altSlice[base+tp] = C.int(A[tp])
			depthSlice[base+tp] = C.int(D[tp])
			logBinConstSlice[base+tp] = C.double(logBinConst[tp])
		}
	}

	// Create batch struct
	batch := C.LikelihoodBatch{
		n_ssm:         C.int(nSSM),
		n_timepoints:  C.int(nTP),
		phi:           e.cPhi,
		alt_counts:    e.cAlt,
		total_depth:   e.cDepth,
		mu_r:          e.cMuR,
		mu_v:          e.cMuV,
		log_bin_const: e.cLogBinConst,
	}

	ret := C.cuda_compute_batch(&batch, e.cOutLLH)
	if ret != 0 {
		return nil, fmt.Errorf("CUDA computation failed")
	}

	// Copy results
	result := make([]float64, nSSM)
	outSlice := unsafe.Slice(e.cOutLLH, e.nSSM)
	for i := 0; i < nSSM; i++ {
		result[i] = float64(outSlice[i])
	}

	return result, nil
}

// Close releases GPU resources
func (e *GPULikelihoodEngine) Close() {
	e.mu.Lock()
	defer e.mu.Unlock()

	// Free C-allocated memory
	if e.cPhi != nil {
		C.free(unsafe.Pointer(e.cPhi))
		e.cPhi = nil
	}
	if e.cAlt != nil {
		C.free(unsafe.Pointer(e.cAlt))
		e.cAlt = nil
	}
	if e.cDepth != nil {
		C.free(unsafe.Pointer(e.cDepth))
		e.cDepth = nil
	}
	if e.cMuR != nil {
		C.free(unsafe.Pointer(e.cMuR))
		e.cMuR = nil
	}
	if e.cMuV != nil {
		C.free(unsafe.Pointer(e.cMuV))
		e.cMuV = nil
	}
	if e.cLogBinConst != nil {
		C.free(unsafe.Pointer(e.cLogBinConst))
		e.cLogBinConst = nil
	}
	if e.cOutLLH != nil {
		C.free(unsafe.Pointer(e.cOutLLH))
		e.cOutLLH = nil
	}
	if e.cHasCNV != nil {
		C.free(unsafe.Pointer(e.cHasCNV))
		e.cHasCNV = nil
	}
	if e.cMajorCN != nil {
		C.free(unsafe.Pointer(e.cMajorCN))
		e.cMajorCN = nil
	}
	if e.cMinorCN != nil {
		C.free(unsafe.Pointer(e.cMinorCN))
		e.cMinorCN = nil
	}

	C.cuda_cleanup()
	e.staticUploaded = false

	// Remove finalizer since we've cleaned up
	runtime.SetFinalizer(e, nil)
}
