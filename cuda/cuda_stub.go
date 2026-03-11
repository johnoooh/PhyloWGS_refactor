//go:build !cuda
// +build !cuda

// Package cuda provides GPU-accelerated likelihood computation for PhyloWGS
// This stub file is used when CUDA is not available
package cuda

import "fmt"

// GPULikelihoodEngine manages GPU resources for batched likelihood computation
type GPULikelihoodEngine struct{}

// SSMDataInterface defines the interface for SSM data access
type SSMDataInterface interface {
	GetA() []int
	GetD() []int
	GetMuR() float64
	GetMuV() float64
	GetLogBinNormConst() []float64
}

// Available returns false when CUDA is not compiled in
func Available() bool {
	return false
}

// NewGPUEngine returns an error when CUDA is not available
func NewGPUEngine(nSSM, nTimepoints int) (*GPULikelihoodEngine, error) {
	return nil, fmt.Errorf("CUDA support not compiled in (build with -tags cuda)")
}

// UploadStaticData is a no-op stub
func (e *GPULikelihoodEngine) UploadStaticData(ssms []SSMDataInterface) error {
	return fmt.Errorf("CUDA not available")
}

// ComputeBatchFast is a no-op stub
func (e *GPULikelihoodEngine) ComputeBatchFast(phis []float64, nSSM, nTimepoints int) ([]float64, error) {
	return nil, fmt.Errorf("CUDA not available")
}

// ComputeBatch is a no-op stub
func (e *GPULikelihoodEngine) ComputeBatch(phis [][]float64, ssms []SSMDataInterface) ([]float64, error) {
	return nil, fmt.Errorf("CUDA not available")
}

// Close is a no-op stub
func (e *GPULikelihoodEngine) Close() {}
