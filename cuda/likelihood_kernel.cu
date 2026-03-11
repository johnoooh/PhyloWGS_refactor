// likelihood_kernel.cu - CUDA kernel for batched SSM log-likelihood computation
// Uses double precision throughout for numerical accuracy in MCMC

#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>

// Block size tuned for modern GPUs - 256 threads for good occupancy
#define BLOCK_SIZE 256
#define WARP_SIZE 32

// Warp-level reduction using shuffle operations (double precision)
__device__ double warpReduceSum(double val) {
    for (int offset = WARP_SIZE / 2; offset > 0; offset /= 2) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

// Block-level reduction using shared memory
__device__ double blockReduceSum(double val) {
    __shared__ double shared[32]; // One slot per warp
    int lane = threadIdx.x % WARP_SIZE;
    int wid = threadIdx.x / WARP_SIZE;

    // Warp-level reduction
    val = warpReduceSum(val);

    // Write reduced value from each warp
    if (lane == 0) {
        shared[wid] = val;
    }
    __syncthreads();

    // Read from shared memory only if that warp existed
    val = (threadIdx.x < blockDim.x / WARP_SIZE) ? shared[lane] : 0.0;

    // Final reduction in first warp
    if (wid == 0) {
        val = warpReduceSum(val);
    }

    return val;
}

// Compute mu with CNV correction
// When SSM has an overlapping CNV, the expected VAF changes based on copy number
// We average over two possibilities: variant on major or minor allele
__device__ double computeMuWithCNV(double p, double muR, int majorCN, int minorCN) {
    int totalCN = majorCN + minorCN;
    if (totalCN <= 0) {
        // Complete deletion - fallback to standard formula
        return (1.0 - p) * muR + p * 0.5;
    }
    
    // Case 1: Variant on major allele
    // nr1 = normal cells (2 copies) + tumor cells with minor allele copies
    // nv1 = tumor cells with major allele copies (carrying variant)
    double nr1 = (1.0 - p) * 2.0 + p * (double)minorCN;
    double nv1 = p * (double)majorCN;
    
    // Case 2: Variant on minor allele  
    double nr2 = (1.0 - p) * 2.0 + p * (double)majorCN;
    double nv2 = p * (double)minorCN;
    
    // Compute mu for each case: mu = (nr*muR + nv*(1-muR)) / (nr + nv)
    double mu = 0.0;
    int n_valid = 0;
    
    if (nv1 > 0.0) {
        double denom1 = nr1 + nv1;
        if (denom1 > 1e-15) {
            double mu1 = (nr1 * muR + nv1 * (1.0 - muR)) / denom1;
            mu += mu1;
            n_valid++;
        }
    }
    
    if (nv2 > 0.0) {
        double denom2 = nr2 + nv2;
        if (denom2 > 1e-15) {
            double mu2 = (nr2 * muR + nv2 * (1.0 - muR)) / denom2;
            mu += mu2;
            n_valid++;
        }
    }
    
    if (n_valid > 0) {
        mu /= (double)n_valid;  // Average over valid cases
    } else {
        // Fallback to standard formula
        mu = (1.0 - p) * muR + p * 0.5;
    }
    
    return mu;
}

// Main likelihood kernel (standard, no CNV)
// Each block handles one SSM, threads within block handle timepoints
// This allows efficient reduction per-SSM
__global__ void computeBatchLikelihood(
    const double* __restrict__ phi,        // [nSSM * nTimepoints] - phi values
    const int* __restrict__ A,             // [nSSM * nTimepoints] - alt counts  
    const int* __restrict__ D,             // [nSSM * nTimepoints] - total depth
    const double* __restrict__ muR,        // [nSSM] - reference mu per SSM
    const double* __restrict__ muV,        // [nSSM] - variant mu per SSM
    const double* __restrict__ logBinConst,// [nSSM * nTimepoints] - log binomial coeffs
    double* __restrict__ outLLH,           // [nSSM] - output per-SSM log-likelihood
    int nSSM,
    int nTimepoints
) {
    int ssm = blockIdx.x;
    if (ssm >= nSSM) return;
    
    // Per-SSM constants
    double muR_ssm = muR[ssm];
    double muV_ssm = muV[ssm];
    
    // Each thread accumulates partial LLH for subset of timepoints
    double partialLLH = 0.0;
    
    // Stride over timepoints
    int base = ssm * nTimepoints;
    for (int tp = threadIdx.x; tp < nTimepoints; tp += blockDim.x) {
        int idx = base + tp;
        
        // Get phi value (already clamped to [0,1] in Go)
        double p = phi[idx];
        
        // Compute mu = (1-p)*muR + p*muV
        double mu = (1.0 - p) * muR_ssm + p * muV_ssm;
        
        // Clamp mu to avoid log(0) and log(1)
        if (mu < 1e-15) mu = 1e-15;
        if (mu > 1.0 - 1e-15) mu = 1.0 - 1e-15;
        
        // Get read counts
        int a = A[idx];
        int d = D[idx];
        
        // Log binomial likelihood: a*log(mu) + (d-a)*log(1-mu)
        // Using double precision log() for accuracy
        double logLik = (double)a * log(mu) + (double)(d - a) * log(1.0 - mu);
        
        // Add precomputed log binomial coefficient
        logLik += logBinConst[idx];
        
        partialLLH += logLik;
    }
    
    // Block-level reduction to sum across all threads
    double totalLLH = blockReduceSum(partialLLH);
    
    // Thread 0 writes result
    if (threadIdx.x == 0) {
        outLLH[ssm] = totalLLH;
    }
}

// CNV-aware likelihood kernel
// Uses copy number information to compute corrected mu
__global__ void computeBatchLikelihoodCNV(
    const double* __restrict__ phi,        // [nSSM * nTimepoints] - phi values
    const int* __restrict__ A,             // [nSSM * nTimepoints] - alt counts  
    const int* __restrict__ D,             // [nSSM * nTimepoints] - total depth
    const double* __restrict__ muR,        // [nSSM] - reference mu per SSM
    const double* __restrict__ muV,        // [nSSM] - variant mu per SSM (unused when CNV)
    const double* __restrict__ logBinConst,// [nSSM * nTimepoints] - log binomial coeffs
    const int* __restrict__ hasCNV,        // [nSSM] - 1 if has CNV, 0 otherwise
    const int* __restrict__ majorCN,       // [nSSM] - major copy number
    const int* __restrict__ minorCN,       // [nSSM] - minor copy number
    double* __restrict__ outLLH,           // [nSSM] - output per-SSM log-likelihood
    int nSSM,
    int nTimepoints
) {
    int ssm = blockIdx.x;
    if (ssm >= nSSM) return;
    
    // Per-SSM constants
    double muR_ssm = muR[ssm];
    double muV_ssm = muV[ssm];
    int has_cnv = hasCNV[ssm];
    int major_cn = majorCN[ssm];
    int minor_cn = minorCN[ssm];
    
    // Each thread accumulates partial LLH for subset of timepoints
    double partialLLH = 0.0;
    
    // Stride over timepoints
    int base = ssm * nTimepoints;
    for (int tp = threadIdx.x; tp < nTimepoints; tp += blockDim.x) {
        int idx = base + tp;
        
        // Get phi value (already clamped to [0,1] in Go)
        double p = phi[idx];
        
        // Compute mu - either standard or CNV-corrected
        double mu;
        if (has_cnv) {
            mu = computeMuWithCNV(p, muR_ssm, major_cn, minor_cn);
        } else {
            mu = (1.0 - p) * muR_ssm + p * muV_ssm;
        }
        
        // Clamp mu to avoid log(0) and log(1)
        if (mu < 1e-15) mu = 1e-15;
        if (mu > 1.0 - 1e-15) mu = 1.0 - 1e-15;
        
        // Get read counts
        int a = A[idx];
        int d = D[idx];
        
        // Log binomial likelihood: a*log(mu) + (d-a)*log(1-mu)
        double logLik = (double)a * log(mu) + (double)(d - a) * log(1.0 - mu);
        
        // Add precomputed log binomial coefficient
        logLik += logBinConst[idx];
        
        partialLLH += logLik;
    }
    
    // Block-level reduction to sum across all threads
    double totalLLH = blockReduceSum(partialLLH);
    
    // Thread 0 writes result
    if (threadIdx.x == 0) {
        outLLH[ssm] = totalLLH;
    }
}

// Alternative kernel for small timepoint counts (< 32)
// Uses one warp per SSM for efficiency when nTimepoints is small
__global__ void computeBatchLikelihoodSmall(
    const double* __restrict__ phi,
    const int* __restrict__ A,
    const int* __restrict__ D,
    const double* __restrict__ muR,
    const double* __restrict__ muV,
    const double* __restrict__ logBinConst,
    double* __restrict__ outLLH,
    int nSSM,
    int nTimepoints
) {
    // Multiple SSMs per block, one warp per SSM
    int warpId = threadIdx.x / WARP_SIZE;
    int lane = threadIdx.x % WARP_SIZE;
    int warpsPerBlock = blockDim.x / WARP_SIZE;
    
    int ssm = blockIdx.x * warpsPerBlock + warpId;
    if (ssm >= nSSM) return;
    
    double muR_ssm = muR[ssm];
    double muV_ssm = muV[ssm];
    
    double partialLLH = 0.0;
    int base = ssm * nTimepoints;
    
    // Each lane handles subset of timepoints
    for (int tp = lane; tp < nTimepoints; tp += WARP_SIZE) {
        int idx = base + tp;
        double p = phi[idx];
        double mu = (1.0 - p) * muR_ssm + p * muV_ssm;
        
        if (mu < 1e-15) mu = 1e-15;
        if (mu > 1.0 - 1e-15) mu = 1.0 - 1e-15;
        
        int a = A[idx];
        int d = D[idx];
        
        double logLik = (double)a * log(mu) + (double)(d - a) * log(1.0 - mu);
        logLik += logBinConst[idx];
        
        partialLLH += logLik;
    }
    
    // Warp reduction
    double totalLLH = warpReduceSum(partialLLH);
    
    if (lane == 0) {
        outLLH[ssm] = totalLLH;
    }
}

// CNV-aware version for small timepoint counts
__global__ void computeBatchLikelihoodSmallCNV(
    const double* __restrict__ phi,
    const int* __restrict__ A,
    const int* __restrict__ D,
    const double* __restrict__ muR,
    const double* __restrict__ muV,
    const double* __restrict__ logBinConst,
    const int* __restrict__ hasCNV,
    const int* __restrict__ majorCN,
    const int* __restrict__ minorCN,
    double* __restrict__ outLLH,
    int nSSM,
    int nTimepoints
) {
    // Multiple SSMs per block, one warp per SSM
    int warpId = threadIdx.x / WARP_SIZE;
    int lane = threadIdx.x % WARP_SIZE;
    int warpsPerBlock = blockDim.x / WARP_SIZE;
    
    int ssm = blockIdx.x * warpsPerBlock + warpId;
    if (ssm >= nSSM) return;
    
    double muR_ssm = muR[ssm];
    double muV_ssm = muV[ssm];
    int has_cnv = hasCNV[ssm];
    int major_cn = majorCN[ssm];
    int minor_cn = minorCN[ssm];
    
    double partialLLH = 0.0;
    int base = ssm * nTimepoints;
    
    // Each lane handles subset of timepoints
    for (int tp = lane; tp < nTimepoints; tp += WARP_SIZE) {
        int idx = base + tp;
        double p = phi[idx];
        
        double mu;
        if (has_cnv) {
            mu = computeMuWithCNV(p, muR_ssm, major_cn, minor_cn);
        } else {
            mu = (1.0 - p) * muR_ssm + p * muV_ssm;
        }
        
        if (mu < 1e-15) mu = 1e-15;
        if (mu > 1.0 - 1e-15) mu = 1.0 - 1e-15;
        
        int a = A[idx];
        int d = D[idx];
        
        double logLik = (double)a * log(mu) + (double)(d - a) * log(1.0 - mu);
        logLik += logBinConst[idx];
        
        partialLLH += logLik;
    }
    
    // Warp reduction
    double totalLLH = warpReduceSum(partialLLH);
    
    if (lane == 0) {
        outLLH[ssm] = totalLLH;
    }
}

// Host-side structures for persistent GPU memory
struct GPUBuffers {
    double* d_phi;
    int* d_A;
    int* d_D;
    double* d_muR;
    double* d_muV;
    double* d_logBinConst;
    double* d_outLLH;
    
    // CNV data buffers
    int* d_hasCNV;
    int* d_majorCN;
    int* d_minorCN;
    int cnvDataUploaded;
    
    int nSSM;
    int nTimepoints;
    int allocated;
    cudaStream_t stream;
};

static GPUBuffers g_buffers = {0};

// Initialize CUDA and allocate persistent buffers
extern "C" int cuda_init(int n_ssm, int n_timepoints) {
    if (g_buffers.allocated) {
        // Already allocated - check if sizes match
        if (g_buffers.nSSM >= n_ssm && g_buffers.nTimepoints >= n_timepoints) {
            return 0; // Can reuse existing buffers
        }
        // Need to reallocate - cleanup first
        if (g_buffers.d_phi) cudaFree(g_buffers.d_phi);
        if (g_buffers.d_A) cudaFree(g_buffers.d_A);
        if (g_buffers.d_D) cudaFree(g_buffers.d_D);
        if (g_buffers.d_muR) cudaFree(g_buffers.d_muR);
        if (g_buffers.d_muV) cudaFree(g_buffers.d_muV);
        if (g_buffers.d_logBinConst) cudaFree(g_buffers.d_logBinConst);
        if (g_buffers.d_outLLH) cudaFree(g_buffers.d_outLLH);
        if (g_buffers.d_hasCNV) cudaFree(g_buffers.d_hasCNV);
        if (g_buffers.d_majorCN) cudaFree(g_buffers.d_majorCN);
        if (g_buffers.d_minorCN) cudaFree(g_buffers.d_minorCN);
        if (g_buffers.stream) cudaStreamDestroy(g_buffers.stream);
        g_buffers.allocated = 0;
    }
    
    // Check CUDA device
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess || deviceCount == 0) {
        fprintf(stderr, "CUDA: No devices found\n");
        return -1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    fprintf(stderr, "CUDA: Using %s (compute %d.%d, %.1f GB)\n", 
            prop.name, prop.major, prop.minor, 
            prop.totalGlobalMem / 1e9);
    
    // Allocate device memory
    size_t sizeSSMxTP = (size_t)n_ssm * n_timepoints;
    
    err = cudaMalloc(&g_buffers.d_phi, sizeSSMxTP * sizeof(double));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_A, sizeSSMxTP * sizeof(int));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_D, sizeSSMxTP * sizeof(int));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_muR, n_ssm * sizeof(double));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_muV, n_ssm * sizeof(double));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_logBinConst, sizeSSMxTP * sizeof(double));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_outLLH, n_ssm * sizeof(double));
    if (err != cudaSuccess) goto error;
    
    // CNV data buffers
    err = cudaMalloc(&g_buffers.d_hasCNV, n_ssm * sizeof(int));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_majorCN, n_ssm * sizeof(int));
    if (err != cudaSuccess) goto error;
    
    err = cudaMalloc(&g_buffers.d_minorCN, n_ssm * sizeof(int));
    if (err != cudaSuccess) goto error;
    
    err = cudaStreamCreate(&g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    g_buffers.cnvDataUploaded = 0;
    
    g_buffers.nSSM = n_ssm;
    g_buffers.nTimepoints = n_timepoints;
    g_buffers.allocated = 1;
    
    fprintf(stderr, "CUDA: Allocated buffers for %d SSMs x %d timepoints (%.1f MB)\n",
            n_ssm, n_timepoints,
            (double)(sizeSSMxTP * (sizeof(double)*3 + sizeof(int)*2) + n_ssm * sizeof(double)*3) / 1e6);
    
    return 0;

error:
    fprintf(stderr, "CUDA: Memory allocation failed: %s\n", cudaGetErrorString(err));
    // Cleanup partial allocations
    if (g_buffers.d_phi) cudaFree(g_buffers.d_phi);
    if (g_buffers.d_A) cudaFree(g_buffers.d_A);
    if (g_buffers.d_D) cudaFree(g_buffers.d_D);
    if (g_buffers.d_muR) cudaFree(g_buffers.d_muR);
    if (g_buffers.d_muV) cudaFree(g_buffers.d_muV);
    if (g_buffers.d_logBinConst) cudaFree(g_buffers.d_logBinConst);
    if (g_buffers.d_outLLH) cudaFree(g_buffers.d_outLLH);
    if (g_buffers.d_hasCNV) cudaFree(g_buffers.d_hasCNV);
    if (g_buffers.d_majorCN) cudaFree(g_buffers.d_majorCN);
    if (g_buffers.d_minorCN) cudaFree(g_buffers.d_minorCN);
    memset(&g_buffers, 0, sizeof(g_buffers));
    return -1;
}

// C-compatible batch structure (defined in header)
typedef struct {
    int n_ssm;
    int n_timepoints;
    double* phi;
    int* alt_counts;
    int* total_depth;
    double* mu_r;
    double* mu_v;
    double* log_bin_const;
} LikelihoodBatch;

// Compute batch likelihoods
extern "C" int cuda_compute_batch(const LikelihoodBatch* batch, double* out_llh) {
    if (!g_buffers.allocated) {
        fprintf(stderr, "CUDA: Not initialized\n");
        return -1;
    }
    
    if (batch->n_ssm > g_buffers.nSSM || batch->n_timepoints > g_buffers.nTimepoints) {
        fprintf(stderr, "CUDA: Batch size exceeds allocated buffers\n");
        return -1;
    }
    
    int nSSM = batch->n_ssm;
    int nTP = batch->n_timepoints;
    size_t sizeSSMxTP = (size_t)nSSM * nTP;
    
    // Copy data to device
    cudaError_t err;
    err = cudaMemcpyAsync(g_buffers.d_phi, batch->phi, 
                          sizeSSMxTP * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    err = cudaMemcpyAsync(g_buffers.d_A, batch->alt_counts,
                          sizeSSMxTP * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    err = cudaMemcpyAsync(g_buffers.d_D, batch->total_depth,
                          sizeSSMxTP * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    err = cudaMemcpyAsync(g_buffers.d_muR, batch->mu_r,
                          nSSM * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    err = cudaMemcpyAsync(g_buffers.d_muV, batch->mu_v,
                          nSSM * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    err = cudaMemcpyAsync(g_buffers.d_logBinConst, batch->log_bin_const,
                          sizeSSMxTP * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    // Choose kernel based on timepoint count
    if (nTP <= 32) {
        // Use warp-based kernel for small timepoint counts
        // Multiple warps per block, one warp per SSM
        int warpsPerBlock = 8;
        int threadsPerBlock = warpsPerBlock * WARP_SIZE;
        int numBlocks = (nSSM + warpsPerBlock - 1) / warpsPerBlock;
        
        computeBatchLikelihoodSmall<<<numBlocks, threadsPerBlock, 0, g_buffers.stream>>>(
            g_buffers.d_phi,
            g_buffers.d_A,
            g_buffers.d_D,
            g_buffers.d_muR,
            g_buffers.d_muV,
            g_buffers.d_logBinConst,
            g_buffers.d_outLLH,
            nSSM,
            nTP
        );
    } else {
        // Use block-based kernel for larger timepoint counts
        // One block per SSM
        computeBatchLikelihood<<<nSSM, BLOCK_SIZE, 0, g_buffers.stream>>>(
            g_buffers.d_phi,
            g_buffers.d_A,
            g_buffers.d_D,
            g_buffers.d_muR,
            g_buffers.d_muV,
            g_buffers.d_logBinConst,
            g_buffers.d_outLLH,
            nSSM,
            nTP
        );
    }
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA: Kernel launch failed: %s\n", cudaGetErrorString(err));
        return -1;
    }
    
    // Copy results back
    err = cudaMemcpyAsync(out_llh, g_buffers.d_outLLH,
                          nSSM * sizeof(double), cudaMemcpyDeviceToHost, g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    // Synchronize
    err = cudaStreamSynchronize(g_buffers.stream);
    if (err != cudaSuccess) goto error;
    
    return 0;

error:
    fprintf(stderr, "CUDA: Operation failed: %s\n", cudaGetErrorString(err));
    return -1;
}

// Upload static data (A, D, muR, muV, logBinConst, CNV data) once
// Only phi changes between calls during MCMC
extern "C" int cuda_upload_static(const LikelihoodBatch* batch) {
    if (!g_buffers.allocated) {
        return -1;
    }
    
    int nSSM = batch->n_ssm;
    int nTP = batch->n_timepoints;
    size_t sizeSSMxTP = (size_t)nSSM * nTP;
    
    cudaError_t err;
    err = cudaMemcpyAsync(g_buffers.d_A, batch->alt_counts,
                          sizeSSMxTP * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    err = cudaMemcpyAsync(g_buffers.d_D, batch->total_depth,
                          sizeSSMxTP * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    err = cudaMemcpyAsync(g_buffers.d_muR, batch->mu_r,
                          nSSM * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    err = cudaMemcpyAsync(g_buffers.d_muV, batch->mu_v,
                          nSSM * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    err = cudaMemcpyAsync(g_buffers.d_logBinConst, batch->log_bin_const,
                          sizeSSMxTP * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    // Upload CNV data if provided
    g_buffers.cnvDataUploaded = 0;
    if (batch->has_cnv != NULL && batch->major_cn != NULL && batch->minor_cn != NULL) {
        err = cudaMemcpyAsync(g_buffers.d_hasCNV, batch->has_cnv,
                              nSSM * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
        if (err != cudaSuccess) return -1;
        
        err = cudaMemcpyAsync(g_buffers.d_majorCN, batch->major_cn,
                              nSSM * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
        if (err != cudaSuccess) return -1;
        
        err = cudaMemcpyAsync(g_buffers.d_minorCN, batch->minor_cn,
                              nSSM * sizeof(int), cudaMemcpyHostToDevice, g_buffers.stream);
        if (err != cudaSuccess) return -1;
        
        g_buffers.cnvDataUploaded = 1;
    }
    
    err = cudaStreamSynchronize(g_buffers.stream);
    return (err == cudaSuccess) ? 0 : -1;
}

// Fast path: only upload phi and compute
// Uses CNV-aware kernels if CNV data was uploaded via cuda_upload_static
extern "C" int cuda_compute_phi_only(const double* phi, int n_ssm, int n_timepoints, double* out_llh) {
    if (!g_buffers.allocated) {
        return -1;
    }
    
    size_t sizeSSMxTP = (size_t)n_ssm * n_timepoints;
    
    cudaError_t err;
    err = cudaMemcpyAsync(g_buffers.d_phi, phi,
                          sizeSSMxTP * sizeof(double), cudaMemcpyHostToDevice, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    // Launch kernel - use CNV-aware version if CNV data is available
    if (n_timepoints <= 32) {
        int warpsPerBlock = 8;
        int threadsPerBlock = warpsPerBlock * WARP_SIZE;
        int numBlocks = (n_ssm + warpsPerBlock - 1) / warpsPerBlock;
        
        if (g_buffers.cnvDataUploaded) {
            computeBatchLikelihoodSmallCNV<<<numBlocks, threadsPerBlock, 0, g_buffers.stream>>>(
                g_buffers.d_phi,
                g_buffers.d_A,
                g_buffers.d_D,
                g_buffers.d_muR,
                g_buffers.d_muV,
                g_buffers.d_logBinConst,
                g_buffers.d_hasCNV,
                g_buffers.d_majorCN,
                g_buffers.d_minorCN,
                g_buffers.d_outLLH,
                n_ssm,
                n_timepoints
            );
        } else {
            computeBatchLikelihoodSmall<<<numBlocks, threadsPerBlock, 0, g_buffers.stream>>>(
                g_buffers.d_phi,
                g_buffers.d_A,
                g_buffers.d_D,
                g_buffers.d_muR,
                g_buffers.d_muV,
                g_buffers.d_logBinConst,
                g_buffers.d_outLLH,
                n_ssm,
                n_timepoints
            );
        }
    } else {
        // Large timepoint count - use block-based kernel
        if (g_buffers.cnvDataUploaded) {
            computeBatchLikelihoodCNV<<<n_ssm, BLOCK_SIZE, 0, g_buffers.stream>>>(
                g_buffers.d_phi,
                g_buffers.d_A,
                g_buffers.d_D,
                g_buffers.d_muR,
                g_buffers.d_muV,
                g_buffers.d_logBinConst,
                g_buffers.d_hasCNV,
                g_buffers.d_majorCN,
                g_buffers.d_minorCN,
                g_buffers.d_outLLH,
                n_ssm,
                n_timepoints
            );
        } else {
            computeBatchLikelihood<<<n_ssm, BLOCK_SIZE, 0, g_buffers.stream>>>(
                g_buffers.d_phi,
                g_buffers.d_A,
                g_buffers.d_D,
                g_buffers.d_muR,
                g_buffers.d_muV,
                g_buffers.d_logBinConst,
                g_buffers.d_outLLH,
                n_ssm,
                n_timepoints
            );
        }
    }
    
    err = cudaGetLastError();
    if (err != cudaSuccess) return -1;
    
    err = cudaMemcpyAsync(out_llh, g_buffers.d_outLLH,
                          n_ssm * sizeof(double), cudaMemcpyDeviceToHost, g_buffers.stream);
    if (err != cudaSuccess) return -1;
    
    err = cudaStreamSynchronize(g_buffers.stream);
    return (err == cudaSuccess) ? 0 : -1;
}

// Cleanup
extern "C" void cuda_cleanup(void) {
    if (!g_buffers.allocated) return;
    
    if (g_buffers.d_phi) cudaFree(g_buffers.d_phi);
    if (g_buffers.d_A) cudaFree(g_buffers.d_A);
    if (g_buffers.d_D) cudaFree(g_buffers.d_D);
    if (g_buffers.d_muR) cudaFree(g_buffers.d_muR);
    if (g_buffers.d_muV) cudaFree(g_buffers.d_muV);
    if (g_buffers.d_logBinConst) cudaFree(g_buffers.d_logBinConst);
    if (g_buffers.d_outLLH) cudaFree(g_buffers.d_outLLH);
    if (g_buffers.d_hasCNV) cudaFree(g_buffers.d_hasCNV);
    if (g_buffers.d_majorCN) cudaFree(g_buffers.d_majorCN);
    if (g_buffers.d_minorCN) cudaFree(g_buffers.d_minorCN);
    if (g_buffers.stream) cudaStreamDestroy(g_buffers.stream);
    
    memset(&g_buffers, 0, sizeof(g_buffers));
    
    fprintf(stderr, "CUDA: Cleanup complete\n");
}

// Check if CUDA is available
extern "C" int cuda_available(void) {
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    return (err == cudaSuccess && deviceCount > 0) ? 1 : 0;
}
