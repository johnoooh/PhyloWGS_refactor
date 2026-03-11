// cuda_bridge.h - C interface for CUDA likelihood computation
// CGo-compatible header (no C++ constructs)

#ifndef CUDA_BRIDGE_H
#define CUDA_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

// Batch data structure for likelihood computation
typedef struct {
    int n_ssm;            // Number of SSMs
    int n_timepoints;     // Number of timepoints per SSM
    double* phi;          // [n_ssm * n_timepoints] - cellular prevalence values
    int*    alt_counts;   // [n_ssm * n_timepoints] - alternate read counts
    int*    total_depth;  // [n_ssm * n_timepoints] - total read depth
    double* mu_r;         // [n_ssm] - reference allele probability per SSM
    double* mu_v;         // [n_ssm] - variant allele probability per SSM
    double* log_bin_const;// [n_ssm * n_timepoints] - precomputed log binomial coefficients
    
    // CNV data (optional - set has_cnv[i]=1 if SSM i has CNV)
    int*    has_cnv;      // [n_ssm] - 1 if SSM has CNV, 0 otherwise
    int*    major_cn;     // [n_ssm] - major copy number (paternal)
    int*    minor_cn;     // [n_ssm] - minor copy number (maternal)
} LikelihoodBatch;

// Initialize CUDA device and allocate persistent GPU buffers.
// Call once at startup with maximum expected batch size.
// Returns 0 on success, -1 on failure.
int cuda_init(int n_ssm, int n_timepoints);

// Compute batch log-likelihoods.
// All fields in batch must be filled.
// out_llh must be pre-allocated with n_ssm doubles.
// Returns 0 on success, -1 on failure.
int cuda_compute_batch(const LikelihoodBatch* batch, double* out_llh);

// Upload static data (A, D, muR, muV, logBinConst) to GPU.
// Call once per MCMC run before iterating.
// Only phi values change between MCMC iterations.
// Returns 0 on success, -1 on failure.
int cuda_upload_static(const LikelihoodBatch* batch);

// Fast path: only upload phi values and compute.
// Requires cuda_upload_static() to have been called first.
// phi must be [n_ssm * n_timepoints] doubles.
// out_llh must be [n_ssm] doubles.
// Returns 0 on success, -1 on failure.
int cuda_compute_phi_only(const double* phi, int n_ssm, int n_timepoints, double* out_llh);

// Free GPU buffers.
// Safe to call multiple times.
void cuda_cleanup(void);

// Check if CUDA is available on this system.
// Returns 1 if available, 0 otherwise.
int cuda_available(void);

#ifdef __cplusplus
}
#endif

#endif // CUDA_BRIDGE_H
