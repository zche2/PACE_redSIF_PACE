# =========================================================================================================
# This script 
# =========================================================================================================

# series of L1B, L2AOP, and L2BGC file name as input
# set params metadata, save metadata
# do retrieval for each pixel and save the result in a new dataset

using PACE_SIF
using NCDatasets, Glob
using Plots    # just for fun XD
using JLD2

# include dependencies - to be updated
include("./set_parameters.jl")
include("./pre_process.jl")
include("./pixel_retrieval.jl")
include("./granule_retrieval.jl")

println("=== Starting script with $(Threads.nthreads()) threads ===")



# ===========================================
# preprocess PACE OCI data
# ==========================================
# input dir
L1B_dir = "/home/zhe2/data/PACE/L1B_V3"
L2AOP_dir = "/home/zhe2/data/PACE/L2_AOP"
L2BGC_dir = "/home/zhe2/data/PACE/L2_BGC"
# interim dir 
interim_dir = "/home/zhe2/data/PACE/interim"

# load one scene - to be updated!
L1B_file_pattern = "PACE_OCI.20250130T202059.L1B.V3.nc"
L2AOP_file_pattern = "PACE_OCI.20250130T202059.L2.OC_AOP.V3_0.nc"
L2BGC_file_pattern = "PACE_OCI.20250130T202059.L2.OC_BGC.V3_0.nc"
# L1B_file_pattern = "PACE_OCI.20250928T2224??.L1B.V3.nc"
# L2AOP_file_pattern = "PACE_OCI.20250928T2224??.L2.OC_AOP.V3_1.nc"
# L2BGC_file_pattern = "PACE_OCI.20250928T2224??.L2.OC_BGC.V3_1.nc"

# file paths
L1B_file = first(glob(L1B_file_pattern, L1B_dir))
L2AOP_file = first(glob(L2AOP_file_pattern, L2AOP_dir))
L2BGC_file = first(glob(L2BGC_file_pattern, L2BGC_dir))

# ===========================================
# Pixelwise Retrieval
# ==========================================
# set parameters
nflh_threshold = 0.05;           # threshold for valid pixels based on NFLH
DecompositionMethod = :SVD;    # "NMF" or "SVD"
if_log = true;                 # whether to do log-SVD for transmittance
n     = 10;
ranks = 15;
nPC   = ranks;
nSIF  = 2;
nIter = 25;
thr_Converge = 1e-6;

λ_min = 620.0
λ_max = 860.0

λ_remove_min = 750.0
λ_remove_max = 749.0

λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 
670.76, 671.99, 673.24, 674.51, 675.73, 676.97, 678.21, 679.45, 
754.3, 779.33, 867.11, 869.61, 872.13]

scale_factor_SIF = 20

configuration = "Configurations: \n" *
          "Wavelength range: $λ_min - $λ_max nm\n" *
          "SNR degradation range: $λ_remove_min - $λ_remove_max nm\n" *
          "SIF scale factor: $scale_factor_SIF\n" *
          "Decomposition method: $DecompositionMethod with if_log=$if_log\n" *
          "NFLH threshold: $nflh_threshold\n" *
          "NMF rank: $rank\n" *
          "Order of polynomials to fit: $n, Number of retrieval PCs: $nPC, SIF PCs: $nSIF\n" *
          "Number of iterations: $nIter, Convergence threshold: $thr_Converge\n"


println(configuration)

params = setup_retrieval_parameters(
    ds_subset_L1B["red_wavelength"][:],
    ds_subset_L1B["red_solar_irradiance"][:],
    λ_min, λ_max, scale_factor_SIF, DecompositionMethod, if_log, n, nPC, nSIF, nIter, thr_Converge,
);
# save the parameters for record-keeping
# ===================
# placeholder
# ===================

println("=== Starting retrieving ===")
output_str = "20250130T202059"  # to be updated
granule_retrieval(
    L1B_file, 
    L2AOP_file, 
    L2BGC_file, 
    interim_dir, 
    output_str,
    params, 
    configuration
)