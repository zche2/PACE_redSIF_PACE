# =========================================================================================================
# This script conduct daily retrieval for a series of PACE OCI granules
# =========================================================================================================

# series of L1B, L2AOP, and L2BGC file name as input
# set params metadata, save metadata
# do retrieval for each pixel and save the result in a new dataset
using PACE_SIF
using Dates
using NCDatasets, Glob
using Plots    # just for fun XD
using JLD2

# include dependencies - to be updated
include("./set_parameters.jl")
include("./pre_process.jl")
include("./pixel_retrieval.jl")
include("./granule_retrieval.jl")

println("=== Starting script with $(Threads.nthreads()) threads ===")

L1B_dir = "/home/zhe2/data/PACE/L1B_V3"
L2AOP_dir = "/home/zhe2/data/PACE/L2_AOP_V3.1"
L2BGC_dir = "/home/zhe2/data/PACE/L2_BGC_V3.1"

# choose and search files of the day
choose_date  = "20250927";
L1B_file_lst = glob("PACE_OCI.$choose_date*T??????.L1B.V3.nc", L1B_dir);
L2AOP_file_lst = glob("PACE_OCI.$choose_date*T??????.L2.OC_AOP.V3_1.nc", L2AOP_dir);
L2BGC_file_lst = glob("PACE_OCI.$choose_date*T??????.L2.OC_BGC.V3_1.nc", L2BGC_dir);
interim_retrieval_dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_realData/$choose_date"

# ============
# Set params
# ============
# open the first L1B to get the wavelength and solar irradiance for setting up retrieval parameters
ds_L1B = Dataset(L1B_file_lst[1])
red_wavelength       = ds_L1B.group["sensor_band_parameters"]["red_wavelength"][:]
red_solar_irradiance = ds_L1B.group["sensor_band_parameters"]["red_solar_irradiance"][:]
close(ds_L1B)

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
    red_wavelength,
    red_solar_irradiance,
    λ_min, λ_max, scale_factor_SIF, DecompositionMethod, if_log, n, nPC, nSIF, nIter, thr_Converge,
);
# save the parameters for record-keeping
params_name = joinpath(interim_retrieval_dir, "retrieval_params_$(choose_date).jld2")
@save params_name params
println("Retrieval parameters saved to: $params_name")

# ===========================================
# Start retrieval for each granule
# ===========================================
println("=== Starting retrieving ===")

# loop through granules
len = length(L1B_file_lst)
println("Number of granules to process: $len")
for i in 1:len
    start_time = now();
    L1B_file = L1B_file_lst[i];
    L2AOP_file = L2AOP_file_lst[i];
    L2BGC_file = L2BGC_file_lst[i];
    output_str = String(split(basename(L1B_file), ".")[2]);

    granule_retrieval(
        L1B_file, 
        L2AOP_file, 
        L2BGC_file, 
        interim_retrieval_dir, 
        output_str,
        params,
        configuration
    );
    end_time = now();
    elapsed = end_time - start_time
    
    println("  Time elapsed: $output_str: $elapsed")
end

println("=== Retrieval completed for all granules ===")