# =========================================================================================================
# Run retrieval for one granule with SIF addition: preprocess L1B+L2 → merge → add SIF+noise → retrieval.
# Mirrors global_fit granule_retrieval but adds SIF and SNR-consistent noise before retrieval.
# Delegates retrieval to Run_batch_full_nc_adding_sif_parallel (like granule_retrieval → Run_batch_full_nc_parallel).
# =========================================================================================================

using Base.Threads
using TOML
using NCDatasets

const _ADDING_SIF_DIR = @__DIR__
const _GLOBAL_FIT_DIR = joinpath(_ADDING_SIF_DIR, "..", "global_fit")
const _BATCH_FIT_DIR = joinpath(_ADDING_SIF_DIR, "..", "batch_fit")

include(joinpath(_GLOBAL_FIT_DIR, "pre_process.jl"))
include(joinpath(_GLOBAL_FIT_DIR, "merge_inputs.jl"))
include(joinpath(_BATCH_FIT_DIR, "Run_batch_full_nc_adding_sif_parallel.jl"))

const L1B_RENAME_DIMS = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
    "number_of_bands" => "red_wavelength",
    "red_bands" => "red_wavelength",
)
const L2_RENAME_DIMS = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
)
const L1B_VARS = [
    "red_wavelength",
    "red_solar_irradiance",
    "watermask",
    "latitude",
    "longitude",
    "rhot_red",
]
const L2AOP_VARS = ["nflh"]
const L2BGC_VARS = ["chlor_a"]

"""
    run_granule_retrieval_adding_sif(
        L1B_path, L2AOP_path, L2BGC_path,
        output_dir, interim_dir, config_path;
        use_threads=true, keep_interim=false, shared_config=false, use_in_memory_merge=true
    )

Preprocess one granule, then run retrieval with SIF+noise addition per pixel (via parallel batch runner).
Skips pixels where nflh is missing. Output includes spectrally resolved rmse.
"""
function run_granule_retrieval_adding_sif(
    L1B_path::AbstractString,
    L2AOP_path::AbstractString,
    L2BGC_path::Union{AbstractString, Nothing},
    output_dir::AbstractString,
    interim_dir::AbstractString,
    config_path::AbstractString;
    use_threads::Bool=true,
    keep_interim::Bool=false,
    shared_config::Bool=false,
    use_in_memory_merge::Bool=true,
)
    mkpath(interim_dir)
    mkpath(output_dir)

    stem = splitext(basename(L1B_path))[1]
    m = match(r"^PACE_OCI\.(.+)\.L1B\.V3$", stem)
    granule_id = m !== nothing ? m.captures[1] : stem
    interim_path = joinpath(interim_dir, "interim_$(granule_id).nc")

    if use_in_memory_merge
        preprocess_and_merge_in_memory(L1B_path, L2AOP_path, L2BGC_path, interim_path)
    else
        subset_L1B = subset_netcdf_dataset(
            L1B_path, L1B_VARS, interim_dir, rename_dims=L1B_RENAME_DIMS, post_process_func=TOA_radiance)
        subset_L2AOP = subset_netcdf_dataset(
            L2AOP_path, L2AOP_VARS, interim_dir, rename_dims=L2_RENAME_DIMS)
        subset_L2BGC = isfile(L2BGC_path) ? subset_netcdf_dataset(
            L2BGC_path, L2BGC_VARS, interim_dir, rename_dims=L2_RENAME_DIMS) : nothing
        merge_global_fit_inputs(subset_L1B, subset_L2AOP, interim_path, L2BGC_subset_path=subset_L2BGC)
    end

    cfg = TOML.parsefile(config_path)
    add_cfg = get(cfg, "adding_sif", Dict{String, Any}())
    batch_cfg = merge(get(cfg, "batch_fit", Dict{String, Any}()), Dict(
        "output_dir" => output_dir,
        "use_threads" => use_threads,
    ))
    cfg["batch_fit"] = batch_cfg
    cfg["pace_observation"] = merge(
        get(cfg, "pace_observation", Dict{String, Any}()),
        Dict("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"),
    )
    cfg["batch_fit"]["pixel_filter_vars"] = get(cfg["batch_fit"], "pixel_filter_vars", ["nflh"])

    retrieval_config = config_path
    if !shared_config
        tmp_config = joinpath(interim_dir, "tmp_config_$(granule_id).toml")
        open(tmp_config, "w") do io
            TOML.print(io, cfg)
        end
        retrieval_config = tmp_config
    end

    n_threads = use_threads ? Threads.nthreads() : 1
    cores = [_build_retrieval_core(retrieval_config) for _ in 1:n_threads]

    output_path = run_orbit_full_nc_adding_sif_parallel(cores, interim_path, retrieval_config)

    if !shared_config
        rm(retrieval_config, force=true)
    end
    if !keep_interim
        rm(interim_path, force=true)
    end

    return output_path
end
