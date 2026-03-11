# =========================================================================================================
# Run retrieval for one granule: preprocess L1B+L2 → merge → batch retrieval.
# Uses batch_fit Run_orbit_full_nc_parallel for the actual retrieval.
# =========================================================================================================

using TOML
using NCDatasets

const _GLOBAL_FIT_DIR = @__DIR__
const _BATCH_FIT_DIR = joinpath(_GLOBAL_FIT_DIR, "..", "batch_fit")

include(joinpath(_GLOBAL_FIT_DIR, "pre_process.jl"))
include(joinpath(_GLOBAL_FIT_DIR, "merge_inputs.jl"))
include(joinpath(_BATCH_FIT_DIR, "Run_batch_full_nc_parallel.jl"))

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

# Only vars needed for merge (Rtoa_red, red_wavelength, latitude, longitude, watermask)
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
    run_granule_retrieval(
        L1B_path, L2AOP_path, L2BGC_path,
        output_dir, interim_dir, config_path;
        use_threads=true, keep_interim=false, shared_config=false, use_in_memory_merge=true
    )

Preprocess one granule (subset → TOA → merge), then run batch retrieval.
- use_threads: if true, enable pixel-level parallelism in batch_fit
- keep_interim: if false, remove interim files after retrieval
- shared_config: if true, use config_path directly (no per-granule tmp file)
- use_in_memory_merge: if true, read L1B/L2 directly and write single interim (no subset files)
"""
function run_granule_retrieval(
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
        # Read L1B/L2 directly, TOA in memory, write single interim (no subset files)
        preprocess_and_merge_in_memory(L1B_path, L2AOP_path, L2BGC_path, interim_path)
        subset_L1B = subset_L2AOP = subset_L2BGC = nothing
    else
        # 1) Subset L1B with TOA conversion
        subset_L1B = subset_netcdf_dataset(
            L1B_path,
            L1B_VARS,
            interim_dir,
            rename_dims=L1B_RENAME_DIMS,
            post_process_func=TOA_radiance,
        )

        # 2) Subset L2AOP
        subset_L2AOP = subset_netcdf_dataset(
            L2AOP_path,
            L2AOP_VARS,
            interim_dir,
            rename_dims=L2_RENAME_DIMS,
        )

        # 3) Subset L2BGC (optional)
        subset_L2BGC = nothing
        if L2BGC_path !== nothing && isfile(L2BGC_path)
            subset_L2BGC = subset_netcdf_dataset(
                L2BGC_path,
                L2BGC_VARS,
                interim_dir,
                rename_dims=L2_RENAME_DIMS,
            )
        end

        # 4) Merge into single interim file
        merge_global_fit_inputs(
            subset_L1B,
            subset_L2AOP,
            interim_path,
            L2BGC_subset_path=subset_L2BGC,
        )
    end

    # 5) Build cores and run batch retrieval
    # When shared_config: config_path already has output_dir, use_threads, pace_observation overrides
    retrieval_config = config_path
    if !shared_config
        cfg = TOML.parsefile(config_path)
        batch_cfg = merge(get(cfg, "batch_fit", Dict{String, Any}()), Dict("output_dir" => output_dir, "use_threads" => use_threads))
        cfg["batch_fit"] = batch_cfg
        cfg["pace_observation"] = merge(
            get(cfg, "pace_observation", Dict{String, Any}()),
            Dict("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"),
        )
        cfg["batch_fit"]["pixel_filter_vars"] = get(cfg["batch_fit"], "pixel_filter_vars", ["nflh"])
        tmp_config = joinpath(interim_dir, "tmp_config_$(granule_id).toml")
        open(tmp_config, "w") do io
            TOML.print(io, cfg)
        end
        retrieval_config = tmp_config
    end

    n_threads = use_threads ? Threads.nthreads() : 1
    cores = [_build_retrieval_core(retrieval_config) for _ in 1:n_threads]

    run_orbit_full_nc_parallel(cores, interim_path, retrieval_config)

    # 6) Cleanup
    if !shared_config
        rm(retrieval_config, force=true)
    end
    if !keep_interim
        if !use_in_memory_merge
            subset_L1B !== nothing && rm(subset_L1B, force=true)
            subset_L2AOP !== nothing && rm(subset_L2AOP, force=true)
            subset_L2BGC !== nothing && rm(subset_L2BGC, force=true)
        end
        rm(interim_path, force=true)
    end

    return output_dir
end
