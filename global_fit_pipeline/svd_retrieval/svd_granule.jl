# One-granule SVD retrieval: preprocess/merge L1B+L2 → interim NetCDF → SVD swath batch.

using NCDatasets

const _SVD_DIR = @__DIR__
const _PIPE_DIR = dirname(_SVD_DIR)
const _REPO_ROOT = dirname(_PIPE_DIR)

include(joinpath(_PIPE_DIR, "julia", "merge_interim.jl"))
include(joinpath(_SVD_DIR, "svd_swath_parallel.jl"))

function run_svd_granule(
    L1B_path::AbstractString,
    L2AOP_path::AbstractString,
    L2BGC_path::Union{AbstractString, Nothing},
    output_dir::AbstractString,
    interim_dir::AbstractString,
    retrieval_cfg::AbstractDict,
    pipeline_config_path::AbstractString;
    use_in_memory_merge::Bool = true,
    pixel_range::Union{Nothing,UnitRange{Int}} = nothing,
    scan_range::Union{Nothing,UnitRange{Int}} = nothing,
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
        error("SVD pipeline currently requires use_in_memory_merge=true (subset path not wired)")
    end

    run_svd_orbit_full_nc_parallel(
        interim_path,
        L1B_path,
        retrieval_cfg,
        pipeline_config_path;
        pixel_range = pixel_range,
        scan_range = scan_range,
    )

    rm(interim_path, force = true)
    return output_dir
end
