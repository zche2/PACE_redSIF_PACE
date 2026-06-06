# Quick smoke tests for rasterize_svd_retrievals.jl
#
# Usage:
#   julia --project=. global_fit_pipeline/rasterize/run_rasterize_tests.jl
#
# Optional env:
#   RASTER_TEST_RET=path/to/interim_*_svd_retrieval_*.nc
#   RASTER_TEST_L2AOP_DIR=path/to/L2_AOP

using Dates
using NCDatasets
using Test

include(joinpath(@__DIR__, "rasterize_svd_retrievals.jl"))

const _DEFAULT_RET = "/home/zhe2/data/PACE/svd_retrieval_output/interim_20250703T211157_svd_retrieval_full_gpu.nc"
const _DEFAULT_L2DIR = "/home/zhe2/data/PACE/L2_AOP_V3.1"

function _test_retrieval_path()
    p = get(ENV, "RASTER_TEST_RET", _DEFAULT_RET)
    isfile(p) || error("Set RASTER_TEST_RET to an existing retrieval granule NC; missing: $p")
    return p
end

function _test_cfg(; exclude_nflh::Bool=true)
    return RasterConfig(
        dirname(_test_retrieval_path()),
        mktempdir(),
        _DEFAULT_L2DIR,
        Date(2025, 7, 3),
        Date(2025, 7, 3),
        1,
        0,
        36,  # coarse grid for speed
        Set{Int16}([Int16(1)]),
        true,
        false,
        false,
        exclude_nflh,
        "nflh",
        Inf,
        Inf,
    )
end

@testset "granule id and L2 path" begin
    ret = _test_retrieval_path()
    gid = _granule_id_from_retrieval_path(ret)
    @test gid == "20250703T211157"
    @test gid isa String
    l2 = _resolve_l2aop_path(_DEFAULT_L2DIR, gid)
    @test l2 !== nothing
    @test isfile(l2)
end

@testset "L2 path coarse time match (last 2 digits)" begin
    @test _granule_id_coarse_prefix("20250703T211157") == "20250703T2111"
    l2 = _resolve_l2aop_path(_DEFAULT_L2DIR, "20250703T211157")
    @test l2 !== nothing
    @test occursin("20250703T2111", basename(l2))
    # exact file may be …211148 while retrieval id is …211157
    @test occursin("PACE_OCI.20250703T211", basename(l2))
end

@testset "read swath fields (pixels, scans)" begin
    ret = _test_retrieval_path()
    NCDatasets.Dataset(ret, "r") do ds
        lat = _read_pixels_scans_2d(ds, "latitude")
        sif = _read_pixels_scans_2d(ds, "sif_radiance_678nm")
        @test lat !== nothing
        @test sif !== nothing
        @test size(lat) == size(sif)
        @test size(lat, 1) == length(ds["source_pixel_index"])
        @test size(lat, 2) == length(ds["source_scan_index"])
    end
end

@testset "nflh mask from L2 AOP" begin
    ret = _test_retrieval_path()
    gid = _granule_id_from_retrieval_path(ret)
    l2 = _resolve_l2aop_path(_DEFAULT_L2DIR, gid)
    NCDatasets.Dataset(ret, "r") do ds
        mask = _nflh_present_mask(ds, l2, "nflh")
        @test mask !== nothing
        lat = _read_pixels_scans_2d(ds, "latitude")
        @test size(mask) == size(lat)
        @test count(mask) > 0
        println("  nflh present: ", count(mask), " / ", length(mask))
    end
end

@testset "accumulate one granule (with nflh filter)" begin
    ret = _test_retrieval_path()
    cfg = _test_cfg(exclude_nflh=true)
    g = Grid(cfg.resolution)
    accumulate_file!(g, ret, cfg)
    n_obs = sum(g.counts)
    @test n_obs > 0
    println("  binned observations: ", n_obs)
end

@testset "accumulate one granule (no nflh filter)" begin
    ret = _test_retrieval_path()
    cfg = _test_cfg(exclude_nflh=false)
    g = Grid(cfg.resolution)
    accumulate_file!(g, ret, cfg)
    n_obs = sum(g.counts)
    @test n_obs > 0
    println("  binned observations (no nflh): ", n_obs)
end

@testset "write mini raster netcdf" begin
    ret = _test_retrieval_path()
    cfg = _test_cfg(exclude_nflh=true)
    g, _ = accumulate_grid([ret], cfg; progress_desc = "test")
    out = mktempdir()
    fpath = save_netcdf(g, cfg.start_date, cfg.end_date, out, cfg)
    @test isfile(fpath)
    NCDatasets.Dataset(fpath, "r") do ds
        @test haskey(ds, "sif_radiance_678nm")
        @test ds.attrib["exclude_missing_nflh"] == "true"
    end
    println("  wrote ", fpath)
end

println("All rasterize tests passed.")
