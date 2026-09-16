# Load one open-ocean column from geoschem_ocean_columns_n500.nc into
# params.scattering_params.rt_aerosols (Gaussian-in-pressure RT_Aerosol path).
#
# Requires `using vSmartMOM` (and thus CoreRT / Scattering / Aerosols) in the caller.

import YAML as _YAML
using NCDatasets
using Distributions
using Printf

const _RI_DB_PATH = joinpath(dirname(dirname(pathof(vSmartMOM))), "data", "refractive_indices_database.yaml")
const _AOD_MIN = 1e-6
const _SIGP_MIN = 20.0  # hPa

"""
    load_ocean_column_aerosols!(params, nc_path, i_sample; yaml_path, ri_db_path)

Read column `i_sample` (1-based) from the n500 ocean-column NetCDF and replace
`params.scattering_params.rt_aerosols` with one `RT_Aerosol` per species listed
under `aerosol_scheme.species` in `yaml_path`.

For each species:
- τ_ref = column sum of layer AOD
- μ = AOD-weighted mean wet radius (μm)
- σ = sigma_g from YAML
- profile = Normal(p₀, σp) with AOD-weighted mean/std of `p_mid`
- nᵣ, nᵢ from the vSmartMOM refractive-index database at λ_ref
"""
function load_ocean_column_aerosols!(
    params,
    nc_path::AbstractString,
    i_sample::Integer;
    yaml_path::AbstractString,
    ri_db_path::AbstractString=_RI_DB_PATH,
    aod_min::Real=_AOD_MIN,
)
    sp = params.scattering_params
    sp === nothing && error("params.scattering_params is nothing; add a scattering: block to the YAML")

    cfg = _YAML.load_file(yaml_path)
    scheme = get(cfg, "aerosol_scheme", Dict())
    species_cfg = get(scheme, "species", Dict())
    isempty(species_cfg) && error("No aerosol_scheme.species in $yaml_path")

    λ_ref = Float64(sp.λ_ref)
    r_max = Float64(sp.r_max)
    FT = typeof(sp.λ_ref)

    ri_db = vSmartMOM.Aerosols.load_refractive_index_database(ri_db_path, Float64)

    ds = NCDataset(nc_path)
    n_samp = ds.dim["sample"]
    n_lev = ds.dim["lev"]
    1 ≤ i_sample ≤ n_samp || error("i_sample=$i_sample out of range 1:$n_samp")
    p_mid = Float64.(ds["p_mid"][:])
    length(p_mid) == n_lev || error("p_mid length $(length(p_mid)) ≠ n_lev=$n_lev")
    lat_i = Float64(ds["lat"][i_sample])
    lon_i = Float64(ds["lon"][i_sample])

    # NetCDF may present (lev, sample) or (sample, lev) depending on the writer.
    function column_profile(varname)
        haskey(ds, varname) || error("Missing var $varname in $nc_path")
        A = Array(ds[varname])
        ndims(A) == 2 || error("$varname should be 2-D, got ndims=$(ndims(A))")
        if size(A) == (n_lev, n_samp)
            return Float64.(A[:, i_sample])
        elseif size(A) == (n_samp, n_lev)
            return Float64.(A[i_sample, :])
        else
            error("$varname shape $(size(A)) incompatible with (lev=$n_lev, sample=$n_samp)")
        end
    end

    rt_list = vSmartMOM.CoreRT.RT_Aerosol{FT}[]
    summaries = String[]

    for (name, sc) in sort(collect(species_cfg); by=first)
        aod_var = String(sc["aod_var"])
        rad_var = String(sc["radius_var"])
        sigma_g = Float64(sc["sigma_g"])
        ri_key = String(sc["refractive_index"])

        aod = column_profile(aod_var)
        rad = column_profile(rad_var)
        length(aod) == n_lev || error("$aod_var column length $(length(aod)) ≠ n_lev=$n_lev")
        τ_ref = sum(aod)
        τ_ref < aod_min && continue

        w = aod ./ τ_ref
        μ = clamp(sum(w .* rad), 1e-3, r_max)
        p₀ = sum(w .* p_mid)
        σp = max(sqrt(sum(w .* (p_mid .- p₀) .^ 2)), _SIGP_MIN)

        n_c = try
            vSmartMOM.Aerosols.get_refractive_index(ri_db, ri_key, λ_ref)
        catch e
            @warn "RI lookup failed; using seasalt-like fallback" species=name ri_key=ri_key exception=e
            1.5 + 1e-8im
        end
        nᵣ = FT(real(n_c))
        nᵢ = FT(abs(imag(n_c)))  # Mie path expects non-negative imag

        size_dist = LogNormal(log(μ), log(sigma_g))
        aero = vSmartMOM.Scattering.Aerosol(size_dist, nᵣ, nᵢ)
        profile = Normal(FT(p₀), FT(σp))
        push!(rt_list, vSmartMOM.CoreRT.RT_Aerosol(aero, FT(τ_ref), profile))
        push!(summaries, @sprintf("%s: τ=%.4f μ=%.3fμm p₀=%.0fhPa σp=%.0f", name, τ_ref, μ, p₀, σp))
    end
    close(ds)

    isempty(rt_list) && error("No species with AOD ≥ $aod_min in column $i_sample")
    sp.rt_aerosols = rt_list

    println("Ocean column $i_sample @ ($(round(lat_i; digits=1))°, $(round(lon_i; digits=1))°): $(length(rt_list)) aerosols")
    for s in summaries
        println("  ", s)
    end
    return params
end
