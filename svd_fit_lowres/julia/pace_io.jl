# Thin path helpers for SVD pipeline data files (absolute or relative to base_dir).

"""
    resolve_data_path(path, base_dir) -> String

If `path` is absolute, return it; otherwise join with `base_dir`.
"""
function resolve_data_path(path::AbstractString, base_dir::AbstractString)
    p = String(path)
    return isabspath(p) ? p : joinpath(String(base_dir), p)
end

"""
    resolve_svd_data_paths(data_cfg) -> NamedTuple

Resolve `[data]` paths used by SVD retrieval:
`summer_nc`, `winter_nc`, `sif_file`, `pace_snr_file` (optional).
"""
function resolve_svd_data_paths(data_cfg::AbstractDict)
    base_dir = String(get(data_cfg, "base_dir", ""))
    summer_nc = get(data_cfg, "summer_nc", nothing)
    winter_nc = get(data_cfg, "winter_nc", nothing)
    sif_file = get(data_cfg, "sif_file", nothing)
    pace_snr_file = get(data_cfg, "pace_snr_file", nothing)

    summer_abs = summer_nc === nothing ? nothing : resolve_data_path(String(summer_nc), base_dir)
    winter_abs = winter_nc === nothing ? nothing : resolve_data_path(String(winter_nc), base_dir)
    sif_path = sif_file === nothing ? nothing : resolve_data_path(String(sif_file), base_dir)
    pace_snr_path = pace_snr_file === nothing || isempty(String(pace_snr_file)) ?
        nothing : resolve_data_path(String(pace_snr_file), base_dir)

    return (
        base_dir = base_dir,
        summer_nc = summer_abs,
        winter_nc = winter_abs,
        sif_path = sif_path,
        pace_snr_path = pace_snr_path,
    )
end
