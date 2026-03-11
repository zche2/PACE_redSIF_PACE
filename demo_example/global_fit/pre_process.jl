# =========================================================================================================
# Pre-processes PACE OCI data:
# 1) rescale packed NetCDF data
# 2) compute TOA radiance from rhot
# 3) subset and merge variables from L1B/L2 groups into flat NetCDF
# =========================================================================================================

using NCDatasets

# Function to rescale packed data
function rescale_data(
    x,
    scale_factor,
    offset,
)
    return (x .- offset) ./ scale_factor
end

# Find variable and return (var, dims) from dataset (searches groups and root)
function _find_var_from_dataset(ds, varname::String)
    groups_to_check = Pair{String, Any}[]
    try
        if hasproperty(ds, :group)
            for group_name in keys(ds.group)
                push!(groups_to_check, group_name => ds.group[group_name])
            end
        end
    catch
        nothing
    end
    if isempty(groups_to_check)
        push!(groups_to_check, "" => ds)
    end
    for (_, group) in groups_to_check
        haskey(group, varname) || continue
        var = group[varname]
        try
            dimnames(var)
        catch
            continue
        end
        return (var=var, dims=collect(String.(dimnames(var))))
    end
    # Also check root
    if haskey(ds, varname)
        var = ds[varname]
        return (var=var, dims=collect(String.(dimnames(var))))
    end
    error("Variable '$varname' not found in dataset")
end

# Read var with optional rescale (scale_factor, add_offset)
function _read_var_rescaled(var)
    data = var[:]
    if haskey(var.attrib, "scale_factor") && haskey(var.attrib, "add_offset")
        data = rescale_data(data, var.attrib["scale_factor"], var.attrib["add_offset"])
    end
    return data
end

# Compute Rtoa from rhot: Rtoa = rhot * E * cos(sza) / π / earth_sun_distance_correction
function TOA_radiance(ds::Dataset)
    rhot_red_var = ds["rhot_red"]
    dims = dimnames(rhot_red_var)
    data_type = eltype(rhot_red_var)

    rhot_red = rhot_red_var[:]
    solar_irradiance = ds["red_solar_irradiance"][:]
    solar_zenith_angle = ds["solar_zenith"][:]
    earth_sun_correction = ds.attrib["earth_sun_distance_correction"]

    solar_irradiance = reshape(solar_irradiance, (1, 1, size(solar_irradiance)...))
    solar_zenith_angle = reshape(solar_zenith_angle, (size(solar_zenith_angle)..., 1))

    Rtoa = rhot_red .* solar_irradiance .* cosd.(solar_zenith_angle) ./ π / earth_sun_correction

    nc_type = Missing <: data_type ? Base.nonmissingtype(data_type) : data_type
    Rtoa_var = defVar(ds, "Rtoa_red", nc_type, dims, fillvalue=nc_type(-9999.0))
    Rtoa_var.attrib["long_name"] = "Top-of-Atmosphere radiance at Red Band"
    Rtoa_var.attrib["units"] = "W/m^2/sr/μm"
    Rtoa_var[:] = Rtoa

    return ds
end

# Subset NetCDF dataset: copy selected vars from groups into flat output with optional rename_dims and post_process
function subset_netcdf_dataset(
    filepath::String,
    selected_vars::Vector{String},
    output_path::String;
    indices::Dict=Dict(),
    prefix_groups::Bool=false,
    rename_dims=Dict(),
    post_process_func::Union{Function, Nothing}=nothing,
)
    ds = Dataset(filepath)
    output_name = "subset_" * basename(filepath)
    output_file = joinpath(output_path, output_name)
    merged_data = Dataset(output_file, "c")

    try
        defined_dims = Set{String}()

        # Collect groups to iterate: child groups (NetCDF4) or root (flat NetCDF)
        groups_to_check = Pair{String, Any}[]
        try
            if hasproperty(ds, :group)
                for group_name in keys(ds.group)
                    push!(groups_to_check, group_name => ds.group[group_name])
                end
            end
        catch
            nothing
        end
        if isempty(groups_to_check)
            push!(groups_to_check, "" => ds)
        end

        for (group_name, group) in groups_to_check
            for var_name in keys(group)
                var_name in selected_vars || continue
                var = group[var_name]
                # Skip subgroups (only process variables)
                try
                    dimnames(var)
                catch
                    continue
                end

                key = prefix_groups && !isempty(group_name) ? "$(group_name)_$(var_name)" : var_name
                var_dims = dimnames(var)

                idx_tuple = []
                subset_dims = []
                subset_sizes = []

                for (i, dim_name) in enumerate(var_dims)
                    if haskey(indices, dim_name)
                        push!(idx_tuple, indices[dim_name])
                        push!(subset_dims, dim_name)
                        push!(subset_sizes, length(indices[dim_name]))
                    else
                        push!(idx_tuple, :)
                        push!(subset_dims, dim_name)
                        push!(subset_sizes, size(var, i))
                    end
                end

                for (dim_name, dim_size) in zip(subset_dims, subset_sizes)
                    output_dim_name = haskey(rename_dims, dim_name) ? rename_dims[dim_name] : dim_name
                    if !(output_dim_name in defined_dims)
                        defDim(merged_data, output_dim_name, dim_size)
                        push!(defined_dims, output_dim_name)
                    end
                end

                var_type = eltype(var)
                nc_type = Missing <: var_type ? Base.nonmissingtype(var_type) : var_type
                fillvalue = get(var.attrib, "_FillValue", nothing)
                renamed_subset_dims = [haskey(rename_dims, d) ? rename_dims[d] : d for d in subset_dims]

                merged_var = if fillvalue !== nothing
                    defVar(
                        merged_data,
                        key,
                        nc_type,
                        tuple(renamed_subset_dims...),
                        fillvalue=fillvalue,
                    )
                else
                    defVar(merged_data, key, nc_type, tuple(renamed_subset_dims...))
                end

                for (attr_name, attr_value) in var.attrib
                    attr_name ∈ ["_FillValue", "missing_value"] && continue
                    try
                        merged_var.attrib[attr_name] = attr_value
                    catch e
                        @warn "Could not copy attribute $attr_name for $key: $e"
                    end
                end

                if haskey(var.attrib, "scale_factor") && haskey(var.attrib, "add_offset")
                    scale_factor = var.attrib["scale_factor"]
                    add_offset = var.attrib["add_offset"]
                    original_data = var[idx_tuple...]
                    rescaled_data = rescale_data(original_data, scale_factor, add_offset)
                    merged_var[:] = rescaled_data
                    merged_var.attrib["rescaled"] = "true"
                else
                    merged_var[:] = var[idx_tuple...]
                end
            end
        end

        for (attr_name, attr_value) in ds.attrib
            merged_data.attrib[attr_name] = attr_value
        end

        if post_process_func !== nothing
            post_process_func(merged_data)
        end
    finally
        close(ds)
        close(merged_data)
    end

    return output_file
end
