# =========================================================================================================
# This script pre-processes PACE OCI data, to:
# 1) get TOA radiance
# 2) merge with the corresponding geolocation data
# =========================================================================================================

using NCDatasets

# functino to rescale data
function rescale_data(
        x,   # input packed data
        scale_factor,  # scale factor
        offset,        # offset
    )
    # Rescale data
    x_rescaled = (x .- offset) ./ scale_factor;
    return x_rescaled
end

# function to calculate Rtoa from radiance and solar irradiance
function TOA_radiance(ds::Dataset)
    # println("Creating TOA radiance variable...")

    # Get the NetCDF variable (not the data)
    rhot_red_var = ds["rhot_red"]
    
    # Get dimension names from the variable
    dims = dimnames(rhot_red_var)
    data_type = eltype(rhot_red_var)
    
    # Load the data
    rhot_red = rhot_red_var[:]
    solar_irradiance = ds["red_solar_irradiance"][:]
    solar_zenith_angle = ds["solar_zenith"][:]
    earth_sun_correction = ds.attrib["earth_sun_distance_correction"]   # unitless

    # Ensure dimensions are compatible for broadcasting
    solar_irradiance = reshape(solar_irradiance, (1, 1, size(solar_irradiance)...))
    solar_zenith_angle = reshape(solar_zenith_angle, (size(solar_zenith_angle)..., 1))

    # Shape check
    # println("  | Radiance shape: ", size(rhot_red))
    # println("  | Solar irradiance shape: ", size(solar_irradiance))
    # println("  | Solar zenith angle shape: ", size(solar_zenith_angle))
    # println("  | Earth-sun distance correction: ", earth_sun_correction)

    # Calculate TOA radiance
    Rtoa = rhot_red .* solar_irradiance .* cosd.(solar_zenith_angle) ./ π / earth_sun_correction
    
    # Get the non-missing type for NetCDF
    nc_type = Missing <: data_type ? Base.nonmissingtype(data_type) : data_type
    
    # Create new variable with proper dimensions
    Rtoa_var = defVar(ds, "Rtoa_red", nc_type, dims, fillvalue=nc_type(-9999.0))
    
    # Add attributes
    Rtoa_var.attrib["long_name"] = "Top-of-Atmosphere radiance at Red Band"
    Rtoa_var.attrib["units"] = "W/m^2/sr/μm"
    
    # Write data
    Rtoa_var[:] = Rtoa
    
    # println("  | TOA radiance variable created with shape: ", size(Rtoa_var[:]))
    
    return ds
end

# subset the data
function subset_netcdf_dataset(
        filepath::String, 
        selected_vars::Vector{String}, 
        output_path::String; 
        indices::Dict=Dict(), 
        prefix_groups::Bool=false,
        rename_dims=Dict(),
        post_process_func::Union{Function, Nothing}=nothing
    ) 
    # Open input dataset 
    ds = Dataset(filepath)

    # Create output filename
    output_name = "subset_" * basename(filepath)
    output_file = joinpath(output_path, output_name)

    # Create output dataset, overwrite if exists
    merged_data = Dataset(output_file, "c")

    try
        # Keep track of defined dimensions
        defined_dims = Set{String}()
        
        for group_name in keys(ds.group)
            group = ds.group[group_name]
            
            for var_name in keys(group)
                # Check if this variable is in the selected list
                if var_name in selected_vars
                    var = group[var_name]
                    
                    # Determine output variable name
                    key = prefix_groups ? "$(group_name)_$(var_name)" : var_name
                    
                    # Get dimensions
                    var_dims = dimnames(var)

                    # println("Variable: $var_name")
                    
                    # Build indexing tuple for subsetting
                    idx_tuple = []
                    subset_dims = []
                    subset_sizes = []
                    
                    for (i, dim_name) in enumerate(var_dims)
                        if haskey(indices, dim_name)
                            # Use provided indices
                            push!(idx_tuple, indices[dim_name])
                            push!(subset_dims, dim_name)
                            push!(subset_sizes, length(indices[dim_name]))
                        else
                            # Use all indices
                            push!(idx_tuple, :)
                            push!(subset_dims, dim_name)
                            push!(subset_sizes, size(var, i))
                        end
                    end
                    
                    # Define dimensions in output dataset if not already defined
                    for (dim_name, dim_size) in zip(subset_dims, subset_sizes)
                        # Use renamed dimension if mapping exists
                        output_dim_name = haskey(rename_dims, dim_name) ? rename_dims[dim_name] : dim_name
                        
                        if !(output_dim_name in defined_dims)
                            defDim(merged_data, output_dim_name, dim_size)
                            push!(defined_dims, output_dim_name)
                        end
                    end
                    
                    # Create variable in output dataset
                    var_type = eltype(var)
                    # println("  | Original type: ", var_type)

                    nc_type = if Missing <: var_type
                        Base.nonmissingtype(var_type)
                    else
                        var_type
                    end
                    # println("  | Output type: ", nc_type)

                    # Create variable
                    fillvalue = haskey(var.attrib, "_FillValue") ? var.attrib["_FillValue"] : nothing
                    renamed_subset_dims = [haskey(rename_dims, d) ? rename_dims[d] : d for d in subset_dims]
                    if fillvalue !== nothing
                        merged_var = defVar(
                            merged_data, 
                            key, 
                            nc_type, 
                            tuple(renamed_subset_dims...),
                            fillvalue = fillvalue
                        )
                    else 
                        merged_var = defVar(
                            merged_data, 
                            key, 
                            nc_type, 
                            tuple(renamed_subset_dims...)
                        )
                    end

                    # Copy attributes, skipping _FillValue and missing_value
                    for (attr_name, attr_value) in var.attrib
                        # println("  | Attribute: $attr_name = $attr_value")
                        if attr_name ∉ ["_FillValue", "missing_value"]
                            try
                                merged_var.attrib[attr_name] = attr_value
                            catch e
                                @warn "Could not copy attribute $attr_name for $key: $e"
                            end
                        end
                    end

                    # rescale data by checking whether scale_factor and add_offset exist
                    if haskey(var.attrib, "scale_factor") && haskey(var.attrib, "add_offset")
                        scale_factor = var.attrib["scale_factor"]
                        add_offset = var.attrib["add_offset"]
                        # println("  | Rescaling data with scale_factor=$scale_factor and add_offset=$add_offset")
                        # Load original data
                        original_data = var[idx_tuple...]
                        # Rescale data
                        rescaled_data = rescale_data(original_data, scale_factor, add_offset)
                        # Write rescaled data to output variable
                        merged_var[:] = rescaled_data
                        # add a flag attribute to indicate that the data has been rescaled
                        merged_var.attrib["rescaled"] = "true"
                    else
                        # otherwise keep the original data (NCDatasets handles fill values automatically)
                        merged_var[:] = var[idx_tuple...]
                    end
                    # println("  ✓ Copied variable: $key ($(subset_sizes))")
                end
            end
        end
        
        # Copy global attributes
        for (attr_name, attr_value) in ds.attrib
            merged_data.attrib[attr_name] = attr_value
        end

        # post process if needed
        if post_process_func !== nothing
            post_process_func(merged_data)
        end
        # println("  ✓ Subset created: $output_file")
        
    finally
        close(ds)
        close(merged_data)
    end

    return output_file
end