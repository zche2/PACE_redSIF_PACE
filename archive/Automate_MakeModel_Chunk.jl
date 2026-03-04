#=
This file generates configuration files, executes the interpolation model, and manages logging.
to run the final, type:
nohup julia --project=.. --threads 4 [Filename] > [log file name] 2>&1 &
=#
using JLD2
using Parameters  # For @with_kw macro (Pkg.add("Parameters"))
using YAML        # For writing YAML config files (Pkg.add("YAML"))
using Base.Filesystem # For mkpath, joinpath, abspath

using Interpolations

using Base.Threads # Essential for threading primitives
using Distributed 

using vSmartMOM, vSmartMOM.Absorption
import vSmartMOM.Absorption: make_interpolation_model_reduced

# --- Define a struct to hold all model parameters ---
@with_kw struct ModelParameters{FT <: AbstractFloat}
    #=
    e.g., ModelParameters{Float64}(deltaTemp_K=4)
    =#

    # --- General Simulation Parameters ---
    simulation_id::String    = "default_run" # Unique ID for the simulation
    
    # Paths relative to the script's execution directory
    results_base_dir::String = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection"

    # --- Species to be considered ---
    gas_species::Vector{String}  = ["H2O", "O2"]
    gas_id::Vector{Int64}        = [1, 7]   

    # --- Temperature and pressure range
    UnevenGrid_Pressure::Bool    = false
    MinPressure_Pa::FT           = 19250.
    MaxPressure_Pa::FT           = 105000.
    deltaPressure_Pa::FT         = 20.
    PressureGrid_dir::String     = "default"      # if customized grid, need to specify the jld2 file used

    UnevenGrid_Temperature::Bool = false
    MinTemp_K::FT                = 150.
    MaxTemp_K::FT                = 330.
    deltaTemp_K::FT              = .5
    TempGrid_dir::String         = "default"      # if customized grid, need to specify the jld2 file used

    # --- Wavenumber of interest ---
    res::FT      = 2
    ν_min::FT    = 11627
    ν_max::FT    = 16129

    # --- Paralle processing ---
    chunk_size::Int64   = 500
end

# function MakeGrid(
#             # Required
#             ScaleFactor::AbstractFloat,
#             # Optional
#             MinVal::AbstractFloat = 0.,
#             MaxVal::AbstractFloat = 0.,
#             delta::AbstractFloat  = 1.,
#             GridDir::String   = "Default",
#             Uneven_grid::bool = false,
#         )
#     # to be updated 🌟
#     if Uneven_grid
#         @load GridDir t_grid???
#     else

#     end

# end

function ChunkData(
    MyData::Union{AbstractRange{<:Real}, Vector{<:AbstractFloat}};
    chunk_size::Int64 = 500,
    )
    # ========= create chunks for multi thread processing =========
    FT        = eltype(collect(MyData))
    chunks    = Vector{Vector{FT}}()
    data_size = length(MyData)

    for i in 1:chunk_size:data_size
        push!(chunks, MyData[i : min(i + chunk_size - 1, data_size)])
    end

    return chunks

end

if abspath(PROGRAM_FILE) == @__FILE__
    
    # ================== specify the model =====================
    println("--- make interpolation model ---")

    # Model Params
    MyModel = ModelParameters{Float64}(
        simulation_id       = "UpperAtm_Jul04_Parallel2_HierarchicalThreads",
        UnevenGrid_Pressure = true,
        PressureGrid_dir    = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Pressure_UpperAtm.jld2",
        chunk_size          = 15,
        gas_species         = ["H2O"],
        gas_id              = [1],
    )
    println(MyModel)

    # ==================== generate grids ======================
    println("--- generate grids ---")
    ν_min  = MyModel.ν_min
    ν_max  = MyModel.ν_max
    ν_grid = ν_min:MyModel.res:ν_max
    println("ν_grid created, length = ", length(ν_grid))

    if MyModel.UnevenGrid_Temperature
        @load MyModel.TempGrid_dir t_grid
    else
        t_grid = MyModel.MinTemp_K : MyModel.deltaTemp_K : MyModel.MaxTemp_K
    end
    println("t_grid created, length = ", length(t_grid))

    if MyModel.UnevenGrid_Pressure
        @load MyModel.PressureGrid_dir p_grid
    else
        p_grid = MyModel.MinPressure_Pa : MyModel.deltaPressure_Pa : MyModel.MaxPressure_Pa
    end
    # convert to hPa 
    p_grid /= 100.
    println("p_grid created, converted to hPa, length = ", length(p_grid))

    # create folder to store
    MyFolder = joinpath(MyModel.results_base_dir, MyModel.simulation_id)
    mkpath(MyFolder)

    # Output file name pattern
    OutPattern = joinpath(MyFolder, MyModel.simulation_id * "_MOLEC.jld2")

    # ======================== threading ==========================
    # Check if threading is enabled
    if Threads.nthreads() == 1
        println("WARNING: Julia not started with multiple threads. Use `julia --threads auto` or `JULIA_NUM_THREADS`.")
        println("Running sequentially for demonstration.")
    else
        println("Running with $(Threads.nthreads()) threads.")
    end

    flush(stdout) 

    # ============ Pre-chunk Pressure data for parallel processing ===============
    chunks = ChunkData(p_grid, chunk_size = MyModel.chunk_size)
    
    @sync begin
    
        for (molec, num) in zip(MyModel.gas_species, MyModel.gas_id)

            Threads.@spawn begin   
                println("\n--- Spawning task for $molec on $(Threads.threadid()) ---")
                flush(stdout)        # <--- IMPORTANT: Force flush stdout
                par  = Absorption.read_hitran(artifact(molec), mol=num, iso=1, ν_min=ν_min, ν_max=ν_max);
                file = replace(OutPattern, r"MOLEC" => molec)

                processed_chunks_list = Vector{Array{Float64, 3}}(undef, length(chunks))
                
                # This iteration runs on different threads    
                @threads for chunk_idx in eachindex(chunks) # <--- INNER PARALLELISM
                    
                    p_chunk = chunks[chunk_idx]
                    println(" Mol '$molec': Thread $(Threads.threadid()) - Processed chunk $(chunk_idx).") 
                    println("    Some values in chunk $(chunk_idx): $(p_chunk[1:5])")
                    flush(stdout) 
                    processed_chunks_list[chunk_idx] = make_interpolation_model_reduced(
                                                                                        par,
                                                                                        Voigt(),
                                                                                        ν_grid,
                                                                                        p_chunk,
                                                                                        t_grid,
                                                                                        wavelength_flag=false,
                                                                                        wing_cutoff=10,
                                                                                    );
                    println("    Thread $(Threads.threadid()) completed!")
                    flush(stdout) 
                end

                # concate the final results
                final_processed_data = hcat(processed_chunks_list...)

                # Perform the interpolation
                itp = interpolate(final_processed_data, BSpline(Cubic(Line(OnGrid()))))

                # Get the molecule and isotope numbers from the HitranTable
                mol = par.mol[1]
                iso = all(x->x==par.iso[1], par.iso) ? par.iso[1] : -1

                # Return the interpolation model with all the proper parameters
                itp_model = InterpolationModel(itp, mol, iso,  ν_grid, p_grid, t_grid)

                println("\nConstructed => Saved: $file") 
                flush(stdout)        # <--- IMPORTANT: Force flush stdout
            end

        end

    end

end



#=
function generate_config_file(params::ModelParameters)
    filedir     = params.results_base_dir
    output_dir  = joinpath(filedir, "$(params.simulation_id)")
    # Ensure the output directory exists
    mkpath(output_dir)
    config_path = joinpath(output_dir, "$(params.simulation_id)_config.yaml")

    # generate grids
    ν_grid = ν_min:res:ν_max

    if MyModel.UnevenGrid_Temperature
        @load MyModel.TempGrid_dir t_grid
    else
        t_grid = MyModel.MinTemp_K : MyModel.deltaTemp_K : MyModel.MaxTemp_K
    end

    if MyModel.UnevenGrid_Pressure
        @load MyModel.PressureGrid_dir p_grid
    else
        p_grid = MyModel.MinPressure_Pa : MyModel.deltaPressure_Pa : MyModel.MaxPressure_Pa
    end
    # convert to hPa 
    p_grid /= 100.

end
=#


