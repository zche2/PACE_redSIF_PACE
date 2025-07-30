#=
This file generates configuration files, executes the interpolation model, and manages logging.
to run the final, type:
nohup julia --project=.. --threads 4 [Filename] > [log file name] 2>&1 &
=#
using JLD2
using Parameters  # For @with_kw macro (Pkg.add("Parameters"))
using YAML        # For writing YAML config files (Pkg.add("YAML"))
using Base.Filesystem # For mkpath, joinpath, abspath

using Base.Threads # Essential for threading primitives
using Distributed 

using vSmartMOM, vSmartMOM.Absorption
import vSmartMOM.Absorption: make_interpolation_model

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
    MinPressure_Pa::FT           = 1.
    MaxPressure_Pa::FT           = 105000.
    deltaPressure_Pa::FT         = 1000.
    PressureGrid_dir::String     = "default"      # if customized grid, need to specify the jld2 file used

    UnevenGrid_Temperature::Bool = false
    MinTemp_K::FT                = 150.
    MaxTemp_K::FT                = 330.
    deltaTemp_K::FT              = 5.
    TempGrid_dir::String         = "default"      # if customized grid, need to specify the jld2 file used

    # --- Wavenumber of interest ---
    res::FT      = 0.1
    ν_min::FT    = 11627
    ν_max::FT    = 16150
    
end

if abspath(PROGRAM_FILE) == @__FILE__
    
    # ================== specify the model =====================
    println("--- make interpolation model ---")

    # Model Params 🌟
    MyModel = ModelParameters{Float64}(
        simulation_id       = "Finer_Wavenumber_grid_Jul11",
        res                 = 0.01
        # deltaTemp_K         = 1.,
        # UnevenGrid_Pressure = true,
        # PressureGrid_dir    = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Pressure_UpperAtm.jld2",
    )
    println(MyModel)

    # 🌟

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

    @sync begin
    
        for (molec, num) in zip(MyModel.gas_species, MyModel.gas_id)

            Threads.@spawn begin   
                println("\n--- Spawning task for $molec on $(Threads.threadid()) ---")
                flush(stdout)        # <--- IMPORTANT: Force flush stdout
                par  = Absorption.read_hitran(artifact(molec), mol=num, iso=1, ν_min=ν_min, ν_max=ν_max);
                file = replace(OutPattern, r"MOLEC" => molec)

                itp_model = make_interpolation_model(
                                    par,
                                    Voigt(),
                                    ν_grid,
                                    p_grid,
                                    t_grid,
                                    wavelength_flag=false,
                                    wing_cutoff=10,
                                );
                save_interpolation_model(itp_model, file)
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


