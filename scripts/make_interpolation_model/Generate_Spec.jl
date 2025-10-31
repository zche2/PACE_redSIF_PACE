using vSmartMOM, vSmartMOM.Absorption

function main()

    res   = 4.          # cm^{-1} - default: 2
    ν_min = 11627
    ν_max = 16129
    ν_grid = ν_min:res:ν_max;

    dT    = 5.       # K - default: 5
    dp    = 50.      # Pa - default: 50.
    T_min = 150.     # K
    T_max = 330.     # K
    p_min = 1.       # Pa
    p_max = 105000.  # Pa - default: 105000
    T_grid = T_min:dT:T_max;
    p_grid = p_min:dp:p_max;
    p_grid /= 100; # convert to hPa

    SaveDir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/"

    println("T, p, ν grid created!\n")


    # hitran model 
    AllMolec  = ["O3"];
    AllNum    = [3];
    
    for (molec, num) in zip(AllMolec, AllNum)
    
        println("Now working on $molec")
        flush(stdout) 
        par = Absorption.read_hitran(artifact(molec), mol=num, iso=1, ν_min=ν_min, ν_max=ν_max);
        filepath = joinpath(SaveDir, "interp_xSection_$(molec)_Coarse_Grid.jld2");

        itp_model = make_interpolation_model(
                            par,
                            Voigt(),
                            ν_grid,
                            p_grid,
                            T_grid,
                            wavelength_flag=false,
                            wing_cutoff=10,
                        );
        save_interpolation_model(itp_model, filepath)
        println("Constructed => Saved: $filepath") 
        flush(stdout)        

    end

end

main()