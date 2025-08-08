function layer_VCD(
                    p,                          # pressure: default-hPa
                    q;                          # specific humidity
                    n_layers :: Int = 72,       # number of layers
                    use_hPa :: Bool = true,
                    g₀ = 9.8196,                
                    Na = 6.0221415e23           # Avogradro's number
                )
	
    # Dry and wet mass
    dryMass = 28.9647e-3  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = 18.01528e-3 / Na  # just H2O
    ratio = dryMass / wetMass 

	vmr_h2o = zeros(Float64, n_layers)
    vcd_dry = zeros(Float64, n_layers)
    vcd_h2o = zeros(Float64, n_layers)
	
    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = use_hPa ? ( (p[i + 1] - p[i]) * 100. ) : (p[i + 1] - p[i])
		
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

	return vcd_dry, vcd_h2o, vmr_h2o
end