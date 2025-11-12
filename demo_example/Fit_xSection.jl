### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d0d1237b-2c32-46ed-9060-c12ba438810f
begin
	import Pkg; 
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");
	# Pkg.develop(path="/home/zhe2/FraLab/PACE_redSIF_PACE")
end

# ╔═╡ f3191894-fa07-4ecb-9253-6b9375fd17df
using JLD2, Interpolations, Revise

# ╔═╡ 9896ec16-ce32-4edd-8542-793a8a14f67c
using Polynomials, ForwardDiff, DiffResults, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ 8c72db65-85a4-4dce-aca6-ab25174e561e
using Plots, PlutoUI, Printf, Dates

# ╔═╡ f88aaa45-088e-456c-8efa-22d3215e7e10
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools

# ╔═╡ c92250bf-7299-4c69-becd-a356015f04bd
using PACE_SIF

# ╔═╡ b68ca8b3-df51-447d-ac8e-830f9d3e4f29
md"""
## Cross Section Fit
Created: 2025-11-09
"""

# ╔═╡ c7d1022a-b2b2-4cbb-bc7f-fc3d6cd5b782
begin
	o2_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2";
	o2_sitp = read_rescale(o2_jld2);

	h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2";
	h2o_sitp = read_rescale(h2o_jld2);

	metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log";

	ν_grid, p_grid_hPa, t_grid = o2_sitp.ranges;
	println(
		"wavenumber: Min=$(ν_grid[1]), Max=$(ν_grid[end]), res=$(ν_grid[2]-ν_grid[1])"
	)
	println(
		"pressure (hPa): Min=$(p_grid_hPa[1]), Max=$(p_grid_hPa[end]), Δp=$(p_grid_hPa[2]-p_grid_hPa[1])"
	)
	println(
		"temperature (K): Min=$(t_grid[1]), Max=$(t_grid[end]), ΔT=$(t_grid[2]-t_grid[1])"
	)
end

# ╔═╡ d989b6f1-f079-4e45-96e1-b86103bfb38e
begin
	# set apparent height (pressure) and temperature for 
	# O2 and water vapor, respectively
	pres_O2  = 1030.;   # hPa
	temp_O2  = 200.;    # K
	pres_H2O = 1030.;   # hPa
	temp_H2O = 250.;    # K
	# obtain the high-resolution cross section
	xSec_O2  = o2_sitp(ν_grid, pres_O2, temp_O2); 
	xSec_H2O = h2o_sitp(ν_grid, pres_H2O, temp_H2O);
	# plot
	p = plot(
		size=(800, 400), 
		# yscale=:log10,
		xlabel="wavenumber",
		ylabel="cross section",
		margin=5Plots.mm
	)
	plot!(p, ν_grid, xSec_O2, label="O₂")
	plot!(p, ν_grid, xSec_H2O, label="H₂O")
end

# ╔═╡ cfad865f-d09b-4d98-b50a-bb6c39542a02
begin
	ds_summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc");
	# estimate total vertical column density
	k = 3220;
	vcd_dry = ds_summer["vcd_dry"][k,:];
	vmr_o2  = 0.21;
	vcd_h2o = ds_summer["vcd_h2o"][k,:];
	vmr_h2o = ds_summer["vmr_h2o_var"][k,:];
	∑vcd_o2  = sum(vcd_dry .* vmr_o2);     # molec/cm2
	∑vcd_h2o = sum(vcd_h2o);    # molec/cm2
	println("Total column density of O₂ (molec/cm²): $∑vcd_o2")
	println("Total column density of H₂O (molec/cm²): $∑vcd_h2o")
	# get optical depth
	τ_o2    = xSec_O2 .* ∑vcd_o2;
	τ_h2o   = xSec_H2O .* ∑vcd_h2o;
	τ       = τ_o2 .+ τ_h2o;
	# transmittance
	T_o2    = exp.( - τ_o2 );
	T_h2o   = exp.( - τ_h2o );
	T_tot   = exp.( - τ );
	# plot
	p₁ = plot(
		size=(800, 400), 
		# yscale=:log10,
		xlabel="wavenumber",
		ylabel="transmittance (µ=1.0)",
		margin=5Plots.mm
	)
	plot!(p₁, ν_grid, T_tot, label="O₂ + H₂O")
	plot!(p₁, ν_grid, T_o2, label="O₂")
	plot!(p₁, ν_grid, T_h2o, label="H₂O")
end

# ╔═╡ 1d3c114f-fdd3-4dd0-b66a-17f13e98185a
begin
	# Spectral response function
	# read data
	filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_RSRs.nc";
	pace = Dataset(filename, "r");

	wavlen = pace["wavelength"][:];
	RSR    = pace["RSR"];
	band   = pace["bands"];

	if_wavenumber = true;

	λ = if_wavenumber ? reverse(1 ./ collect(ν_grid) .* 1E7) : (1 ./ collect(ν_grid) .* 1E7);

	ind₁   = findall( λ[1] .< wavlen .< λ[end]);
	ind₂   = findall( λ[1] .< band   .< λ[end]);
	λ_msr  = wavlen[ind₁];
	MyRSR  = RSR[ind₁, ind₂];

	MyKernel = KernelInstrument(
		band[ind₂],
		λ_msr,
		MyRSR,
		λ,
		ν_grid.step.hi
	)
end

# ╔═╡ 3c626557-4612-4b18-8a9e-58bfe739b8a4
md"""
> ##### note: constructing this `Kernel Instrument` takes time, need precompile.
"""

# ╔═╡ 28a29145-c387-4220-b569-761d3d797885
begin
	# reverse the matrix if originally using wavenumber
	T_for_anly = reverse(T_tot);
	# convolution
	T_conv = conv_matx(MyKernel, T_for_anly);

	# plot
	p₂ = plot(
		size=(800, 400), 
		xlabel="wavelength [nm]",
		ylabel="transmittance (µ=1.0)",
		margin=5Plots.mm
	)
	plot!(p₂, MyKernel.wvlen_out, T_for_anly, label="O₂ + H₂O (before)", alpha=.3);
	plot!(p₂, MyKernel.band, T_conv, label="O₂ + H₂O (after convolution)", lw=2.)
end

# ╔═╡ cd8a53cc-7bdd-4a0d-a790-26b21785926b
md"""
> ##### Assemble this into a forward model, given "apparent" pressure, temperature, and vertical column density, calculate convolved transmittance.
"""

# ╔═╡ b4ee4292-9aa2-4e6d-9863-25c375c1acbb
@benchmark o2_sitp($ν_grid, $pres_O2, $temp_O2)

# ╔═╡ 9ee985fd-72fd-4ab5-bbf5-6b9aa71b9349
@benchmark o2_sitp($ν_grid[10000:end-10000], $pres_O2, $temp_O2)

# ╔═╡ 4c6dadd6-7b56-45c9-a948-1f06f214f273
"""
    compute_transmittance(P_O2, VCD_O2, T_O2, P_H2O, VCD_H2O, T_H2O, MyKernel, o2_sitp, h2o_sitp, ν_grid)

	Calculate convolved transmittance from O2 and H2O atmospheric parameters.
	
	# Arguments
	- `P_O2`: O2 pressure [hPa]
	- `VCD_O2`: O2 vertical column density [molec/cm²]
	- `T_O2`: O2 temperature [K]
	- `P_H2O`: H2O pressure [hPa]
	- `VCD_H2O`: H2O vertical column density [molec/cm²]
	- `T_H2O`: H2O temperature [K]
	- `MyKernel`: Instrument convolution kernel matrix
	- `o2_sitp`: O2 cross-section interpolator (ν, P, T)
	- `h2o_sitp`: H2O cross-section interpolator (ν, P, T)
	- `ν_grid`: Wavenumber grid [cm⁻¹]
	
	# Returns
	- `T_conv`: Convolved transmittance spectrum in wavelength space
"""
function compute_transmittance(
    P_O2,
    VCD_O2,
    T_O2,
    P_H2O,
    VCD_H2O,
    T_H2O,
    MyKernel,
    o2_sitp,
    h2o_sitp,
)
    
    # Get high-resolution cross-sections
    σ_O2 = o2_sitp(MyKernel.ν_grid, P_O2, T_O2)   # O2 cross-section [cm²/molec]
    σ_H2O = h2o_sitp(MyKernel.ν_grid, P_H2O, T_H2O)  # H2O cross-section
    
    # Calculate optical depths
    τ_O2 = σ_O2 .* VCD_O2
    τ_H2O = σ_H2O .* VCD_H2O
    τ_total = τ_O2 .+ τ_H2O
    
    # High-resolution transmittance
    T_highres = exp.(-τ_total)
    T_highres_wvlen = reverse(T_highres)
    
    # Convolve to instrument resolution
    T_conv = MyKernel.RSR_out * T_highres_wvlen
    
    return T_conv
end

# ╔═╡ 5ff92177-17e9-439e-ad8c-38217a878a83
begin
	# test
	T_conv₁ = compute_transmittance(
		pres_O2, ∑vcd_o2, temp_O2,
		pres_H2O, ∑vcd_h2o, temp_H2O,
		MyKernel, o2_sitp, h2o_sitp
	);

end

# ╔═╡ 130e11f7-fdc3-44fa-be84-c6f94f7ffa3a
begin
	vcd_dry_all = ds_summer["vcd_dry"][:,:];
	vcd_h2o_all = ds_summer["vcd_h2o"][:,:];
	vmr_h2o_all = ds_summer["vmr_h2o_var"][:,:];
	∑vcd_o2_all  = sum(vcd_dry_all * vmr_o2, dims=2);  
	∑vcd_h2o_all = sum(vcd_h2o_all, dims=2);   
	@show maximum(∑vcd_o2_all), minimum(∑vcd_o2_all);
	@show maximum(∑vcd_h2o_all), minimum(∑vcd_h2o_all);
end

# ╔═╡ 9a1216d4-9ae4-4fed-b796-13f8c8222efc
md"""
### Interactive Transmittance

P(O2): $(@bind pressures_o2 Slider(0.01:10.0:1040.01, default=1030, show_value=true)) hPa

T(O2): $(@bind temperatures_o2 Slider(150.0:5.0:330.0, default=273.15, show_value=true)) K

VCD(O2): $(@bind totVCD_o2 Slider(2.3:0.005:4.6, default=4, show_value=true)) E24

P(H2O): $(@bind pressures_h2o Slider(0.01:10.0:1040.01, default=1030, show_value=true)) hPa

T(H2O): $(@bind temperatures_h2o Slider(150.0:5.0:330.0, default=273.15, show_value=true)) K

log(VCD-H2O): $(@bind totVCD_h2o Slider(20.:.1:23.3, default=21.0, show_value=true))

"""

# ╔═╡ e41487a4-56e9-47b0-8ec6-a5f97d97e46c
begin
	# pressure and temperature
	# fixed params: VCD
	# # high res
	# τ₂      = o2_sitp(ν_grid, pressures_o2, temperatures_o2) .* ∑vcd_o2 .+ h2o_sitp(ν_grid, pressures_h2o, temperatures_h2o) .* ∑vcd_h2o;

	# convolve
	T_conv₂ = compute_transmittance(
		pressures_o2, totVCD_o2*1e24, temperatures_o2,
		pressures_h2o, 10^totVCD_h2o, temperatures_h2o,
		MyKernel, o2_sitp, h2o_sitp
	);

	# plot
	MyTitle₃ = @sprintf("P_O2=%.0f hPa | T_O2=%.0fK | P_H2O=%.0f hPa | T_H2O=%.0fK", pressures_o2, temperatures_o2, pressures_h2o, temperatures_h2o) * @sprintf("\n VCD_O2=%.1g | VCD_H2O=%.1g", totVCD_o2*1e24, 10^totVCD_h2o)
	p₃ = plot(
		size=(800, 350), 
		xlabel="wavelength [nm]",
		ylabel="transmittance (µ=1.0)",
		title = MyTitle₃,
		margin=5Plots.mm
	)
	# reference
	plot!(p₃, MyKernel.band, T_conv, label="O₂ + H₂O (just as ref.)", lw=1.)
	plot!(p₃, MyKernel.band, T_conv₂, label="O₂ + H₂O (after convolution, function)", lw=2.)
end

# ╔═╡ 72eb78e2-8c07-4198-bd8c-2102850776fc
md"""
> `compute_transmittance` as a forward model, p, T, and VCD as state vector. Given a real transmittance spectrum, how accurate can it be?
"""

# ╔═╡ 53051ea6-2500-47b0-b6e1-1d56a3c659dd
function forward_model(
		x
	)
	return compute_transmittance(
		x[1], 10^(x[3]), x[2],
		x[4], 10^(x[6]), x[5],
		MyKernel, o2_sitp, h2o_sitp
	)
end

# ╔═╡ 17375ba6-b754-4185-8324-a83b58507a1d
mean(log10.(∑vcd_o2_all))

# ╔═╡ ee23ddfa-3feb-45e4-81b4-3e3cf1c0f04c
mean(log10.(∑vcd_h2o_all))

# ╔═╡ d9ed49f2-4052-4ff6-b476-e8296cd2573b
begin
	# priori
	xₐ = [pres_O2; temp_O2; log10(∑vcd_o2); pres_H2O; temp_H2O; log10(∑vcd_h2o)];
	# Sa
	Sₐ      = I(6) .* 0.;
	Sₐ[1,1] = 1e3; Sₐ[4,4] = 1e3;
	Sₐ[2,2] = 200; Sₐ[5,5] = 200;
	Sₐ[3,3] = var(log10.(∑vcd_o2_all)); Sₐ[6,6] = var(log10.(∑vcd_h2o_all)); # log-transformed

	# Se
	n  = length(MyKernel.band);
	Se = I(n) .* 1e-6;

	# Jacobian
	K, val = PACE_SIF.Retrieval.Jacobian(xₐ, forward_model, n);

	# show Jacobian
	plot(
		MyKernel.band, 
		K,
		size=(800, 350), 
		label=["P_O2" "T_O2" "log(VCD_O2)" "P_H2O" "T_H2O" "log(VCD_H2O)"],
		title="Jacobian @ priori",
		margin=5Plots.mm
	)
end

# ╔═╡ 3fdbb03f-867f-469c-a8be-12bb7ef779ec
begin
	# get a gain matrix
	Se_inv = inv(Se);
	Sa_inv = inv(Sₐ);
	G = inv( K' * Se_inv * K + Sa_inv )K' * Se_inv;
	# AK
	A = G * K
	
	#heatmap
	function symlog(x)
	    sign(x) * log10(abs(x))
	end
	
	symlogA = symlog.(A);
	
	heatmap(symlogA,
	    xlabel = "Column",
	    ylabel = "Row",
	    title = "Averaging Kernel (symlog scale)",
	    color = :RdBu,
	    clims = (-maximum(abs.(symlogA)), maximum(abs.(symlogA)))
	)
end

# ╔═╡ ac775343-96af-444b-9cc7-36ad99b44178
begin
	# use a transmittance
	j = 1220
	y_obs = ds_summer["transmittance"][j,:];
	
	# start_time
    start_time = time();
	γ_init = 10.0;
	γ⁺     = 10.0;          
	γ⁻     = 2.0;
	nIter  = 20;
	Kₙ, y  = PACE_SIF.Retrieval.Jacobian(xₐ, forward_model, n);;
	x      = xₐ;
	iter_label   = 0;
	thr_Converge = 1e-6;
	max_runtime  = 20;
    
    # Initialize
    len    = length(y_obs);
    γ      = γ_init;

    global RMSE₀ = 1e20;
    RMSE₁ = root_mean_square(y_obs, y);
    ΔRMSE = RMSE₁ - RMSE₀;
    
    # Iteration loop
    while ( abs(ΔRMSE) > thr_Converge ) && ( iter_label < nIter )
		global iter_label, Kₙ, x, y, Sa_inv, Se_inv, γ, RMSE, RMSE₁, ΔRMSE, start_time
		
        # check runtime
        elapsed_time = time() - start_time;
        if elapsed_time > max_runtime
            @warn "LM_Iteration! exceeded max runtime" iteration=iter_label time=elapsed
            break
        end
        
        iter_label += 1;

        # Update x
        Δy = Kₙ' * Se_inv * (y_obs .- y) - Sa_inv * (x .- xₐ)
        Gₙ = (1 + γ) * Sa_inv + Kₙ' * Se_inv * Kₙ
        Δx = inv(Gₙ) * Δy
        x_trial = x .+ Δx
        
        # Evaluate at trial point
        K_trial, y_trial = PACE_SIF.Retrieval.Jacobian(x_trial, forward_model, n)
        RMSE_trial       = root_mean_square(y_obs, y_trial)
        println("RMSE try #$iter_label | $RMSE_trial")
		
        # Check if step improves the fit
        if RMSE_trial < RMSE₁
            # Accept step
            x   = x_trial
            y   = y_trial
            Kₙ  = K_trial
            # Update RMSE
            RMSE₀ = RMSE₁
            RMSE₁ = RMSE_trial
            ΔRMSE = RMSE₁ - RMSE₀
            # Decrease damping (success)
            γ /= γ⁻
            
        else
            # Reject, increase damping
            # iter_label -= 1
            γ *= γ⁺
        end
    end
end

# ╔═╡ a87ac53c-32dc-4e2a-85f0-3cd529a40a23
begin
	println(@sprintf("P_O2=%.0f hPa | T_O2=%.0fK | P_H2O=%.0f hPa | T_H2O=%.0fK", x[1], x[2], x[4], x[5]) * @sprintf("\n VCD_O2=%.1g | VCD_H2O=%.1g", 10^(x[3]), 10^(x[6])))
	
	plot(MyKernel.band, y_obs, size=(800, 350), label="true val", lw=1.)
	# plot!(MyKernel.band, val, label="priori", lw=2)
	plot!(MyKernel.band, y, label="retrieved", lw=2, linestyle=:dash)
end

# ╔═╡ 3442d690-092c-4a9c-b677-f463de97614f
plot(
	MyKernel.band, y_obs .- y, size=(800, 350), label="residual", lw=1.,
	xticks = (620:20:860, string.(620:20:860))
)

# ╔═╡ Cell order:
# ╠═d0d1237b-2c32-46ed-9060-c12ba438810f
# ╠═f3191894-fa07-4ecb-9253-6b9375fd17df
# ╠═9896ec16-ce32-4edd-8542-793a8a14f67c
# ╠═8c72db65-85a4-4dce-aca6-ab25174e561e
# ╠═f88aaa45-088e-456c-8efa-22d3215e7e10
# ╠═c92250bf-7299-4c69-becd-a356015f04bd
# ╟─b68ca8b3-df51-447d-ac8e-830f9d3e4f29
# ╠═c7d1022a-b2b2-4cbb-bc7f-fc3d6cd5b782
# ╠═d989b6f1-f079-4e45-96e1-b86103bfb38e
# ╟─cfad865f-d09b-4d98-b50a-bb6c39542a02
# ╠═1d3c114f-fdd3-4dd0-b66a-17f13e98185a
# ╟─3c626557-4612-4b18-8a9e-58bfe739b8a4
# ╠═28a29145-c387-4220-b569-761d3d797885
# ╟─cd8a53cc-7bdd-4a0d-a790-26b21785926b
# ╟─b4ee4292-9aa2-4e6d-9863-25c375c1acbb
# ╟─9ee985fd-72fd-4ab5-bbf5-6b9aa71b9349
# ╠═4c6dadd6-7b56-45c9-a948-1f06f214f273
# ╠═5ff92177-17e9-439e-ad8c-38217a878a83
# ╠═130e11f7-fdc3-44fa-be84-c6f94f7ffa3a
# ╟─9a1216d4-9ae4-4fed-b796-13f8c8222efc
# ╟─e41487a4-56e9-47b0-8ec6-a5f97d97e46c
# ╟─72eb78e2-8c07-4198-bd8c-2102850776fc
# ╠═53051ea6-2500-47b0-b6e1-1d56a3c659dd
# ╠═17375ba6-b754-4185-8324-a83b58507a1d
# ╠═ee23ddfa-3feb-45e4-81b4-3e3cf1c0f04c
# ╠═d9ed49f2-4052-4ff6-b476-e8296cd2573b
# ╟─3fdbb03f-867f-469c-a8be-12bb7ef779ec
# ╠═ac775343-96af-444b-9cc7-36ad99b44178
# ╠═a87ac53c-32dc-4e2a-85f0-3cd529a40a23
# ╟─3442d690-092c-4a9c-b677-f463de97614f
