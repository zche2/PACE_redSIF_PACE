# %% import
import Pkg
Pkg.activate("/home/zhe2/FraLab/ESE156.jl")

# %% instantiate - download all dependencies for the very first time
# see explanations: https://pkgdocs.julialang.org/v1/environments/
Pkg.instantiate()

# %% dependencies
using Markdown, InteractiveUtils, Plots
using PlutoUI
using ESE156         # This package loads most of our utilities
using vSmartMOM,  vSmartMOM.Absorption

# %% CO2 absorption
co2_par = Absorption.read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6214.4, ν_max=6214.8);

# ╔═╡ e03f3b09-4a22-4708-ac08-7532420c22a1
line_voigt   = make_hitran_model(co2_par, Voigt(), architecture=CPU())

# ╔═╡ 4ad454c3-eb72-4a31-82e9-951df1dfbebf
# Specify our wavenumber grid
#ν_min=6214.4, ν_max=6214.8
ν = 6213:0.001:6216.5;

# ╔═╡ e7f5682f-90e4-453f-a7b9-cb03796bfcf0
plotly()

# ╔═╡ f17f8699-45dc-490a-a1d8-99aad3e4a016
@bind p Slider(1.0:5:6500.0, default=1000.0)

# ╔═╡ 04f69164-b312-467b-9911-6d648b771efb
@bind T Slider(180.0:5:550.0, default=290.0)

# ╔═╡ 664161fa-d2dc-48ec-8f09-989957306582
σ = absorption_cross_section(line_voigt, ν, p    , T);

# ╔═╡ bb7f28b3-a54e-4378-8f36-2c5546a9857d
plot(ν, σ)

# ╔═╡ a72715f1-f580-4816-bf9a-410818af4bb3
@show p, T