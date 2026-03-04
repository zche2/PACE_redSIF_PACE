### Make cross section model out of gridded pressure and temperature
> #####  About the Model
The most up-to-date script should be `Automate_MakeModel.jl`

Save the model in jld2 files, see directory: `../data/MyProjects/PACE_redSIF_PACE/interp_xSection`

with meta data of temperature, pressure, and wavenumber grids used

Should have a .log file documenting the configuration

if not, don't trust and use it in successional analysis :\

p.s., something has messed up, .log files in this folder should not be distribute to the data folder!

> ##### File Used in Subsequent Analysis

- Absorptance of O$_2$ and H$_2$O stored in `../data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01`