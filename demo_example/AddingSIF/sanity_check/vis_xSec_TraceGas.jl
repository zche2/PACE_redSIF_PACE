using DelimitedFiles
using Plots

# Path to the O3 cross section data file
o3_file = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/GOME-2_FM3_O3_Temp_cross-section_V5_0.txt"

# Read the O3 data: assuming two columns (wavelength [nm], cross section [cm^2/molecule])
o3_data = readdlm(o3_file, skipstart=37)

# column 2 is wavelength in nm
wavelength_o3 = o3_data[:, 2]
# column 9 cross section in cm^2/molecule at 273 k
xSec_o3_273 = o3_data[:, 9]
# column 11 cross section in cm^2/molecule at 293 k
xSec_o3_293 = o3_data[:, 11]

# set wavelength bounds
wavelength_min = 640
wavelength_max = 756
wavelength_o3_clip = wavelength_o3[wavelength_o3 .>= wavelength_min .&& wavelength_o3 .<= wavelength_max]
xSec_o3_273 = xSec_o3_273[wavelength_o3 .>= wavelength_min .&& wavelength_o3 .<= wavelength_max]
xSec_o3_293 = xSec_o3_293[wavelength_o3 .>= wavelength_min .&& wavelength_o3 .<= wavelength_max]

# Plot O3 cross-section spectrum
plot(
    size=(800, 400),
    xlabel="Wavelength [nm]",
    ylabel="O₃ Cross Section [cm²/molecule]",
    title="O₃ Absorption Cross Section",
    legend = true,
    yscale = :log10,
    linewidth = 2,
    dpi = 300,
    grid = true,
)
plot!(wavelength_o3_clip, xSec_o3_273, label="273 K")
plot!(wavelength_o3_clip, xSec_o3_293, label="293 K")



# file 2: O3 cross section from SerdyuchenkoGorshelev5digits_latest.dat
o3_file_2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SerdyuchenkoGorshelev5digits_latest.dat"
o3_data_2 = readdlm(o3_file_2, skipstart=45)
wavelength_o3_2 = o3_data_2[:, 1]
xSec_o3_2 = o3_data_2[:, 6]

# set wavelength bounds
wavelength_min = 640
wavelength_max = 756
wavelength_o3_clip_2 = wavelength_o3_2[wavelength_o3_2 .>= wavelength_min .&& wavelength_o3_2 .<= wavelength_max]
xSec_o3_2 = xSec_o3_2[wavelength_o3_2 .>= wavelength_min .&& wavelength_o3_2 .<= wavelength_max]

# Plot O3 cross-section spectrum
plot(
    size=(800, 400),
    xlabel="Wavelength [nm]",
    ylabel="O₃ Cross Section [cm²/molecule]",
    title="O₃ Absorption Cross Section",
    legend = false,
    yscale = :log10,
    linewidth = 2,
    dpi = 300,
    grid = true,
)
plot!(wavelength_o3_clip_2, xSec_o3_2, label="high res 0.01 nm")