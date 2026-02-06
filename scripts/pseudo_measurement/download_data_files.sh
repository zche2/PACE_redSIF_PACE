#!/bin/bash
# Script to download required data files for Create_pseudo_v2_visualize_xSecFit_standalone.jl
# Run this from /Users/cfranken/data directory
# Usage: bash download_data_files.sh <username>@<server>

if [ $# -eq 0 ]; then
    echo "Usage: $0 <ssh_server>"
    echo "Example: $0 zhe2@servername.edu"
    exit 1
fi

SERVER=$1
REMOTE_BASE="/home/zhe2/data/MyProjects/PACE_redSIF_PACE"
LOCAL_BASE="/Users/cfranken/data"

echo "Downloading required data files from $SERVER"
echo "Local directory: $LOCAL_BASE"
echo ""

# Create necessary directories
echo "Creating local directory structure..."
mkdir -p "$LOCAL_BASE/reference_spectra"
mkdir -p "$LOCAL_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01"
mkdir -p "$LOCAL_BASE/PACE_OCI"

echo ""
echo "====================================================================="
echo "File 1/6: Winter transmittance (optional but recommended)"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc" "$LOCAL_BASE/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "File 2/6: SNR lookup table"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt" "$LOCAL_BASE/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "File 3/6: SIF singular vectors"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/reference_spectra/SIF_singular_vector.jld2" "$LOCAL_BASE/reference_spectra/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "File 4/6: O2 cross section LUT"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2" "$LOCAL_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "File 5/6: H2O cross section LUT"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2" "$LOCAL_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "File 6/6: Instrument kernel"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/KernelInstrument.jld2" "$LOCAL_BASE/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "Optional: Metadata log file"
echo "====================================================================="
scp "$SERVER:$REMOTE_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log" "$LOCAL_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/" || echo "Failed - continuing..."

echo ""
echo "====================================================================="
echo "Download complete!"
echo "====================================================================="
echo ""
echo "Checking downloaded files:"
echo ""

# Check what we have
files=(
    "$LOCAL_BASE/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
    "$LOCAL_BASE/transmittance_winter_FineWvResModel_FullRange_Aug01.nc"
    "$LOCAL_BASE/sample_granule_20240830T131442_new_chl.nc"
    "$LOCAL_BASE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt"
    "$LOCAL_BASE/reference_spectra/SIF_singular_vector.jld2"
    "$LOCAL_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2"
    "$LOCAL_BASE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2"
    "$LOCAL_BASE/KernelInstrument.jld2"
)

for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        size=$(du -h "$file" | cut -f1)
        echo "✓ $size  $(basename $file)"
    else
        echo "✗ MISSING: $(basename $file)"
    fi
done

echo ""
echo "You can now run the standalone script!"
