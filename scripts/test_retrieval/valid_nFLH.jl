
using NCDatasets
# ===========================================
# path to OCI data
# ===========================================
granule_name = "sample_granule_20240830T131442_new_chl"
path_oci = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/$(granule_name).nc"

nflh_threshold = 0.2  # Threshold for valid NFLH

println("Loading OCI data from $path_oci ...")

# ===========================================

# load data
ds = Dataset(path_oci)
# find pixels with valid NFLH
nflh = ds["nflh"][:]
valid_indices = findall(coalesce.(nflh .> 0.2, false))