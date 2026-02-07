# =========================================================================================================
# This script helps to download PACE OCI L1B V3 data
# Direct data access from: https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-1B/PACE-OCI/
# =========================================================================================================

import earthaccess

# Authenticate with Earthdata
auth = earthaccess.login(strategy="netrc")  # Uses .netrc for credentials

# Specify the date (YYYYMMDD format) and output directory
start_date = "20250927"  # Start date
end_date   = "20250928"    # End date
output_dir = "/home/zhe2/data/PACE/L1B_V3/"
output_dir_aop = "/home/zhe2/data/PACE/L2_AOP_V3.1/"
output_dir_bgc = "/home/zhe2/data/PACE/L2_BGC_V3.1/"

# Search for PACE OCI L1B V3 files
results = earthaccess.search_data(
    short_name="PACE_OCI_L1B_SCI",
    version="3",
    temporal=(f"{start_date[:4]}-{start_date[4:6]}-{start_date[6:8]}", 
              f"{end_date[:4]}-{end_date[4:6]}-{end_date[6:8]}")
)

# download all results that match the search criteria
earthaccess.download(results, local_path=output_dir)

# search for PACE OCI V3.1 L2 AOP and BGC
results_aop = earthaccess.search_data(
    short_name="PACE_OCI_L2_AOP",
    version="3.1",
    temporal=(f"{start_date[:4]}-{start_date[4:6]}-{start_date[6:8]}", 
              f"{end_date[:4]}-{end_date[4:6]}-{end_date[6:8]}")
)
results_bgc = earthaccess.search_data(
    short_name="PACE_OCI_L2_BGC",
    version="3.1",
    temporal=(f"{start_date[:4]}-{start_date[4:6]}-{start_date[6:8]}", 
              f"{end_date[:4]}-{end_date[4:6]}-{end_date[6:8]}")
)
# download
earthaccess.download(results_aop, local_path=output_dir_aop)
earthaccess.download(results_bgc, local_path=output_dir_bgc)
