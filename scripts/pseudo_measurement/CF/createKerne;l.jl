filename = "/Users/cfranken/data/PACE_OCI_RSRs.nc";
pace = Dataset(filename, "r");

wavlen = pace["wavelength"][:];
RSR    = pace["RSR"];
band   = pace["bands"];