# red SIF retrieval for OCI (PACE)
try atm. transmission calculation for PACE red SIF retrieval

Known issue:

- diagonal prior σ on each PC coefficient: OCI is now setting the prior σ using "singular value", i.e., σ_k = S[k] / √n_profiles; TROPOMI is now using "loading variance", i.e., σ_k = S[k] × std(Vt[k,:])
