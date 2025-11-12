### Descriptions of files

##### <b>Forward models</b> 
- v1: Some preliminary try-out. Have some evaluation of **numerical efficiency**
- v2: 
    Start from simple linear model with reflectance only, add transmittance (scaled to maximum value=1)

    correction to the Gauss-Newton iteration.
    - v2.1:
        Add Gaussian-shaped SIF (bug detected: Wrong Legendre Polynomials)
    - v2.2:
        2025-Aug-26: Try "baseline fitting": retrieval outside of SIF range (similar idea as nFLH)

        2025-Aug-28: Correction to the scaling of transmittance spectrum (Now use baseline band and force it to 1)
    - v2.3: 
        2025-Aug-29: Add a gaussian-shaped chlorophyll absorption feature (wasn't manage to do that, and given a known SIF shape, it is no more necessary)
- v3:
    Allow fitting one-way and two-way transmittance separately

    Adopt the full spectral shape of SIF (interpolated to OCI res.)

    Update priori covariance terms for PCs

- v4:
    - v4.1: Non-negative matrix factorization fit 
    - v4.2: Add chlorophyll absorption (bad)

- v5: Cross section fit (see demo)

##### <b>Retrieval</b>
Apply the up-to-date forward model
