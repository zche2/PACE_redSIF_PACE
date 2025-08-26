### Descriptions of files
- v1: Some preliminary try-out. Have some evaluation of **numerical efficiency**
- v2: 
    Start from simple linear model with reflectance only, add transmittance (scaled to maximum value=1)
    correction to the Gauss-Newton iteration.
    - v2.1:
        Add Gaussian-shaped SIF (bug detected: Wrong Legendre Polynomials)
    - v2.2:
        Try "baseline fitting": retrieval outside of SIF range (similar idea as nFLH)
