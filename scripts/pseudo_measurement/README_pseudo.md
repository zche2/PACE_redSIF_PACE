### Descriptions of files

##### <b>Create pseudo measurements</b> 
- v1:

    <b> Retrieval results in JLD2 files:</b>
    v1. use two different T1 and T2 

    v1.1 only choose one transmittance for constructing pseudo measurements.

    v1.2 only choose one transmittance for constructing pseudo measurements, nPoly=6, rank=20
    
    ----- 

    above do not see any improvement, SIF consistently being underestimated by ~25%

    - v1.3 
        nPoly=10, rank=15

        Change penalty function for T₁ and T₂ conversion: 

        smooth_x = 10. / (1 + exp( -x[px.nPoly+px.nPC+2]) ) + 1.  
    - v1.4
        order=4, nPoly=10, rank=15, new penalty function, no considering viewing geometry.

- v3: 
    Enable fitting more PCs of SIF shape

    - v3.0 nSIF = 1 as reference retrieval
        v3.0.1 zero noise: no measurement noise when creating pseudo data

    - v3.1 nSIF = 2

        v3.1.1 Set SNR degradation bands  
        v3.1.2 Wavelength 670-780nm (try to avoid Water band as much as possible).  
        v3.1.3 Sₐ updated.  
        
    - v3.2 Try log SVD to transmittance spectra and apply the new penalty function.  
        v3.2.1 log SVD with updated Sₐ.  
        v3.2.1_ReducedNoise reduce standard deviation by a factor of 10 (equivalent to pixel aggregation).  
        v3.2.1_ZeroNoise just zero noise :\  


    - v3.3 NMF to optical thickness (add if_log to NMF decomposition), other params are kept the same as v3.1.3