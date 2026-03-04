### Descriptions of files

##### <b>Create pseudo measurements</b> 
- v1:

    <b> Retrieval results in JLD2 files:</b>
    - v1. use two different T1 and T2 

    - v1.1 only choose one transmittance for constructing pseudo measurements.

    - v1.2 only choose one transmittance for constructing pseudo measurements,nPoly=6, rank=20
    
    ----- 

    above do not see any improvement, SIF consistently being underestimated by ~25%

    - v1.3 
        nPoly=10, rank=15  
        Change penalty function for T₁ and T₂ conversion:   
        smooth_x = 10. / (1 + exp( -x[px.nPoly+px.nPC+2]) ) + 1.  

    -----

    - v2.1 - v2.3 | cross section fits (p, T, vcd in state vector): are time-consuming. T1 and T2 are essentially unrelated
    
    - v2.4 | incorporate viewing geometry, true sza and vza are known and remain the same both when constructing pseudo-measurements and doing the retrieval. => since µ-scaling is applied after convolution, some non-linearities might be introduced. Given that the ensemble of transmittance is constructed with varying viewing angles, should not be re-accounted for ...

