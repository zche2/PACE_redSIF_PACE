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