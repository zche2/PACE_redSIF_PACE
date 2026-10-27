# Usage: python nflh_compare.py <outdir>  (needs <outdir>/svd_retrievals.csv; writes nflh_mock.csv, compare_all.csv)
import numpy as np, netCDF4 as nc, pandas as pd, sys
S=sys.argv[1]; D="/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction/output_rayleigh_scaling/"
F={30.0:"rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sunlight_on_AerColumn_100_grid_surface.nc",
   60.5:"rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sza60p5_sunlight_on_AerColumn_100_grid_surface.nc"}
import os
sb=np.genfromtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),"rt_solar_oci.csv"),delimiter=",")  # vSmartMOM solar, OCI-convolved
sb=sb[~np.isnan(sb).any(1)]
rows=[]
for sza,f in F.items():
    d=nc.Dataset(D+f); lam=d["wavelength"][:].data
    Rrs=np.transpose(d["Rrs"][:].data); sw=np.transpose(d["sif_waterleaving"][:].data); Ed=d["Ed"][:].data
    sif=d["sif_678"][:].data; vza=d["vza"][:].data; rs=d["rayleigh_scale"][:].data; ae=d["aerosol_scale"][:].data
    dims=d["Rrs"].dimensions
    if sza==30.0: print("Rrs dims",dims, "shape",Rrs.shape)
    ib=lambda L: int(np.argmin(abs(lam-L)))
    ip=ib(678.2); F0=np.interp(lam[ip],sb[:,0],sb[:,1])
    if sza==30.0:
        print(f"peak band {lam[ip]:.2f} nm, F0(vSmartMOM, OCI-conv) = {F0:.1f} W m-2 um-1")
        for L in (630,667,748,800): print(f"  band {L} -> {lam[ib(L)]:.2f} nm, water-leaving SIF/peak = {sw[ib(L),-1]/sw[ip,-1]:.3e}")
    for tag,(a,b) in {"std667_748":(667.0,748.0),"wide630_800":(630.0,800.0)}.items():
        i1,i3=ib(a),ib(b); l1,l2,l3=lam[i1],lam[ip],lam[i3]
        # Rrs layout per dims: (band, sif, view, rayleigh, aerosol) in Julia == reversed in netCDF4
        R=Rrs
        nflh=F0*(R[ip]-(R[i1]+(R[i3]-R[i1])*(l2-l1)/(l3-l1)))
        S_=sw
        lh_true=S_[ip]-(S_[i1]+(S_[i3]-S_[i1])*(l2-l1)/(l3-l1))  # line height of water-leaving SIF itself
        for k in range(len(sif)):
            for iv in range(len(vza)):
                for ir in range(len(rs)):
                    for ia in range(len(ae)):
                        rows.append(dict(sza=sza,aerosol_scale=ae[ia],rayleigh_scale=rs[ir],vza=vza[iv],sif_678=sif[k],
                            baseline=tag,nflh=nflh[k,iv,ir,ia],sif_wl_678=S_[ip,k],sif_lh_true=lh_true[k],F0=F0))
    d.close()
n=pd.DataFrame(rows); n.to_csv(S+"/nflh_mock.csv",index=False)
v=pd.read_csv(S+"/svd_retrievals.csv")
print("SVD status counts:\n",v.groupby(["config","solar"]).status.value_counts().unstack(fill_value=0))
key=["sza","aerosol_scale","rayleigh_scale","vza","sif_678"]
for c in ("sza","aerosol_scale","rayleigh_scale","vza","sif_678"): v[c]=v[c].round(6); n[c]=n[c].round(6)
w=n.pivot_table(index=key,columns="baseline",values="nflh").reset_index()
w=w.merge(n[n.baseline=="std667_748"][key+["sif_wl_678","sif_lh_true"]],on=key)
for (cfg,sol),g in v.groupby(["config","solar"]):
    w=w.merge(g[key+["sif_ret"]].rename(columns={"sif_ret":f"{cfg}|{sol}"}),on=key)
w.to_csv(S+"/compare_all.csv",index=False)
