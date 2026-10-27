import numpy as np, pandas as pd, sys
S=sys.argv[1]; w=pd.read_csv(S+"/compare_all.csv")
meth=["std667_748","wide630_800","nPC15_npoly3|L1B","nPC15_npoly3|RT","nPC15_npoly3_O2A|L1B","nPC15_npoly3_O2A|RT"]
w["cos"]=np.cos(np.radians(w.sza))
out=[]
for k,g in w.groupby(["sza","aerosol_scale","rayleigh_scale","vza"]):
    r=dict(zip(["sza","aer","ray","vza"],k)); x=g.sif_wl_678.values
    for m in meth:
        y=g[m].values*(g.cos.values if m.startswith(("std","wide")) else 1)
        p=np.polyfit(x,y,1); r[m+"_slope"]=p[0]; r[m+"_icpt"]=p[1]
    out.append(r)
o=pd.DataFrame(out); o.to_csv(S+"/slopes.csv",index=False)
pd.set_option("display.width",250); pd.set_option("display.precision",3)
print("nFLH*cosSZA slope vs water-leaving SIF(678) (intercept 0); no VZA dependence -> vza=15 rows")
print(o[o.vza==15][["sza","aer","ray","std667_748_slope","wide630_800_slope"]].to_string(index=False))
for m in meth[2:]:
    print("\n",m,": slope / intercept of retrieved TOA SIF(678.2) vs water-leaving SIF(678)")
    t=o.pivot_table(index=["sza","aer","ray"],columns="vza",values=[m+"_slope",m+"_icpt"])
    print(t.to_string())
