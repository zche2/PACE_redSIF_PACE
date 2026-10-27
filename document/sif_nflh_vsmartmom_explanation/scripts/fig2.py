import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, sys
from matplotlib.lines import Line2D
S=sys.argv[1]; o=pd.read_csv(S+"/slopes.csv"); w=pd.read_csv(S+"/compare_all.csv")
M=[("std667_748","nFLH·cosθs, baseline 667/748",'C0','-','o'),("wide630_800","nFLH·cosθs, baseline 630/800",'C1','-','o'),
   ("nPC15_npoly3|RT","SVD npoly3 640–756, RT solar",'C2','-','o'),("nPC15_npoly3|L1B","SVD npoly3 640–756, L1B solar",'olive','--','^'),
   ("nPC15_npoly3_O2A|RT","SVD npoly3_O2A 640–780, RT solar",'C3','-','o'),("nPC15_npoly3_O2A|L1B","SVD npoly3_O2A 640–780, L1B solar",'C4','--','^')]
isn=lambda m: m.startswith(("std","wide"))
cosf=lambda g,m: np.cos(np.radians(g.sza.values)) if isn(m) else 1.0
lab=[f"aer×{a:g}\nray×{r:g}" for a in (0,1,5) for r in (0,0.1,1)]
plt.rcParams.update({"font.size":11})
# ---------- Figure 1: summary vs scenario ----------
fig,ax=plt.subplots(3,2,figsize=(15,13),sharex=True)
for j,sza in enumerate((30.0,60.5)):
    for m,l,c,ls,mk in M:
        for vza,al,lw in ((15,1,1.8),(60,.45,1.0)):
            if isn(m) and vza==60: continue
            g=o[(o.sza==sza)&(o.vza==vza)].sort_values(["aer","ray"])
            ax[0,j].plot(range(9),g[m+"_slope"],ls=ls,marker=mk if vza==15 else 's',color=c,ms=5,alpha=al,lw=lw)
            if not isn(m): ax[1,j].plot(range(9),g[m+"_icpt"],ls=ls,marker=mk if vza==15 else 's',color=c,ms=5,alpha=al,lw=lw)
        g=w[(w.sza==sza)&(w.vza==30)&np.isclose(w.sif_678,0.3)].sort_values(["aerosol_scale","rayleigh_scale"])
        ax[2,j].plot(range(9),g[m]*cosf(g,m)/g.sif_wl_678,ls=ls,marker=mk,color=c,ms=5)
    ax[0,j].set_title(f"SZA = {sza}°",fontsize=13); ax[0,j].axhline(1,color='k',lw=.6); ax[1,j].axhline(0,color='k',lw=.6); ax[2,j].axhline(1,color='k',lw=.6)
    ax[1,j].set_ylim(-0.17,0.2); ax[2,j].set_xticks(range(9)); ax[2,j].set_xticklabels(lab,fontsize=9)
    for i,t in enumerate(("(a) slope","(b) SVD intercept","(c) ratio at SIF₆₇₈=0.3, VZA 30")): ax[i,j].text(.01,.97,t,transform=ax[i,j].transAxes,va="top",fontsize=10,bbox=dict(fc="w",ec="none",alpha=.8))
ax[0,0].set_ylabel("slope vs water-leaving SIF(678)"); ax[1,0].set_ylabel("intercept [mW m⁻² nm⁻¹ sr⁻¹]"); ax[2,0].set_ylabel("estimate / water-leaving SIF")
h=[Line2D([],[],color=c,ls=ls,marker=mk,label=l) for m,l,c,ls,mk in M]+[Line2D([],[],color='gray',marker='o',label='VZA 15 (solid)'),Line2D([],[],color='gray',marker='s',alpha=.45,label='VZA 60 (faded, panels a–b)')]
fig.legend(handles=h,loc="lower center",ncol=4,fontsize=10,frameon=False)
fig.suptitle("Mock nFLH (surface Rrs, F0 = 1453.2 W m⁻² µm⁻¹, vSmartMOM) vs SVD on noise-free TOA L_toa",fontsize=13)
fig.tight_layout(rect=(0,.06,1,.97)); fig.savefig(S+"/nflh_vs_svd.png",dpi=130); plt.close(fig)
# ---------- Figures 2/3: scatter + residual, one per SZA ----------
sc=[(a,r) for a in (0,1,5) for r in (0,0.1,1)]
fam={0:plt.get_cmap("Blues"),1:plt.get_cmap("Greens"),5:plt.get_cmap("Reds")}; shade={0:.45,0.1:.7,1:.95}
cm=lambda k: fam[sc[k][0]](shade[sc[k][1]])
vm={15:'o',30:'^',60:'s'}
for sza in (30.0,60.5):
    MS=[("std667_748|raw","nFLH (no cosθs)\nbaseline 667/748"),("wide630_800|raw","nFLH (no cosθs)\nbaseline 630/800")]+[(m,l.replace(", ","\n")) for m,l,*_ in M]
    fig,ax=plt.subplots(2,8,figsize=(33,11),gridspec_kw=dict(height_ratios=[1.15,1]))
    g0=w[w.sza==sza]
    for j,(mm,l) in enumerate(MS):
        raw=mm.endswith("|raw"); m=mm.split("|raw")[0]
        for k,(a,r) in enumerate(sc):
            for vza,mk in vm.items():
                if isn(m) and vza!=15: continue
                g=g0[np.isclose(g0.aerosol_scale,a)&np.isclose(g0.rayleigh_scale,r)&(g0.vza==vza)]
                y=g[m].values*(1.0 if raw else cosf(g,m)); x=g.sif_wl_678.values
                ax[0,j].scatter(x,y,s=34,marker=mk,facecolor=cm(k),edgecolor="k",linewidth=.3,alpha=.9)
                ax[1,j].plot(x,y-x,marker=mk,color=cm(k),ms=5,lw=.6,alpha=.85)
        top=2.3 if (raw and sza>45) else (1.35 if raw else 1.15); lim=[-0.1,top]; ax[0,j].plot(lim,lim,'k',lw=.7); ax[0,j].set_xlim(-0.02,1.05); ax[0,j].set_ylim(*lim); ax[0,j].grid(alpha=.3); ax[1,j].grid(alpha=.3)
        ax[1,j].axhline(0,color='k',lw=.7); ax[1,j].set_ylim(-0.32,1.25 if (raw and sza>45) else 0.62)
        if raw: [ax[i,j].set_facecolor('#fff6e6') for i in (0,1)]
        ax[0,j].set_title(l.replace('nFLH·cosθs','nFLH × cosθs'),fontsize=11); ax[1,j].set_xlabel("water-leaving SIF(678) [mW m⁻² nm⁻¹ sr⁻¹]",fontsize=9)
    ax[0,0].set_ylabel("estimate (nFLH × cosθs, or SVD SIF₆₇₈)"); ax[1,0].set_ylabel("estimate − truth")
    h=[Line2D([],[],ls='',marker='o',color=cm(k),label=f"aer×{a:g}, ray×{r:g}") for k,(a,r) in enumerate(sc)]+[Line2D([],[],ls='',marker=mk,color='k',label=f"VZA {v}") for v,mk in vm.items()]
    fig.legend(handles=h,loc="lower center",ncol=12,fontsize=10,frameon=False)
    fig.suptitle(f"SZA = {sza}°: estimate vs water-leaving SIF (nFLH from surface Rrs has no VZA dependence → VZA 15 only; shaded = raw nFLH, note different y-range)",fontsize=14)
    fig.tight_layout(rect=(0,.07,1,.96)); fig.savefig(S+f"/scatter_sza{str(sza).replace('.','p')}.png",dpi=110); plt.close(fig)
