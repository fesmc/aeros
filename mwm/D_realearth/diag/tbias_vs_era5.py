"""Zonal-mean T / RH / q bias of an aeros monthly climatology vs ERA5 on the model sigma levels.
usage: python3 tbias_vs_era5.py aeros_monthly.nc [zonal_dump.nc]"""
import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/"
def era(v, sl=False):
    f = nc.Dataset(E+("monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc" if sl else
                      "monthly-pressure-levels/era5_monthly-pressure-levels_%s_1991-2020_clim.nc")%v)
    x = f[v][:].astype(float); lat=f["latitude"][:]
    x = x.mean(axis=-1)                    # zonal mean
    if lat[0] > lat[-1]: lat=lat[::-1]; x = x[..., ::-1]
    return lat, (f["pressure_level"][:] if not sl else None), x
A = nc.Dataset(sys.argv[1]); latA=A["lat"][:]; sig=A["lev"][:]; tA=A["t"][:]; rhA=A["rh"][:]  # (month,lev,lat)
o = np.argsort(latA); latA=latA[o]; tA=tA[:,:,o]; rhA=rhA[:,:,o]
latE, pl, tE = era("t"); _,_,qE = era("q"); _,_,rE = era("r"); _,_,spE = era("sp", True)
if pl[0] < pl[-1]: pass
else: pl=pl[::-1]; tE=tE[:,::-1]; qE=qE[:,::-1]; rE=rE[:,::-1]
def onlat(x, lat): return np.array([np.interp(latA, lat, xi) for xi in x.reshape(-1, x.shape[-1])]).reshape(x.shape[:-1]+(len(latA),))
tE=onlat(tE,latE); qE=onlat(qE,latE); rE=onlat(rE,latE); spE=onlat(spE,latE)/100.
def tosig(x):  # ERA5 (month,plev,lat) -> (month,sig,lat) at p = sig*sp(month,lat), log-p interp
    out=np.zeros((12,len(sig),len(latA)))
    for m in range(12):
        for j in range(len(latA)):
            out[m,:,j]=np.interp(np.log(sig*spE[m,j]), np.log(pl), x[m,:,j])
    return out
tEs, qEs, rEs = tosig(tE), tosig(qE), tosig(rE)
# aeros q from rh and t (mixed-phase-ish: liquid sat over 273, ice below 250, blend) - approximate
def qsat(T,p):
    esl=611.2*np.exp(17.67*(T-273.15)/(T-29.65)); esi=611.2*np.exp(22.46*(T-273.15)/(T-0.53))
    w=np.clip((T-250.)/23.,0,1); es=w*esl+(1-w)*esi; return 0.622*es/(p*100.-0.378*es)
pA = sig[None,:,None]*spE[:,None,:]
qA = rhA/100.*qsat(tA,pA)
seas={"ANN":range(12),"DJF":[11,0,1],"JJA":[5,6,7]}
bands=[("60-75S",-75,-60),("30-45S",-45,-30),("15S-15N",-15,15),("30-45N",30,45),("45-60N",45,60),("60-75N",60,75)]
w=np.cos(np.deg2rad(latA))
def band(x,lo,hi):
    m=(latA>=lo)&(latA<=hi); return (x[...,m]*w[m]).sum(-1)/w[m].sum()
for nm,ms in seas.items():
    ms=list(ms)
    print(f"\n== {nm}: T bias aeros-ERA5 [K] per sigma level (rows) x band (cols); then RH aeros/ERA5 [%]; then q aeros/ERA5 [g/kg]")
    print("sigma   "+"".join(f"{b[0]:>11s}" for b in bands))
    dT=tA[ms].mean(0)-tEs[ms].mean(0)
    for k in range(len(sig)):
        print(f"{sig[k]:.3f}  "+"".join(f"{band(dT[k],lo,hi):11.1f}" for _,lo,hi in bands))
    print("RH")
    for k in range(len(sig)):
        print(f"{sig[k]:.3f}  "+"".join(f"{band(rhA[ms].mean(0)[k],lo,hi):5.0f}/{band(rEs[ms].mean(0)[k],lo,hi):<5.0f}" for _,lo,hi in bands))
    print("q")
    for k in range(len(sig)):
        print(f"{sig[k]:.3f}  "+"".join(f"{band(qA[ms].mean(0)[k],lo,hi)*1e3:5.2f}/{band(qEs[ms].mean(0)[k],lo,hi)*1e3:<5.2f}" for _,lo,hi in bands))
    # lapse rate dT/dlnp between adjacent levels [K per ln p]
    print("dT/dlnp (K per unit ln p), aeros/ERA5, between levels k,k+1")
    for k in range(len(sig)-1):
        dl=np.log(sig[k+1]/sig[k])
        ga=(tA[ms].mean(0)[k+1]-tA[ms].mean(0)[k])/dl; ge=(tEs[ms].mean(0)[k+1]-tEs[ms].mean(0)[k])/dl
        print(f"{sig[k]:.2f}-{sig[k+1]:.2f} "+"".join(f"{band(ga,lo,hi):5.1f}/{band(ge,lo,hi):<5.1f}" for _,lo,hi in bands))
if len(sys.argv)>2:
    Z=nc.Dataset(sys.argv[2]); latZ=Z["lat"][:]; o=np.argsort(latZ)
    print("\n== per-process heating [K/day], 2nd-half time mean, zonal mean (dump), plus omega [hPa/day]")
    for v in ["q_cnv","q_cnd","q_rad","q_surf","q_net","omega","t","q","rh","pfull"]:
        x=Z[v][:]; x = x if x.shape[0]==len(latZ) else x.T   # (lat,lev)
        x=x[o]; print(v); print("sigma   "+"".join(f"{b[0]:>11s}" for b in bands))
        for k in range(len(sig)):
            print(f"{sig[k]:.3f}  "+"".join(f"{band(x[:,k],lo,hi):11.2f}" for _,lo,hi in bands))
