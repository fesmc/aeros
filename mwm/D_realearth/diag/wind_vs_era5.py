import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/"
def era(v, sl=False):
    f = nc.Dataset(E+("monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc" if sl else "monthly-pressure-levels/era5_monthly-pressure-levels_%s_1991-2020_clim.nc")%v)
    x=f[v][:].astype(float); lat=f["latitude"][:]; 
    if lat[0]>lat[-1]: lat=lat[::-1]; x=x[...,::-1,:]
    return lat,(f["pressure_level"][:] if not sl else None),x
latE,pl,uE=era("u"); _,_,vE=era("v"); _,_,spE=era("sp",True); _,_,lsm=era("lsm",True); _,_,t2m=era("t2m",True); _,_,skt=era("skt",True)
spE=spE/100.
A=nc.Dataset(sys.argv[1]); latA=A["lat"][:]; sig=A["lev"][:]; o=np.argsort(latA); latA=latA[o]
uA=A["u"][:,:,o]; vA=A["v"][:,:,o]; vlow=A["vlow"][:][:,o,:].mean(-1); ks=len(sig)-1
# ERA5 |V| at p = sig_low*sp per point, then zonal mean; and zonal-mean u
def at_sig(x, s):
    out=np.zeros(x.shape[:1]+x.shape[2:])
    for m in range(12):
        for j in range(x.shape[2]):
            for i in range(x.shape[3]):
                out[m,j,i]=np.interp(np.log(s*spE[m,j,i]), np.log(pl[::-1]) if pl[0]>pl[-1] else np.log(pl), (x[m,::-1,j,i] if pl[0]>pl[-1] else x[m,:,j,i]))
    return out
uEs=at_sig(uE,sig[ks]); vEs=at_sig(vE,sig[ks]); spd=np.sqrt(uEs**2+vEs**2)
ocean=(lsm[0]<0.5)
def zm(x,mask=None):
    if mask is None: return x.mean(-1)
    return np.array([[np.mean(x[m,j][mask[j]]) if mask[j].any() else np.nan for j in range(x.shape[1])] for m in range(12)])
w=np.cos(np.deg2rad(latA))
bands=[("45-60S",-60,-45),("30-45S",-45,-30),("15-30S",-30,-15),("15S-15N",-15,15),("15-30N",15,30),("30-45N",30,45),("45-60N",45,60)]
def band(x,lat,lo,hi):
    xi=np.array([np.interp(latA,lat,r) for r in x]); m=(latA>=lo)&(latA<=hi); return np.nansum(xi[:,m]*w[m],-1).mean()/w[m].sum()
print("annual, lowest level (sigma %.3f): aeros zonal-mean u / ERA5 u ; aeros zonal-mean |V|(vlow) / ERA5 |V| of monthly means ; ERA5 t2m-skt over ocean"%sig[ks])
print("band      uA/uE      vlowA/|V|E   t2m-skt(oc)")
for nm,lo,hi in bands:
    print(f"{nm:8s} {band(uA[:,ks,:],latA,lo,hi):5.1f}/{band(zm(uEs),latE,lo,hi):<5.1f}  {band(vlow,latA,lo,hi):5.1f}/{band(zm(spd),latE,lo,hi):<5.1f}  {band(zm(t2m-skt,ocean),latE,lo,hi):6.2f}")
