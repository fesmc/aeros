"""NH 40-70N (zonal, and land-band TOA): JJA RH / cloud-fraction profiles vs ERA5, TOA reflected SW over land."""
import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/"
def era(v, sl=False):
    f = nc.Dataset(E+("monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc" if sl else
                      "monthly-pressure-levels/era5_monthly-pressure-levels_%s_1991-2020_clim.nc")%v)
    x=f[v][:].astype(float); lat=f["latitude"][:]; lon=f["longitude"][:]
    if lat[0]>lat[-1]: lat=lat[::-1]; x=x[...,::-1,:]
    return lat, lon, (f["pressure_level"][:] if not sl else None), x
A=nc.Dataset(sys.argv[1]); la=A["lat"][:]; lo=A["lon"][:]; o=np.argsort(la); la=la[o]; sig=A["lev"][:]
lat,lon,_,lsm=era("lsm",True); land=lsm[0]>0.5
w=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,len(lon)))
L=np.array([[lsm[0][np.argmin(abs(lat-y)),np.argmin(abs(lon-x%360))]>0.5 for x in lo] for y in la])
wa=np.cos(np.deg2rad(la))[:,None]*np.ones((1,len(lo)))
months=[5,6,7]  # JJA
for name,bandE,bandA in [("land",(lat[:,None]>=40)&(lat[:,None]<=70)&land,(la[:,None]>=40)&(la[:,None]<=70)&L),
                         ("ocean",(lat[:,None]>=40)&(lat[:,None]<=70)&~land,(la[:,None]>=40)&(la[:,None]<=70)&~L)]:
    def mE(x): return np.mean([(x[k]*w*bandE).sum()/(w*bandE).sum() for k in months])
    def mA(v): x=A[v][:][:,o,:]; return np.mean([(x[k]*wa*bandA).sum()/(wa*bandA).sum() for k in months])
    tisr=era("tisr",True)[3]/86400; tsr=era("tsr",True)[3]/86400; tsrc=era("tsrc",True)[3]/86400
    ttr=era("ttr",True)[3]/86400; ttrc=era("ttrc",True)[3]/86400; tcc=era("tcc",True)[3]
    print(f"JJA 40-70N {name}: TOA refl SW model {mA('swup'):.0f} / ERA5 {mE(tisr-tsr):.0f} (clear {mE(tisr-tsrc):.0f}, cloud {mE(tsrc-tsr):.0f});  OLR {mA('olr'):.0f}/{mE(-ttr):.0f} (clear {mA('olr_cs'):.0f}/{mE(-ttrc):.0f});  cover {mA('cover'):.2f}/{mE(tcc):.2f}")
# zonal profiles 40-70N (model zonal fields; ERA5 zonal), JJA
_,_,pl,rE=era("r"); _,_,_,spE=era("sp",True)
try: _,_,_,ccE=era("cc")
except Exception: ccE=None
if pl[0]>pl[-1]: pl=pl[::-1]; rE=rE[:,::-1]; ccE=None if ccE is None else ccE[:,::-1]
def zb(x,lt): 
    wz=np.cos(np.deg2rad(lt)); b=(lt>=40)&(lt<=70); return (x[...,b]*wz[b]).sum(-1)/wz[b].sum()
rEz=zb(rE.mean(-1),lat); spz=zb(spE.mean(-1),lat)/100.; ccz=None if ccE is None else zb(ccE.mean(-1),lat)
rA=zb(A["rh"][:][:,:,o],la); cfA=zb(A["cf"][:][:,:,o],la)
print("JJA 40-70N zonal:  sigma   RH model / ERA5   cf model / ERA5(cc)")
for k in range(len(sig)):
    p=sig[k]*np.mean(spz[months])
    re=np.mean([np.interp(np.log(p),np.log(pl),rEz[m]) for m in months]); ra=np.mean([rA[m,k] for m in months])
    ce=None if ccz is None else np.mean([np.interp(np.log(p),np.log(pl),ccz[m]) for m in months]); ca=np.mean([cfA[m,k] for m in months])
    print(f"  {sig[k]:5.3f} ({p:4.0f} hPa)  {ra*100 if ra<=1.5 else ra:5.0f}/{re:<5.0f}   {ca:5.2f}/{'--' if ce is None else f'{ce:.2f}'}")
