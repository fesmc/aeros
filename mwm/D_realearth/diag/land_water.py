"""NH 40-70N land: P, E, P-E, SW down all-sky vs clear-sky, albedo; model vs ERA5."""
import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc"
def era(v):
    f=nc.Dataset(E%v); x=f[v][:].astype(float); lat=f["latitude"][:]; lon=f["longitude"][:]
    if lat[0]>lat[-1]: lat=lat[::-1]; x=x[:,::-1,:]
    return lat,lon,x
lat,lon,lsm=era("lsm"); land=lsm[0]>0.5; w=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,len(lon)))
band=(lat[:,None]>=40)&(lat[:,None]<=70)&land
def m(x): return np.array([(x[k]*w*band).sum()/(w*band).sum() for k in range(12)])
Ev={v:m(era(v)[2]) for v in ["tp","e","ssrd","ssrdc","ssr","slhf"]}
Ev["tp"]*=1000; Ev["e"]*=-1000; Ev["ssrd"]/=86400; Ev["ssrdc"]/=86400; Ev["ssr"]/=86400
A=nc.Dataset(sys.argv[1]); la=A["lat"][:]; lo=A["lon"][:]; o=np.argsort(la); la=la[o]
L=np.array([[lsm[0][np.argmin(abs(lat-y)),np.argmin(abs(lon-x%360))]>0.5 for x in lo] for y in la])
wa=np.cos(np.deg2rad(la))[:,None]*np.ones((1,len(lo))); b=(la[:,None]>=40)&(la[:,None]<=70)&L
def ma(v,s=1.0): x=A[v][:][:,o,:]*s; return np.array([(x[k]*wa*b).sum()/(wa*b).sum() for k in range(12)])
P=ma("precip2d"); Ee=ma("evap"); swn=ma("swn"); swup=ma("swup")
if P.max()<1e-2: P*=86400
if Ee.max()<1e-2: Ee*=86400
swd=swn+swup
print("mon   P mm/d      E mm/d      P-E        SWdown(all)   SWnet       albedo   | ERA5 SWdown clr, cloud effect on SWd")
for k in range(12):
    print(f"{k+1:3d}  {P[k]:5.2f}/{Ev['tp'][k]:<5.2f}  {Ee[k]:5.2f}/{Ev['e'][k]:<5.2f}  {P[k]-Ee[k]:5.2f}/{Ev['tp'][k]-Ev['e'][k]:<5.2f}  {swd[k]:5.0f}/{Ev['ssrd'][k]:<5.0f}  {swn[k]:5.0f}/{Ev['ssr'][k]:<5.0f}  {swup[k]/max(swd[k],1):5.2f}/{1-Ev['ssr'][k]/max(Ev['ssrd'][k],1):<5.2f} | {Ev['ssrdc'][k]:5.0f}  {Ev['ssrd'][k]-Ev['ssrdc'][k]:5.0f}")
print("annual P-E model %.2f ERA5 %.2f mm/d"%((P-Ee).mean(),(Ev['tp']-Ev['e']).mean()))
