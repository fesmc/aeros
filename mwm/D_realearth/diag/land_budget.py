import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc"
def era(v):
    f=nc.Dataset(E%v); x=f[v][:].astype(float); lat=f["latitude"][:]; lon=f["longitude"][:]
    if lat[0]>lat[-1]: lat=lat[::-1]; x=x[:,::-1,:]
    return lat,lon,x
lat,lon,lsm=era("lsm"); land=lsm[0]>0.5; w=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,len(lon)))
band=(lat[:,None]>=40)&(lat[:,None]<=70)&land
def m(x): return np.array([(x[k]*w*band).sum()/(w*band).sum() for k in range(12)])
Ev={v:m(era(v)[2]) for v in ["strd","ssr","sshf","slhf","tcc","sd","t2m","skt"]}
for v in ["strd","ssr","sshf","slhf"]: Ev[v]/=86400.
A=nc.Dataset(sys.argv[1]); la=A["lat"][:]; lo=A["lon"][:]; o=np.argsort(la); la=la[o]
L=np.array([[lsm[0][np.argmin(abs(lat-y)),np.argmin(abs(lon-x%360))]>0.5 for x in lo] for y in la])
wa=np.cos(np.deg2rad(la))[:,None]*np.ones((1,len(lo))); b=(la[:,None]>=40)&(la[:,None]<=70)&L
def ma(v,s=1.0): x=A[v][:][:,o,:]*s; return np.array([(x[k]*wa*b).sum()/(wa*b).sum() for k in range(12)])
Av={"lwd":ma("lwd"),"swn":ma("swn"),"shf":ma("shf"),"lh":ma("evap",2.5e6/86400.),"cover":ma("cover"),"t2m":ma("t2m"),"skin":ma("skin")}
print("NH 40-70N land, monthly: model / ERA5   (fluxes W/m2: LWdown, SWnet(sfc), SH up, LH up; cover)")
print("mon   LWd          SWn          SH           LH           cover        skin-air")
for k in range(12):
    print(f"{k+1:3d}  {Av['lwd'][k]:5.0f}/{Ev['strd'][k]:<5.0f}  {Av['swn'][k]:5.0f}/{Ev['ssr'][k]:<5.0f}  {Av['shf'][k]:5.0f}/{-Ev['sshf'][k]:<5.0f}  {Av['lh'][k]:5.0f}/{-Ev['slhf'][k]:<5.0f}  {Av['cover'][k]:5.2f}/{Ev['tcc'][k]:<5.2f}  {Av['skin'][k]-Av['t2m'][k]:5.1f}/{Ev['skt'][k]-Ev['t2m'][k]:<5.1f}")
