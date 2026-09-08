"""NH 40-70N land seasonal cycle: aeros lowest-level T (t2m field) and skin vs ERA5 t2m and skt."""
import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc"
def era(v):
    f=nc.Dataset(E%v); x=f[v][:].astype(float); lat=f["latitude"][:]; lon=f["longitude"][:]
    if lat[0]>lat[-1]: lat=lat[::-1]; x=x[:,::-1,:]
    return lat,lon,x
lat,lon,lsm=era("lsm"); _,_,t2m=era("t2m"); _,_,skt=era("skt"); _,_,sd=era("sd")
land=lsm[0]>0.5; w=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,len(lon)))
def m(x,mask): return np.array([(x[k]*w*mask).sum()/(w*mask).sum() for k in range(12)])
band=(lat[:,None]>=40)&(lat[:,None]<=70)&land
print("month      ERA5 t2m   ERA5 skt |", "  ".join(f"{n.split('/')[-1][:12]:>12s}" for n in sys.argv[1:]))
rows=[]
for f in sys.argv[1:]:
    A=nc.Dataset(f); la=A["lat"][:]; lo=A["lon"][:]; o=np.argsort(la); la=la[o]
    L=np.array([[lsm[0][np.argmin(abs(lat-y)),np.argmin(abs(lon-x%360))]>0.5 for x in lo] for y in la])
    wa=np.cos(np.deg2rad(la))[:,None]*np.ones((1,len(lo))); b=(la[:,None]>=40)&(la[:,None]<=70)&L
    ta=A["t2m"][:][:,o,:]; sk=A["skin"][:][:,o,:]
    rows.append((np.array([(ta[k]*wa*b).sum()/(wa*b).sum() for k in range(12)]), np.array([(sk[k]*wa*b).sum()/(wa*b).sum() for k in range(12)])))
e2=m(t2m,band); es=m(skt,band)
for k in range(12):
    print(f"{k+1:5d}   {e2[k]:9.1f} {es[k]:9.1f} |", "  ".join(f"{r[0][k]:6.1f}/{r[1][k]:<5.1f}" for r in rows))
print("Jul-Jan  ERA5 t2m %.1f skt %.1f |"%(e2[6]-e2[0], es[6]-es[0]), "  ".join(f"{r[0][6]-r[0][0]:6.1f}/{r[1][6]-r[1][0]:<5.1f}" for r in rows))
print("(model columns: lowest-level T / skin)")
