"""Zonal-mean omega [hPa/day] from continuity on the zonal-mean v: model (monthly v on sigma) vs ERA5 (v on p)."""
import sys, numpy as np, netCDF4 as nc
a=6.371e6; g=9.81
E="/Users/alrobi001/data/era5/"
f=nc.Dataset(E+"monthly-pressure-levels/era5_monthly-pressure-levels_v_1991-2020_clim.nc"); latE=f["latitude"][:]; pl=f["pressure_level"][:].astype(float); vE=f["v"][:].astype(float).mean(-1).mean(0)  # (plev,lat) annual
if latE[0]>latE[-1]: latE=latE[::-1]; vE=vE[:,::-1]
o=np.argsort(pl); pl=pl[o]; vE=vE[o]
def omega_from_v(v, p, lat):   # v(p,lat), p increasing [hPa]; omega(p,lat) = -(1/(a cos)) d/dphi [cos * int_0^p v dp]
    phi=np.deg2rad(lat); c=np.cos(phi)
    V=np.zeros_like(v); 
    for k in range(1,len(p)): V[k]=V[k-1]+0.5*(v[k]+v[k-1])*(p[k]-p[k-1])   # hPa m/s
    om=-np.gradient(V*c[None,:], phi, axis=1)/(a*c[None,:])                   # hPa/s
    return om*86400
# remove the (small) mass-weighted mean of v so the streamfunction closes at the surface
def clean(v,p): w=np.gradient(p); return v-(v*w[:,None]).sum(0)/w.sum()
omE=omega_from_v(clean(vE,pl),pl,latE)
A=nc.Dataset(sys.argv[1]); latA=A["lat"][:]; sig=A["lev"][:]; o=np.argsort(latA); latA=latA[o]; vA=A["v"][:].mean(0)[:,o]   # (lev,lat)
pA=sig*1000.
omA=omega_from_v(clean(vA,pA),pA,latA)
print("annual zonal-mean omega [hPa/day] from [v] continuity; model rows at its sigma levels, ERA5 interpolated to the same p")
bands=[("60-75S",-75,-60),("45-60S",-60,-45),("30-45S",-45,-30),("15-30S",-30,-15),("15S-15N",-15,15),("15-30N",15,30),("30-45N",30,45),("45-60N",45,60),("60-75N",60,75)]
print("sigma   "+"".join(f"{b[0]:>12s}" for b in bands))
w=np.cos(np.deg2rad(latA))
for k in range(len(sig)):
    e=np.array([np.interp(pA[k],pl,omE[:,j]) for j in range(len(latE))]); e=np.interp(latA,latE,e)
    row=[]
    for _,lo,hi in bands:
        m=(latA>=lo)&(latA<=hi); row.append(f"{(omA[k][m]*w[m]).sum()/w[m].sum():5.1f}/{(e[m]*w[m]).sum()/w[m].sum():<5.1f}")
    print(f"{sig[k]:.3f}  "+" ".join(row))
