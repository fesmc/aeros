import sys, numpy as np, netCDF4 as nc
E="/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_%s_1991-2020_clim.nc"
def era(v):
    f=nc.Dataset(E%v); x=f[v][:].astype(float); lat=f["latitude"][:]; lon=f["longitude"][:]
    if lat[0]>lat[-1]: lat=lat[::-1]; x=x[:,::-1,:]
    return lat,lon,x
lat,lon,lsm=era("lsm"); _,_,slhf=era("slhf"); _,_,sshf=era("sshf")
# ERA5 accumulations J/m2 per day -> W/m2 (sign: positive downward in ERA5)
slhf=-slhf/86400.; sshf=-sshf/86400.
w=np.cos(np.deg2rad(lat))[None,:,None]
def mean(x,m): return (x*w*m).sum()/(w*m*np.ones_like(x)).sum()
land=(lsm[0]>0.5)[None]; oc=~land
print("ERA5   LH ocean/land %.1f/%.1f  SH ocean/land %.1f/%.1f"%(mean(slhf,oc),mean(slhf,land),mean(sshf,oc),mean(sshf,land)))
for f in sys.argv[1:]:
    A=nc.Dataset(f); la=A["lat"][:]; lo=A["lon"][:]; o=np.argsort(la)
    ev=A["evap"][:][:,o,:]*2.5e6/86400.; sh=A["shf"][:][:,o,:]; la=la[o]
    # land mask on the model grid: nearest ERA5 lsm
    L=np.array([[lsm[0][np.argmin(abs(lat-y)),np.argmin(abs(lon-x%360))]>0.5 for x in lo] for y in la])[None]
    wa=np.cos(np.deg2rad(la))[None,:,None]
    def m(x,mask): return (x*wa*mask).sum()/(wa*mask*np.ones_like(x)).sum()
    print("%-45s LH ocean/land %.1f/%.1f  SH ocean/land %.1f/%.1f"%(f.split('/')[-1],m(ev,~L),m(ev,L),m(sh,~L),m(sh,L)))
