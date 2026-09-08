import sys, numpy as np, netCDF4 as nc
tgt=np.loadtxt("mwm/D_realearth/diag/lid_target_era5.txt")
for n in sys.argv[1:]:
    A=nc.Dataset(f"mwm/D_realearth/diag/{n}.nc"); lat=A["lat"][:]; d=A["t"][:][:,0,:]-tgt; w=np.cos(np.deg2rad(lat)); u=A["u"][:].mean(0)
    b=lambda x,lo,hi: (lambda m:(x[m]*w[m]).sum()/w[m].sum())((lat>=lo)&(lat<=hi))
    print(f"{n:9s} lid-target: 60-90S {b(d.mean(0),-90,-60):5.1f} 30-60S {b(d.mean(0),-60,-30):5.1f} trop {b(d.mean(0),-15,15):5.1f} 30-60N {b(d.mean(0),30,60):5.1f} 60-90N {b(d.mean(0),60,90):5.1f} | max [u] {abs(u).max():5.1f}")
