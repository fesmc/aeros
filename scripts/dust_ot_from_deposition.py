"""Dust optical thickness climatology from a CLIMBER-X dust-deposition map.

CLIMBER-X (src/atm/dust.f90) relates the atmospheric dust load L to deposition D by
    D = (c_dust_dry / h_dust + c_dust_wet * P) * L,      tau_ot = c_dust_mec * L,
with c_dust_dry 2.5e-3 m/s, h_dust = c_dhs_1 2000 m, c_dust_wet 0.075 m2/kg (Bauer 2010),
c_dust_mec 950 m2/kg (nml/atm_par.nml). Inverting that with the ERA5 precipitation
climatology turns its prescribed deposition maps (input/fake_dust_*.nc, kg/m2/s, monthly,
5x5) into the monthly optical thickness aeros_radiation reads (dust_file).

usage: python3 scripts/dust_ot_from_deposition.py fake_dust_Lambert2015_present.nc input/dust_ot_lambert2015_present.nc
"""
import sys, numpy as np, netCDF4 as nc
C_DRY, H_DUST, C_WET, C_MEC = 2.5e-3, 2000.0, 0.075, 950.0
ERA5_TP = "/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_tp_1991-2020_clim.nc"
src, dst = sys.argv[1], sys.argv[2]
d = nc.Dataset(src); dep = np.maximum(0.0, np.ma.filled(d["dust_dep"][:], 0.0)); lat = d["lat"][:]; lon = d["lon"][:]
g = nc.Dataset(ERA5_TP); tp = g["tp"][:]; la = g["latitude"][:]; lo = g["longitude"][:]
prc = np.array([[[tp[m, np.argmin(abs(la - y)), np.argmin(abs(lo - x % 360))] for x in lon] for y in lat]
                for m in range(12)]) * 1000.0 / 86400.0          # m/day -> kg/m2/s
ot = dep / (C_DRY / H_DUST + C_WET * prc) * C_MEC
o = nc.Dataset(dst, "w"); o.createDimension("lon", len(lon)); o.createDimension("lat", len(lat)); o.createDimension("time", 12)
for n, v in (("lon", lon), ("lat", lat), ("time", np.arange(1, 13))):
    x = o.createVariable(n, "f8", (n,)); x[:] = v
x = o.createVariable("dust_ot", "f8", ("time", "lat", "lon")); x[:] = ot
x.long_name = "dust optical thickness"; x.units = "1"
o.source = f"{src}: deposition -> load by the CLIMBER-X dust.f90 balance with ERA5 1991-2020 tp, x c_dust_mec"
o.close()
w = np.cos(np.deg2rad(lat))[None, :, None] * np.ones_like(ot)
print(f"{dst}: global mean ot {np.sum(ot*w)/np.sum(w):.3f}, max {ot.max():.2f}")
