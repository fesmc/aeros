#!/usr/bin/env python3
"""Summarize and plot an aeros Held-Suarez diagnostics file."""
import sys
import numpy as np
from netCDF4 import Dataset


def load(path):
    d = Dataset(path)
    out = {k: np.array(d[k][:]) for k in
           ("u_zm", "v_zm", "temp_zm", "uv_eddy", "vt_eddy", "eke",
            "temp_var", "pressure", "ps_zm", "u_sfc", "lat")}
    out["nsample"] = d.getncattr("nsample")
    out["t0"] = d.getncattr("time_start")
    out["t1"] = d.getncattr("time_end")
    d.close()
    # ncio writes (lat, lev) as (lev, lat) in numpy order; detect and fix.
    nlat = out["lat"].size
    for k in ("u_zm", "v_zm", "temp_zm", "uv_eddy", "vt_eddy", "eke",
              "temp_var", "pressure"):
        if out[k].shape[0] != nlat:
            out[k] = out[k].T
    return out


def summarize(d, label):
    lat, u, p = d["lat"], d["u_zm"], d["pressure"]
    j, k = np.unravel_index(np.argmax(u), u.shape)
    nh = lat > 0
    sh = lat < 0
    jn, kn = np.unravel_index(np.argmax(np.where(nh[:, None], u, -1e30)), u.shape)
    js, ks = np.unravel_index(np.argmax(np.where(sh[:, None], u, -1e30)), u.shape)

    vt, uv, eke = d["vt_eddy"], d["uv_eddy"], d["eke"]
    jv, kv = np.unravel_index(np.argmax(vt), vt.shape)
    ju, ku = np.unravel_index(np.argmax(uv), uv.shape)

    usfc = d["u_sfc"]
    jw = np.argmax(usfc)
    # tropical easterlies
    trop = np.abs(lat) < 15
    T = d["temp_zm"]
    jeq = np.argmin(np.abs(lat))

    # hemispheric symmetry of [u]
    asym = np.sqrt(np.sum((u - u[::-1, :]) ** 2) / np.sum(u ** 2))

    print(f"\n=== {label} ===")
    print(f"  window: days {d['t0']:.0f}-{d['t1']:.0f}, {d['nsample']} samples")
    print(f"  jet max [u]            {u[j,k]:7.2f} m/s  at {lat[j]:6.1f} deg, {p[j,k]/100:6.1f} hPa")
    print(f"    NH jet               {u[jn,kn]:7.2f} m/s  at {lat[jn]:6.1f} deg, {p[jn,kn]/100:6.1f} hPa")
    print(f"    SH jet               {u[js,ks]:7.2f} m/s  at {lat[js]:6.1f} deg, {p[js,ks]/100:6.1f} hPa")
    print(f"  surface [u] max        {usfc[jw]:7.2f} m/s  at {lat[jw]:6.1f} deg")
    print(f"  surface [u] tropics    {usfc[trop].min():7.2f} m/s  (min, should be easterly)")
    print(f"  max [v'T']             {vt[jv,kv]:7.2f} K m/s at {lat[jv]:6.1f} deg, {p[jv,kv]/100:6.1f} hPa")
    print(f"  max [u'v']             {uv[ju,ku]:7.2f} m2/s2 at {lat[ju]:6.1f} deg, {p[ju,ku]/100:6.1f} hPa")
    print(f"  max EKE                {eke.max():7.2f} m2/s2")
    print(f"  max [T'^2]             {d['temp_var'].max():7.2f} K2")
    print(f"  surface [T] equator    {T[jeq,-1]:7.2f} K")
    print(f"  surface [T] pole       {T[0,-1]:7.2f} K   {T[-1,-1]:7.2f} K")
    print(f"  [T] equator-pole       {T[jeq,-1]-0.5*(T[0,-1]+T[-1,-1]):7.2f} K")
    print(f"  hemispheric asymmetry  {asym:9.4f}")
    # sign checks
    nh_mid = (lat > 20) & (lat < 70)
    sh_mid = (lat < -20) & (lat > -70)
    print(f"  [v'T'] NH midlat mean  {vt[nh_mid].mean():7.3f} (poleward => >0)")
    print(f"  [v'T'] SH midlat mean  {vt[sh_mid].mean():7.3f} (poleward => <0)")
    print(f"  [u'v'] NH midlat mean  {uv[nh_mid].mean():7.3f} (into jet => >0)")
    print(f"  [u'v'] SH midlat mean  {uv[sh_mid].mean():7.3f} (into jet => <0)")
    return dict(u=u, lat=lat, p=p, vt=vt, uv=uv, eke=eke, T=T, usfc=usfc)


def plot(runs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(runs)
    fig, axes = plt.subplots(n, 4, figsize=(17, 4.0 * n), squeeze=False)
    for r, (label, d) in enumerate(runs):
        lat = d["lat"]
        P = d["p"] / 100.0
        LAT = np.broadcast_to(lat[:, None], P.shape)
        panels = [
            ("[u]  m s$^{-1}$", d["u"], np.arange(-30, 41, 5), "RdBu_r"),
            ("[T]  K", d["T"], np.arange(180, 321, 10), "inferno"),
            ("[v'T']  K m s$^{-1}$", d["vt"], np.arange(-24, 25, 3), "RdBu_r"),
            ("[u'v']  m$^2$ s$^{-2}$", d["uv"], np.arange(-60, 61, 10), "RdBu_r"),
        ]
        for c, (title, F, levels, cmap) in enumerate(panels):
            ax = axes[r][c]
            cf = ax.contourf(LAT, P, F, levels=levels, cmap=cmap, extend="both")
            ax.contour(LAT, P, F, levels=levels, colors="k", linewidths=0.3)
            ax.invert_yaxis()
            ax.set_xlim(-90, 90)
            ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
            ax.set_title(f"{label}   {title}", fontsize=10)
            if c == 0:
                ax.set_ylabel("pressure [hPa]")
            ax.set_xlabel("latitude")
            fig.colorbar(cf, ax=ax, pad=0.02)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    runs = []
    for arg in sys.argv[1:]:
        label, path = arg.split("=", 1)
        d = load(path)
        runs.append((label, summarize(d, label)))
    if runs:
        plot(runs, "output/held_suarez.png")
