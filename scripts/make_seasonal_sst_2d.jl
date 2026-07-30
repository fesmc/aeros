# Build a full-2D seasonal-SST boundary condition from the ERA5 1991-2020 monthly
# climatology, for the coupled land runs (real geography needs a real 2D ocean).
#
# ERA5 skin temperature (skt) is masked to ocean (land-sea mask lsm < 0.5); land
# points are filled with the zonal-mean ocean SST at their latitude so the field
# is defined everywhere (model land points ignore SST -- they carry a prognostic
# soil temperature -- but coastal model ocean points then interpolate from clean
# ocean values, not land skin temps). The result stays on the ERA5 grid; the
# model regrids each month at run time via aeros_bcinput_read_field(itime=month).
#
#   julia scripts/make_seasonal_sst_2d.jl [out.nc]
#
# Writes input/sst_seasonal_2d.nc: sst(longitude, latitude, month) [K], laid out
# like the ERA5 source (coords named longitude/latitude for the bcinput reader).

using NCDatasets
using Statistics
using Printf

const ERA5_SL = "/Users/alrobi001/data/era5/monthly-single-levels"
slc(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")
out_nc = length(ARGS) >= 1 ? ARGS[1] : "input/sst_seasonal_2d.nc"

dskt = NCDataset(slc("skt")); dlsm = NCDataset(slc("lsm"))
skt = Array(dskt["skt"][:, :, :])          # (lon, lat, month)
lsm = Array(dlsm["lsm"][:, :, :])
lon = Array(dskt["longitude"][:]); lat = Array(dskt["latitude"][:])
close(dskt); close(dlsm)

nlon, nlat, nmon = size(skt)
@assert nmon == 12
ocean = mean(lsm; dims = 3)[:, :, 1] .< 0.5   # (lon, lat)

sst = similar(skt)
for m in 1:nmon
    for j in 1:nlat
        oc = @view ocean[:, j]
        colf = skt[:, j, m]
        if any(oc)
            zm = mean(colf[oc])              # zonal-mean ocean SST at this lat
            colf[.!oc] .= zm                 # fill land points
        end                                   # (no ocean at this lat -> keep skt)
        sst[:, j, m] = colf
    end
end

isfile(out_nc) && rm(out_nc)
NCDataset(out_nc, "c") do ds
    defDim(ds, "longitude", nlon); defDim(ds, "latitude", nlat); defDim(ds, "month", nmon)
    vlo = defVar(ds, "longitude", Float64, ("longitude",)); vlo[:] = lon
    vlo.attrib["units"] = "degrees_east"
    vla = defVar(ds, "latitude",  Float64, ("latitude",));  vla[:] = lat
    vla.attrib["units"] = "degrees_north"
    vmo = defVar(ds, "month", Float64, ("month",)); vmo[:] = collect(1.0:nmon)
    vs  = defVar(ds, "sst", Float64, ("longitude", "latitude", "month"))
    vs[:, :, :] = sst
    vs.attrib["units"] = "K"
    vs.attrib["long_name"] = "ocean-masked, land-filled ERA5 skt climatology"
end

@printf("wrote %s : sst(%d lon, %d lat, %d month)\n", out_nc, nlon, nlat, nmon)
# sanity: land-fill removed land extremes -> ocean-like range at a NH-land latitude
j = argmin(abs.(lat .- 60.0))
@printf("  lat %+.0f : Jan SST range over lon %.1f-%.1f K (ocean-filled)\n",
        lat[j], minimum(sst[:, j, 1]), maximum(sst[:, j, 1]))
