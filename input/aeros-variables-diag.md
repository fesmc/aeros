# aeros diagnostic variables

Zonal-mean, time-mean statistics accumulated over an averaging window by
`src/aeros_diagnostics.f90`, written once at the end of a run. These are the
fields the Held & Suarez (1994) benchmark is judged on.

Primes are departures from the instantaneous zonal mean, so the eddy fluxes are
TOTAL — transient plus stationary. There is no `time` dimension: one file is
one averaging window, and its bounds are written as attributes.

| id | variable          | dimensions       | units       | long_name                                          |
|----|-------------------|------------------|-------------|----------------------------------------------------|
|  1 | u_zm              | lat, lev         | m s-1       | Zonal-mean zonal wind                              |
|  2 | v_zm              | lat, lev         | m s-1       | Zonal-mean meridional wind                         |
|  3 | temp_zm           | lat, lev         | K           | Zonal-mean temperature                             |
|  4 | uv_eddy           | lat, lev         | m2 s-2      | Eddy momentum flux                                 |
|  5 | vt_eddy           | lat, lev         | K m s-1     | Eddy heat flux                                     |
|  6 | eke               | lat, lev         | m2 s-2      | Eddy kinetic energy                                |
|  7 | temp_var          | lat, lev         | K2          | Eddy temperature variance                          |
|  8 | pressure          | lat, lev         | Pa          | Full-level pressure at the mean surface pressure   |
|  9 | ps_zm             | lat              | Pa          | Zonal-mean surface pressure                        |
| 10 | u_sfc             | lat              | m s-1       | Zonal-mean lowest-level zonal wind                 |
