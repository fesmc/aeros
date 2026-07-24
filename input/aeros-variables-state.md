# aeros state variables

Instantaneous model state, written on the Gaussian grid at each output time.
Consumed by `aeros_write_state` in `src/aeros_io.f90` via fesm-utils'
`variable_io`, following yelmo's convention.

Editing this table changes what a run writes. Adding a row also needs a `case`
in `aeros_write_var_state`, which fails loudly on a name it does not know
rather than silently skipping it — a variable in the table but not the
dispatcher is a typo, not a request.

Dimensions are the netCDF dimension names created by `aeros_write_init`;
`time` is appended automatically and must not appear here.

| id | variable          | dimensions       | units       | long_name                                          |
|----|-------------------|------------------|-------------|----------------------------------------------------|
|  1 | ps                | lon, lat         | Pa          | Surface pressure                                   |
|  2 | u                 | lon, lat, lev    | m s-1       | Zonal wind                                         |
|  3 | v                 | lon, lat, lev    | m s-1       | Meridional wind                                    |
|  4 | temp              | lon, lat, lev    | K           | Temperature                                        |
|  5 | qv                | lon, lat, lev    | kg kg-1     | Specific humidity                                  |
