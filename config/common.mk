# Shared build configuration for aeros (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). This file references
# variables those provide: FFLAGS / FFLAGS_OPENMP (compiler) and INC_NC / LIB_NC
# (machine or auto-detected netCDF).
#
# aeros has three external dependencies, all supplied by the fesm-utils
# checkout at the repo root:
#
#   fesm-utils  precision, nml, ncio, variable_io, coords (projections and
#               SCRIP conservative remapping, for the polar nests at M3/M4),
#               timestepping, timer.
#   SHTns       spherical-harmonic transforms -- the spectral core's kernel.
#               fesm-utils builds it under SHTns/shtns-{serial,omp}/.
#   FFTW        SHTns' longitudinal FFT. Never called directly by aeros, but it
#               must be linked, and after SHTns (static archives resolve
#               left-to-right, so a dependency must follow its dependent).
#
# No linear solver (lis) and no LAPACK are wired: the semi-implicit solve is one
# DENSE nlev x nlev system per spectral degree (src/dynamics/aeros_semiimp.f90),
# which is a forty-line pivoted LU and ~0.9 Mflop per timestep at T42L20.

# --- fesm-utils
# OpenMP is the DEFAULT build (openmp ?= 1). Unlike chion, aeros' OpenMP results
# are NOT bit-identical to serial: the spectral transforms reduce over latitude
# rows, so a threaded reduction reorders floating-point additions. Compare
# serial and OpenMP output with a tolerance, never with `cmp`.
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

# --- SHTns (spherical-harmonic transforms), built by fesm-utils' build.py.
# `shtns.f03`, the Fortran 2003 ISO_C_BINDING interface, lives in include/ and is
# `include`d verbatim by src/aeros_spectral.f90 -- there is no Fortran module
# wrapper in fesm-utils, so aeros supplies its own (docs/design.md section 8).
SHTNSROOT = fesm-utils/SHTns/shtns-serial
INC_SHTNS = -I${SHTNSROOT}/include
LIB_SHTNS = -L${SHTNSROOT}/lib -lshtns

# --- FFTW (pulled in by SHTns; -lm because FFTW references libm)
FFTWROOT = fesm-utils/fftw/fftw-serial
INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib -lfftw3 -lm

# OpenMP build (make openmp=1, the aeros default): swap the serial dependency
# builds for their OpenMP variants. $(FFLAGS_OPENMP) is appended to FFLAGS here.
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    SHTNSROOT = fesm-utils/SHTns/shtns-omp
    INC_SHTNS = -I${SHTNSROOT}/include
    LIB_SHTNS = -L${SHTNSROOT}/lib -lshtns_omp

    FFTWROOT = fesm-utils/fftw/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

    FFLAGS += $(FFLAGS_OPENMP)
endif

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

# SHTns precedes FFTW (it calls FFTW); both precede fesm-utils only by
# convention -- they share no symbols. SHTns is C++-linked internally, so
# -lstdc++ is required when the link is driven by the Fortran compiler.
# fesmc/insol (Laskar 2004 insolation). Prebuilt static lib, no dependencies;
# INSOLROOT is a symlink to the insol checkout (gitignored, like fesm-utils).
INSOLROOT = insol
INC_INSOL = -I${INSOLROOT}/libinsol/include
LIB_INSOL = -L${INSOLROOT}/libinsol/include -linsol

LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_SHTNS) $(LIB_FFTW) $(LIB_INSOL) -lstdc++ $(LFLAGS_EXTRA)
