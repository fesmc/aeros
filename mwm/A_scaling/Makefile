# Makefile — aeros M0a scaling benchmark (Option A)
#
# Re-runnable. Assumes the Intel oneAPI compilers are already sourced
# (build.sh does that for you). Uses ifx (LLVM-based); ifort/gfortran will NOT
# link the OpenMP SHTns build.

# Force ifx: the environment may export FC=f77, and make's implicit default is
# also f77 — neither can link the OpenMP SHTns build. Override on the make line
# to use a different compiler (not recommended).
FC      := ifx

# --- external OpenMP libraries (verified paths) ------------------------------
SH  ?= /albedo/work/projects/p_forclima/robinson/models/fesm-utils/SHTns/shtns-omp
FFT ?= /albedo/work/projects/p_forclima/robinson/models/fesm-utils/fftw/fftw-omp

# Directory holding libiomp5.so — rpath'd in so the binary runs on a compute
# node without sourcing setvars (verified with `find` under the oneapi tree).
INTEL_RT ?= /albedo/soft/sw/spack-sw/intel-oneapi-compilers/2024.1.0-fdumqva/compiler/2024.1/lib

# MKL runtime dir (for variant 2). Defaults to $(MKLROOT)/lib/intel64, which
# build.sh sets by sourcing the MKL setvars; rpath'd so the binary resolves
# libmkl_*.so on a compute node without sourcing.
MKL_RT ?= $(MKLROOT)/lib/intel64

# --- Variant 1 (spectral_shtns): SHTns transforms, one 2-D field at a time ---
FFLAGS  ?= -qopenmp -O3 -xHost -fp-model=precise -I$(SH)/include
LDFLAGS ?= -L$(SH)/lib -lshtns_omp -L$(FFT)/lib -lfftw3_omp -lfftw3 -lm \
           -Wl,-rpath,$(INTEL_RT)

# --- Variant 2 (batched): custom batched Legendre = sequential MKL dgemm per
#     wavenumber (OpenMP over m) + SERIAL FFTW (per-thread execution). No SHTns.
#     -qmkl=sequential so OUR OpenMP-over-m loop owns the parallelism, not MKL.
FFLAGS_B  ?= -qopenmp -O3 -xHost -fp-model=precise -qmkl=sequential -I$(FFT)/include
LDFLAGS_B ?= -qmkl=sequential -L$(FFT)/lib -lfftw3 -lm \
             -Wl,-rpath,$(INTEL_RT) -Wl,-rpath,$(MKL_RT)

SRC    = src/aeros_bench.f90
BIN    = bin/aeros_bench.x
SRC_B  = src/aeros_bench_batched.f90
BIN_B  = bin/aeros_bench_batched.x

.PHONY: all clean

all: $(BIN) $(BIN_B)

$(BIN): $(SRC) | bin
	$(FC) $(FFLAGS) $(SRC) $(LDFLAGS) -o $(BIN)

$(BIN_B): $(SRC_B) | bin
	$(FC) $(FFLAGS_B) $(SRC_B) $(LDFLAGS_B) -o $(BIN_B)

bin:
	mkdir -p bin

clean:
	rm -f $(BIN) $(BIN_B) *.mod src/*.mod
