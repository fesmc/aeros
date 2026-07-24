#!/usr/bin/env bash
# build.sh — source the Intel oneAPI compilers and build bin/aeros_bench.x
#
# Usage:  ./build.sh          (build)
#         ./build.sh clean    (clean then build)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# Intel MKL first (sets MKLROOT for variant 2's -qmkl=sequential), then the ifx
# compilers setvars so ifx 2024 stays the active compiler on PATH.
# SHTns-omp cannot be linked with gfortran/ifort.
# The Intel setvars scripts reference unbound vars, so relax `set -u` around
# them. MKL's setvars marks SETVARS_COMPLETED; clear it so the compilers setvars
# still runs (else it "skips re-execution" and returns non-zero under set -e).
set +u
source /albedo/soft/sw/spack-sw/intel-oneapi-mkl/2022.1.0-akthm3n/setvars.sh intel64
unset SETVARS_COMPLETED
source /albedo/soft/sw/spack-sw/intel-oneapi-compilers/2024.1.0-fdumqva/setvars.sh
set -u

if [[ "${1:-}" == "clean" ]]; then
    make clean
fi

make -j1
echo
echo "Built: $HERE/bin/aeros_bench.x            (variant 1: spectral_shtns)"
echo "Built: $HERE/bin/aeros_bench_batched.x    (variant 2: batched DGEMM)"
echo
echo "ldd check variant 1 (no Intel lib should be 'not found'):"
ldd bin/aeros_bench.x | grep -Ei 'iomp5|shtns|fftw|mkl|not found' || true
echo
echo "ldd check variant 2 (no Intel lib should be 'not found'):"
ldd bin/aeros_bench_batched.x | grep -Ei 'iomp5|fftw|mkl|not found' || true
