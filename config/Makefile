.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug  ?= 0
openmp ?= 1

# NOTE: there is no `precision=` switch. aeros is double precision
# throughout, fixed in src/aeros_defs.f90.
#
# That reverses docs/design.md section 4 ("Float32 throughout the core") and it
# was decided by measurement: SHTns has no single-precision CPU path, so an sp
# core converts up and back at every transform, and transforms are 86-97% of a
# dynamics step. Measured, sp was ~17% SLOWER than dp at every production
# configuration. See docs/m0a_results.md section 5.
#
# Single precision survives only as wp_ext in aeros_defs -- the kind used at
# the coupling boundary, where an external model hands over sp fields.

# PATH options
#
# NOTE: openmp=0 and openmp=1 share these directories, so switching openmp does
# need a `make clean`. Only the dependency include dirs are swapped.
srcdir  = src
dyndir  = src/dynamics
physdir = src/physics
objdir  = libaeros/include
bindir  = libaeros/bin
drvdir  = drivers
testdir = tests

## COMPILER CONFIGURATION ##
# Assembled by `configme` from a compiler fragment, a machine fragment and an
# auto-detected netCDF block. The fragment defines the compiler (FC, FFLAGS,
# FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and any
# host-specific overrides (e.g. LFLAGS_EXTRA).
#
# Note the compiler fragment's FFLAGS_BASE references $(objdir), which is why
# the PATH options above must be defined before this point.

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC and LFLAGS_EXTRA set there. Defines the fesm-utils,
# SHTns and FFTW paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

# Include flags shared by every compile in this project. Collected once here
# rather than repeated per rule: aeros links three external dependencies
# (fesm-utils, SHTns, FFTW) and the list will not get shorter.
INCFLAGS = $(INC_NC) $(INC_FESMUTILS) $(INC_SHTNS) $(INC_FFTW) $(INC_INSOL)

###############################################
##
## List of aeros rules and source files
##
###############################################

include config/Makefile_aeros.mk

###############################################
##
## Compilation of complete programs
##
###############################################

check-hash:
	./check_githash > git_aeros.txt

# Static library compilation
aeros-static: $(aeros_base) $(aeros_dynamics) $(aeros_physics) $(aeros_core)
	ar rc $(objdir)/libaeros.a $(aeros_base) $(aeros_dynamics) $(aeros_physics) $(aeros_core)
	ranlib $(objdir)/libaeros.a
	./check_githash > git_aeros.txt
	@echo " "
	@echo "    $(objdir)/libaeros.a is ready."
	@echo " "

## Drivers ####################################
#
# Executables intended to be staged and run by `runme` (see .runme/info.json).
# Tests live under tests/ and are built by the `tests` target below.

run : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/aeros_run.x $(drvdir)/aeros_run.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    aeros_run.x is ready."

# M0a. Both need openmp=1 to mean anything.
#
# bench_sht  : which threading mode (settled -- see aeros_sht_pool_class)
# bench_m0a  : what a timestep costs, i.e. the number docs/design.md section 11
#              admits is currently an extrapolation
#
# bench_m0a is .F90 so the machine and compiler configme resolved can be
# compiled in, and land in the CSV alongside the timings. A benchmark result
# without its toolchain attached is not a result.
bench : bench-sht bench-m0a

bench-sht : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/bench_sht.x $(drvdir)/bench_sht.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    bench_sht.x is ready."

bench-m0a : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) \
			-DAEROS_MACHINE='"$(MACHINE)"' -DAEROS_COMPILER='"$(COMPILER)"' \
			-o $(bindir)/bench_m0a.x $(drvdir)/bench_m0a.F90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    bench_m0a.x is ready."

held-suarez : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/held_suarez.x $(drvdir)/held_suarez.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    held_suarez.x is ready."

## Acceptance tests ###########################

defs : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_defs.x $(testdir)/test_defs.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_defs.x is ready."

spectral : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_spectral.x $(testdir)/test_spectral.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_spectral.x is ready."

vertical : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_vertical.x $(testdir)/test_vertical.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_vertical.x is ready."

vordiv : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_vordiv.x $(testdir)/test_vordiv.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_vordiv.x is ready."

tendency : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_tendency.x $(testdir)/test_tendency.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_tendency.x is ready."

grid : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_grid.x $(testdir)/test_grid.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_grid.x is ready."

semiimp : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_semiimp.x $(testdir)/test_semiimp.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_semiimp.x is ready."

timestep : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_timestep.x $(testdir)/test_timestep.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_timestep.x is ready."

correction : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_correction.x $(testdir)/test_correction.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_correction.x is ready."

transport : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_transport.x $(testdir)/test_transport.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_transport.x is ready."

condensation : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_condensation.x $(testdir)/test_condensation.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_condensation.x is ready."

convection : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_convection.x $(testdir)/test_convection.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_convection.x is ready."

moist-run : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_moist_run.x $(testdir)/test_moist_run.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_moist_run.x is ready."

radiation : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_radiation.x $(testdir)/test_radiation.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_radiation.x is ready."

surface : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_surface.x $(testdir)/test_surface.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_surface.x is ready."

ocean : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_ocean.x $(testdir)/test_ocean.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_ocean.x is ready."

seaice : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_seaice.x $(testdir)/test_seaice.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_seaice.x is ready."

rce : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_rce.x $(testdir)/test_rce.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_rce.x is ready."

held-suarez-test : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_held_suarez.x $(testdir)/test_held_suarez.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_held_suarez.x is ready."

ecckd : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_ecckd.x $(testdir)/test_ecckd.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_ecckd.x is ready."

restart : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_restart.x $(testdir)/test_restart.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_restart.x is ready."

topography : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_topography.x $(testdir)/test_topography.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_topography.x is ready."

land : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_land.x $(testdir)/test_land.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_land.x is ready."

cloud-prog : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_cloud_prog.x $(testdir)/test_cloud_prog.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_cloud_prog.x is ready."

# All acceptance tests
tests : defs grid spectral vertical vordiv tendency semiimp timestep correction transport condensation convection moist-run radiation surface ocean rce held-suarez-test ecckd restart topography land cloud-prog seaice

# Driver executables
drivers : run bench held-suarez

validate-era5 : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/validate_era5.x $(drvdir)/validate_era5.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    validate_era5.x is ready."

rce-long : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/rce_long.x $(drvdir)/rce_long.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    rce_long.x is ready."

probe-insol : aeros-static
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/probe_insol.x $(drvdir)/probe_insol.f90 \
		-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
	@echo "    probe_insol.x is ready."

probe-lw : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/probe_lw.x $(drvdir)/probe_lw.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    probe_lw.x is ready."

all : drivers tests

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make aeros-static : compiles libaeros.a, the aeros static library."
	@echo " make run          : compiles aeros_run.x, the standalone driver."
	@echo " make bench        : compiles the M0a benchmarks (bench_sht, bench_m0a)."
	@echo " make held-suarez  : compiles held_suarez.x, the M1 validation benchmark."
	@echo " make tests        : compiles all acceptance tests."
	@echo " make all          : drivers + tests."
	@echo " make clean        : cleans object files"
	@echo ""
	@echo " Options: debug=1 (bounds/traps), debug=2 (profiling),"
	@echo "          openmp=0 (serial; omp is the default)"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
	rm -f ./git_aeros.txt
