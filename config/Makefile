.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug     ?= 0
openmp    ?= 1
precision ?= sp

# Working precision of the aeros state and public interfaces (docs/design.md
# section 4, "Float32 throughout the core"). `sp` is the production setting.
#
# `dp` is not merely a debugging aid here. SHTns -- the spherical-harmonic
# transform library, and the single largest cost centre in the model (25-40% of
# runtime at T42L19, docs/design.md section 4) -- has a DOUBLE-PRECISION C
# interface: `real(c_double)` spatial arrays and `complex(c_double_complex)`
# spectral arrays. Single precision exists in SHTns only on GPU. So an sp build
# cannot pass its state across the transform boundary without a copy-convert,
# which a dp build does not need (see src/aeros_spectral.f90).
#
# Both builds are therefore real configurations to be measured against each
# other at M0a, not a production setting plus a validation toggle:
#   sp: half the memory traffic and SIMD width in the physics and the
#       semi-implicit solve, paid for with a convert on every transform.
#   dp: zero-copy transforms, twice the bandwidth everywhere else.
# Which wins is an empirical question about a machine, and M0a answers it.
ifeq ($(precision), dp)
	psuffix = -dp
else ifeq ($(precision), sp)
	psuffix =
else
$(error precision must be 'sp' or 'dp', got '$(precision)')
endif

# PATH options
#
# The sp build keeps the unsuffixed paths, so the artifact locations named by
# .configme/manifest.toml and .runme/info.json are unchanged and the dp build
# is purely additive. The two builds therefore coexist without `make clean`.
#
# NOTE: openmp=0 and openmp=1 share these directories, so switching openmp does
# need a `make clean`. Only the dependency include dirs are swapped.
srcdir  = src
dyndir  = src/dynamics
physdir = src/physics
objdir  = libaeros/include$(psuffix)
bindir  = libaeros/bin$(psuffix)
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

# Select the working precision. src/aeros_defs.F90 is the only preprocessed
# source; its .F90 extension makes both gfortran and ifort preprocess it
# automatically, so no -cpp/-fpp flag is needed and no configme compiler
# fragment has to know about this switch.
ifeq ($(precision), dp)
	DFLAGS += -DAEROS_DP
endif

# Include flags shared by every compile in this project. Collected once here
# rather than repeated per rule: aeros links three external dependencies
# (fesm-utils, SHTns, FFTW) and the list will not get shorter.
INCFLAGS = $(INC_NC) $(INC_FESMUTILS) $(INC_SHTNS) $(INC_FFTW)

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

## Acceptance tests ###########################

defs : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_defs.x $(testdir)/test_defs.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_defs.x is ready."

spectral : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_spectral.x $(testdir)/test_spectral.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_spectral.x is ready."

grid : aeros-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -o $(bindir)/test_grid.x $(testdir)/test_grid.f90 \
			-L${CURDIR}/$(objdir) -laeros $(LFLAGS)
		@echo "    test_grid.x is ready."

# All acceptance tests
tests : defs grid spectral

# Driver executables
drivers : run

all : drivers tests

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make aeros-static : compiles libaeros.a, the aeros static library."
	@echo " make run          : compiles aeros_run.x, the standalone driver."
	@echo " make tests        : compiles all acceptance tests."
	@echo " make all          : drivers + tests."
	@echo " make clean        : cleans object files"
	@echo ""
	@echo " Options: debug=1 (bounds/traps), debug=2 (profiling),"
	@echo "          openmp=0 (serial; omp is default), precision=dp (sp is default)"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
	rm -f ./git_aeros.txt
