###############################################
##
## Rules for individual aeros files
##
###############################################
#
# One explicit rule per object, with hand-written dependencies. There is no
# auto-dependency generation and no wildcard, matching yelmo and chion. When
# adding a source file, add its rule here AND add it to the appropriate object
# list at the bottom.
#
# Build order matters: aeros_defs must come first, since every other module
# uses it. Object lists are ordered so that `ar` receives dependencies before
# dependents.
#
# $(INCFLAGS) (defined in config/Makefile) carries the netCDF, fesm-utils,
# SHTns and FFTW include paths. Every rule takes all four rather than a minimal
# subset: the subset would be a maintenance trap for no measurable gain, since
# a -I to an unused directory costs nothing.

## aeros base #################################

# aeros_defs is the ONLY preprocessed source. The .F90 extension makes both
# gfortran and ifort preprocess it without a -cpp/-fpp flag; the Makefile
# passes -DAEROS_DP when precision=dp.
$(objdir)/aeros_defs.o: $(srcdir)/aeros_defs.F90
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# The SHTns wrapper. `include 'shtns.f03'` resolves through $(INC_SHTNS),
# so this object is the one that would fail first if fesm-utils' SHTns build
# were missing -- which makes it a useful canary, and a reason to keep it
# early in the build order.
$(objdir)/aeros_spectral.o: $(srcdir)/aeros_spectral.f90 \
							$(objdir)/aeros_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

$(objdir)/aeros_grid.o: $(srcdir)/aeros_grid.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

$(objdir)/aeros_state.o: $(srcdir)/aeros_state.f90 \
							$(objdir)/aeros_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

## aeros dynamics #############################
#
# The spectral primitive-equation core (docs/design.md sections 3.2, 4).
# Empty at M0; filled at M1, validated against Held-Suarez.

## aeros physics ##############################
#
# Column physics: radiation, condensation, surface tiles (sections 5, 6).
# Empty at M0; filled at M2.

## aeros core #################################
#
# IO and the public facade. These sit above everything else and must be
# archived last.

$(objdir)/aeros_io.o: $(srcdir)/aeros_io.f90 \
							$(objdir)/aeros_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

$(objdir)/aeros.o: $(srcdir)/aeros.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_grid.o $(objdir)/aeros_state.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

###############################################
##
## Object lists
##
###############################################

aeros_base =     $(objdir)/aeros_defs.o \
                 $(objdir)/aeros_spectral.o \
                 $(objdir)/aeros_grid.o \
                 $(objdir)/aeros_state.o

# M1
aeros_dynamics =

# M2
aeros_physics =

aeros_core =     $(objdir)/aeros_io.o \
                 $(objdir)/aeros.o
