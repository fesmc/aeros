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

# No source in src/ is preprocessed: there is no build-time precision switch
# (see config/Makefile). Only drivers/bench_m0a.F90 is .F90, to compile in the
# machine and compiler names.
$(objdir)/aeros_defs.o: $(srcdir)/aeros_defs.f90
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
# The spectral primitive-equation core (docs/design.md sections 3.2, 4),
# validated against Held-Suarez.

# M1.1: hybrid sigma-pressure coordinate and the hydrostatic relation on it.
$(objdir)/aeros_vertical.o: $(dyndir)/aeros_vertical.f90 \
							$(objdir)/aeros_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M1.2: vorticity/divergence <-> (u,v). Sits on the SHTns wrapper's vector
# transforms, so it must be archived after aeros_spectral.
$(objdir)/aeros_vordiv.o: $(dyndir)/aeros_vordiv.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M1.3: primitive-equation right-hand sides by the transform method. Sits on
# everything above it -- the vertical coordinate, the vor/div mapping and the
# transform pool -- so it is archived last among the dynamics.
$(objdir)/aeros_tendency.o: $(dyndir)/aeros_tendency.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_vertical.o $(objdir)/aeros_vordiv.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M1.4: the semi-implicit gravity-wave solve. Needs the vertical coordinate
# (the reference state and its alphas), the truncation (one factorization per
# degree) and the tendency type it consumes.
$(objdir)/aeros_semiimp.o: $(dyndir)/aeros_semiimp.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_vertical.o $(objdir)/aeros_tendency.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M2.2: the additive tendency correction (sections 3.7, 3.8). Consumes the
# assembled tendency type and nothing else in the dynamics -- it is applied
# spectrally, after aeros_tendency_spectral, precisely so that it does not have
# to know about the grid-space machinery.
$(objdir)/aeros_correction.o: $(dyndir)/aeros_correction.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_tendency.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M2.3: prognostic humidity, transported on the grid by a positive-definite
# finite-volume scheme. Off the spectral core entirely -- it needs only the
# vertical coordinate (to diagnose layer masses) and the grid geometry.
$(objdir)/aeros_moisture.o: $(dyndir)/aeros_moisture.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M1.4: the leapfrog itself -- time levels, RAW filter, hyperdiffusion. Sits on
# top of every other dynamics module and is the only place a prognostic
# variable changes value.
$(objdir)/aeros_timestep.o: $(dyndir)/aeros_timestep.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_state.o $(objdir)/aeros_vertical.o \
							$(objdir)/aeros_vordiv.o $(objdir)/aeros_tendency.o \
							$(objdir)/aeros_semiimp.o $(objdir)/aeros_held_suarez.o \
							$(objdir)/aeros_correction.o $(objdir)/aeros_moisture.o \
							$(objdir)/aeros_held_suarez.o $(objdir)/aeros_condensation.o \
							$(objdir)/aeros_convection.o $(objdir)/aeros_surface.o \
							$(objdir)/aeros_radiation.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

## aeros physics ##############################
#
# Column physics: radiation, condensation, surface tiles (sections 5, 6).
# Empty at M0; the Held-Suarez benchmark forcing lands at M1.5 and the real
# parameterizations at M2.

# M1.5: Held & Suarez (1994) idealized forcing, the M1 validation benchmark.
# Consumes the tendency's grid-space work arrays, so it is compiled after the
# dynamics and archived between them and the core.
$(objdir)/aeros_held_suarez.o: $(physdir)/aeros_held_suarez.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_vertical.o $(objdir)/aeros_tendency.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M2.3b: large-scale condensation. A column process -- it needs the vertical
# coordinate for the layer pressures and nothing spectral -- applied at the
# grid seam like the Held-Suarez forcing.
$(objdir)/aeros_condensation.o: $(physdir)/aeros_condensation.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M2.3d: moist convective adjustment. A column process; reuses q_sat from
# condensation, so it is compiled after it.
$(objdir)/aeros_convection.o: $(physdir)/aeros_convection.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o \
							$(objdir)/aeros_condensation.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M2.4: radiation (section 5). Ported SESAM broadband LW band kernel on the
# resolved column. A column process like the others; needs the vertical
# coordinate for the layer pressures and the grid for its geometry.
$(objdir)/aeros_radiation.o: $(physdir)/aeros_radiation.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o \
							$(objdir)/aeros_grid.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# M2.4c: surface energy/moisture budget (section 6.1). Prescribed-SST
# aquaplanet with bulk turbulent fluxes; reuses q_sat from condensation.
$(objdir)/aeros_surface.o: $(physdir)/aeros_surface.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o \
							$(objdir)/aeros_condensation.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

## aeros core #################################
#
# IO and the public facade. These sit above everything else and must be
# archived last.

# Global conservation integrals. In src/ rather than src/dynamics/ because the
# budgets are a property of the model state, not of the integrator, and M2 adds
# water and latent energy to them. Depends on the vertical coordinate for the
# layer masses, so it is archived after the dynamics.
$(objdir)/aeros_budget.o: $(srcdir)/aeros_budget.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# Zonal-mean and eddy statistics, accumulated online. In src/ for the same
# reason as aeros_budget: a property of the model state, not of the integrator.
$(objdir)/aeros_diagnostics.o: $(srcdir)/aeros_diagnostics.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

# Output, driven by the variable tables under input/. Depends on the
# diagnostics type it writes.
$(objdir)/aeros_io.o: $(srcdir)/aeros_io.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_vertical.o \
							$(objdir)/aeros_diagnostics.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INCFLAGS) -c -o $@ $<

$(objdir)/aeros.o: $(srcdir)/aeros.f90 \
							$(objdir)/aeros_defs.o $(objdir)/aeros_spectral.o \
							$(objdir)/aeros_grid.o $(objdir)/aeros_state.o \
							$(objdir)/aeros_vertical.o $(objdir)/aeros_timestep.o
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

aeros_dynamics = $(objdir)/aeros_vertical.o \
                 $(objdir)/aeros_vordiv.o \
                 $(objdir)/aeros_tendency.o \
                 $(objdir)/aeros_correction.o \
                 $(objdir)/aeros_moisture.o \
                 $(objdir)/aeros_semiimp.o

aeros_physics =  $(objdir)/aeros_held_suarez.o \
                 $(objdir)/aeros_condensation.o \
                 $(objdir)/aeros_convection.o \
                 $(objdir)/aeros_radiation.o \
                 $(objdir)/aeros_surface.o

# aeros_timestep sits above the physics as well as the dynamics -- it is what
# calls the forcing -- so it is archived after both.
aeros_core =     $(objdir)/aeros_timestep.o \
                 $(objdir)/aeros_budget.o \
                 $(objdir)/aeros_diagnostics.o \
                 $(objdir)/aeros_io.o \
                 $(objdir)/aeros.o
