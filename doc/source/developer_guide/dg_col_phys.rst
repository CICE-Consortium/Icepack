:tocdepth: 3

.. _dev_colphys:

Icepack Column Physics
========================

File List
------------------------------------

The column physics source code contains the following files

| **columnphysics/**   the column physics code
|    **icepack_aerosol.F90**       handles most work associated with the aerosol tracers
|    **icepack_age.F90**           handles most work associated with the age tracer
|    **icepack_algae.F90**         biogeochemistry
|    **icepack_atmo.F90**          stability-based parameterization for calculation of turbulent ice–atmosphere fluxes
|    **icepack_brine.F90**         evolves the brine height tracer
|    **icepack_firstyear.F90**     handles most work associated with the first-year ice area tracer
|    **icepack_flux.F90**          fluxes needed/produced by the model
|    **icepack_fsd.F90**           supports floe size distribution
|    **icepack_intfc.F90**         interface routines for linking Icepack with a host sea ice model
|    **icepack_itd.F90**           utilities for managing ice thickness distribution
|    **icepack_kinds.F90**         basic definitions of reals, integers, etc.
|    **icepack_mechred.F90**       mechanical redistribution (ridging)
|    **icepack_meltpond_cesm.F90** CESM melt pond parameterization
|    **icepack_meltpond_lvl.F90**  level-ice melt pond parameterization
|    **icepack_meltpond_topo.F90** topo melt pond parameterization
|    **icepack_mushy_physics.F90** physics routines for mushy thermodynamics
|    **icepack_ocean.F90**         mixed layer ocean model
|    **icepack_orbital.F90**       orbital parameters for Delta-Eddington shortwave parameterization
|    **icepack_parameters.F90**    basic model parameters including physical and numerical constants requried for column package
|    **icepack_shortwave.F90**     shortwave and albedo parameterizations
|    **icepack_therm_0layer.F90**  zero-layer thermodynamics of :cite:`Semtner76`
|    **icepack_therm_bl99.F90**    multilayer thermodynamics of :cite:`Bitz99`
|    **icepack_therm_itd.F90**     thermodynamic changes mostly related to ice thickness distribution
|    **icepack_therm_mushy.F90**   mushy-theory thermodynamics of :cite:`Turner13`
|    **icepack_therm_shared.F90**  code shared by all thermodynamics parameterizations
|    **icepack_therm_vertical.F90**  vertical growth rates and fluxes
|    **icepack_tracers.F90**       tracer information
|    **icepack_warnings.F90**      utilities for writing warning and error messages
|    **icepack_wavefracspec.F90**  wave impact on sea ice
|    **icepack_zbgc.F90**          driver for ice biogeochemistry and brine tracer motion
|    **icepack_zbgc_shared.F90**   parameters and shared code for biogeochemistry and brine height
|    **icepack_zsalinity.F90**     vertical salinity parameterization of :cite:`Jeffery11`


Coding Standard
------------------------------------

The column physics is a library that solves the sea ice column physics on a 
timestep by timestep and gridpoint by gridpoint basis.  It consists of Fortran routines with 
input and output arguments.  The model state is saved in the host model.  There is no 
communication between gridcells so the underlying implementation
supports no parallelization.  It however can be called in parallel by a driver
that is running on multiple pes with a decomposed grid.

The column physics does not have a time manager.  Calendaring is expected to be
dealt with by the host model.  The column physics does not read any forcing data,
that is passed into the column physics though interfaces.  In fact, 
there are no direct IO capabilities in the column physics.  That is to say, the
column physics does not open files to read or write.  The column physics is able to write 
data via several different routines that specifically have a fortran unit number as an input
argument.  In addition, there is a warning and abort package (see section :ref:`aborts`) that
provides the column package with the ability to store log output.  That output can
be queried by the host model or it can be written directly via a specific routine.
The warning package also provides access to an abort flag that the host model can
query after each call to check for successful completion of the column physics package.

All column physics public interfaces and public data are defined in the **icepack_intfc.F90**
file (see section :ref:`calling`).  Internal column physics settings should all be accessible through interfaces.
The internal constants, parameters, and tracer settings have init (set), query (get), and
write interfaces that provides access to internal column physics settings.  The host model
should not have to use "use" statements to access any of the column physics data outside
of what is provided through the icepack_intfc module.  
The public column physics interfaces use optional arguments where it makes sense and
there is an ongoing effort to make more of the interfaces support keyword=value arguments
for clarity and backwards compatibility.

