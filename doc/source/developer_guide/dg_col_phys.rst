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
argument.  In addition, there is a warning package (see section :ref:`warning`) that
provides the column package with the ability to store log output.  That output can
be queried by the host model or it can be written directly via a specific routine.
The warning package also provides access to an abort flag that the host model can
query after each call to check for successful completion of the column physics package.

All column physics public interfaces and public data are defined in the **icepack_intfc.F90**
file.  Internal column physics settings should all be accessible through interfaces.
The internal constants, parameters, and tracer settings have init (set), query (get), and
write interfaces that provides access to internal column physics settings.  The host model
should not have to use "use" statements to access any of the column physics data outside
of what is provided through the icepack_intfc module.  
The public column physics interfaces use optional arguments where it makes sense and
there is an ongoing effort to make more of the interfaces support keyword=value arguments
for clarity and backwards compatibility.


Using Icepack
------------------------------------

In this section, the various public icepack interfaces will be defined and 
how to use them will be described.

.. dev_intfc:

Interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Column physics data and subroutines are made public through the **icepack_intfc.F90**
file.  That file contains the entire list of data and subroutines needed to
initialize, setup, and run the column physics package.  That file points
to other modules within the column physics where the interfaces are located.

Within **icepack_intfc.F90**, internal icepack kinds are defined via the
icepack_kinds module::

      use icepack_kinds, only: icepack_char_len  => char_len
      use icepack_kinds, only: icepack_char_len_long  => char_len_long
      use icepack_kinds, only: icepack_log_kind  => log_kind
      use icepack_kinds, only: icepack_int_kind  => int_kind
      use icepack_kinds, only: icepack_real_kind => real_kind
      use icepack_kinds, only: icepack_dbl_kind  => dbl_kind
      use icepack_kinds, only: icepack_r16_kind  => r16_kind

icepack_tracers defines a handful of parameters that provide information
about maximum array sizes for static dimensioning::

      use icepack_tracers,   only: icepack_max_nbtrcr => max_nbtrcr
      use icepack_tracers,   only: icepack_max_algae  => max_algae
      use icepack_tracers,   only: icepack_max_dic    => max_dic
      use icepack_tracers,   only: icepack_max_doc    => max_doc
      use icepack_tracers,   only: icepack_max_don    => max_don
      use icepack_tracers,   only: icepack_max_fe     => max_fe
      use icepack_tracers,   only: icepack_max_aero   => max_aero
      use icepack_tracers,   only: icepack_nmodal1    => nmodal1
      use icepack_tracers,   only: icepack_nmodal2    => nmodal2
      use icepack_parameters, only: icepack_nspint     => nspint

icepack_parameters provides init, query, write, and recompute methods to
define constant values and model parameters.  These constants have defaults 
that the caller can query or reset::

      use icepack_parameters, only: icepack_init_parameters
      use icepack_parameters, only: icepack_query_parameters
      use icepack_parameters, only: icepack_write_parameters
      use icepack_parameters, only: icepack_recompute_constants

icepack_tracers provides init, query, and write methods to
define various tracer sizes, flags, indices, and numbers.  The
tracers have some defaults that the caller can query or reset::

      use icepack_tracers, only: icepack_compute_tracers
      use icepack_tracers, only: icepack_query_tracer_sizes
      use icepack_tracers, only: icepack_write_tracer_sizes
      use icepack_tracers, only: icepack_init_tracer_flags
      use icepack_tracers, only: icepack_query_tracer_flags
      use icepack_tracers, only: icepack_write_tracer_flags
      use icepack_tracers, only: icepack_init_tracer_indices
      use icepack_tracers, only: icepack_query_tracer_indices
      use icepack_tracers, only: icepack_write_tracer_indices
      use icepack_tracers, only: icepack_init_tracer_numbers
      use icepack_tracers, only: icepack_query_tracer_numbers
      use icepack_tracers, only: icepack_write_tracer_numbers

icepack_itd provides three public interfaces to compute the ice
thickness distribution::

      use icepack_itd, only: icepack_init_itd
      use icepack_itd, only: icepack_init_itd_hist
      use icepack_itd, only: icepack_aggregate

icepack_mechred contains two public interfaces to compute ridging
and ice strength::

      use icepack_mechred, only: icepack_step_ridge
      use icepack_mechred, only: icepack_ice_strength

icepack_shortwave provides a routine to initialize the radiation
computation and an routine to update the radiation computation::

      use icepack_shortwave, only: icepack_prep_radiation
      use icepack_shortwave, only: icepack_step_radiation

icepack_brine address brine and zsalinity computations::

      use icepack_brine, only: icepack_init_hbrine
      use icepack_brine, only: icepack_init_zsalinity

icepack_zbgc contains several public interfaces to support initialization
and computation for the skeletal layer bgc and zbgc options::

      use icepack_zbgc , only: icepack_init_bgc
      use icepack_zbgc , only: icepack_init_zbgc
      use icepack_zbgc , only: icepack_biogeochemistry
      use icepack_zbgc , only: icepack_init_OceanConcArray
      use icepack_zbgc , only: icepack_init_ocean_conc

There are a couple of routines to support computation of an atmosphere
and ocean interaction::

      use icepack_atmo , only: icepack_atm_boundary
      use icepack_ocean, only: icepack_ocn_mixed_layer

icepack_step_therm1 and icepack_step_therm2 compute the ice
thermodynamics in two steps::

      use icepack_therm_vertical, only: icepack_step_therm1
      use icepack_therm_itd     , only: icepack_step_therm2

icepack_therm_shared provides several methods to compute different
internal terms::

      use icepack_therm_shared  , only: icepack_ice_temperature
      use icepack_therm_shared  , only: icepack_snow_temperature
      use icepack_therm_shared  , only: icepack_liquidus_temperature
      use icepack_therm_shared  , only: icepack_sea_freezing_temperature
      use icepack_therm_shared  , only: icepack_enthalpy_snow
      use icepack_therm_shared  , only: icepack_init_thermo
      use icepack_therm_shared  , only: icepack_init_trcr

icepack_orbital provides a routine to set orbital parameters needed
for some albedo computations::

      use icepack_orbital , only: icepack_init_orbit

icepack_warnings provides several methods for getting, writing,
and clearing messages.  There is also a function that returns
a logical flag indicating whether the column physics has aborted::

      use icepack_warnings, only: icepack_warnings_clear
      use icepack_warnings, only: icepack_warnings_getall
      use icepack_warnings, only: icepack_warnings_print
      use icepack_warnings, only: icepack_warnings_flush
      use icepack_warnings, only: icepack_warnings_aborted

icepack_configure is a standalone icepack method that should always be called
first::

      public :: icepack_configure


Calling Sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The calling sequence required to setup and run the column physics is generally
described below.  Several steps may be needed to be taken by the host between
icepack calls in order to support the icepack interfaces.  
The icepack driver and the CICE model provide working examples
of how to do this in practice.  The sample below does not include bgc.

start driver

* call *icepack_configure*

initialize driver and read in driver namelist

* call *icepack_init_parameters*
* call *icepack_init_tracers_*
* call *icepack_init_trcr*
* call *icepack_init_thermo*
* call *icepack_init_itd*
* call *icepack_init_itd_hist*
* call *icepack_step_radiation*
* call *icepack_init_zsalinity*
* call *icepack_init_hbrine*
* call *icepack_aggregate*

loop over timesteps
loop over gridcells

* call *icepack_prep_radiation*
* call *icepack_step_therm1*
* call *icepack_step_therm2*
* call *icepack_aggregate*
* call *icepack_step_ridge*
* call *icepack_step_radiation*
* call *icepack_atm_boundary*
* call *icepack_ocn_mixed_layer*

end loop over gridcells
end loop over timesteps

end driver

.. _warning:

The Warning Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Icepack has no IO capabilities.  It does not have direct knowledge of
any input or output files.  However, it can write output through specific
interfaces that pass in a fortran file unit number.  There are several 
methods in icepack that support writing data to a file this way including
the various *icepack_write_* interfaces.

Separately, the icepack warning package is where icepack stores internal output and
error messages not directly set in the various write routines.  The warning package
also contains an *icepack_warnings_aborted* function that will be set to true 
if icepack detects an abort.  In that case, icepack will return to the driver.
As a general rule, after each call to icepack, the driver should call::

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

to flush (print and clear) the icepack warning buffer and to check whether icepack 
aborted.  If icepack aborts, it's actually up to the driver to cleanly shut the
model down.

Alternatively, *icepack_warnings_getall* provides the saved icepack messages to
the driver via an array of strings in the argument list.  This allows the driver
to reformat that output as needed.  *icepack_warnings_print*
writes out the messages but does not clear them, and *icepack_warnings_clear* zeros
out the icepack warning messages.
