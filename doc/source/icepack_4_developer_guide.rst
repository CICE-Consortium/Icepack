:tocdepth: 4

******************
Developers Guide
******************

The Icepack model consists of three different parts, the column physics
code, the icepack driver, and the scripts.  Development of each of these
pieces will be described below separately.


Subroutine calls and other linkages into Icepack from the host model should only
need to access the **icepack\_intfc\*.F90** interface modules within the 
``columnphysics/`` directory.  
The Icepack driver in the ``configuration/driver/`` directory is based on the CICE
model and provides an example of the sea ice host model capabilities needed for inclusion
of Icepack.  In particular, host models will need to include code equivalent to that
in the modules **icedrv\_\*_column.F90**.  Calls into the Icepack interface routines
are primarily from **icedrv\_step\_mod.F90** but there are others (search the driver code
for ``intfc``).

Guiding principles for the creation of Icepack include the following: 
CHECK THAT THESE ARE TRUE

- The column physics modules shall be independent of all sea ice model infrastructural
  elements that may vary from model to model.  Examples include input/output, timers,
  references to CPUs or computational tasks, initialization other than that necessary for
  strictly physical reasons, and anything related to a horizontal grid.
- The column physics modules shall not call or reference any routines or code that 
  reside outside of the **columnphysics/** directory.
- Any capabilities required by a host sea ice model (e.g. calendar variables, tracer 
  flags, diagnostics) shall be implemented in the driver and passed into or out of the 
  column physics modules via array arguments.


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
|    **icepack_constants.F90**     physical and numerical constants required for column package
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
|    **icepack_parameters.F90**    basic model parameters
|    **icepack_shortwave.F90**     shortwave and albedo parameterizations
|    **icepack_therm_0layer.F90**  zero-layer thermodynamics of :cite:`Semtner76`
|    **icepack_therm_bl99.F90**    multilayer thermodynamics of :cite:`BL99`
|    **icepack_therm_itd.F90**     thermodynamic changes mostly related to ice thickness distribution
|    **icepack_therm_mushy.F90**   mushy-theory thermodynamics of :cite:`THB13`
|    **icepack_therm_shared.F90**  code shared by all thermodynamics parameterizations
|    **icepack_therm_vertical.F90**  vertical growth rates and fluxes
|    **icepack_tracers.F90**       tracer information
|    **icepack_warnings.F90**      utilities for writing warning and error messages
|    **icepack_zbgc.F90**          driver for ice biogeochemistry and brine tracer motion
|    **icepack_zbgc_shared.F90**   parameters and shared code for biogeochemistry and brine height
|    **icepack_zsalinity.F90**     vertical salinity parameterization of :cite:`JHE11`


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
that is passed into the column physics thru interfaces.  In fact, 
there are no direct IO capabilities in the column physics.  That is to say, the
column physics does not open files to read or write.  The column physics is able to write 
data via several different routines that specifically have a fortran unit number as an input
argument.  In addition, there is a warning package (see section :ref:`warning`) that
provides the column package with the ability to store log output.  That output can
be queried by the host model or it can be written directly via a specific routine.
The warning package also provides access to an abort flag that the host model can
query after each call to check for successful completion of the column physics package.

All column physics public interfaces and public data are defined in the **icepack_intfc.F90**
file.  Internal column physics settings should all be accessible thru interfaces.
The internal constants, parameters, and tracer settings have init (set), query (get), and
write interfaces that provides access to internal column physics settings.  The host model
should not have to use "use" statements to access any of the column physics data outside
of what is provided thru the icepack_intfc module.  
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

Column physics data and subroutines are made public thru the **icepack_intfc.F90**
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

icepack_tracers defines a handful of parameters constants that provide information
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
      use icepack_constants, only: icepack_nspint     => nspint

icepack_constants provides a list of static parameter constants::

      use icepack_constants, only: c0 => c0

icepack_constants provides init, query, write, and recompute methods to
define constant values.  These constants have defaults that the caller
can query or reset::

      use icepack_constants, only: icepack_init_constants
      use icepack_constants, only: icepack_query_constants
      use icepack_constants, only: icepack_write_constants
      use icepack_constants, only: icepack_recompute_constants

icepack_parameters provides init, query, and write methods to
define model parameters.  These parameters have defaults that the caller
can query or reset::

      use icepack_parameters, only: icepack_init_parameters
      use icepack_parameters, only: icepack_query_parameters
      use icepack_parameters, only: icepack_write_parameters

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

* call *icepack_init_constants*
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
any input or output files.  However, it can write output thru specific
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



.. _dev_driver:

Driver Implementation
========================

The icepack driver is Fortran source code and exists to test the column physics
in a stand-alone mode for some simple column configurations.

File List
-------------------

The icepack driver consists of the following files 

|  **configuration/driver/**       driver for testing Icepack in stand-alone mode
|        **icedrv_MAIN.F90**        main program
|        **icedrv_InitMod.F90**     routines for initializing a run
|        **icedrv_RunMod.F90**      main driver routines for time stepping
|        **icedrv_arrays_column.F90**    essential arrays to describe the state of the ice
|        **icedrv_calendar.F90**    keeps track of what time it is
|        **icedrv_constants.F90**   physical and numerical constants and parameters
|        **icedrv_diagnostics.F90** miscellaneous diagnostic and debugging routines
|        **icedrv_diagnostics_bgc.F90**  diagnostic routines for biogeochemistry
|        **icedrv_domain_size.F90** domain sizes
|        **icedrv_flux.F90**        fluxes needed/produced by the model
|        **icedrv_forcing.F90**     routines to read and interpolate forcing data for stand-alone model runs
|        **icedrv_forcing_bgc.F90** routines to read and interpolate forcing data for bgc stand-alone model runs
|        **icedrv_init.F90**        general initialization routines
|        **icedrv_init_column.F90** initialization routines specific to the column physics
|        **icedrv_restart.F90**     driver for reading/writing restart files
|        **icedrv_restart_column.F90**  (CHECK: RENAME bgc) restart routines specific to the column physics
|        **icedrv_restart_shared.F90**  code shared by all restart options
|        **icedrv_state.F90**       essential arrays to describe the state of the ice
|        **icedrv_step.F90**        routines for time stepping the major code components
|        **icedrv_system.F90**      overall system management calls

Overview
------------

The icepack driver exists to test the column physics.  At the present time, it is hardwired
to run 4 different gridcells on one processor with the same forcing used for all gridcells.  
There is no MPI and no threading built into the icepack driver.  There is limited IO capabilities,
no history files, and no netcdf restart files.  The model generally runs very quickly.

There are a few different forcings available.


.. _dev_scripts:

Scripts Implementation
========================

The scripts are the third part of the icepack package.  They support setting up
cases, building, and running the icepack stand-alone model.

File List
--------------

The directory structure under configure/scripts is as follows.

| **configuration/scripts/**
|        **Makefile**              primary makefile
|        **icepack.batch.csh**     creates batch scripts for particular machines
|        **icepack.build**         compiles the code
|        **icepack.launch.csh**    creates script logic that runs the executable
|        **icepack.run.setup.csh** sets up the run scripts
|        **icepack.run.suite.csh** sets up the test suite
|        **icepack.settings**      defines environment, model configuration and run settings
|        **icepack.test.setup.csh**   creates configurations for testing the model
|        **icepack_decomp.csh**    defines the grid size
|        **icepack_in**            namelist input data
|        **machines/**             machine specific files to set env and Macros
|        **makdep.c**              determines module dependencies
|        **options/**              other namelist configurations available from the icepack.create.case command line
|        **parse_namelist.sh**     replaces namelist with command-line configuration
|        **parse_namelist_from_settings.sh**   replaces namelist with values from icepack.settings
|        **parse_settings.sh**     replaces settings with command-line configuration
|        **tests/**                scripts for configuring and running basic tests

:: _dev_strategy:

Strategy
-----------

The icepack scripts are implemented such that everything is resolved after
**icepack.create.case** is called.  This is done by both copying specific files
into the case directory and running scripts as part of the **icepack.create.case**
command line to setup various files.

**icepack.create.case** drives the case setup.  It is written in csh.  All supporting
scripts are relatively simple csh or sh scripts.

The file **icepack.settings** specifies a set of env defaults for the case.  The file
**icepack_in** defines the namelist input for the icepack driver.

:: _dev_options:

Preset Case Options
---------------------


``icepack.create.case -s`` option allows the user to choose some predetermined icepack
settings and namelist.  Those options are defined in **configurations/scripts/options/**
and the files are prefixed by either set_env, set_nml, or test_nml.  When **icepack.create.case**
is executed, the appropriate files are read from **configurations/scripts/options/**
and the **icepack.settings** and/or **icepack_in** files are updated in the case directory
based on the values in those files.

The filename suffix determines the name of the -s option.  So, for instance, 

  ``icepack.create.case -s diag1,debug,bgcISPOL``

will search for option files with suffixes of diag1, debug, and bgcISPOL and then
apply those settings.  

**parse_namelist.sh**, **parse_settings.sh**, and **parse_namelist_from_settings.sh** 
are the three scripts that modify **icepack_in** and **icepack.settings**.

To add new options, just add new files to the **configurations/scripts/options/** directory
with appropriate names and syntax.  The set_nml file syntax is the same as namelist
syntax and the set_env files are consistent with csh setenv syntax.  See other files for
examples of the syntax.

:: _dev_machines:

Machines
-----------

Machine specific information is contained in **configuration/scripts/machines**.  That
directory contains a Macros file and an env file for each supported machine.
One other files will need to be
changed to support a port, that is **configuration/scripts/icepack.batch.csh**.
To port to a new machine, see :ref:`porting`.  

:: _dev_testing:

Test scripts
-------------

Under **configuration/scripts/tests** are several files including the scripts to 
setup the smoke and restart tests (**test_smoke.script**, **test_restart.script*).
A baseline test script (**baseline.script**) is also there to setup the regression
and comparison testing.  That directory also contains the preset test suites 
(ie. **base_suite.ts**) and a file that supports post-processing on the model
output (**timeseries.csh**).  

There is a subdirectory, **configuration/scripts/tests/CTest**, that supports the
CTest scripts.  These scripts allow test reporting to CDash.

To add a new test, a file associated with that test will need to be added to the
**configuration/scripts/tests** directory similar to **test_smoke.script** 
and **test_restart.script**.  In addition, some new options files in 
**configuration/scripts/options** may need to be added similar to **test_nml.restart1**,
**test_nml.restart2**, and **set_nml.restart**.  

.. _addtrcr:

Adding tracers
====================

We require that any changes made to the code be implemented in such a way that they can
be "turned off" through namelist flags.  In most cases, code run with such changes should 
be bit-for-bit identical with the unmodified code.  Occasionally, non-bit-for-bit changes
are necessary, e.g. associated with an unavoidable change in the order of operations. In
these cases, changes should be made in stages to isolate the non-bit-for-bit changes, 
so that those that should be bit-for-bit can be tested separately.

Tracers added to Icepack will also require extensive modifications to the host
sea ice model, including initialization on the horizontal grid, namelist flags 
and restart capabilities.  Modifications to the Icepack driver should reflect
the modifications needed in the host model but are not expected to match completely.
We recommend that the logical namelist variable
``tr_[tracer]`` be used for all calls involving the new tracer outside of
**ice\_[tracer].F90**, in case other users do not want to use that
tracer.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, level-ice area and all of the
biogeochemistry tracers in this release are area-weighted tracers. In
the absence of sources and sinks, the total mass of a volume-weighted
tracer such as aerosol (kg) is conserved under transport in horizontal
and thickness space (the mass in a given grid cell will change), whereas
the aerosol concentration (kg/m) is unchanged following the motion, and
in particular, the concentration is unchanged when there is surface or
basal melting. The proper units for a volume-weighted mass tracer in the
tracer array are kg/m.

In several places in the code, tracer computations must be performed on
the conserved "tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies, and code must be included to account for new
dependencies (e.g., :math:`a_{lvl}` and :math:`a_{pnd}` in
**ice\_itd.F90** and **ice\_mechred.F90**).

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **icedrv\_domain\_size.F90**: increase ``max_ntrcr`` (can also add option
   to **icepack.settings** and **icepack.build**)

#. **icedrv\_state.F90**: declare ``nt_[tracer]`` and ``tr_[tracer]``

#. **icepack\_[tracer].F90**: create initialization, physics routines

#. **ice\_drv\_init.F90**: (some of this may be done in **ice\_[tracer].F90**
   instead)

   -  add new module and ``tr_[tracer]`` to list of used modules and
      variables

   -  add logical namelist variable ``tr_[tracer]``

   -  initialize namelist variable

   -  print namelist variable to diagnostic output file

   -  increment number of tracers in use based on namelist input (``ntrcr``)

   -  define tracer types (``trcr_depend`` = 0 for ice area tracers, 1 for
      ice volume, 2 for snow volume, 2+``nt_``[tracer] for dependence on
      other tracers)

#. **icepack\_itd.F90**, **icepack\_mechred.F90**: Account for new dependencies
   if needed.

#. **icedrv\_InitMod.F90**: initialize tracer (includes reading restart
   file)

#. **icedrv\_RunMod.F90**, **icedrv\_step\_mod.F90**:

   -  call routine to write tracer restart data

   -  call physics routines in **icepack\_[tracer].F90** (often called from
      **icedrv\_step\_mod.F90**)

#. **icedrv\_restart.F90**: define restart variables

#. **icepack\_in**: add namelist variables to *tracer\_nml* and
   *icefields\_nml*

#. If strict conservation is necessary, add diagnostics as noted for
   topo ponds in Section :ref:`ponds`.
