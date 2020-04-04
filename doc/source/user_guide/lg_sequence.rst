:tocdepth: 3

.. _sequence_and_interface:

Sequencing and Interfaces
---------------------------

Access to Interfaces
~~~~~~~~~~~~~~~~~~~~~~

Icepack public parameters and interfaces as accessed via a single module in
icepack, **icepack\_intfc.F90**.  The standard syntax to gain access to Icepack
parameters and interfaces is through Fortran90 use.  For example::

      use icepack_intfc, only: icepack_warnings_flush
      use icepack_intfc, only: icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_tracer_indices
      use icepack_intfc, only: icepack_configure

The full suite of public parameters and interfaces is documented in :ref:`icepack_intfc.F90`.

Interfaces Module
~~~~~~~~~~~~~~~~~~~~~~

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
      use icepack_kinds, only: icepack_int8_kind => int8_kind
      use icepack_kinds, only: icepack_real_kind => real_kind
      use icepack_kinds, only: icepack_dbl_kind  => dbl_kind
      use icepack_kinds, only: icepack_r16_kind  => r16_kind

icepack_tracers defines a handful of parameters that provide information
about maximum array sizes for static dimensioning::

      use icepack_tracers,   only: icepack_max_iso    => max_iso
      use icepack_tracers,   only: icepack_max_nbtrcr => max_nbtrcr
      use icepack_tracers,   only: icepack_max_algae  => max_algae
      use icepack_tracers,   only: icepack_max_dic    => max_dic
      use icepack_tracers,   only: icepack_max_doc    => max_doc
      use icepack_tracers,   only: icepack_max_don    => max_don
      use icepack_tracers,   only: icepack_max_fe     => max_fe
      use icepack_tracers,   only: icepack_max_aero   => max_aero
      use icepack_tracers,   only: icepack_nmodal1    => nmodal1
      use icepack_tracers,   only: icepack_nmodal2    => nmodal2
      use icepack_parameters,only: icepack_nspint     => nspint

icepack_parameters provides init, query, write, and recompute methods to
define constant values and model parameters.  These constants have defaults 
that the caller can query or reset::

      use icepack_parameters, only: icepack_init_parameters
      use icepack_parameters, only: icepack_query_parameters
      use icepack_parameters, only: icepack_write_parameters
      use icepack_parameters, only: icepack_recompute_constants

icepack_parameters also provides a set of constants::

      use icepack_parameters, only: c0, c1, c1p5, c2, c3, c4, c5, c6, c8
      use icepack_parameters, only: c10, c15, c16, c20, c25, c100, c1000
      use icepack_parameters, only: p001, p01, p1, p2, p4, p5, p6, p05
      use icepack_parameters, only: p15, p25, p75, p333, p666

icepack_tracers provides init, query, and write methods to
define various tracer sizes, flags, and indices.  The
tracers have some defaults that the caller can query or reset::

      use icepack_tracers, only: icepack_compute_tracers
      use icepack_tracers, only: icepack_init_tracer_flags
      use icepack_tracers, only: icepack_query_tracer_flags
      use icepack_tracers, only: icepack_write_tracer_flags
      use icepack_tracers, only: icepack_init_tracer_indices
      use icepack_tracers, only: icepack_query_tracer_indices
      use icepack_tracers, only: icepack_write_tracer_indices
      use icepack_tracers, only: icepack_init_tracer_sizes
      use icepack_tracers, only: icepack_query_tracer_sizes
      use icepack_tracers, only: icepack_write_tracer_sizes

icepack_itd provides three public interfaces to compute the ice
thickness distribution::

      use icepack_itd, only: icepack_init_itd
      use icepack_itd, only: icepack_init_itd_hist
      use icepack_itd, only: icepack_aggregate

icepack_fsd provides three public interfaces to compute the floe
size distribution::

      use icepack_fsd, only: icepack_init_fsd_bounds
      use icepack_fsd, only: icepack_init_fsd
      use icepack_fsd, only: icepack_cleanup_fsd

icepack_mechred contains two public interfaces to compute ridging
and ice strength::

      use icepack_mechred, only: icepack_step_ridge
      use icepack_mechred, only: icepack_ice_strength

icepack_wavefracspec provides two public interface to compute the
impact of waves on sea ice::

      use icepack_wavefracspec, only: icepack_init_wave
      use icepack_wavefracspec, only: icepack_step_wavefracture

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
      use icepack_zbgc , only: icepack_init_ocean_bio
      use icepack_zbgc , only: icepack_load_ocean_bio_array

There are a couple of routines to support computation of an atmosphere
and ocean interaction::

      use icepack_atmo , only: icepack_atm_boundary
      use icepack_ocean, only: icepack_ocn_mixed_layer

icepack_orbital provides methods to set and query orbital parameters::

      use icepack_orbital       , only: icepack_init_orbit
      use icepack_orbital       , only: icepack_query_orbit

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

icepack_mushy_physics provides three public interfaces to compute various
functions::

      use icepack_mushy_physics , only: icepack_mushy_density_brine
      use icepack_mushy_physics , only: icepack_mushy_liquid_fraction
      use icepack_mushy_physics , only: icepack_mushy_temperature_mush

icepack_warnings provides several methods for getting, writing,
and clearing messages.  There is also a function that returns
a logical flag indicating whether the column physics has aborted::

      use icepack_warnings, only: icepack_warnings_clear
      use icepack_warnings, only: icepack_warnings_print
      use icepack_warnings, only: icepack_warnings_flush
      use icepack_warnings, only: icepack_warnings_aborted

**icepack_configure** is a standalone icepack method that should always be called
first::

      public :: icepack_configure

The actual interfaces are documented in :ref:`docintfc`


.. _callingseq:

Calling Sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The calling sequence required to setup and run the column physics is generally
described below.  Several steps may be needed to be taken by the host between
icepack calls in order to support the icepack interfaces.  
The icepack driver and the CICE model provide working examples
of how to do this in practice.  The sample below does not include bgc::

  start driver

    call *icepack_configure*

  initialize driver and read in driver namelist

    call *icepack_init_parameters*
    call *icepack_init_tracers_*
    call *icepack_init_trcr*
    call *icepack_init_thermo*
    call *icepack_init_itd*
    call *icepack_init_itd_hist*
    loop over gridcells
      call *icepack_step_radiation*
      call *icepack_init_zsalinity*
    end loop over gridcells
    call *icepack_init_hbrine*
    loop over gridcells
       call *icepack_aggregate*
    end loop over gridcells

    loop over timesteps
      loop over gridcells
        call *icepack_prep_radiation*
        call *icepack_step_therm1*
        call *icepack_step_therm2*
        call *icepack_aggregate*
        call *icepack_step_ridge*
        call *icepack_step_radiation*
        call *icepack_atm_boundary*
        call *icepack_ocn_mixed_layer*
      end loop over gridcells
    end loop over timesteps

  end driver

