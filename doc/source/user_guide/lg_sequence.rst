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

The full suite of public parameters and interfaces is documented in :ref:`icepack_intfc.F90`

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

