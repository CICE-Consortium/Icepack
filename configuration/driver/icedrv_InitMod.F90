!=======================================================================
!
!  This module contains the Icepack initialization routine that sets model
!  parameters and initializes the grid and state variables.
!
!  authors Elizabeth C. Hunke, LANL

      module icedrv_InitMod

      use icedrv_kinds
      use icedrv_constants, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_flags
      use icepack_intfc, only: icepack_write_tracer_flags, icepack_write_tracer_indices
      use icepack_intfc, only: icepack_write_tracer_sizes
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private
      public :: icedrv_initialize

!=======================================================================

      contains

!=======================================================================


!  Initialize Icepack

      subroutine icedrv_initialize

      use icedrv_arrays_column, only: hin_max, c_hi_range
      use icedrv_arrays_column, only: floe_rad_l, floe_rad_c, &
          floe_binwidth, c_fsd_range
      use icedrv_calendar, only: dt, time, istep, istep1, &
          init_calendar, calendar
      use icepack_intfc, only: icepack_init_itd, icepack_init_itd_hist
      use icepack_intfc, only: icepack_init_fsd_bounds
      use icepack_intfc, only: icepack_warnings_flush
      use icedrv_domain_size, only: ncat, nfsd
!     use icedrv_diagnostics, only: icedrv_diagnostics_debug
      use icedrv_flux, only: init_coupler_flux, init_history_therm, &
          init_flux_atm_ocn
      use icedrv_forcing, only: init_forcing, get_forcing, get_wave_spec
      use icedrv_forcing_bgc, only: get_forcing_bgc, faero_default, init_forcing_bgc 
      use icedrv_restart_shared, only: restart
      use icedrv_init, only: input_data, init_state, init_grid2, init_fsd
      use icedrv_init_column, only: init_thermo_vertical, init_shortwave, init_zbgc
      use icepack_intfc, only: icepack_configure

      logical (kind=log_kind) :: &
         skl_bgc, &    ! from icepack
         z_tracers, &  ! from icepack
         tr_aero, &    ! from icepack
         tr_zaero, &   ! from icepack
         tr_fsd, wave_spec

      character(len=*), parameter :: subname='(icedrv_initialize)'

      call icepack_configure()  ! initialize icepack
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call input_data           ! namelist variables
      call init_zbgc            ! vertical biogeochemistry namelist

      ! generate some output
      call icepack_write_tracer_flags(nu_diag)
      call icepack_write_tracer_sizes(nu_diag)
      call icepack_write_tracer_indices(nu_diag)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call init_grid2           ! grid variables
      call init_calendar        ! initialize some calendar stuff
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      call icepack_init_itd(ncat=ncat, hin_max=hin_max)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted(subname)) then
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      call icepack_init_itd_hist(ncat=ncat, c_hi_range=c_hi_range, hin_max=hin_max) ! output

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted(subname)) then
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (tr_fsd) call icepack_init_fsd_bounds(   &
         nfsd=nfsd,                   &  ! floe size distribution
         floe_rad_l=floe_rad_l,       &  ! fsd size lower bound in m (radius)
         floe_rad_c=floe_rad_c,       &  ! fsd size bin centre in m (radius)
         floe_binwidth=floe_binwidth, &  ! fsd size bin width in m (radius)
         c_fsd_range=c_fsd_range)        ! string for history output
      call init_fsd

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted(subname)) then
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      call calendar(time)       ! determine the initial date

      call init_state           ! initialize the ice state
      call init_restart         ! initialize restart variables
      call init_history_therm   ! initialize thermo history variables

      if (restart) &
         call init_shortwave    ! initialize radiative transfer

      istep  = istep  + 1    ! update time step counters
      istep1 = istep1 + 1
      time = time + dt       ! determine the time and date
      call calendar(time)    ! at the end of the first timestep

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------
      call icepack_query_parameters(skl_bgc_out=skl_bgc)
      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_query_parameters(wave_spec_out=wave_spec)
      call icepack_query_tracer_flags(tr_aero_out=tr_aero)
      call icepack_query_tracer_flags(tr_zaero_out=tr_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call init_forcing      ! initialize forcing (standalone)     
      if (skl_bgc .or. z_tracers) call init_forcing_bgc !cn
      if (tr_fsd .and. wave_spec) call get_wave_spec ! wave spectrum in ice
      call get_forcing(istep1)       ! get forcing from data arrays

      ! aerosols
      ! if (tr_aero)  call faero_data                   ! data file
      ! if (tr_zaero) call fzaero_data                  ! data file (gx1)
      if (tr_aero .or. tr_zaero)  call faero_default    ! default values
      if (skl_bgc .or. z_tracers) call get_forcing_bgc  ! biogeochemistry

      if (.not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

      call init_flux_atm_ocn    ! initialize atmosphere, ocean fluxes

      end subroutine icedrv_initialize

!=======================================================================

      subroutine init_restart

      use icedrv_calendar, only: time, calendar
      use icepack_intfc, only: icepack_aggregate
      use icedrv_domain_size, only: ncat, max_ntrcr, nx
      use icedrv_init, only: ice_ic
      use icedrv_init, only: tmask
      use icedrv_init_column, only: init_hbrine, init_bgc
      use icedrv_restart, only: restartfile
      use icedrv_restart_shared, only: restart
      use icedrv_restart_bgc, only: read_restart_bgc
      use icedrv_state ! almost everything

      integer(kind=int_kind) :: &
         i                            ! horizontal indices

      logical (kind=log_kind) :: &
         skl_bgc, &    ! from icepack
         z_tracers, &  ! from icepack
         solve_zsal, & ! from icepack
         tr_brine, &   ! from icepack
         tr_fsd        ! from icepack

      character(len=*), parameter :: subname='(init_restart)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(skl_bgc_out=skl_bgc)
      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_query_parameters(solve_zsal_out=solve_zsal)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------

      if (tr_brine .or. skl_bgc) then ! brine height tracer
        call init_hbrine
      endif

      if (restart) then
         call restartfile (ice_ic)
         call calendar (time)
      endif

      if (solve_zsal .or. skl_bgc .or. z_tracers) then
        if (tr_fsd) then
            write (nu_diag,*) 'FSD implementation incomplete for use with BGC'
            call icedrv_system_abort(string=subname,file=__FILE__,line=__LINE__)
         endif
         call init_bgc
         if (restart) call read_restart_bgc ! complete BGC initialization
      endif

      !-----------------------------------------------------------------
      ! aggregate tracers
      !-----------------------------------------------------------------
      do i = 1, nx
         if (tmask(i)) &
         call icepack_aggregate(ncat=ncat,          &
                                aicen=aicen(i,:),   &
                                vicen=vicen(i,:),   &
                                vsnon=vsnon(i,:),   &
                                trcrn=trcrn(i,:,:), &
                                aice=aice (i),      &
                                vice=vice (i),      &
                                vsno=vsno (i),      &
                                trcr=trcr (i,:),    &
                                aice0=aice0(i),     &
                                ntrcr=max_ntrcr,    &
                                trcr_depend=trcr_depend, &
                                trcr_base=trcr_base,     &
                                n_trcr_strata=n_trcr_strata, &
                                nt_strata=nt_strata)
      enddo
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      end subroutine init_restart

!=======================================================================

      end module icedrv_InitMod

!=======================================================================
