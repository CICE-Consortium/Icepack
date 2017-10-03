!=======================================================================
!
!  This module contains the Icepack initialization routine that sets model
!  parameters and initializes the grid and state variables.
!
!  authors Elizabeth C. Hunke, LANL

      module icepack_drv_InitMod

      use icepack_kinds_mod

      implicit none
      private
      public :: icepack_initialize
      save

!=======================================================================

      contains

!=======================================================================


!  Initialize Icepack

      subroutine icepack_initialize

      use icepack_drv_arrays_column, only: hin_max, c_hi_range, zfswin, trcrn_sw, &
          ocean_bio_all, ice_bio_net, snow_bio_net
      use icepack_drv_calendar, only: dt, dt_dyn, time, istep, istep1, write_ic, &
          init_calendar, calendar
      use icepack_intfc, only: icepack_init_itd, icepack_init_itd_hist, &
          icepack_clear_warnings, icepack_print_warnings
      use icepack_drv_domain_size, only: ncat
      use icepack_drv_constants, only: nu_diag
      use icepack_drv_flux, only: init_coupler_flux, init_history_therm, &
          init_history_dyn, init_flux_atm_ocn
      use icepack_drv_forcing, only: init_forcing, get_forcing
      use icepack_drv_forcing_bgc, only: get_forcing_bgc, faero_default !, get_atm_bgc, &
!          faero_data, faero_optics
      use icepack_drv_restart_shared, only: restart
      use icepack_drv_init, only: input_data, init_state, init_grid2
      use icepack_drv_init_column, only: init_thermo_vertical, init_shortwave, init_zbgc
      use icepack_intfc_tracers, only: tr_aero, tr_zaero
      use icepack_intfc_shared, only: skl_bgc, z_tracers

      logical(kind=log_kind) :: l_stop
      character(char_len) :: stop_label

      call input_data           ! namelist variables
      call init_zbgc            ! vertical biogeochemistry namelist

      call init_grid2           ! grid variables

      call init_calendar        ! initialize some calendar stuff

      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      call icepack_init_itd(ncat, hin_max, l_stop, stop_label) ! ice thickness distribution
      if (l_stop) then
         write(nu_diag,*) trim(stop_label)
         stop
      endif

      call icepack_clear_warnings()
      call icepack_init_itd_hist(ncat, hin_max, c_hi_range) ! output
      call icepack_print_warnings(nu_diag)

      call calendar(time)       ! determine the initial date

      call init_state           ! initialize the ice state
      call init_restart         ! initialize restart variables
      call init_history_therm   ! initialize thermo history variables

!      if (tr_aero .or. tr_zaero) call faero_optics !initialize aerosol optical 
                                                   !property tables

      istep  = istep  + 1    ! update time step counters
      istep1 = istep1 + 1
      time = time + dt       ! determine the time and date
      call calendar(time)    ! at the end of the first timestep

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

      call init_forcing      ! initialize forcing (standalone)
!?      call init_coupler_flux ! complete forcing initialization
      call get_forcing(istep1)       ! get forcing from data arrays

!      call get_forcing_atmo     ! atmospheric forcing from data
!      call get_forcing_ocn(dt)  ! ocean forcing from data

      ! aerosols
      ! if (tr_aero)  call faero_data                   ! data file
      ! if (tr_zaero) call fzaero_data                  ! data file (gx1)
      if (tr_aero .or. tr_zaero)  call faero_default    ! default values

      if (skl_bgc .or. z_tracers) call get_forcing_bgc  ! biogeochemistry
!      if (z_tracers) call get_atm_bgc                   ! biogeochemistry

      call init_shortwave    ! initialize radiative transfer using current swdn

      call init_flux_atm_ocn    ! initialize atmosphere, ocean fluxes

      end subroutine icepack_initialize

!=======================================================================

      subroutine init_restart

      use icepack_drv_arrays_column, only: dhsn
      use icepack_drv_calendar, only: time, calendar
      use icepack_constants, only: c0
      use icepack_intfc, only: icepack_aggregate
      use icepack_drv_domain_size, only: ncat, max_ntrcr, n_aero, nx
      use icepack_drv_flux, only: sss
      use icepack_drv_init, only: ice_ic, tmask
      use icepack_drv_init_column, only: init_hbrine, init_bgc
      use icepack_drv_restart, only: restartfile, read_restart_hbrine
      use icepack_drv_restart_shared, only: restart
      use icepack_drv_state ! almost everything
      use icepack_intfc_tracers, only: tr_iage, tr_FY, tr_lvl, nt_alvl, nt_vlvl, &
          tr_pond_cesm, nt_apnd, nt_hpnd, tr_pond_lvl, nt_ipnd, &
          tr_pond_topo, tr_aero, tr_brine, nt_iage, nt_FY, nt_aero
      use icepack_intfc_shared, only: skl_bgc, z_tracers, solve_zsal

      integer(kind=int_kind) :: &
         i                            ! horizontal indices

      if (restart) then
         call restartfile (ice_ic)
      endif      



      !in CICE, the following line:
      !if (tr_brine .or. skl_bgc) call init_hbrine ! brine height tracer
      !is called like this:
      if (tr_brine .or. skl_bgc) then ! brine height tracer
        call init_hbrine
        !if (tr_brine .and. restart_hbrine) call read_restart_hbrine
        if (tr_brine .and. restart) call read_restart_hbrine
      endif



      !the bgc restarts are contained in this subroutine
      !if (solve_zsal .or. skl_bgc .or. z_tracers) call init_bgc ! biogeochemistry

      !-----------------------------------------------------------------
      ! aggregate tracers
      !-----------------------------------------------------------------
      do i = 1, nx
         if (tmask(i)) &
         call icepack_aggregate (ncat,               &
                                aicen(i,:),  &
                                trcrn(i,:,:),&
                                vicen(i,:),  &
                                vsnon(i,:),  &
                                aice (i),  &
                                trcr (i,:),  &
                                vice (i),  &
                                vsno (i),  &
                                aice0(i),  &
                                max_ntrcr,          &
                                trcr_depend,        &
                                trcr_base,          &
                                n_trcr_strata,      &
                                nt_strata)
      enddo
      end subroutine init_restart

!=======================================================================

      end module icepack_drv_InitMod

!=======================================================================
