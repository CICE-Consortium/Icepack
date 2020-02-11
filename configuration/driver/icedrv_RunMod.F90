!=======================================================================
!
!  Main driver for time stepping of Icepack
!
!  authors Elizabeth C. Hunke, LANL

      module icedrv_RunMod

      use icedrv_kinds
      use icedrv_constants, only: c0, c1, nu_diag
      use icepack_intfc, only: icepack_warnings_flush
      use icepack_intfc, only: icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_sizes
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private
      public :: icedrv_run, ice_step

!=======================================================================

      contains

!=======================================================================
!
!  This is the main driver routine for advancing CICE forward in time.
!
!  author Elizabeth C. Hunke, LANL

      subroutine icedrv_run

      use icedrv_calendar, only: istep, istep1, time, dt, stop_now, calendar
      use icedrv_forcing, only: get_forcing, get_wave_spec
      use icedrv_forcing_bgc, only: faero_default, get_forcing_bgc
      use icedrv_flux, only: init_flux_atm_ocn

      logical (kind=log_kind) :: skl_bgc, z_tracers, tr_aero, tr_zaero, &
                                 wave_spec, tr_fsd

      character(len=*), parameter :: subname='(icedrv_run)'

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_zaero_out=tr_zaero, &
                                      tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      timeLoop: do

         call ice_step

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date

         call calendar(time)    ! at the end of the timestep

         if (stop_now >= 1) exit timeLoop

         call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers,&
                                       wave_spec_out=wave_spec)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

         if (tr_fsd .and. wave_spec) call get_wave_spec ! wave spectrum in ice
         call get_forcing(istep1)  ! get forcing from data arrays

         ! aerosols
         if (tr_aero .or. tr_zaero)  call faero_default    ! default values

         if (skl_bgc .or. z_tracers) call get_forcing_bgc  ! biogeochemistry

         call init_flux_atm_ocn ! initialize atmosphere, ocean fluxes

      enddo timeLoop

      end subroutine icedrv_run

!=======================================================================
!
!  Calls drivers for physics components, some initialization, and output
!
!  author Elizabeth C. Hunke, LANL

      subroutine ice_step

      use icedrv_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep
      use icedrv_diagnostics, only: runtime_diags, init_mass_diags
!     use icedrv_diagnostics, only: icedrv_diagnostics_debug
      use icedrv_diagnostics_bgc, only: hbrine_diags, zsal_diags, bgc_diags
      use icedrv_flux, only: init_history_therm, init_history_bgc, &
          daidtt, daidtd, dvidtt, dvidtd, dagedtt, dagedtd, init_history_dyn
      use icedrv_restart, only: dumpfile, final_restart
      use icedrv_restart_bgc, only: write_restart_bgc
      use icedrv_step, only: prep_radiation, step_therm1, step_therm2, &
          update_state, step_dyn_ridge, step_radiation, &
          biogeochemistry, step_dyn_wave

      integer (kind=int_kind) :: &
         k               ! dynamics supercycling index

      logical (kind=log_kind) :: &
         calc_Tsfc, skl_bgc, solve_zsal, z_tracers, tr_brine, &  ! from icepack
         tr_fsd, wave_spec

      real (kind=dbl_kind) :: &
         offset          ! d(age)/dt time offset

      character(len=*), parameter :: subname='(ice_step)'

!      call icedrv_diagnostics_debug ('beginning time step')

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_query_parameters(solve_zsal_out=solve_zsal, & 
                                    calc_Tsfc_out=calc_Tsfc, &
                                    wave_spec_out=wave_spec)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine,tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

      call init_mass_diags   ! diagnostics per timestep
      call init_history_therm
      call init_history_bgc

      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------
      
      if (calc_Tsfc) call prep_radiation ()

!      call icedrv_diagnostics_debug ('post prep_radiation')

      !-----------------------------------------------------------------
      ! thermodynamics and biogeochemistry
      !-----------------------------------------------------------------

      call step_therm1     (dt) ! vertical thermodynamics
      call biogeochemistry (dt) ! biogeochemistry
      call step_therm2     (dt) ! ice thickness distribution thermo

      ! clean up, update tendency diagnostics
      offset = dt
      call update_state (dt, daidtt, dvidtt, dagedtt, offset)

!      call icedrv_diagnostics_debug ('post thermo')

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------
      
      call init_history_dyn
      
      ! wave fracture of the floe size distribution
      ! note this is called outside of the dynamics subcycling loop
      if (tr_fsd .and. wave_spec) call step_dyn_wave(dt)

      do k = 1, ndtd
        
        ! ridging
        call step_dyn_ridge (dt_dyn, ndtd)
        
        ! clean up, update tendency diagnostics
        offset = c0
        call update_state (dt_dyn, daidtd, dvidtd, dagedtd, offset)
        
      enddo

!      call icedrv_diagnostics_debug ('post dynamics')
      
      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------
      
      call step_radiation (dt)

      !-----------------------------------------------------------------
      ! get ready for coupling and the next time step
      !-----------------------------------------------------------------
      
      call coupling_prep

!      call icedrv_diagnostics_debug ('post step_rad, cpl')

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------
      
      if (mod(istep,diagfreq) == 0) then
         call runtime_diags(dt)       ! log file
         if (solve_zsal)              call zsal_diags
         if (skl_bgc .or. z_tracers)  call bgc_diags
         if (tr_brine)                call hbrine_diags
      endif
      
      if (write_restart == 1) then
         call dumpfile     ! core variables for restarting
         if (solve_zsal .or. skl_bgc .or. z_tracers) &
            call write_restart_bgc         ! biogeochemistry
         call final_restart
      endif
      
    end subroutine ice_step
    
!=======================================================================
!
! Prepare for coupling
!
! authors: Elizabeth C. Hunke, LANL

      subroutine coupling_prep

      use icedrv_arrays_column, only: alvdfn, alidfn, alvdrn, alidrn, &
          albicen, albsnon, albpndn, apeffn, fzsal_g, fzsal, snowfracn
      use icedrv_calendar, only: dt
      use icedrv_domain_size, only: ncat, nx
      use icedrv_flux, only: alvdf, alidf, alvdr, alidr, albice, albsno, &
          albpnd, apeff_ai, fpond, fresh, l_mpond_fresh, &
          alvdf_ai, alidf_ai, alvdr_ai, alidr_ai, fhocn_ai, &
          fresh_ai, fsalt_ai, fsalt, &
          fswthru_ai, fhocn, fswthru, scale_factor, snowfrac, &
          swvdr, swidr, swvdf, swidf, &
          frzmlt_init, frzmlt, &
          fzsal_ai, fzsal_g_ai, flux_bio, flux_bio_ai
      use icedrv_forcing, only: oceanmixed_ice
      use icedrv_state, only: aicen
      use icedrv_step, only: ocean_mixed_layer

      ! local variables

      integer (kind=int_kind) :: & 
         n           , & ! thickness category index
         i           , & ! horizontal index
         k           , & ! tracer index
         nbtrcr

      real (kind=dbl_kind) :: &
         netsw, &        ! flag for shortwave radiation presence
         rhofresh, &     !
         puny            !

      character(len=*), parameter :: subname='(coupling_prep)'

      !-----------------------------------------------------------------
      ! Save current value of frzmlt for diagnostics.
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         call icepack_query_parameters(puny_out=puny, rhofresh_out=rhofresh)
         call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

         do i = 1, nx
            frzmlt_init  (i) = frzmlt(i)
         enddo

         if (oceanmixed_ice) &
         call ocean_mixed_layer (dt) ! ocean surface fluxes and sst

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do i = 1, nx
            alvdf(i) = c0
            alidf(i) = c0
            alvdr(i) = c0
            alidr(i) = c0

            albice(i) = c0
            albsno(i) = c0
            albpnd(i) = c0
            apeff_ai(i) = c0
            snowfrac(i) = c0
         enddo
         do n = 1, ncat
         do i = 1, nx
            if (aicen(i,n) > puny) then
                  
            alvdf(i) = alvdf(i) + alvdfn(i,n)*aicen(i,n)
            alidf(i) = alidf(i) + alidfn(i,n)*aicen(i,n)
            alvdr(i) = alvdr(i) + alvdrn(i,n)*aicen(i,n)
            alidr(i) = alidr(i) + alidrn(i,n)*aicen(i,n)

            netsw = swvdr(i) + swidr(i) + swvdf(i) + swidf(i)
            if (netsw > puny) then ! sun above horizon
               albice(i) = albice(i) + albicen(i,n)*aicen(i,n)
               albsno(i) = albsno(i) + albsnon(i,n)*aicen(i,n)
               albpnd(i) = albpnd(i) + albpndn(i,n)*aicen(i,n)
            endif

            apeff_ai(i) = apeff_ai(i) + apeffn(i,n)*aicen(i,n) ! for history
            snowfrac(i) = snowfrac(i) + snowfracn(i,n)*aicen(i,n) ! for history
               
            endif ! aicen > puny
         enddo
         enddo

         do i = 1, nx

      !-----------------------------------------------------------------
      ! reduce fresh by fpond for coupling
      !-----------------------------------------------------------------

            if (l_mpond_fresh) then
               fpond(i) = fpond(i) * rhofresh/dt
               fresh(i) = fresh(i) - fpond(i)
            endif

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i) = alvdf  (i)
            alidf_ai  (i) = alidf  (i)
            alvdr_ai  (i) = alvdr  (i)
            alidr_ai  (i) = alidr  (i)
            fresh_ai  (i) = fresh  (i)
            fsalt_ai  (i) = fsalt  (i)
            fhocn_ai  (i) = fhocn  (i)
            fswthru_ai(i) = fswthru(i)
            fzsal_ai  (i) = fzsal  (i) 
            fzsal_g_ai(i) = fzsal_g(i)  

            if (nbtrcr > 0) then
            do k = 1, nbtrcr
               flux_bio_ai  (i,k) = flux_bio  (i,k)
            enddo
            endif

      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i) = &
                       swvdr(i)*(c1 - alvdr_ai(i)) &
                     + swvdf(i)*(c1 - alvdf_ai(i)) &
                     + swidr(i)*(c1 - alidr_ai(i)) &
                     + swidf(i)*(c1 - alidf_ai(i))
         enddo

      end subroutine coupling_prep

!=======================================================================

      end module icedrv_RunMod

!=======================================================================
