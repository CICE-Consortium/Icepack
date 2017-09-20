!=======================================================================
!
!  Main driver for time stepping of Icepack
!
!  authors Elizabeth C. Hunke, LANL

      module icepack_drv_RunMod

      use icepack_kinds_mod

      implicit none
      private
      public :: icepack_run, ice_step
      save

!=======================================================================

      contains

!=======================================================================
!
!  This is the main driver routine for advancing CICE forward in time.
!
!  author Elizabeth C. Hunke, LANL

      subroutine icepack_run

      use icepack_drv_calendar, only: istep, istep1, time, dt, stop_now, calendar
      use icepack_drv_forcing, only: get_forcing
!      use icepack_drv_forcing_bgc, only: get_forcing_bgc, get_atm_bgc, fzaero_data, & 
!          faero_default
      use icepack_drv_flux, only: init_flux_atm_ocn
      use icepack_intfc_tracers, only: tr_aero, tr_zaero
      use icepack_intfc_shared, only: skl_bgc, z_tracers

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      timeLoop: do

         call ice_step

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date

         call calendar(time)    ! at the end of the timestep

         if (stop_now >= 1) exit timeLoop

          call get_forcing(istep1)  ! get forcing from data arrays
!         call get_forcing_atmo     ! atmospheric forcing from data
!         call get_forcing_ocn(dt)  ! ocean forcing from data

         ! aerosols
!         if (tr_aero .or. tr_zaero)  call faero_default    ! default values

!         if (skl_bgc .or. z_tracers) call get_forcing_bgc  ! biogeochemistry
!         if (z_tracers) call get_atm_bgc                   ! biogeochemistry

         call init_flux_atm_ocn ! initialize atmosphere, ocean fluxes

      enddo timeLoop

      end subroutine icepack_run

!=======================================================================
!
!  Calls drivers for physics components, some initialization, and output
!
!  author Elizabeth C. Hunke, LANL

      subroutine ice_step

      use icepack_drv_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep
      use icepack_drv_constants, only: c0
      use icepack_drv_diagnostics, only: runtime_diags, init_mass_diags
      use icepack_drv_diagnostics_bgc, only: hbrine_diags, zsal_diags, bgc_diags
      use icepack_drv_domain_size, only: nslyr
      use icepack_drv_flux, only: scale_factor, init_history_therm, init_history_bgc, &
          daidtt, daidtd, dvidtt, dvidtd, dagedtt, dagedtd, init_history_dyn
      use icepack_drv_restart, only: dumpfile, final_restart
!      use icepack_drv_restart_column, only: write_restart_age, write_restart_FY, &
!          write_restart_lvl, write_restart_pond_cesm, write_restart_pond_lvl, &
!          write_restart_pond_topo, write_restart_aero, &
!          write_restart_bgc, write_restart_hbrine
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: tr_iage, tr_FY, tr_lvl, &
          tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_brine, tr_aero
      use icepack_drv_step_mod, only: prep_radiation, step_therm1, step_therm2, &
          update_state, step_dyn_ridge, step_radiation, &
          biogeochemistry
      use icepack_intfc_shared, only: calc_Tsfc, skl_bgc, solve_zsal, z_tracers

      integer (kind=int_kind) :: &
         k               ! dynamics supercycling index

      real (kind=dbl_kind) :: &
         offset          ! d(age)/dt time offset

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

      call init_mass_diags   ! diagnostics per timestep
      call init_history_therm
      call init_history_bgc

      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------
      
      if (calc_Tsfc) call prep_radiation (dt)
      
      !-----------------------------------------------------------------
      ! thermodynamics and biogeochemistry
      !-----------------------------------------------------------------

      call step_therm1     (dt) ! vertical thermodynamics
      call biogeochemistry (dt) ! biogeochemistry
      call step_therm2     (dt) ! ice thickness distribution thermo
      
      ! clean up, update tendency diagnostics
      offset = dt
      call update_state (dt, daidtt, dvidtt, dagedtt, offset)
      
      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------
      
      call init_history_dyn
      
      do k = 1, ndtd
        
        ! ridging
        call step_dyn_ridge (dt_dyn, ndtd)
        
        ! clean up, update tendency diagnostics
        offset = c0
        call update_state (dt_dyn, daidtd, dvidtd, dagedtd, offset)
        
      enddo
      
      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------
      
      call step_radiation (dt)
      
      !-----------------------------------------------------------------
      ! get ready for coupling and the next time step
      !-----------------------------------------------------------------
      
      call coupling_prep
      
      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------
      
      if (mod(istep,diagfreq) == 0) then
        call runtime_diags(dt)          ! log file
        if (solve_zsal) call zsal_diags(dt)  
        if (skl_bgc .or. z_tracers)  call bgc_diags (dt)
        if (tr_brine) call hbrine_diags(dt)
      endif
      
      if (write_restart == 1) then
        call dumpfile     ! core variables for restarting
        !            if (tr_iage)      call write_restart_age
        !            if (tr_FY)        call write_restart_FY
        !            if (tr_lvl)       call write_restart_lvl
        !            if (tr_pond_cesm) call write_restart_pond_cesm
        !            if (tr_pond_lvl)  call write_restart_pond_lvl
        !            if (tr_pond_topo) call write_restart_pond_topo
        !            if (tr_aero)      call write_restart_aero
        !            if (solve_zsal .or. skl_bgc .or. z_tracers) &
        !                              call write_restart_bgc 
        !            if (tr_brine)     call write_restart_hbrine
        !            if (kdyn == 2)    call write_restart_eap
        call final_restart
      endif
      
    end subroutine ice_step
    
!=======================================================================
!
! Prepare for coupling
!
! authors: Elizabeth C. Hunke, LANL

      subroutine coupling_prep

      use icepack_drv_arrays_column, only: alvdfn, alidfn, alvdrn, alidrn, &
          albicen, albsnon, albpndn, apeffn, fzsal_g, fzsal, snowfracn
      use icepack_drv_calendar, only: dt
      use icepack_intfc_shared, only: calc_Tsfc, oceanmixed_ice, max_aero
      use icepack_intfc_tracers, only: nbtrcr
      use icepack_drv_constants, only: c0, c1, puny, rhofresh
      use icepack_drv_domain_size, only: ncat, nx
      use icepack_drv_flux, only: alvdf, alidf, alvdr, alidr, albice, albsno, &
          albpnd, apeff_ai, coszen, fpond, fresh, l_mpond_fresh, &
          alvdf_ai, alidf_ai, alvdr_ai, alidr_ai, fhocn_ai, &
          fresh_ai, fsalt_ai, fsalt, &
          fswthru_ai, fhocn, fswthru, scale_factor, snowfrac, &
          swvdr, swidr, swvdf, swidf, Tf, Tair, Qa, strairxT, strairyt, &
          fsens, flat, fswabs, flwout, evap, Tref, Qref, &
          fsurfn_f, flatn_f, frzmlt_init, frzmlt, &
          faero_ocn, fzsal_ai, fzsal_g_ai, flux_bio, flux_bio_ai
      use icepack_drv_state, only: aicen, aice, aice_init
      use icepack_drv_step_mod, only: ocean_mixed_layer

      ! local variables

      integer (kind=int_kind) :: & 
         n           , & ! thickness category index
         i           , & ! horizontal index
         k               ! tracer index

      real (kind=dbl_kind) :: &
         netsw           ! flag for shortwave radiation presence

      !-----------------------------------------------------------------
      ! Save current value of frzmlt for diagnostics.
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

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
                  
            alvdf(i) = alvdf(i) &
               + alvdfn(i,n)*aicen(i,n)
            alidf(i) = alidf(i) &
               + alidfn(i,n)*aicen(i,n)
            alvdr(i) = alvdr(i) &
               + alvdrn(i,n)*aicen(i,n)
            alidr(i) = alidr(i) &
               + alidrn(i,n)*aicen(i,n)

            netsw = swvdr(i) + swidr(i) &
                  + swvdf(i) + swidf(i)
            if (netsw > puny) then ! sun above horizon
            albice(i) = albice(i) &
               + albicen(i,n)*aicen(i,n)
            albsno(i) = albsno(i) &
               + albsnon(i,n)*aicen(i,n)
            albpnd(i) = albpnd(i) &
               + albpndn(i,n)*aicen(i,n)
            endif

            apeff_ai(i) = apeff_ai(i) &       ! for history
               + apeffn(i,n)*aicen(i,n)
            snowfrac(i) = snowfrac(i) &       ! for history
               + snowfracn(i,n)*aicen(i,n)
               
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

      end module icepack_drv_RunMod

!=======================================================================
