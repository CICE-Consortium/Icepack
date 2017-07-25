!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_drv_diagnostics

      use icepack_kinds_mod
      use icepack_drv_constants, only: c0, nu_diag, nu_diag_out
      use icepack_drv_calendar, only: diagfreq, istep1, istep
      use icepack_intfc_shared, only: max_aero
      use icepack_drv_domain_size, only: nx

      implicit none
      private
      public :: runtime_diags, diagnostic_abort, init_mass_diags

      save

      ! diagnostic output file
      character (len=char_len), public :: diag_file

      ! point print data

      logical (kind=log_kind), public :: &
         print_points         ! if true, print point data

      integer (kind=int_kind), parameter, public :: &
         npnt = 2             ! total number of points to be printed

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(nx) :: &
         pdhi             , & ! change in mean ice thickness (m)
         pdhs             , & ! change in mean snow thickness (m)
         pde                  ! change in ice and snow energy (W m-2)

!=======================================================================

      contains

!=======================================================================

! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW

      subroutine runtime_diags (dt)

      use icepack_intfc_shared, only: calc_Tsfc, ktherm
      use icepack_drv_constants, only: c1, c1000, c2, p001, p5, puny, rhoi, rhos, rhow, &
          rhofresh, Tffresh, Lfresh, Lvap, ice_ref_salinity, &
          m2_to_km2, awtvdr, awtidr, awtvdf, awtidf
      use icepack_drv_domain_size, only: ncat, n_aero
      use icepack_drv_flux, only: alvdr, alidr, alvdf, alidf, evap, fsnow, frazil, &
          fswabs, fswthru, flw, flwout, fsens, fsurf, flat, frzmlt_init, frain, fpond, &
          coszen, fhocn_ai, fsalt_ai, fresh_ai, frazil_diag, &
          update_ocn_f, Tair, Qa, fsw, fcondtop, meltt, meltb, meltl, snoice, &
          dsnow, congel, sst, sss, Tf, fhocn, &
          swvdr, swvdf, swidr, swidf, &
          alvdr_init, alvdf_init, alidr_init, alidf_init, faero_atm, faero_ocn
      use icepack_drv_state   ! everything
      use icepack_intfc_tracers ! everything
#ifdef CCSMCOUPLED
      use icepack_drv_prescribed_mod, only: prescribed_ice
#endif

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         n

      ! fields at diagnostic points
      real (kind=dbl_kind) :: & 
          pTair, &
          pfsnow, pfrain, &
          paice, hiavg, hsavg, hbravg, psalt, pTsfc, &
          pevap, &
          pfhocn

      real (kind=dbl_kind), dimension (nx) :: &
         work1, work2

      character (len=char_len), dimension(nx) :: nx_names

      nx_names(1) = ' icefree'
      nx_names(2) = '    slab'
      nx_names(3) = 'full ITD'
      nx_names(4) = '    land'

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         call total_energy (work1)
         call total_salt   (work2)

         do n = 1, nx
           pTair = Tair(n) - Tffresh ! air temperature
           pfsnow = fsnow(n)*dt/rhos ! snowfall
           pfrain = frain(n)*dt/rhow ! rainfall

           paice = aice(n)           ! ice area           
           hiavg = c0                ! avg snow/ice thickness
           hsavg = c0
           hbravg = c0               ! avg brine thickness
           psalt = c0 
           if (paice /= c0) then
             hiavg = vice(n)/paice
             hsavg = vsno(n)/paice
             if (tr_brine) hbravg = trcr(n,nt_fbri)* hiavg
           endif
           if (vice(n) /= c0) psalt = work2(n)/vice(n)
           pTsfc = trcr(n,nt_Tsfc)   ! ice/snow sfc temperature
           pevap = evap(n)*dt/rhoi   ! sublimation/condensation
           pdhi(n) = vice(n) - pdhi(n)  ! ice thickness change
           pdhs(n) = vsno(n) - pdhs(n)  ! snow thickness change
           pde(n) =-(work1(n)- pde(n))/dt ! ice/snow energy change
           pfhocn = -fhocn(n)        ! ocean heat used by ice
           
           !-----------------------------------------------------------------
           ! start spewing
           !-----------------------------------------------------------------
          
           write(nu_diag_out+n-1,899) nx_names(n)
           
           write(nu_diag_out+n-1,*) '                         '
           write(nu_diag_out+n-1,*) '----------atm----------'
           write(nu_diag_out+n-1,900) 'air temperature (C)    = ',pTair  ! air temperature
           write(nu_diag_out+n-1,900) 'specific humidity      = ',Qa(n)  ! specific humidity
           write(nu_diag_out+n-1,900) 'snowfall (m)           = ',pfsnow ! snowfall
           write(nu_diag_out+n-1,900) 'rainfall (m)           = ',pfrain ! rainfall
           if (.not.calc_Tsfc) then
             write(nu_diag_out+n-1,900) 'total surface heat flux= ', fsurf(n)  ! total sfc heat flux
             write(nu_diag_out+n-1,900) 'top sfc conductive flux= ',fcondtop(n)! top sfc cond flux
             write(nu_diag_out+n-1,900) 'latent heat flx        = ',flat(n)    ! latent heat flux
           else
             write(nu_diag_out+n-1,900) 'shortwave radiation sum= ',fsw(n)  ! shortwave radiation
             write(nu_diag_out+n-1,900) 'longwave radiation     = ',flw(n)  ! longwave radiation
           endif
           write(nu_diag_out+n-1,*) '----------ice----------'
           write(nu_diag_out+n-1,900) 'area fraction          = ',aice(n)  ! ice area    
           write(nu_diag_out+n-1,900) 'avg ice thickness (m)  = ',hiavg  ! avg snow/ice thickness
           write(nu_diag_out+n-1,900) 'avg snow depth (m)     = ',hsavg  ! avg snow/ice thickness
           write(nu_diag_out+n-1,900) 'avg salinity (ppt)     = ',psalt  ! avg salinity 
           write(nu_diag_out+n-1,900) 'avg brine thickness (m)= ',hbravg   ! avg brine thickness
           
           if (calc_Tsfc) then
             write(nu_diag_out+n-1,900) 'surface temperature(C) = ',pTsfc ! ice/snow sfc temperature
             write(nu_diag_out+n-1,900) 'absorbed shortwave flx = ',fswabs(n)! absorbed solar flux
             write(nu_diag_out+n-1,900) 'outward longwave flx   = ',flwout(n)! outward longwave flux
             write(nu_diag_out+n-1,900) 'sensible heat flx      = ',fsens(n) ! sensible heat flux
             write(nu_diag_out+n-1,900) 'latent heat flx        = ',flat(n)  ! latent heat flux
           endif
           write(nu_diag_out+n-1,900) 'subl/cond (m ice)      = ',pevap   ! sublimation/condensation
           write(nu_diag_out+n-1,900) 'top melt (m)           = ',meltt(n)! top melt
           write(nu_diag_out+n-1,900) 'bottom melt (m)        = ',meltb(n)! bottom melt
           write(nu_diag_out+n-1,900) 'lateral melt (m)       = ',meltl(n)! lateral melt
           write(nu_diag_out+n-1,900) 'new ice (m)            = ',frazil(n) ! frazil ice
           write(nu_diag_out+n-1,900) 'congelation (m)        = ',congel(n) ! congelation ice
           write(nu_diag_out+n-1,900) 'snow-ice (m)           = ',snoice(n) ! snow ice
           write(nu_diag_out+n-1,900) 'snow change (m)        = ',dsnow(n)  ! snow change
           write(nu_diag_out+n-1,900) 'effective dhi (m)      = ',pdhi(n)   ! ice thickness change
           write(nu_diag_out+n-1,900) 'effective dhs (m)      = ',pdhs(n)   ! snow thickness change
           write(nu_diag_out+n-1,900) 'intnl enrgy chng(W/m^2)= ',pde (n)   ! ice/snow energy change
           write(nu_diag_out+n-1,*) '----------ocn----------'
           write(nu_diag_out+n-1,900) 'sst (C)                = ',sst(n)  ! sea surface temperature
           write(nu_diag_out+n-1,900) 'sss (ppt)              = ',sss(n)  ! sea surface salinity
           write(nu_diag_out+n-1,900) 'freezing temp (C)      = ',Tf(n)   ! freezing temperature
           write(nu_diag_out+n-1,900) 'heat used (W/m^2)      = ',pfhocn  ! ocean heat used by ice
           
         end do
       endif                    ! print_points
       !  799 format (27x,a24)
       !  800 format (a25,2x,f24.17)
       !  801 format (a25,2x,1pe24.17)
       !  899 format (27x,a24,2x,a24)
899    format (43x,a24)
900    format (a25,2x,f24.17)
       !  901 format (a25,2x,1pe24.17,2x,1pe24.17)
       !  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
       !  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine runtime_diags

!=======================================================================

! Computes global combined ice and snow mass sum
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_mass_diags

      use icepack_drv_domain_size, only: nx
      use icepack_drv_state, only: vice, vsno

      integer (kind=int_kind) :: n, i

      real (kind=dbl_kind), dimension (nx) :: &
         work1

      call total_energy (work1)
      if (print_points) then
         do n = 1, nx
            i = n ! NOTE this must be the same as above (FIX THIS)
               pdhi(n) = vice(i)
               pdhs(n) = vsno(i)
               pde (n) = work1(i)
         enddo  ! npnt
      endif                     ! print_points

      end subroutine init_mass_diags

!=======================================================================

! Computes total energy of ice and snow in a grid cell.
!
! authors: E. C. Hunke, LANL

      subroutine total_energy (work)

      use icepack_drv_domain_size, only: ncat, nilyr, nslyr, nx
      use icepack_drv_state, only: vicen, vsnon, trcrn
      use icepack_intfc_tracers, only: nt_qice, nt_qsno

      real (kind=dbl_kind), dimension (nx),  &
         intent(out) :: &
         work      ! total energy

      ! local variables

      integer (kind=int_kind) :: &
        i, k, n

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         work(:) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do i = 1, nx
                  work(i) = work(i) &
                                 + trcrn(i,nt_qice+k-1,n) &
                                 * vicen(i,n) / real(nilyr,kind=dbl_kind)
               enddo            ! i
            enddo               ! k

            do k = 1, nslyr
               do i = 1, nx
                  work(i) = work(i) &
                                 + trcrn(i,nt_qsno+k-1,n) &
                                 * vsnon(i,n) / real(nslyr,kind=dbl_kind)
               enddo            ! i
            enddo               ! k
         enddo                  ! n

      end subroutine total_energy

!=======================================================================

! Computes bulk salinity of ice and snow in a grid cell.
! author: E. C. Hunke, LANL

      subroutine total_salt (work)

      use icepack_drv_domain_size, only: ncat, nilyr, nslyr, nx
      use icepack_drv_state, only: vicen, trcrn
      use icepack_intfc_tracers, only: nt_sice

      real (kind=dbl_kind), dimension (nx),  &
         intent(out) :: &
         work      ! total salt

      ! local variables

      integer (kind=int_kind) :: &
        i, k, n

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         work(:) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do i = 1, nx
                  work(i) = work(i) &
                                 + trcrn(i,nt_sice+k-1,n) &
                                 * vicen(i,n) / real(nilyr,kind=dbl_kind)
               enddo            ! i
            enddo               ! k
         enddo                  ! n

      end subroutine total_salt

!=======================================================================

! prints error information prior to aborting

      subroutine diagnostic_abort(istop, istep1, stop_label)

      use icepack_drv_constants, only: nu_diag
      use icepack_drv_state, only: aice

      integer (kind=int_kind), intent(in) :: &
         istop       , & ! indices of grid cell where model aborts
         istep1          ! time step number

      character (char_len), intent(in) :: stop_label

      ! local variables

      write (nu_diag,*) 'istep1, i, aice =', &
                         istep1, istop, aice(istop)
      write (nu_diag,*) stop_label
      stop

      end subroutine diagnostic_abort

!=======================================================================

      end module icepack_drv_diagnostics

!=======================================================================
