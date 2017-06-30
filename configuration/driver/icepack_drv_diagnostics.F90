!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_drv_diagnostics

      use icepack_kinds_mod
      use icepack_drv_constants, only: c0, nu_diag
      use icepack_drv_calendar, only: diagfreq, istep1, istep
      use icepack_intfc_shared, only: max_aero

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
      real (kind=dbl_kind), dimension(npnt) :: &
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
      use icepack_drv_domain_size, only: ncat, n_aero, nx
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
         i, n

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         paice, pTair, pQa, pfsnow, pfrain, pfsw, pflw, & 
         pTsfc, pevap, pfswabs, pflwout, pflat, pfsens, &
         pfsurf, pfcondtop, psst, psss, pTf, hiavg, hsavg, hbravg, &
         pfhocn, psalt, &
         pmeltt, pmeltb, pmeltl, psnoice, pdsnow, pfrazil, pcongel

      real (kind=dbl_kind), dimension (nx) :: &
         work1, work2

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         call total_energy (work1)
         call total_salt   (work2)

         do n = 1, npnt
               i = n + 1  ! NOTE this prints the slab and full-ITD cells

               pTair(n) = Tair(i) - Tffresh ! air temperature
               pQa(n) = Qa(i)               ! specific humidity
               pfsnow(n) = fsnow(i)*dt/rhos ! snowfall
               pfrain(n) = frain(i)*dt/rhow ! rainfall
               pfsw(n) = fsw(i)             ! shortwave radiation
               pflw(n) = flw(i)             ! longwave radiation
               paice(n) = aice(i)           ! ice area
               
               hiavg(n) = c0                       ! avg snow/ice thickness
               hsavg(n) = c0
               hbravg(n) = c0                      ! avg brine thickness
               if (paice(n) /= c0) then
                  hiavg(n) = vice(i)/paice(n)
                  hsavg(n) = vsno(i)/paice(n)
                  if (tr_brine) hbravg(n) = trcr(i,nt_fbri)* hiavg(n)
               endif
               if (vice(i) /= c0) psalt(n) = work2(i)/vice(i)
               pTsfc(n) = trcr(i,nt_Tsfc)   ! ice/snow sfc temperature
               pevap(n) = evap(i)*dt/rhoi   ! sublimation/condensation
               pfswabs(n) = fswabs(i)       ! absorbed solar flux
               pflwout(n) = flwout(i)       ! outward longwave flux
               pflat(n) = flat(i)           ! latent heat flux
               pfsens(n) = fsens(i)         ! sensible heat flux
               pfsurf(n) = fsurf(i)         ! total sfc heat flux
               pfcondtop(n) = fcondtop(i)   ! top sfc cond flux
               pmeltt(n) = meltt(i)         ! top melt
               pmeltb(n) = meltb(i)         ! bottom melt
               pmeltl(n) = meltl(i)         ! lateral melt
               psnoice(n) = snoice(i)       ! snow ice
               pdsnow(n) = dsnow(i)         ! snow change
               pfrazil(n) = frazil(i)       ! frazil ice
               pcongel(n) = congel(i)       ! congelation ice
               pdhi(n) = vice(i) - pdhi(n)  ! ice thickness change
               pdhs(n) = vsno(i) - pdhs(n)  ! snow thickness change
               pde(n) =-(work1(i)- pde(n))/dt ! ice/snow energy change 
               psst(n) = sst(i)             ! sea surface temperature
               psss(n) = sss(i)             ! sea surface salinity
               pTf(n) = Tf(i)               ! freezing temperature
               pfhocn(n) = -fhocn(i)        ! ocean heat used by ice

         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

        write(nu_diag,899) '  slab',' full ITD'

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

       if (print_points) then

        write(nu_diag,*) '                         '
        write(nu_diag,*) '----------atm----------'
        write(nu_diag,900) 'air temperature (C)    = ',pTair(1),pTair(2)
        write(nu_diag,900) 'specific humidity      = ',pQa(1),pQa(2)
        write(nu_diag,900) 'snowfall (m)           = ',pfsnow(1), &
                                                       pfsnow(2)
        write(nu_diag,900) 'rainfall (m)           = ',pfrain(1), &
                                                       pfrain(2)
        if (.not.calc_Tsfc) then
           write(nu_diag,900) 'total surface heat flux= ',pfsurf(1),pfsurf(2)
           write(nu_diag,900) 'top sfc conductive flux= ',pfcondtop(1), &
                                                          pfcondtop(2)
           write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        else
           write(nu_diag,900) 'shortwave radiation sum= ',pfsw(1),pfsw(2)
           write(nu_diag,900) 'longwave radiation     = ',pflw(1),pflw(2)
        endif
        write(nu_diag,*) '----------ice----------'
        write(nu_diag,900) 'area fraction          = ',paice(1),paice(2)
        write(nu_diag,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
        write(nu_diag,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
        write(nu_diag,900) 'avg salinity (ppt)     = ',psalt(1),psalt(2)
        write(nu_diag,900) 'avg brine thickness (m)= ',hbravg(1),hbravg(2)

        if (calc_Tsfc) then
           write(nu_diag,900) 'surface temperature(C) = ',pTsfc(1),pTsfc(2)
           write(nu_diag,900) 'absorbed shortwave flx = ',pfswabs(1), &
                                                          pfswabs(2)
           write(nu_diag,900) 'outward longwave flx   = ',pflwout(1), &
                                                          pflwout(2)
           write(nu_diag,900) 'sensible heat flx      = ',pfsens(1), &
                                                          pfsens(2)
           write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        endif
        write(nu_diag,900) 'subl/cond (m ice)      = ',pevap(1),pevap(2)
        write(nu_diag,900) 'top melt (m)           = ',pmeltt(1) &
                                                      ,pmeltt(2)
        write(nu_diag,900) 'bottom melt (m)        = ',pmeltb(1) &
                                                      ,pmeltb(2)
        write(nu_diag,900) 'lateral melt (m)       = ',pmeltl(1) &
                                                      ,pmeltl(2)
        write(nu_diag,900) 'new ice (m)            = ',pfrazil(1), &
                                                       pfrazil(2)
        write(nu_diag,900) 'congelation (m)        = ',pcongel(1), &
                                                       pcongel(2)
        write(nu_diag,900) 'snow-ice (m)           = ',psnoice(1), &
                                                       psnoice(2)
        write(nu_diag,900) 'snow change (m)        = ',pdsnow(1), &
                                                       pdsnow(2)
        write(nu_diag,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
        write(nu_diag,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
        write(nu_diag,900) 'intnl enrgy chng(W/m^2)= ',pde (1),pde (2)
        write(nu_diag,*) '----------ocn----------'
        write(nu_diag,900) 'sst (C)                = ',psst(1),psst(2)
        write(nu_diag,900) 'sss (ppt)              = ',psss(1),psss(2)
        write(nu_diag,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
        write(nu_diag,900) 'heat used (W/m^2)      = ',pfhocn(1), &
                                                       pfhocn(2)

       endif                    ! print_points

  799 format (27x,a24)
  800 format (a25,2x,f24.17)
  801 format (a25,2x,1pe24.17)
  899 format (27x,a24,2x,a24)
  900 format (a25,2x,f24.17,2x,f24.17)
  901 format (a25,2x,1pe24.17,2x,1pe24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

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
         do n = 1, npnt
            i = n+1 ! NOTE this must be the same as above (FIX THIS)
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
