!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL

      module icedrv_diagnostics

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, nu_diag_out
      use icedrv_domain_size, only: nx
      use icepack_intfc, only: c0
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private
      public :: runtime_diags, &
                init_mass_diags, &
                icedrv_diagnostics_debug, &
                print_state

      ! diagnostic output file
      character (len=char_len), public :: diag_file

      ! point print data

      logical (kind=log_kind), public :: &
         print_points         ! if true, print point data

      integer (kind=int_kind), parameter, public :: &
         npnt = 2             ! total number of points to be printed

      character (len=char_len), dimension(nx), public :: nx_names

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

      use icedrv_arrays_column, only: floe_rad_c
      use icedrv_domain_size, only: ncat, nfsd
      use icedrv_flux, only: evap, fsnow, frazil
      use icedrv_flux, only: fswabs, flw, flwout, fsens, fsurf, flat
      use icedrv_flux, only: frain
      use icedrv_flux, only: Tair, Qa, fsw, fcondtop
      use icedrv_flux, only: meltt, meltb, meltl, snoice
      use icedrv_flux, only: dsnow, congel, sst, sss, Tf, fhocn
      use icedrv_state, only: aice, vice, vsno, trcr, trcrn, aicen

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         n, nc, k

      logical (kind=log_kind) :: &
         calc_Tsfc, tr_fsd

      ! fields at diagnostic points
      real (kind=dbl_kind) :: & 
         pTair, pfsnow, pfrain, &
         paice, hiavg, hsavg, hbravg, psalt, pTsfc, &
         pevap, pfhocn, fsdavg

      real (kind=dbl_kind), dimension (nx) :: &
         work1, work2

      real (kind=dbl_kind) :: &
         Tffresh, rhos, rhow, rhoi

      logical (kind=log_kind) :: tr_brine
      integer (kind=int_kind) :: nt_fbri, nt_Tsfc, nt_fsd

      character(len=*), parameter :: subname='(runtime_diags)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine,tr_fsd_out=tr_fsd)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, nt_Tsfc_out=nt_Tsfc,&
                                        nt_fsd_out=nt_fsd)
      call icepack_query_parameters(Tffresh_out=Tffresh, rhos_out=rhos, &
           rhow_out=rhow, rhoi_out=rhoi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
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
        fsdavg = c0               ! FSD rep radius 
        hsavg = c0
        hbravg = c0               ! avg brine thickness
        psalt = c0 
        if (paice /= c0) then
          hiavg = vice(n)/paice
          hsavg = vsno(n)/paice
          if (tr_brine) hbravg = trcr(n,nt_fbri)* hiavg
          if (tr_fsd) then
              do nc = 1, ncat
              do k = 1, nfsd
                  fsdavg  = fsdavg &
                          + trcrn(n,nt_fsd+k-1,nc) * floe_rad_c(k) &
                          * aicen(n,nc) / paice
              end do
              end do
          end if

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
        write(nu_diag_out+n-1,900) 'air temperature (C)    = ',pTair
        write(nu_diag_out+n-1,900) 'specific humidity      = ',Qa(n)
        write(nu_diag_out+n-1,900) 'snowfall (m)           = ',pfsnow
        write(nu_diag_out+n-1,900) 'rainfall (m)           = ',pfrain
        if (.not.calc_Tsfc) then
          write(nu_diag_out+n-1,900) 'total surface heat flux= ', fsurf(n)
          write(nu_diag_out+n-1,900) 'top sfc conductive flux= ',fcondtop(n)
          write(nu_diag_out+n-1,900) 'latent heat flux       = ',flat(n)
        else
          write(nu_diag_out+n-1,900) 'shortwave radiation sum= ',fsw(n)
          write(nu_diag_out+n-1,900) 'longwave radiation     = ',flw(n)
        endif
        write(nu_diag_out+n-1,*) '----------ice----------'
        write(nu_diag_out+n-1,900) 'area fraction          = ',aice(n)! ice area
        write(nu_diag_out+n-1,900) 'avg ice thickness (m)  = ',hiavg
        write(nu_diag_out+n-1,900) 'avg snow depth (m)     = ',hsavg
        write(nu_diag_out+n-1,900) 'avg salinity (ppt)     = ',psalt
        write(nu_diag_out+n-1,900) 'avg brine thickness (m)= ',hbravg
        if (tr_fsd) &
        write(nu_diag_out+n-1,900) 'avg fsd rep radius (m) = ',fsdavg

       
        if (calc_Tsfc) then
          write(nu_diag_out+n-1,900) 'surface temperature(C) = ',pTsfc ! ice/snow
          write(nu_diag_out+n-1,900) 'absorbed shortwave flx = ',fswabs(n)
          write(nu_diag_out+n-1,900) 'outward longwave flx   = ',flwout(n)
          write(nu_diag_out+n-1,900) 'sensible heat flx      = ',fsens(n)
          write(nu_diag_out+n-1,900) 'latent heat flx        = ',flat(n)
        endif
        write(nu_diag_out+n-1,900) 'subl/cond (m ice)      = ',pevap   ! sublimation/condensation
        write(nu_diag_out+n-1,900) 'top melt (m)           = ',meltt(n)
        write(nu_diag_out+n-1,900) 'bottom melt (m)        = ',meltb(n)
        write(nu_diag_out+n-1,900) 'lateral melt (m)       = ',meltl(n)
        write(nu_diag_out+n-1,900) 'new ice (m)            = ',frazil(n) ! frazil
        write(nu_diag_out+n-1,900) 'congelation (m)        = ',congel(n)
        write(nu_diag_out+n-1,900) 'snow-ice (m)           = ',snoice(n)
        write(nu_diag_out+n-1,900) 'snow change (m)        = ',dsnow(n)
        write(nu_diag_out+n-1,900) 'effective dhi (m)      = ',pdhi(n)   ! ice thickness change
        write(nu_diag_out+n-1,900) 'effective dhs (m)      = ',pdhs(n)   ! snow thickness change
        write(nu_diag_out+n-1,900) 'intnl enrgy chng(W/m^2)= ',pde (n)   ! ice/snow energy change
        write(nu_diag_out+n-1,*) '----------ocn----------'
        write(nu_diag_out+n-1,900) 'sst (C)                = ',sst(n)  ! sea surface temperature
        write(nu_diag_out+n-1,900) 'sss (ppt)              = ',sss(n)  ! sea surface salinity
        write(nu_diag_out+n-1,900) 'freezing temp (C)      = ',Tf(n)   ! freezing temperature
        write(nu_diag_out+n-1,900) 'heat used (W/m^2)      = ',pfhocn  ! ocean heat used by ice
        
      end do
899   format (43x,a24)
900   format (a25,2x,f24.17)

      end subroutine runtime_diags

!=======================================================================

! Computes global combined ice and snow mass sum
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_mass_diags

      use icedrv_domain_size, only: nx
      use icedrv_state, only: vice, vsno

      integer (kind=int_kind) :: i

      real (kind=dbl_kind), dimension (nx) :: work1

      character(len=*), parameter :: subname='(init_mass_diags)'

      call total_energy (work1)
      do i = 1, nx
         pdhi(i) = vice (i)
         pdhs(i) = vsno (i)
         pde (i) = work1(i)
      enddo

      end subroutine init_mass_diags

!=======================================================================

! Computes total energy of ice and snow in a grid cell.
!
! authors: E. C. Hunke, LANL

      subroutine total_energy (work)

      use icedrv_domain_size, only: ncat, nilyr, nslyr, nx
      use icedrv_state, only: vicen, vsnon, trcrn

      real (kind=dbl_kind), dimension (nx), intent(out) :: &
         work      ! total energy

      ! local variables

      integer (kind=int_kind) :: &
        i, k, n

      integer (kind=int_kind) :: nt_qice, nt_qsno

      character(len=*), parameter :: subname='(total_energy)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_tracer_indices(nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

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

      use icedrv_domain_size, only: ncat, nilyr, nx
      use icedrv_state, only: vicen, trcrn

      real (kind=dbl_kind), dimension (nx),  &
         intent(out) :: &
         work      ! total salt

      ! local variables

      integer (kind=int_kind) :: &
        i, k, n

      integer (kind=int_kind) :: nt_sice

      character(len=*), parameter :: subname='(total_salt)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_tracer_indices(nt_sice_out=nt_sice)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

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
!
! Wrapper for the print_state debugging routine.
! Useful for debugging in the main driver (see ice.F_debug)
!
! author Elizabeth C. Hunke, LANL
!
      subroutine icedrv_diagnostics_debug (plabeld)

      use icedrv_calendar, only: istep1

      character (*), intent(in) :: plabeld

      character(len=*), parameter :: &
         subname='(icedrv_diagnostics_debug)'

      ! printing info for routine print_state

      integer (kind=int_kind), parameter :: &
         check_step = 1439, & ! begin printing at istep1=check_step
         ip = 3               ! i index

      if (istep1 >= check_step) then
         call print_state(plabeld,ip)
      endif

      end subroutine icedrv_diagnostics_debug

!=======================================================================

! This routine is useful for debugging.
! Calls to it should be inserted in the form (after thermo, for example)
!     plabel = 'post thermo'
!     if (istep1 >= check_step) call print_state(plabel,ip)
! 'use ice_diagnostics' may need to be inserted also
! author: Elizabeth C. Hunke, LANL

      subroutine print_state(plabel,i)

      use icedrv_calendar,  only: istep1, time
      use icedrv_domain_size, only: ncat, nilyr, nslyr, nfsd
      use icedrv_state, only: aice0, aicen, vicen, vsnon, uvel, vvel, trcrn
      use icedrv_flux, only: uatm, vatm, potT, Tair, Qa, flw, frain, fsnow
      use icedrv_flux, only: fsens, flat, evap, flwout
      use icedrv_flux, only: swvdr, swvdf, swidr, swidf, rhoa
      use icedrv_flux, only: frzmlt, sst, sss, Tf, Tref, Qref, Uref
      use icedrv_flux, only: uocn, vocn
      use icedrv_flux, only: fsw, fswabs, fswint_ai, fswthru, scale_factor
      use icedrv_flux, only: alvdr_ai, alvdf_ai, alidf_ai, alidr_ai

      character (*), intent(in) :: plabel

      integer (kind=int_kind), intent(in) :: & 
          i              ! horizontal index

      ! local variables

      real (kind=dbl_kind) :: &
          eidebug, esdebug, &
          qi, qs, Tsnow, &
          puny, Lfresh, cp_ice, &
          rhoi, rhos

      integer (kind=int_kind) :: n, k

      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_fsd

      logical (kind=log_kind) :: tr_fsd

      character(len=*), parameter :: subname='(print_state)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd) 
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno,nt_fsd_out=nt_fsd)
      call icepack_query_parameters(puny_out=puny, Lfresh_out=Lfresh, cp_ice_out=cp_ice, &
           rhoi_out=rhoi, rhos_out=rhos)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! write diagnostics
      !-----------------------------------------------------------------

      write(nu_diag,*) trim(plabel)
      write(nu_diag,*) 'istep1, i, time', &
                        istep1, i, time
      write(nu_diag,*) ' '
      write(nu_diag,*) 'aice0', aice0(i)
      do n = 1, ncat
         write(nu_diag,*) ' '
         write(nu_diag,*) 'n =',n
         write(nu_diag,*) 'aicen', aicen(i,n)
         write(nu_diag,*) 'vicen', vicen(i,n)
         write(nu_diag,*) 'vsnon', vsnon(i,n)
         if (aicen(i,n) > puny) then
            write(nu_diag,*) 'hin', vicen(i,n)/aicen(i,n)
            write(nu_diag,*) 'hsn', vsnon(i,n)/aicen(i,n)
         endif
         write(nu_diag,*) 'Tsfcn',trcrn(i,nt_Tsfc,n)
         if (tr_fsd) write(nu_diag,*) 'afsdn',trcrn(i,nt_fsd,n) ! fsd cat 1
         write(nu_diag,*) ' '
      enddo                     ! n

      eidebug = c0
      do n = 1,ncat
         do k = 1,nilyr
            qi = trcrn(i,nt_qice+k-1,n)
            write(nu_diag,*) 'qice, cat ',n,' layer ',k, qi
            eidebug = eidebug + qi
            if (aicen(i,n) > puny) then
               write(nu_diag,*)  'qi/rhoi', qi/rhoi
            endif
         enddo
         write(nu_diag,*) ' '
      enddo
      write(nu_diag,*) 'qice(i)',eidebug
      write(nu_diag,*) ' '

      esdebug = c0
      do n = 1,ncat
         if (vsnon(i,n) > puny) then
            do k = 1,nslyr
               qs = trcrn(i,nt_qsno+k-1,n)
               write(nu_diag,*) 'qsnow, cat ',n,' layer ',k, qs
               esdebug = esdebug + qs
               Tsnow = (Lfresh + qs/rhos) / cp_ice
               write(nu_diag,*) 'qs/rhos', qs/rhos
               write(nu_diag,*) 'Tsnow', Tsnow
            enddo
            write(nu_diag,*) ' '
         endif
      enddo
      write(nu_diag,*) 'qsnow(i)',esdebug
      write(nu_diag,*) ' '

      write(nu_diag,*) 'uvel(i)',uvel(i)
      write(nu_diag,*) 'vvel(i)',vvel(i)

      write(nu_diag,*) ' '
      write(nu_diag,*) 'atm states and fluxes'
      write(nu_diag,*) '            uatm    = ',uatm (i)
      write(nu_diag,*) '            vatm    = ',vatm (i)
      write(nu_diag,*) '            potT    = ',potT (i)
      write(nu_diag,*) '            Tair    = ',Tair (i)
      write(nu_diag,*) '            Qa      = ',Qa   (i)
      write(nu_diag,*) '            rhoa    = ',rhoa (i)
      write(nu_diag,*) '            swvdr   = ',swvdr(i)
      write(nu_diag,*) '            swvdf   = ',swvdf(i)
      write(nu_diag,*) '            swidr   = ',swidr(i)
      write(nu_diag,*) '            swidf   = ',swidf(i)
      write(nu_diag,*) '            flw     = ',flw  (i)
      write(nu_diag,*) '            frain   = ',frain(i)
      write(nu_diag,*) '            fsnow   = ',fsnow(i)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'ocn states and fluxes'
      write(nu_diag,*) '            frzmlt  = ',frzmlt (i)
      write(nu_diag,*) '            sst     = ',sst    (i)
      write(nu_diag,*) '            sss     = ',sss    (i)
      write(nu_diag,*) '            Tf      = ',Tf     (i)
      write(nu_diag,*) '            uocn    = ',uocn   (i)
      write(nu_diag,*) '            vocn    = ',vocn   (i)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'srf states and fluxes'
      write(nu_diag,*) '            Tref    = ',Tref  (i)
      write(nu_diag,*) '            Qref    = ',Qref  (i)
      write(nu_diag,*) '            Uref    = ',Uref  (i)
      write(nu_diag,*) '            fsens   = ',fsens (i)
      write(nu_diag,*) '            flat    = ',flat  (i)
      write(nu_diag,*) '            evap    = ',evap  (i)
      write(nu_diag,*) '            flwout  = ',flwout(i)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'shortwave'
      write(nu_diag,*) '            fsw          = ',fsw         (i)
      write(nu_diag,*) '            fswabs       = ',fswabs      (i)
      write(nu_diag,*) '            fswint_ai    = ',fswint_ai   (i)
      write(nu_diag,*) '            fswthru      = ',fswthru     (i)
      write(nu_diag,*) '            scale_factor = ',scale_factor(i)
      write(nu_diag,*) '            alvdr        = ',alvdr_ai    (i)
      write(nu_diag,*) '            alvdf        = ',alvdf_ai    (i)
      write(nu_diag,*) '            alidr        = ',alidr_ai    (i)
      write(nu_diag,*) '            alidf        = ',alidf_ai    (i)
      write(nu_diag,*) ' '

      call icepack_warnings_flush(nu_diag)

      end subroutine print_state

!=======================================================================

      end module icedrv_diagnostics

!=======================================================================
