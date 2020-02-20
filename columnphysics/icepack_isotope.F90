!=======================================================================
!
!BOP
!
! !MODULE: icepack_isotope - Isotope tracer within sea ice
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors: David Bailey, NCAR
!          Jiang Zhu, UW-Madison
!
! 2014: Added i2x evaporation flux
!       Added fractionation options
!       Fixed bugs
!
! !INTERFACE:
!
      module icepack_isotope
!
! !USES:
!
      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, p001, p1, p5, puny
!
!EOP
!
      implicit none

      character(len=5), parameter ::    &
         frac = 'cfrac'   ! fractionation coefficient calculation method
                          !  cfrac, constant fractionation
                          !  nfrac, nonfractionation
                          !  gfrac, growth-rate dependent for H2_18O

      ! Species indicies - public so thay can be seen by water_tracers
      integer, parameter, public  :: ispundef = 0    ! Undefined
      integer, parameter, public  :: isph2o   = 1    ! H2O    ! "regular" water
      integer, parameter, public  :: isph216o = 2    ! H216O  ! H216O, nearly the same as "regular" water
      integer, parameter, public  :: isphdo   = 3    ! HDO
      integer, parameter, public  :: isph218o = 4    ! H218O

      integer, parameter, public  :: pwtspec = 4    ! number of water species (h2o,hdo,h218o,h216o)

!=======================================================================

      contains

!=======================================================================

!BOP
!
! !ROUTINE: update_isotope
!
! !DESCRIPTION:
!
!  Increase isotope in ice or snow surface due to deposition
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine update_isotope (dt,                  &
                                nilyr,    nslyr,      &
                                n_iso,                &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,   evap,       &
                                fsnow,    Tsfc,      &
                                Qref_iso,             &
                                isosno, isoice,       &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                fiso_atm, fiso_evapn, &
                                fiso_ocnn, HDO_ocn,   &
                                H2_16O_ocn, H2_18O_ocn)
!
! !USES:
!
!     use water_isotopes, only: wiso_alpi
      use icepack_parameters, only: ktherm, rhoi, rhos, Tffresh
      use icepack_tracers, only: nt_iso
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nilyr, nslyr, n_iso

      real (kind=dbl_kind), intent(in) :: &
         dt                     ! time step

      real (kind=dbl_kind), intent(in) :: &
         Tsfc,     &
         meltt,    &
         melts,    &
         meltb,    &
         congel,   &            ! congelation ice growth    (m/step)
         snoice,   &            ! ice thickness increase    (m/step)
         evap,     &            ! surface evaporation
         fsnow,    &            ! snowfall       (kg/m^2/s of water)
         vicen,    &            ! volume per unit area of ice    (m)
         vsnon,    &
         aicen,    &
         aice_old, &
         vice_old, &
         vsno_old, &
         HDO_ocn,    &
         H2_16O_ocn, &
         H2_18O_ocn 

      real (kind=dbl_kind), dimension(:), intent(in) ::  &
         fiso_atm,      &       ! isotopic snowfall (kg/m^2/s of water)
         Qref_iso              ! isotope reference humidity

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         fiso_ocnn,     &       ! isotopic freshwater (kg/m^2/s)
         fiso_evapn             ! evaporative water flux (kg/m^2/s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         isosno, isoice         ! mass of isotopes  (kg)

!
!  local variables
!
      real (kind=dbl_kind) :: &
         dzssl,     &
         dzint,     &
         dzssli,    &
         dzinti,    &
         evaps,     &           ! evaporation over snow     (m/step)
         evapi,     &           ! evaporation over ice      (m/step)
         dhs_snoice,&           ! snow thickness reduction  (m/step)
         hi,        &
         hs

      real (kind=dbl_kind), dimension(n_iso) :: &
        isotot, isotot0         ! for diagnostics 

      real (kind=dbl_kind) :: &
         dzssl_new, &
         dzint_new, &
         dzssli_new, &
         dzinti_new, &
         dznew

      real (kind=dbl_kind) :: &
         hs_old, hi_old, hslyr_old, hilyr_old, dhs, dhi,        &
         hslyr, hilyr, sloss1, sloss2,                          &
         TsfK,      &           ! snow/ice temperature (K)
         alphai,    &           ! ice/vapour fractionation coefficient
         ratio,     &           ! isotopic ratio
         work,      &           ! temporary variable
         alpha

      integer (kind=int_kind) :: k

! These need to be the same as in the DE code. Put in a common place?
      real (kind=dbl_kind), parameter :: &
        hi_ssl = .050_dbl_kind, &
        hs_ssl = .040_dbl_kind

! initialize

      hs_old=vsno_old/aice_old
      hi_old=vice_old/aice_old
      hslyr_old=hs_old/real(nslyr,kind=dbl_kind)
      hilyr_old=hi_old/real(nilyr,kind=dbl_kind)

      dzssl=min(hslyr_old/c2,hs_ssl)
      dzint=hs_old-dzssl
      dzssli=min(hilyr_old/c2,hi_ssl)
      dzinti=hi_old-dzssli

      if (aicen > puny) then
         hs = vsnon/aicen
         hi = vicen/aicen
      elseif (aice_old > puny) then
         hs = vsnon/aice_old
         hi = vicen/aice_old
      endif

      if (ktherm == 2) then
         dhs_snoice = snoice
      else
         dhs_snoice = snoice*rhoi/rhos
      endif

!     if (hs > puny) then
!        evaps = evap*dt/rhos
!     else
!        evapi = evap*dt/rhoi
!     endif
      evaps = hs-(hs_old-melts-dhs_snoice+&
          fsnow/rhos*dt)
      evapi = hi-(hi_old-meltt-meltb+congel+snoice)

! condensation of vapor onto snow and ice

      TsfK = Tsfc + Tffresh

      if (evaps > c0) then   ! condensation to snow
         do k = 1, n_iso         
            ratio = c1   ! ratio between 18O(HDO) and 16O in humidity
            alphai = c1  ! fractionation coefficient
            if (frac.ne.'nfrac' .and. Qref_iso(2)>puny) &
               ratio = Qref_iso(k)/Qref_iso(2)
            if (Qref_iso(2) <= puny) ratio = c0
            if (frac.ne.'nfrac' .and. k==1) alphai = wiso_alpi(3,TsfK)
            if (frac.ne.'nfrac' .and. k==2) alphai = wiso_alpi(2,TsfK)
            if (frac.ne.'nfrac' .and. k==3) alphai = wiso_alpi(4,TsfK)
            work = alphai*ratio*rhos*evaps*aicen
            fiso_evapn(k) = fiso_evapn(k)+work/dt
            isosno(k,1) = isosno(k,1)+work
         enddo
         dzssl = dzssl+evaps
      endif

      if (evapi > c0) then   ! condensation to ice
         do k = 1, n_iso         
            ratio = c1 ! ratio between 18O(HDO) and 16O in ref humidity
            alphai = c1  ! fractionation coefficient
            if (frac.ne.'nfrac' .and. Qref_iso(2)>puny) &
               ratio = Qref_iso(k)/Qref_iso(2)
            if (Qref_iso(2) <= puny) ratio = c0
            if (frac.ne.'nfrac' .and. k==1) alphai = wiso_alpi(3,TsfK)
            if (frac.ne.'nfrac' .and. k==2) alphai = wiso_alpi(2,TsfK)
            if (frac.ne.'nfrac' .and. k==3) alphai = wiso_alpi(4,TsfK)
            work = alphai*ratio*rhoi*evapi*aicen
            fiso_evapn(k) = fiso_evapn(k)+work/dt
            isoice(k,1) = isoice(k,1)+work
         enddo
         dzssli = dzssli+evapi
      endif

!     basal ice growth and isotope uptake

      if (congel > c0) then
         alpha = isoice_alpha(congel/dt,'HDO',frac)
         work = alpha*HDO_ocn*rhoi*congel*aicen
         isoice(1,2) = isoice(1,2)+work
         fiso_ocnn(1) = fiso_ocnn(1)-work/dt

         alpha = isoice_alpha(congel/dt,'H2_16O',frac)
         work = alpha*H2_16O_ocn*rhoi*congel*aicen
         isoice(2,2) = isoice(2,2)+work
         fiso_ocnn(2) = fiso_ocnn(2)-work/dt

         alpha = isoice_alpha(congel/dt,'H2_18O',frac)
         work = alpha*H2_18O_ocn*rhoi*congel*aicen
         isoice(3,2) = isoice(3,2)+work
         fiso_ocnn(3) = fiso_ocnn(3)-work/dt

         dzinti = dzinti+congel
      endif

! sublimation of snow and ice

      if (evaps < c0) then   ! snow sublimation (no fractionation)
         do k = 1, n_iso         
            !ratio = c1 ! ratio between 18O(HDO) and 16O in snow ssl
            !if (isosno(2,1) > puny)      &
            !   ratio = isosno(k,1)/isosno(2,1)
            !if (ratio > c5) ratio = c1   !! remove latter?
            !work = ratio*rhos*evaps*aicen
            !fiso_evapn(k) = fiso_evapn(k)+work/dt
               
            sloss1 = c0
            sloss2 = c0
            if (dzssl > puny)                &
               sloss1 = isosno(k,1)*         &
                  min(-evaps,dzssl)/dzssl
            if (isosno(k,1) >= sloss1) then
               isosno(k,1) = isosno(k,1)-sloss1
            else
               sloss1 = isosno(k,1)
               isosno(k,1) = c0
            endif
            if (dzint > puny)                &
               sloss2 = isosno(k,2)*         &
                  max(-evaps-dzssl,c0)/dzint
            if (isosno(k,2) >= sloss2) then
               isosno(k,2) = isosno(k,2)-sloss2
            else
               sloss2 = isosno(k,2)
               isosno(k,2) = c0
            endif
!           if (isosno(k,2) < c0) then
!              write(nu_diag,*) 'Neg isosno(k,2)',isosno(k,2),sloss2
!           endif
            fiso_evapn(k) = fiso_evapn(k)-(sloss1+sloss2)/dt
         enddo

         dzint = dzint+min(dzssl+evaps,c0)
         dzssl = max(dzssl+evaps,c0)
         if ( dzssl <= c0) then  ! ssl goes away
            fiso_evapn(:) = fiso_evapn(:)-isosno(:,1)/dt
            isosno(:,1) = c0
            dzssl = c0
         endif
         if (dzint <= c0) then   ! int goes away
            fiso_evapn(:) = fiso_evapn(:)-isosno(:,2)/dt
            isosno(:,2) = c0
            dzint = c0
         endif
      endif

      if (evapi < c0) then   ! ice sublimation (no fractionation)
         do k = 1, n_iso         
            !!ratio = c1 ! ratio between 18O(HDO) and 16O in ice ssl
            !!if (isoice(2,1) > puny)      &
            !!   ratio = isoice(k,1)/isoice(2,1)
            !!if (ratio > c5) ratio = c1   ! remove latter?
            !!work = ratio*rhoi*evapi*aicen
            !!fiso_evapn(k) = fiso_evapn(k)+work/dt

            sloss1 = c0
            sloss2 = c0
            if (dzssli > puny)               &
               sloss1 = isoice(k,1)*         &
                  min(-evapi,dzssli)/dzssli
            if (isoice(k,1) >= sloss1) then
               isoice(k,1) = isoice(k,1)-sloss1
            else
               sloss1 = isoice(k,1)
               isoice(k,1) = c0
            endif
            if (dzinti > puny)               &
               sloss2 = isoice(k,2)*         &
                  max(-evapi-dzssli,c0)/dzinti
            if (isoice(k,2) >= sloss2) then
               isoice(k,2) = isoice(k,2)-sloss2
            else
               sloss2 = isoice(k,2)
               isoice(k,2) = c0
            endif
            fiso_evapn(k) = fiso_evapn(k)-(sloss1+sloss2)/dt
         enddo

         dzinti = dzinti+min(dzssli+evapi,c0)
         dzssli = max(dzssli+evapi,c0)
         if ( dzssli <= c0) then ! ssl goes away
            fiso_evapn(:) = fiso_evapn(:)-isoice(:,1)/dt
            isoice(:,1) = c0
            dzssli = c0
         endif
         if (dzinti <= c0) then  ! int goes away
            fiso_evapn(:) = fiso_evapn(:)-isoice(:,2)/dt
            isoice(:,2) = c0
            dzinti = c0
         endif
      endif

!     surface snow melt

      if (melts > c0) then
         do k=1,n_iso
            sloss1=c0
            sloss2=c0
            if (dzssl > puny)         &
             sloss1 = isosno(k,1)     &
                   *min(melts,dzssl)/dzssl
            if (isosno(k,1) >= sloss1) then
               isosno(k,1) = isosno(k,1)-sloss1
            else
               sloss1 = isosno(k,1)
               isosno(k,1) = c0
            endif
            if (dzint > puny)             &
               sloss2=isosno(k,2)         &
                     *max(melts-dzssl,c0)/dzint
            if (isosno(k,2) >= sloss2) then
               isosno(k,2) = isosno(k,2)-sloss2
            else
               sloss2 = isosno(k,2)
               isosno(k,2) = c0
            endif
!           if (isosno(k,2) < c0) then
!               write(nu_diag,*) 'Neg isosno(k,2)',isosno(k,2),sloss2
!           endif
            fiso_ocnn(k)=fiso_ocnn(k)+(sloss1+sloss2)/dt
         enddo  ! n_iso

         dzint=dzint+min(dzssl-melts,c0)
         dzssl=max(dzssl-melts,c0)
         if ( dzssl <= c0) then ! ssl melts away
          fiso_ocnn(:)= fiso_ocnn(:)+isosno(:,1)/dt
          isosno(:,1) = c0
          dzssl = c0
         endif
         if (dzint <= c0) then  ! int melts away
          fiso_ocnn(:) = fiso_ocnn(:)+isosno(:,2)/dt
          isosno(:,2) = c0
          dzint = c0
         endif
      endif

!     surface ice melt
        if (meltt > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzssli > puny)    &
           sloss1=isoice(k,1)   &
                 *min(meltt,dzssli)/dzssli
          if (isoice(k,1) >= sloss1) then
             isoice(k,1) = isoice(k,1)-sloss1
          else
             sloss1 = isoice(k,1)
             isoice(k,1) = c0
          endif
          if (dzinti > puny)    &
           sloss2=isoice(k,2)   &
                 *max(meltt-dzssli,c0)/dzinti
          if (isoice(k,2) >= sloss2) then
             isoice(k,2) = isoice(k,2)-sloss2
          else
             sloss2 = isoice(k,2)
             isoice(k,2) = c0
          endif
          fiso_ocnn(k)=fiso_ocnn(k)+(sloss1+sloss2)/dt
         enddo

         dzinti=dzinti+min(dzssli-meltt,c0)
         dzssli=max(dzssli-meltt,c0)
         if (dzssli <= c0) then   ! ssl ice melts away
          fiso_ocnn(:) = fiso_ocnn(:)+isoice(:,1)
          isoice(:,1) = c0
          dzssli = c0
         endif
         if (dzinti <= c0) then   ! int ice melts away
          fiso_ocnn(:) = fiso_ocnn(:)+isoice(:,2)
          isoice(:,2) = c0
          dzinti = c0
         endif
        endif

!      basal ice melt.  Assume all isotopes lost in basal melt

        if (meltb > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzssli > puny)  &
           sloss1=max(meltb-dzinti,c0)  &
                 *isoice(k,1)/dzssli
          if (isoice(k,1) >= sloss1) then
             isoice(k,1) = isoice(k,1)-sloss1
          else
             sloss1 = isoice(k,1)
             isoice(k,1) = c0
          endif
          if (dzinti > puny)  &
           sloss2=min(meltb,dzinti)  &
                 *isoice(k,2)/dzinti
          if (isoice(k,2) >= sloss2) then
             isoice(k,2) = isoice(k,2)-sloss2
          else
             sloss2 = isoice(k,2)
             isoice(k,2) = c0
          endif
          fiso_ocnn(k)=fiso_ocnn(k)+(sloss1+sloss2)/dt
         enddo
 
         dzssli = dzssli+min(dzinti-meltb, c0)
         dzinti = max(dzinti-meltb, c0)           
         if (dzssli <= c0) then   ! ssl ice melts away
          fiso_ocnn(:) = fiso_ocnn(:)+isoice(:,1)
          isoice(:,1) = c0
          dzssli = c0
         endif
         if (dzinti <= c0) then   ! int ice melts away
          fiso_ocnn(:) = fiso_ocnn(:)+isoice(:,2)
          isoice(:,2) = c0
          dzinti = c0
         endif
        endif

!     snowfall and isotope deposition

         if (fsnow > c0) then
           isosno(:,1) = isosno(:,1)+               &
                             fiso_atm(:)*aicen*dt
           dzssl = dzssl+fsnow/rhos*dt
         endif

!     snoice formation

        if (dhs_snoice > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzint > puny)                         &
           sloss2 = min(dhs_snoice,dzint)      &
                  *isosno(k,2)/dzint
          if (isosno(k,2) >= sloss2) then
             isosno(k,2) = isosno(k,2)-sloss2
          else
             sloss2 = isosno(k,2)
             isosno(k,2) = c0
          endif
!         if (isosno(k,2) < c0) then
!             write(nu_diag,*) 'Snow-ice isosno(k,2)',isosno(k,2),sloss2
!         endif
          if (dzssl > puny)                         &
           sloss1 = max(dhs_snoice-dzint,c0)   &
                  *isosno(k,1)/dzssl
          if (isosno(k,1) >= sloss1) then
             isosno(k,1) = isosno(k,1)-sloss1
          else
             sloss1 = isosno(k,1)
             isosno(k,1) = c0
          endif
          isoice(k,1) = isoice(k,1)+sloss1+sloss2
         enddo

         dzssl=dzssl-max(dhs_snoice-dzint,c0)
         dzint=max(dzint-dhs_snoice,c0)
         dzssli=dzssli+snoice
         if ( dzssl <= c0) then ! ssl goes away
          fiso_ocnn(:)= fiso_ocnn(:)+isosno(:,1)/dt
          isosno(:,1) = c0
          dzssl = c0
         endif
         if (dzint <= c0) then  ! int goes away
          fiso_ocnn(:) = fiso_ocnn(:)+isosno(:,2)/dt
          isosno(:,2) = c0
          dzint = c0
         endif
        endif

!     redistribute isotope within vertical layers

         hslyr = hs/real(nslyr,kind=dbl_kind)
         hilyr = hi/real(nilyr,kind=dbl_kind)
         dzssl_new = min(hslyr/c2,hs_ssl)       ! new ssl for snow
         dzint_new = hs-dzssl_new
         dzssli_new = min(hilyr/c2,hi_ssl)      ! new ssl for ice
         dzinti_new = hi-dzssli_new

         do k=1,n_iso

          dznew=min(dzssl_new-dzssl,c0)
          sloss1=c0
          if (dzssl > puny) &
           sloss1=dznew*isosno(k,1)/dzssl ! not neccesarily a loss term
          dznew=max(dzssl_new-dzssl,c0)
          if (dzint > puny) &
           sloss1=sloss1+isosno(k,2)*dznew/dzint ! not really a loss term
          isosno(k,1) =isosno(k,1)+sloss1 
          isosno(k,2) =isosno(k,2)-sloss1
!         if (isosno(k,2) < c0) then
!             write(nu_diag,*) 'redistribute isosno(k,2)',isosno(k,2),sloss1
!             write(nu_diag,*) 'dzssl_new,dzssl,dzint',dzssl_new,dzssl,dzint
!         endif

          sloss2=c0
          dznew=min(dzssli_new-dzssli,c0)
          if (dzssli > puny) & 
           sloss2=dznew*isoice(k,1)/dzssli
          dznew=max(dzssli_new-dzssli,c0)
          if (dzinti > puny) & 
           sloss2=sloss2+isoice(k,2)*dznew/dzinti
          isoice(k,1) = isoice(k,1)+sloss2
          isoice(k,2) = isoice(k,2)-sloss2

          isotot(k)=isosno(k,2)+isosno(k,1)      &
           +isoice(k,2)+isoice(k,1)
!         if ( (isotot(k)-isotot0(k))                 &
!            - fiso_atm(k)*aicen*dt                &
!            - fiso_evapn(k)*dt                         &
!            + fiso_ocnn(k)*dt > 1e-6) then
!           write(nu_diag,*) 'isotope tracer:      ',k
!           write(nu_diag,*) 'isotot-isotot0     ',isotot(k)-isotot0(k) 
!           write(nu_diag,*) 'fiso_atm-fiso_ocnn          ',fiso_atm(k)*aicen*dt &
!                                                        + fiso_evapn(k)*dt &
!                                                        - fiso_ocnn(k)*dt
!         endif

         enddo          ! n_iso

!     reload tracers

!     scale fiso_ocnn. It will be re-scaled by aicen latter in merge_fluxes
      if (aicen > puny) then
         fiso_ocnn(:) = fiso_ocnn(:)/aicen
         fiso_evapn(:) = fiso_evapn(:)/aicen
      endif

      end subroutine update_isotope

!=======================================================================
      function isoice_alpha(growth_rate, sp, frac)

! !DESCRIPTION:
!
! calculate the fractionation coefficient for sea-ice formation
!
! !REVISION HISTORY:
!
! authors: Jiang Zhu, UW-Madison 
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::   &
         growth_rate                     ! sea-ice formation rate (m/s)
      character(*), intent(in) ::   &
         sp,frac                         ! species: H2_16O, H2_18O, HDO
                                         ! calculation methods:
                                         !  cfrac, constant fractionation
                                         !  nfrac, nonfractionation
                                         !  gfrac, growth-rate dependent
      real (kind=dbl_kind) ::   &
         isoice_alpha                    ! return fractionation

      if (frac == 'nfrac') isoice_alpha = c1
      if (sp == 'H2_16O')  isoice_alpha = c1

      ! Lehmann and Siegenthaler, 1991
      !--------------------------------------------------
      if (frac == 'cfrac' .and. sp == 'HDO')            &
         isoice_alpha = 1.02120_dbl_kind
      if (frac == 'cfrac' .and. sp == 'H2_18O')         &
         isoice_alpha = 1.00291_dbl_kind
         
      ! Eq.9, Toyota et al., 2013
      ! For HDO, 7.2852 = 0.2120/0.00291
      !--------------------------------------------------
      if (frac == 'gfrac' .and. sp == 'HDO')                        &
         isoice_alpha = c1+7.2852_dbl_kind*1.2280E-3_dbl_kind+      &
            0.7311E-3_dbl_kind*exp(-growth_rate/8.0100E8_dbl_kind)+ &
            0.8441E-3_dbl_kind*exp(-growth_rate/0.7800E6_dbl_kind)
      if (frac == 'gfrac' .and. sp == 'H2_18O')                     &
         isoice_alpha = c1+1.2280E-3_dbl_kind+                      &
            0.7311E-3_dbl_kind*exp(-growth_rate/8.0100E8_dbl_kind)+ &
            0.8441E-3_dbl_kind*exp(-growth_rate/0.7800E6_dbl_kind)
      return

      end function isoice_alpha

!=======================================================================
      function wiso_alpi(isp,tk)

!-----------------------------------------------------------------------
! Purpose: return ice/vapour fractionation from loop-up tables
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
      integer , intent(in)        :: isp  ! species indes
      real(kind=dbl_kind), intent(in)        :: tk   ! temperature (k)
      real(kind=dbl_kind) :: wiso_alpi               ! return fractionation

!From Merlivat & Nief,1967 for HDO, and Majoube, 1971b for H218O:
      real(kind=dbl_kind), parameter, dimension(pwtspec) :: &  ! ice/vapour
         alpai = (/ 0._dbl_kind, 0._dbl_kind, 16289._dbl_kind,   0._dbl_kind         /), &
         alpbi = (/ 0._dbl_kind, 0._dbl_kind, 0._dbl_kind,       11.839_dbl_kind     /), &
         alpci = (/ 0._dbl_kind, 0._dbl_kind, -9.45e-2_dbl_kind, -28.224e-3_dbl_kind /)

!-----------------------------------------------------------------------
      if (isp == isph2o) then
        wiso_alpi = c1
        return
      end if

      wiso_alpi = exp(alpai(isp)/tk**2 + alpbi(isp)/tk + alpci(isp))

#ifdef NOFRAC
      wiso_alpi = c1
#endif
!
      return
      end function wiso_alpi


!=======================================================================

      end module icepack_isotope

!=======================================================================
