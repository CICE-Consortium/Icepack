!=======================================================================

! sealvl meltpond parameterization
!
! This meltpond parameterization partitions pond volume into depth
! and area based on an assumed, subcategory elevation distribution
! (hypsometry) which permits pond bases to sit below sea level.
! Pond upper surfaces can be below sea level (e.g., when ponds initially
! form in the spring or after mechanical redistribution), but drainage
! processes will not lower pond surfaces below sea level (with the
! default parameter values). Currently, the only physical impacts
! of pond water are in the delta-Eddington radiation scheme.
!
! The sealvl meltpond parameterization was inspired by the level pond
! parameterization of Elizabeth Hunke, David Hebert, and Olivier
! Lecomte. Wherever possible, the code matches the level pond scheme
! (e.g., the pond lid refreezing calculations).
!
! authors David Clemens-Sewall (NCAR/NOAA)

      module icepack_meltpond_sealvl
      
      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c10, p01, p5, puny
      use icepack_parameters, only: viscosity_dyn, rhoi, rhos, rhow 
      use icepack_parameters, only: Timelt, Tffresh, Lfresh, rhofresh
      use icepack_parameters, only: gravit, depressT, rhofresh, kice
      use icepack_parameters, only: rhosi, pndaspect, use_smliq_pnd
      use icepack_parameters, only: ktherm, frzpnd, dpscale, hi_min
      use icepack_parameters, only: pndhyps, pndfrbd, pndhead, apnd_sl
      use icepack_tracers,    only: nilyr
      use icepack_warnings,   only: warnstr, icepack_warnings_add
      use icepack_warnings,   only: icepack_warnings_setabort
      use icepack_warnings,   only: icepack_warnings_aborted

      implicit none

      private 
      public ::   icepack_init_sealvlpnd,   &
                  compute_ponds_sealvl,   &
                  pond_hypsometry,        &
                  pond_height

!=======================================================================

      contains

!=======================================================================

      subroutine icepack_init_sealvlpnd

      use icepack_parameters, only: hpmin, hp0, pndmacr

      ! Set parameters for sealvl pond parameterization
      pndhyps = 'sealevel'
      pndfrbd = 'category'
      pndhead = 'hyps'
      pndmacr = 'head'

      ! Disable hp0 shortwave parameterization
      hp0 = hpmin

      end subroutine icepack_init_sealvlpnd

!=======================================================================

      subroutine compute_ponds_sealvl( dt,                    &
                                       meltt,  melts,  frain, &
                                       Tair,   fsurfn, Tsfcn, &
                                       dhs,    ffrac,         &
                                       aicen,  vicen,  vsnon, &
                                       qicen,  sicen,         &
                                       apnd,   hpnd,  ipnd,   &
                                       meltsliqn, frpndn,     &
                                       ilpndn, flpndn)

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step (s)

      real (kind=dbl_kind), intent(in) :: &
         Tsfcn, &    ! surface temperature (C)
         meltt, &    ! top melt rate (m/s)
         melts, &    ! snow melt rate (m/s)
         frain, &    ! rainfall rate (kg/m2/s)
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &    ! ice area fraction
         vicen, &    ! ice volume (m)
         vsnon, &    ! snow volume (m)
         meltsliqn   ! liquid contribution to meltponds in dt (kg/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         apnd, hpnd, ipnd, & ! pond tracers
         frpndn, &   ! pond drainage rate due to freeboard constraint (m/step)
         ilpndn, &   ! pond loss/gain due to ice lid (m/step)
         flpndn      ! pond flushing rate due to ice permeability (m/s)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         qicen, &    ! ice layer enthalpy (J m-3)
         sicen       ! salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         dhs         ! depth difference for snow on sea ice and pond ice

      real (kind=dbl_kind), intent(out) :: &
         ffrac       ! fraction of fsurfn over pond used to melt ipond

      ! local temporary variables

      real (kind=dbl_kind) :: &
         dhpond, &   ! change in hpond (m)
         dvn_temp    ! local variable for change in volume due to rfrac

      real (kind=dbl_kind), dimension (nilyr) :: &
         Tmlt        ! melting temperature (C)

      real (kind=dbl_kind) :: &
         hi      , & ! ice thickness (m)
         hs      , & ! snow depth (m)
         dTs     , & ! surface temperature diff for freeze-up (C)
         Tp      , & ! pond freezing temperature (C)
         Ts      , & ! surface air temperature (C)
         vpondn  , & ! pond volume per category area (m)
         dvpondn , & ! change in pond volume per category area (m)
         hlid    , & ! refrozen lid thickness
         dhlid   , & ! change in refrozen lid thickness
         bdt     , & ! 2 kice dT dt / (rhoi Lfresh)
         hpsurf  , & ! height of pond surface above mean ice base (m)
         draft, pressure_head, perm, drain ! for permeability

      real (kind=dbl_kind), parameter :: &
         Td   = c2 , & ! temperature difference for freeze-up (C)
         rexp = p01    ! pond contraction scaling

      character(len=*),parameter :: subname='(compute_ponds_sealvl)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      vpondn = hpnd * apnd
      ffrac = c0

      !-----------------------------------------------------------------
      ! Identify grid cells where ponds can be
      !-----------------------------------------------------------------

      if (aicen > puny) then

         hi = vicen/aicen
         hs = vsnon/aicen

         if (hi < hi_min) then

            !-----------------------------------------------------------
            ! Remove ponds on thin ice
            !-----------------------------------------------------------
            frpndn = vpondn
            apnd = c0
            hpnd = c0
            vpondn = c0
            hlid = c0

         else

            !-----------------------------------------------------------
            ! update pond volume
            !-----------------------------------------------------------
            ! add melt water
            if (use_smliq_pnd) then
               dvpondn = (meltt*rhoi + meltsliqn)/rhofresh
            else
               dvpondn = (meltt*rhoi + melts*rhos + frain*dt)/rhofresh
            endif
            dvn_temp = dvpondn

            ! shrink pond volume under freezing conditions
            if (trim(frzpnd) == 'cesm') then
               Tp = Timelt - Td
               dTs = max(Tp - Tsfcn,c0)
               dvpondn = dvpondn - vpondn * (c1 - exp(rexp*dTs/Tp))

            else
               ! trim(frzpnd) == 'hlid' Stefan approximation
               ! assumes pond is fresh (freezing temperature = 0 C)
               ! and ice grows from existing pond ice
               hlid = ipnd
               if (dvpondn == c0) then ! freeze pond
                  Ts = Tair - Tffresh
                  if (Ts < c0) then
                     ! if (Ts < -c2) then ! as in meltpond_cesm
                     bdt = -c2*Ts*kice*dt/(rhoi*Lfresh)
                     dhlid = p5*sqrt(bdt)         ! open water freezing
                     if (hlid > dhlid) dhlid = p5*bdt/hlid ! extant ice
                     dhlid = min(dhlid, hpnd*rhofresh/rhoi)
                     hlid = hlid + dhlid
                  else
                     dhlid = c0 ! to account for surface inversions
                  endif
               else ! convert refrozen pond ice back to water
                  dhlid = max(fsurfn*dt / (rhoi*Lfresh), c0) ! > 0
                  dhlid = -min(dhlid, hlid) ! < 0
                  hlid = max(hlid + dhlid, c0)
                  if (hs - dhs < puny) then ! pond ice is snow-free
                     ! fraction of fsurfn over pond used to melt ipond
                     ffrac = c1 
                     if (fsurfn > puny) &
                        ffrac = min(-dhlid*rhoi*Lfresh/(dt*fsurfn), c1)
                  endif
               endif
               dvpondn = dvpondn - dhlid*apnd*rhoi/rhofresh
            endif

            ! Track lost/gained meltwater per unit category area from 
            ! pond lid freezing/melting. Note sign flip relative to dvn
            ilpndn = dvn_temp - dvpondn
            
            !-----------------------------------------------------------
            ! update pond area and depth
            !-----------------------------------------------------------
            call pond_hypsometry(hpnd, apnd, dvpond=dvpondn, hin=hi)

            ! limit pond depth to maintain nonnegative freeboard
            if (trim(pndfrbd) == 'floor') then
               dhpond = ((rhow-rhoi)*hi - rhos*hs)/rhofresh - hpnd
            elseif (trim(pndfrbd) == 'category') then
               dhpond = ((rhow-rhoi)*hi-rhos*hs)/(rhofresh*apnd) &
                  - hpnd
            else
               call icepack_warnings_add(subname// &
                  " invalid pndfrbd option" )
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               if (icepack_warnings_aborted(subname)) return
            endif
            dhpond = min(dhpond, c0) ! strictly drainage
            frpndn = - dhpond * apnd
            call pond_hypsometry(hpnd, apnd, dhpond=dhpond, hin=hi)
            
            ! clean up empty ponds. Note, this implies that if ponds 
            ! fully drain or freeze, the lid ice also ceases to exist
            if (hpnd <= c0 .or. apnd <= c0) then
               apnd = c0
               hpnd = c0
               hlid = c0
            endif

            !-----------------------------------------------------------
            ! drainage due to permeability (flushing)
            ! setting dpscale = 0 turns this off
            ! NOTE this uses the initial salinity and melting T profiles
            !-----------------------------------------------------------

            if (ktherm /= 2 .and. hpnd > c0 .and. dpscale > puny) then
               draft = (rhos*hs + rhoi*hi + rhofresh*hpnd*apnd)/rhow
               call pond_height(apnd, hpnd, hi, hpsurf)
               pressure_head = gravit * rhow * max(hpsurf - draft, c0)
               Tmlt(:) = -sicen(:) * depressT
               call brine_permeability(qicen, &
                    sicen, Tmlt, perm)
               if (icepack_warnings_aborted(subname)) return
               drain = perm*pressure_head*dt/(viscosity_dyn*hi)*dpscale
               dhpond = -min(drain, hpnd)
               flpndn = -dhpond * apnd               
               call pond_hypsometry(hpnd, apnd, dhpond=dhpond, hin=hi)
            endif
         endif ! hi < hi_min

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------
         if (trim(frzpnd) == 'hlid') ipnd = hlid

      endif

      end subroutine compute_ponds_sealvl

!=======================================================================

! determine the liquid fraction of brine in the ice and the permeability

      subroutine brine_permeability(qicen, salin, Tmlt, perm)

      use icepack_therm_shared, only: calculate_Tin_from_qin

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         qicen, &  ! enthalpy for each ice layer (J m-3)
         salin, &  ! salinity (ppt)
         Tmlt      ! melting temperature (C)

      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability (m^2)

      ! local variables

      real (kind=dbl_kind) ::   &
         Sbr       ! brine salinity

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature (C)
         phi       ! liquid fraction

      integer (kind=int_kind) :: k

      character(len=*),parameter :: subname='(brine_permeability)'

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Tin(k) = calculate_Tin_from_qin(qicen(k),Tmlt(k))
      enddo

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Sbr = c1/(1.e-3_dbl_kind - depressT/Tin(k)) ! Notz thesis eq 3.6
         phi(k) = salin(k)/Sbr ! liquid fraction
         if (phi(k) < 0.05) phi(k) = c0 ! impermeable
      enddo

      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-8_dbl_kind * (minval(phi))**3

      end subroutine brine_permeability

!=======================================================================

! compute the changes in pond area and depth

      subroutine pond_hypsometry(hpnd, apnd, dhpond, dvpond, hin)

      real (kind=dbl_kind), intent(inout) :: &
         hpnd, &   ! pond depth of ponded area tracer
         apnd      ! pond fractional area of category tracer
         
      real (kind=dbl_kind), intent(in), optional :: &
         dvpond, & ! incoming change in pond volume per category area
         dhpond, & ! incoming change in pond depth
         hin       ! category ice thickness
      
      ! local variables
      
      real (kind=dbl_kind) :: &
         dv, &     ! local variable for change in pond volume
         vp        ! local variable for pond volume per category area
      
      character(len=*),parameter :: subname='(pond_hypsometry)'

      ! Behavior is undefined if dhpond and dvpond are both present
      if (present(dhpond) .and. present(dvpond)) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname// &
            " dhpond and dvpond cannot both be input" )
         return
      endif
      ! or both absent
      if ((.not. present(dhpond)) .and. (.not. present(dvpond))) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname// &
            " dhpond and dvpond cannot both be absent" )
         return
      endif
      ! Check that category ice thickness is present if needed
      if ((trim(pndhyps) == 'sealevel') .and. (.not. present(hin))) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname// &
            " hin needed for sealevel ponds")
         return
      endif
      
      ! Get the change in volume
      if (present(dvpond)) then
         dv = dvpond
      else ! compute change in volume from change in depth
         dv = dhpond * apnd
      endif
      ! compute initial pond volume from area and depth
      vp = apnd * hpnd
      vp = vp + dv ! update pond volume
      ! Compute pond area assuming that apond*pndaspect = hpond
      if (vp <= puny) then
         apnd = c0
         hpnd = c0
      else
         call calc_pndaspect(hin)
         apnd = sqrt(vp/pndaspect)
         ! preserve pond volume if pond fills all available area
         hpnd = c0
         if (apnd < c1) then
            hpnd = apnd * pndaspect
         else
            apnd = 1 ! ponds fill entire category
            hpnd = vp ! conserve volume
         endif
      endif
      
      end subroutine pond_hypsometry

!=======================================================================

! compute the height of pond upper surface above mean base of ice

      subroutine pond_height(apond, hpnd, hin, hpsurf)

         real (kind=dbl_kind), intent(in) :: &
            hin  , & ! category mean ice thickness
            apond , & ! pond area fraction of the category
            hpnd     ! mean pond depth (m)
   
         real (kind=dbl_kind), intent(out) :: &
            hpsurf   ! height of pond surface above base of the ice (m)
         
         character(len=*),parameter :: subname='(pond_height)'
   
         if (trim(pndhead) == 'perched') then
            hpsurf = hin + hpnd
         elseif (trim(pndhead) == 'hyps') then
            if ((trim(pndhyps) == 'fixed') .or. &
               (trim(pndhyps) == 'sealevel')) then
               ! Applying a fixed aspect ratio to the ponds implicitly 
               ! assumes that the hypsometric curve has a constant slope
               ! of double the aspect ratio.
               ! If ponds occupy lowest elevations first. 
               call calc_pndaspect(hin)
               if (apond < c1) then
                  hpsurf = hin - pndaspect + c2*pndaspect*apond
               else ! ponds cover all available area
                  hpsurf = hin + hpnd
               endif
            else
               call icepack_warnings_add(subname// &
                  " unsupported pndhyps option" )
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               if (icepack_warnings_aborted(subname)) return
            endif
         else
            call icepack_warnings_add(subname// &
               " invalid pndhead option" )
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            if (icepack_warnings_aborted(subname)) return
         endif
   
      end subroutine pond_height

!=======================================================================

      subroutine calc_pndaspect(hin)

         real (kind=dbl_kind), intent(in), optional :: &
            hin   ! category mean ice thickness (m)

         character(len=*),parameter :: subname='(sealevel_pndaspect)'

         ! Compute the pond aspect ratio for sea level ponds
         if (trim(pndhyps) == 'sealevel') then
            if (.not. present(hin)) then
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname// &
                  " hin needed for sealevel ponds")
               return
            endif
            pndaspect = hin*(rhow - rhosi) / &
               (rhofresh*apnd_sl**c2 - c2*rhow*apnd_sl + rhow)
         endif ! Otherwise do nothing to pond aspect

      end subroutine calc_pndaspect

   end module icepack_meltpond_sealvl
!=======================================================================