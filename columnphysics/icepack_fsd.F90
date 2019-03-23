!=========================================================================
!
! This module contains the subroutines required to define
! a floe size distribution tracer for sea ice
!
! Variable naming convention
! for k = 1, nfsd and n = 1, ncat
!    afsdn(k,n) = trcrn(:,:,nt_nfsd+k-1,n,:)
!    afsd (k) is per thickness category or averaged over n
!    afs is the associated scalar value for (k,n)
!
! authors: liuxy
!          C. M. Bitz, UW
!          Lettie Roach, NIWA
!
! 2015: liuxy modified from ice_fsd module
! 2016: CMB rewrote a lot of it
! 2016: LR made some modifications
! 2019: ECH ported to Icepack
 
      module icepack_fsd

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c3, c4, p01, p1, p5, puny
      use icepack_parameters, only: pi, floeshape, wave_spec, bignum, gravit, rhoi
      use icepack_tracers, only: nt_fsd, tr_fsd
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: icepack_init_fsd_bounds, icepack_init_fsd, icepack_cleanup_fsd, &
         fsd_lateral_growth, fsd_add_new_ice, fsd_weld_thermo

!      logical (kind=log_kind), public :: & 
!         write_diag_wave     ! if .true., write the lat/lons from find_wave to history file

! move to namelist? domain_size?
      integer(kind=int_kind) ::  &
         nwavefreq = 25        ! number of frequencies in wave spectrum

      real (kind=dbl_kind), public :: &
         c_weld = p01     ! constant of proportionality for welding
                          ! see documentation for details

!      logical (kind=log_kind), public :: &
!         rdc_frzmlt      ! if true, only (1-oo) of frzmlt can melt (lat and bot)

      real(kind=dbl_kind), dimension(:), allocatable ::  &
         floe_rad_h,         & ! fsd size higher bound in m (radius)
         floe_area_l,        & ! fsd area at lower bound (m^2)
         floe_area_h,        & ! fsd area at higher bound (m^2)
         floe_area_c,        & ! fsd area at bin centre (m^2)
         floe_area_binwidth, & ! floe area bin width (m^2)
         area_scaled_l,      & ! area bins scaled so they begin at zero
         area_scaled_h,      & ! and no binwidth is greater than 1
         area_scaled_c,      & ! (dimensionless)
         area_scaled_binwidth

      integer(kind=int_kind), dimension(:,:), allocatable, public ::  &
         iweld

!=======================================================================

      contains

!=======================================================================

!  Initialize ice fsd bounds (call whether or not restarting)
!  Define the bounds, midpoints and widths of floe size
!  categories in area and radius
!
!  Note that radius widths cannot be larger than twice previous
!  category width or floe welding will not have an effect
!
!  Note also that the bound of the lowest floe size category is used
!  to define the lead region width and the domain spacing for wave breaking
!
      subroutine icepack_init_fsd_bounds(nfsd, &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth, &  ! fsd size bin width in m (radius)
         c_fsd_range)      ! string for history output

      integer (kind=int_kind), intent(in) :: &
         nfsd              ! number of floe size categories

      real(kind=dbl_kind), dimension(:), intent(inout) ::  &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth     ! fsd size bin width in m (radius)

      character (len=35), intent(out) :: &
           c_fsd_range(nfsd) ! string for history output

      ! local variables

      integer (kind=int_kind) :: n, m, k
      integer (kind=int_kind) :: ierr

      real (kind=dbl_kind) :: test

      real (kind=dbl_kind), dimension (nfsd+1) :: &
         area_lims, area_lims_scaled

      real (kind=dbl_kind), dimension (0:nfsd) :: &
         floe_rad
                                              
      real (kind=dbl_kind), dimension(:), allocatable :: &
         lims

      character(len=8) :: c_fsd1,c_fsd2
      character(len=2) :: c_nf
      character(len=*), parameter :: subname='(init_fsd_bounds)'

      if (nfsd.eq.24) then

         allocate(lims(24+1))

         lims = (/ 6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                   5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                   3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                   9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                   3.35434988e+03,   4.55051413e+03,   6.17323164e+03,   8.37461170e+03, &
                   1.13610059e+04,   1.54123510e+04,   2.09084095e+04,   2.83643675e+04, &
                   3.84791270e+04 /)
        
      elseif (nfsd.eq.16) then

         allocate(lims(16+1))

         lims = (/ 6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                   5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                   3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                   9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                   3.35434988e+03 /)
        
      else if (nfsd.eq.12) then

         allocate(lims(12+1))

         lims = (/ 6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                   5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                   3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                   9.45812834e+02 /)
 
      else if (nfsd.eq.1) then ! default case

         allocate(lims(1+1))

         lims = (/ 6.65000000e-02,   3.0e+02 /)

      else

         call icepack_warnings_add(subname//&
            ' floe size categories not defined for nfsd')
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__) 
         return

      end if

      allocate(                                                   &
         floe_rad_h          (nfsd), & ! fsd size higher bound in m (radius)
         floe_area_l         (nfsd), & ! fsd area at lower bound (m^2)
         floe_area_h         (nfsd), & ! fsd area at higher bound (m^2)
         floe_area_c         (nfsd), & ! fsd area at bin centre (m^2)
         floe_area_binwidth  (nfsd), & ! floe area bin width (m^2)
         area_scaled_l       (nfsd), & ! area bins scaled so they begin at zero
         area_scaled_h       (nfsd), & ! and no binwidth is greater than 1
         area_scaled_c       (nfsd), & ! (dimensionless)
         area_scaled_binwidth(nfsd), & !
         iweld         (nfsd, nfsd), & ! fsd categories that can weld
         stat=ierr)
      if (ierr/=0) then
         call icepack_warnings_add(subname//' Out of Memory fsd')
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__) 
         return
      endif

      floe_rad_l = lims(1:nfsd  )
      floe_rad_h = lims(2:nfsd+1)
      floe_rad_c = (floe_rad_h+floe_rad_l)/c2

      floe_area_l = c4*floeshape*floe_rad_l**c2
      floe_area_c = c4*floeshape*floe_rad_c**c2
      floe_area_h = c4*floeshape*floe_rad_h**c2

      floe_binwidth = floe_rad_h - floe_rad_l

      ! scaled area for welding
      floe_area_binwidth = floe_area_h - floe_area_l
      area_lims(1:nfsd) = floe_area_l(1:nfsd)
      area_lims(nfsd+1) = floe_area_h(nfsd)
      area_lims_scaled = (area_lims - area_lims(1))/MAXVAL(floe_area_binwidth)

      area_scaled_h = area_lims_scaled(2:nfsd+1)
      area_scaled_l = area_lims_scaled(1:nfsd  )
      area_scaled_c = (area_scaled_h + area_scaled_l) / c2
      area_scaled_binwidth = area_scaled_h - area_scaled_l

      ! floe size categories that can combine during welding
      iweld(:,:) = -999
      do n = 1, nfsd
      do m = 1, nfsd
         test = area_scaled_h(n) - area_scaled_c(m)
         do k = 1, nfsd
            if ((test >= area_scaled_l(k)) .and. (test < area_scaled_h(k))) then
               iweld(n,m) = k + 1
            end if
         end do
      end do
      end do

      if (allocated(lims)) deallocate(lims)

      ! write fsd bounds
      floe_rad(0) = floe_rad_l(1)
      do n = 1, nfsd
         floe_rad(n) = floe_rad_h(n)
      enddo

         write(warnstr,*) ' '
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname
         call icepack_warnings_add(warnstr)
         write(warnstr,*) 'floe_rad(n-1) < fsd Cat n < floe_rad(n)'
         call icepack_warnings_add(warnstr)
         do n = 1, nfsd
            write(warnstr,*) floe_rad(n-1),' < fsd Cat ',n, ' < ',floe_rad(n)
            call icepack_warnings_add(warnstr)
            ! Write integer n to character string
            write (c_nf, '(i2)') n    

            ! Write floe_rad to character string
            write (c_fsd1, '(f6.3)') floe_rad(n-1)
            write (c_fsd2, '(f6.3)') floe_rad(n)

            ! Save character string to write to history file
            c_fsd_range(n)=c_fsd1//'m < fsd Cat '//c_nf//' < '//c_fsd2//'m'
         enddo

         write(warnstr,*) ' '
         call icepack_warnings_add(warnstr)

       
      end subroutine icepack_init_fsd_bounds

!=======================================================================
!
!  Initialize the FSD 
!  When growing from no-ice conditions, initialize to zero
!  Otherwise initalize with a power law, following Perovich (2014)

      subroutine icepack_init_fsd(nfsd, ice_ic, &
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth, &  ! fsd size bin width in m (radius)
         afsd)             ! floe size distribution tracer

      integer(kind=int_kind), intent(in) :: &
         nfsd

      character(len=char_len_long), intent(in) :: &
         ice_ic           ! method of ice cover initialization

      real(kind=dbl_kind), dimension(:), intent(inout) ::  &
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth     ! fsd size bin width in m (radius)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         afsd              ! floe size tracer: fraction distribution of floes

      ! local variables

      real (kind=dbl_kind) :: alpha, totfrac

      integer (kind=int_kind) :: k

      real  (kind=dbl_kind), dimension (nfsd) :: &
         num_fsd           ! number distribution of floes

      if (trim(ice_ic) == 'none') then

         afsd(:) = c0

      else            ! Perovich (2014)
 
         ! fraction of ice in each floe size and thickness category
         ! same for ALL cells (even where no ice) initially
         alpha = 2.1_dbl_kind
         totfrac = c0                                   ! total fraction of floes 
         do k = 1, nfsd
            num_fsd(k) = (2*floe_rad_c(k))**(-alpha-c1) ! number distribution of floes
            afsd   (k) = num_fsd(k)*floe_area_c(k)*floe_binwidth(k) ! fraction distribution of floes
            totfrac = totfrac + afsd(k)
         enddo
         afsd = afsd/totfrac                    ! normalize

      endif ! ice_ic

      end subroutine icepack_init_fsd

!=======================================================================
!
!  Clean up small values and renormalize
!
!  author:  Elizabeth Hunke, LANL

      subroutine icepack_cleanup_fsd (ncat, nfsd, afsdn)

      integer (kind=int_kind), intent(in) :: &
         ncat           , & ! number of thickness categories
         nfsd               ! number of floe size categories

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         afsdn              ! floe size distribution tracer

      ! local variables

      integer (kind=int_kind) :: &
         n, k               ! thickness and floe size category indices

      real (kind=dbl_kind) :: &
         tot

      if (tr_fsd) then

      do n = 1, ncat
         do k = 1, nfsd
            if (afsdn(k,n) < puny) afsdn(k,n) = c0
         enddo

         tot = sum(afsdn(:,n))
         if (tot > puny) then
            do k = 1, nfsd
               afsdn(k,n) = afsdn(k,n) / tot ! normalize
            enddo
         else
            afsdn(:,n) = c0
            !afsdn(1,n) = c1                  ! default to smallest floe size
            !do k = 2, nfsd
            !   afsdn(k,n) = c0
            !enddo
         endif
      enddo ! ncat

      endif ! tr_fsd

      end subroutine icepack_cleanup_fsd

!=======================================================================
! 
!  Given the joint ice thickness and floe size distribution, calculate
!  the lead region and the total lateral surface area following Horvat
!  and Tziperman (2015)
!
! author: Lettie Roach, NIWA/VUW

      subroutine partition_area (ncat,       nfsd,      &
                                 floe_rad_c, aice,      &
                                 aicen,      vicen,     &
                                 afsdn,      lead_area, &
                                 latsurf_area)

      integer (kind=int_kind), intent(in) :: &
         ncat           , & ! number of thickness categories
         nfsd               ! number of floe size categories

      real (kind=dbl_kind), dimension(:), intent(in) ::  &
         floe_rad_c         ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), intent(in) :: &
         aice               ! ice concentration

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen          , & ! fractional area of ice
         vicen              ! volume per unit area of ice

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         afsdn              ! floe size distribution tracer

      real (kind=dbl_kind), intent(out) :: &
         lead_area      , & ! the fractional area of the lead region
         latsurf_area       ! the area covered by lateral surface of floes

      ! local variables

      integer (kind=int_kind) :: &
         n              , & ! thickness category index
         k                  ! floe size index

      real (kind=dbl_kind) :: &
         width_leadreg, &   ! width of lead region
         thickness          ! actual thickness of ice in thickness cat

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      lead_area    = c0
      latsurf_area = c0

      ! Set the width of the lead region to be the smallest
      ! floe size category, as per Horvat & Tziperman (2015)
      width_leadreg = floe_rad_c(1)

      ! Only calculate these areas where there is some ice
      if (aice > puny) then

         ! lead area = sum of areas of annuli around all floes
         do n = 1, ncat
            do k = 1, nfsd
               lead_area = lead_area + aicen(n) * afsdn(k,n) &
                         * (c2*width_leadreg    /floe_rad_c(k)     &
                             + width_leadreg**c2/floe_rad_c(k)**2)
            enddo ! k
         enddo    ! n

         ! cannot be greater than the open water fraction
         lead_area=MIN(lead_area, c1-aice)

                ! sanity checks
!                if (lead_area.gt.c1) stop 'lead_area not frac!'
!                if (lead_area.ne.lead_area) stop 'lead_a NaN'

                if (lead_area.lt.c0) then
                        if (lead_area.lt.(c0-puny)) then
!                                stop 'lead_area lt0 in partition_area'
                        else
                                lead_area=MAX(lead_area,c0)
                        end if
                end if

         ! Total fractional lateral surface area in each grid (per unit ocean area)
         do n = 1, ncat
            thickness = c0
            if (aicen(n) > c0) thickness = vicen(n)/aicen(n)

            do k = 1, nfsd
               latsurf_area = latsurf_area &
                   + afsdn(k,n) * aicen(n) * c2 * thickness/floe_rad_c(k)
            end do
         end do

                ! check
                if (latsurf_area.lt.c0) stop 'negative latsurf_area'
                if (latsurf_area.ne.latsurf_area) stop &
                   'latsurf_area NaN'

      end if ! aice

      end subroutine partition_area

!=======================================================================

      subroutine fsd_lateral_growth (ncat,      nfsd,         &
                                     dt,        aice,         &
                                     aicen,     vicen,        &
                                     vi0new,                  &
                                     frazil,    floe_rad_c,   &
                                     afsdn,                   &
                                     lead_area, latsurf_area, &
                                     G_radial,  d_an_latg,    &
                                     tot_latg)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nfsd      ! number of floe size categories

      real (kind=dbl_kind), intent(in) :: &
         dt    , & ! time step (s)
         aice      ! total concentration of ice

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         afsdn              ! floe size distribution tracer

      real (kind=dbl_kind), intent(inout) :: &
         vi0new, & ! volume of new ice added to cat 1
         frazil    ! frazil ice growth        (m/step-->cm/day)

      ! floe size distribution
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         floe_rad_c    ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), dimension(ncat), intent(out) :: &
!      real (kind=dbl_kind), dimension(:), intent(out) :: &
         d_an_latg

      real (kind=dbl_kind), intent(out) :: &
         G_radial    , & ! lateral melt rate (m/s)
         tot_latg        ! total fsd lateral growth in open water

!!      real (kind=dbl_kind), dimension(nfsd), intent(out) :: &
!      real (kind=dbl_kind), dimension(:), intent(out) :: &
!         d_afsd_latg, d_afsd_newi

      ! local variables

      integer (kind=int_kind) :: &
         n, k             ! ice category indices

      real (kind=dbl_kind) :: &
         vi0new_lat       ! volume of new ice added laterally to fsd

      real (kind=dbl_kind), intent(out) :: &
         lead_area      , & ! the fractional area of the lead region
         latsurf_area       ! the area covered by lateral surface of floes

      character(len=*),parameter :: subname='(fsd_lateral_growth)'

      lead_area    = c0
      latsurf_area = c0
      G_radial     = c0
      tot_latg     = c0
      d_an_latg    = c0

      ! partition volume into lateral growth and frazil
      call partition_area (ncat,       nfsd,      &
                           floe_rad_c, aice,      &
                           aicen,      vicen,     &
                           afsdn,      lead_area, &
                           latsurf_area)

      vi0new_lat = c0
      if (latsurf_area > puny) then
         vi0new_lat = vi0new * lead_area / (c1 + aice/latsurf_area)
      end if

!      if (vi0new_lat < c0) print*,'ERROR vi0new_lat < 0', vi0new_lat

      ! for history/diagnostics
      frazil = vi0new - vi0new_lat

      ! lateral growth increment
      if (vi0new_lat > puny) then
         G_radial = vi0new_lat/dt
         do n = 1, ncat
            if (aicen(n) > puny) then
               afsdn(:,n) = afsdn(:,n)/SUM(afsdn(:,n)) ! in case of 10e-10 errors
            end if

            do k = 1, nfsd
               d_an_latg(n) = d_an_latg(n) &
                            + c2*aicen(n)*afsdn(k,n)*G_radial*dt/floe_rad_c(k)
            end do
         end do ! n
      endif ! vi0new_lat > 0

      ! Use remaining ice volume as in standard model,
      ! but ice cannot grow into the area that has grown laterally

      vi0new = vi0new - vi0new_lat
      tot_latg = SUM(d_an_latg(:))

      end subroutine fsd_lateral_growth

!=======================================================================

      ! lateral growth of existing ice and growth of new ice in category 1

      subroutine fsd_add_new_ice (ncat, n,    nfsd,          &
                                  dt,         ai0new,        &
                                  d_an_latg,  d_an_newi,   &
                                  floe_rad_c, floe_binwidth, &
                                  G_radial,   area2,         &
                                  ice_wave_sig_ht,           &
                                  d_afsd_latg,               &
                                  d_afsd_newi,               &
                                  afsdn,      aicen_init,    &
                                  aicen,      trcrn)

      integer (kind=int_kind), intent(in) :: &
         n     , & ! thickness category number
         ncat  , & ! number of thickness categories
         nfsd      ! number of floe size categories

      real (kind=dbl_kind), intent(in) :: &
         dt           , & ! time step (s)
         ai0new       , & ! area of new ice added to cat 1
         G_radial     , & ! lateral melt rate (m/s)
         ice_wave_sig_ht  ! wave significant height in ice

!      real (kind=dbl_kind), dimension(nwavefreq), intent(in)  :: &
!      real (kind=dbl_kind), dimension(:), intent(in)  :: &
      real (kind=dbl_kind), dimension(nwavefreq)  :: &
         wave_spectrum

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         d_an_latg, d_an_newi

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen_init     , & ! fractional area of ice
         floe_rad_c     , & ! fsd size bin centre in m (radius)
         floe_binwidth      ! fsd size bin width in m (radius)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         afsdn     ! floe size distribution tracer (originally areal_mfstd_init)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         area2 , & ! area after lateral growth and before new ice formation
         aicen     ! concentration of ice

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         d_afsd_latg    , & ! change in floe size distribution (area)
         d_afsd_newi        ! due to fsd lateral growth

      ! local variables

!      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
!         d_afsdn_latg, d_afsdn_newi

      integer (kind=int_kind) :: &
         k             ! floe size category index

      real (kind=dbl_kind), dimension (nfsd,ncat) :: &
         afsdn_latg    ! after lateral growth

      real (kind=dbl_kind), dimension (nfsd) :: &
         df_flx, &     ! finite differences for G_r*tilda(L)
         afsd_ni       ! areal mFSTD after new ice added

      real (kind=dbl_kind), dimension(nfsd+1) :: &
         f_flx         !

      integer (kind=int_kind) :: &
         new_size      ! index for floe size of new ice

      character(len=*),parameter :: subname='(fsd_add_new_ice)'

! for now
      wave_spectrum(:) = c0

      afsdn_latg(:,n) = afsdn(:,n)  ! default

      if (d_an_latg(n) > puny) then ! lateral growth

         area2(n) = aicen_init(n) + d_an_latg(n) ! area after lateral growth, before new ice forms

         df_flx(:) = c0 ! NB could stay zero if all in largest FS cat
         f_flx (:) = c0
         do k = 2, nfsd
            f_flx(k) = G_radial * afsdn(k-1,n) / floe_binwidth(k-1)
         end do
         do k = 1, nfsd
            df_flx(k) = f_flx(k+1) - f_flx(k)
         end do

!         if (abs(sum(df_flx)) > puny) print*,'fsd_add_new ERROR df_flx /= 0'

         afsdn_latg(:,n) = c0
         do k = 1, nfsd
            afsdn_latg(k,n) = afsdn(k,n) &
                            + dt * (-df_flx(k) + c2 * G_radial * afsdn(k,n) &
                            * (c1/floe_rad_c(k) - SUM(afsdn(:,n)/floe_rad_c(:))) )
         end do

         ! combining the next 2 lines changes the answers
         afsdn_latg(:,n) = afsdn_latg(:,n) / SUM(afsdn_latg(:,n)) ! just in case (may be errors < 1e-11)
         !trcrn(nt_fsd:nt_fsd+nfsd-1,n) = afsdn_latg(:,n)

!         d_afsdn_latg(:,n) = afsdn_latg(:,n) - afsdn(:,n)

!         if (abs(sum(afsdn(:,n))-c1) > puny) print*,'fsd_add_new WARNING afsdn not normal'
!         if (abs(sum(afsdn_latg(:,n))-c1) > puny) print*,'fsd_add_new WARNING afsdn_latg not normal'
!         if (any(afsdn_latg < c0)) print*,'fsd_add_new ERROR afsdn_latg < 0'

      else ! no lateral growth, no change in floe size

         afsdn_latg(:,n) = afsdn(:,n)

      end if ! lat growth

      afsd_ni(:) = c0

      new_size = nfsd
      if (n == 1) then
         ! add new frazil ice to smallest thickness
         if (d_an_newi(n) > puny) then

!            if (d_an_newi(n) > aicen(n)) print*,'fsd_add_new ERROR d_an_newi > aicen'

            if (SUM(afsdn_latg(:,n)) > puny) then ! fsd lateral growth occurred

               if (wave_spec) then
                  if (ice_wave_sig_ht > puny) &
                     call wave_dep_growth(nfsd, wave_spectrum(:), new_size)

                  ! grow in new_size category
                  afsd_ni(new_size) = (afsdn_latg(new_size,n)*area2(n) + ai0new) &
                                                          / (area2(n) + ai0new)
                  do k = 1, new_size-1  ! diminish other floe cats accordingly
                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n) + ai0new)
                  end do
                  do k = new_size+1, nfsd  ! diminish other floe cats accordingly
                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n) + ai0new)
                  end do

! add this option later, maybe
!               else if () then ! grow in largest category
!                  afsd_ni(nfsd) =  (afsdn_latg(nfsd,n)*area2(n) + ai0new) &
!                                                    / (area2(n) + ai0new)
!                  do k = 1, nfsd-1  ! diminish other floe cats accordingly
!                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n)+ai0new)
!                  enddo

               else ! grow in smallest floe size category
                  afsd_ni(1) = (afsdn_latg(1,n)*area2(n) + ai0new) &
                                             / (area2(n) + ai0new)
                  do k = 2, nfsd  ! diminish other floe cats accordingly
                     afsd_ni(k) = afsdn_latg(k,n)*area2(n) / (area2(n)+ai0new)
                  enddo
               end if ! new_fs_option

            else ! entirely new ice or not

               if (wave_spec) then
                  if (ice_wave_sig_ht > puny) &
                     call wave_dep_growth(nfsd, wave_spectrum(:), new_size)
                  afsd_ni(new_size) = c1
               else
                  afsd_ni(1) = c1
!               elseif () then ! grow in largest category
!                  afsd_ni(nfsd) = c1
               endif      ! wave forcing

            endif ! entirely new ice or not

!            if (abs(sum(afsd_ni)-c1) > puny) print*,'fsd_add_new afsd_ni not normal',abs(sum(afsd_ni))
!            if (any(afsd_ni < c0)) print*,'fsd_add_new afsd_ni < 0', afsd_ni

            afsd_ni(:) = afsd_ni(:) / SUM(afsd_ni(:))
            trcrn(nt_fsd:nt_fsd+nfsd-1,n) = afsd_ni(:)

            ! for history/diagnostics
!            d_afsdn_newi(:,n) = afsd_ni(:) - afsdn_latg(:,n)

         endif ! d_an_newi > puny
      endif    ! n = 1

      ! history/diagnostics
      do k = 1, nfsd
         d_afsd_latg(k) = d_afsd_latg(k) - aicen_init(n)*afsdn(k,n) &
                + (aicen_init(n) + d_an_latg(n))*afsdn_latg(k,n)
         d_afsd_newi(k) = d_afsd_newi(k) + aicen(n)*trcrn(nt_fsd+k-1,n) &
                - (aicen_init(n) + d_an_latg(n))*afsdn_latg(k,n)
      enddo    ! k

      end subroutine fsd_add_new_ice

!=======================================================================
!
! Given a wave spectrum, calculate size of new floes based on tensile failire
! See Shen & Ackley (2004), Roach, Smith & Dean (2018) for further details
! Author: Lettie Roach (NIWA) 2018

      subroutine wave_dep_growth (nfsd, local_wave_spec, new_size)

      integer (kind=int_kind), intent(in) :: &
         nfsd              ! number of floe size categories

      real (kind=dbl_kind), dimension(nwavefreq), intent(in) :: &
         local_wave_spec ! e(f), dimension set in ice_forcing

      integer (kind=int_kind), intent(out) :: &
         new_size ! index of floe size category in which new floes will growh

      ! local variables
      real (kind=dbl_kind), parameter :: &
         tensile_param = 0.167_dbl_kind

      real (kind=dbl_kind)  :: &
         mom0,   & ! zeroth moment of the spectrum (m)
         h_sig,  & ! significant wave height (m)
         w_amp,  & ! wave amplitude (m)
         f_peak, & ! peak frequency (s^-1)
         w_peak, & ! wavelength from peak freqency (m)
         r_max     ! radius

      integer (kind=int_kind) :: k

      ! ----- begin forcing --- these values will come from a coupler or forcing data
      real (kind=dbl_kind), dimension(nwavefreq) :: &
         wavefreq, & ! wave frequencies
         dwavefreq   ! wave frequency bin widths

      ! hardwired for wave coupling with our version of Wavewatch
      ! from Wavewatch, set as f(n+1) = C*f(n) where C is a constant set by the user, typically ~ 1.1.
      ! these freq are for C = 1.1
      wavefreq = (/0.04118,     0.045298,    0.0498278,   0.05481058,  0.06029164,  0.06632081, &
                   0.07295289,  0.08024818,  0.08827299,  0.09710029,  0.10681032,  0.11749136, &
                   0.1292405,   0.14216454,  0.15638101,  0.17201911,  0.18922101,  0.20814312, &
                   0.22895744,  0.25185317,  0.27703848,  0.30474234,  0.33521661,  0.36873826, &
                   0.40561208/)

      ! boundaries of bin n are at f(n)*sqrt(1/C) and f(n)*sqrt(C)
      dwavefreq(:) = wavefreq(:)*(SQRT(1.1_dbl_kind) - SQRT(c1/1.1_dbl_kind))
      ! ----- end forcing ---

      mom0 = SUM(local_wave_spec*dwavefreq)                   ! zeroth moment
      h_sig = c4*SQRT(mom0)                                   ! sig wave height
      w_amp = h_sig/c2                                        ! sig wave amplitude
      f_peak = wavefreq(MAXLOC(local_wave_spec, DIM=1))       ! peak frequency
      if (f_peak > puny) w_peak = gravit / (c2*pi*f_peak**c2) ! wavelength from peak freq

      ! tensile failure
      if (w_amp > puny) then
         r_max = SQRT(c2*tensile_param*w_peak**c2/(pi**c3*w_amp*gravit*rhoi))/c2
      else
         r_max = bignum
      end if

      new_size = nfsd
      do k = nfsd-1, 1, -1
         if (r_max <= floe_rad_h(k)) new_size = k
      end do

      end subroutine wave_dep_growth

!=======================================================================

      subroutine fsd_weld_thermo (ncat,  nfsd,   &
                                  dt,    frzmlt, &
                                  aicen, trcrn,  &
                                  d_afsd_weld)

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nfsd         ! number of floe size categories

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step (s)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen        ! ice concentration

      real (kind=dbl_kind), intent(in) :: &
         frzmlt       ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         d_afsd_weld  ! change in fsd due to welding

      real (kind=dbl_kind), parameter :: &
         aminweld = p1 ! minimum ice concentration likely to weld

      ! local variables

      integer (kind=int_kind) :: &
        nt        , & ! time step index
        n         , & ! thickness category index
        k, kx, ky     ! floe size category indices

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         afsdn        ! floe size distribution tracer

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         d_afsdn_weld ! change in afsdn due to welding

      real (kind=dbl_kind), dimension(nfsd) :: &
         afsd_init , & ! initial values
         afsd_tmp  , & ! work array
         coag          !

      real(kind=dbl_kind) :: &
         afsd_sum  , & ! work array
         subdt     , & ! subcycling time step for stability
         darea     , & ! total area lost due to welding
         darea_nfsd, & ! area lost from category nfsd
         stability     ! to satisfy stability condition for Smol. eqn

      integer(kind=int_kind) :: &
         ndt_weld        ! number of sub-timesteps required for stability

      afsdn  (:,:) = c0
      afsd_init(:) = c0
      coag     (:) = c0
      afsd_sum     = c0
      darea        = c0
      darea_nfsd   = c0
      stability    = c0
      ndt_weld     = 1
      subdt        = dt

      do n = 1, ncat

         d_afsd_weld (:)   = c0
         d_afsdn_weld(:,n) = c0
         afsdn(:,n) = trcrn(nt_fsd:nt_fsd+nfsd-1,n)
         call icepack_cleanup_fsd (ncat, nfsd, afsdn)


         ! If there is some ice in the lower (nfsd-1) categories
         ! and there is freezing potential
         if ((frzmlt > puny) .and. &               ! freezing potential
             (aicen(n) > aminweld) .and. &         ! low concentrations area unlikely to weld
             (SUM(afsdn(1:nfsd-1,n)) > puny)) then ! some ice in nfsd-1 categories

            ! time step limitations for welding
            stability = dt * c_weld * aicen(n) * area_scaled_h(nfsd)
            ndt_weld = NINT(stability+p5) ! add 0.5 to round up number of subcycles
            subdt = dt/FLOAT(ndt_weld)    ! subcycling time step

            afsd_init(:) = afsdn(:,n)     ! save initial values
            afsd_tmp (:) = afsd_init(:)   ! work array

!            ! sanity checks
!            afsd_sum = SUM(afsd_init)
!            if (ABS(afsd_sum - c1) > puny) then
!                print*,'afsd_sum',afsd_sum
!                stop 'not 1 b4 weld'
!            endif
!            if (ANY(afsd_init < c0-puny)) stop 'negative mFSTD b4 weld'
!            if (ANY(afsd_init > c1+puny)) stop 'mFSTD>1 b4 weld'

            darea_nfsd = c0
            do nt = 1, ndt_weld
               do kx = 1, nfsd
               coag(kx) = c0
               do ky = 1, kx
                  k = iweld(kx,ky)
                  coag(kx) = coag(kx)  &
                           + area_scaled_c(ky) * afsd_tmp(ky) * aicen(n) &
                           * (SUM(afsd_tmp(k:nfsd)) &
                             + (afsd_tmp(k-1)/area_scaled_binwidth(k-1)) &
                             * (area_scaled_h(k-1) - area_scaled_h(kx) + area_scaled_c(ky)))
               end do ! ky
               end do ! kx

               afsd_tmp(1) = afsd_tmp(1) - subdt*c_weld*coag(1)
               do k = 2, nfsd
                  afsd_tmp(k) = afsd_tmp(k) - subdt*c_weld*(coag(k) - coag(k-1))
               enddo

               ! sanity checks
               if (ANY(afsd_tmp < c0-puny)) then
                        print *, 'afsd_init',afsd_init
                        print *, 'coag',coag
                        print *, 'afsd_tmp ',afsd_tmp
                        print *, 'WARNING negative mFSTD weld, l'
               end if
               if (ANY(afsd_tmp < c0-puny)) stop &
			'negative mFSTD weld, l'
               if (ANY(afsd_tmp > c1+puny)) stop &
			' mFSTD> 1 weld, l'
               if (ANY(dt*c_weld*coag < -puny)) stop &
		 'not positive'

               ! update
               darea_nfsd = darea_nfsd + subdt*c_weld*coag(nfsd)

            end do ! time

            ! ignore loss in largest cat
            afsd_tmp(nfsd) = afsd_tmp(nfsd) + darea_nfsd

            afsd_sum = c0
            do k = 1, nfsd
               afsd_sum = afsd_sum + afsd_tmp(k)
            enddo
            darea = SUM(afsd_init) - afsd_sum

            call icepack_cleanup_fsd (ncat, nfsd, afsdn)

            do k = 1, nfsd
               afsdn(k,n) = afsd_tmp(k)/afsd_sum ! in case of small numerical errors
               trcrn(nt_fsd+k-1,n) = max(afsdn(k,n), c0)
               trcrn(nt_fsd+k-1,n) = afsdn(k,n)
               ! history/diagnostics
               d_afsdn_weld(k,n) = afsdn(k,n) - afsd_init(k)
            enddo


            ! more sanity checks
            if (darea < -puny) stop 'area gain'
            if (ABS(darea) > puny) stop 'area change after correction'
            if (ANY(afsdn(:,n) < -puny)) stop 'neg, weld'
            if (afsdn(1,n) > afsd_init(1)+puny) then
                        !print *, afsdn(:,n)
			!print *, afsd_init(:)
			stop 'gain in smallest cat'
            end if
         endif ! try to weld
      enddo ! n

      ! history/diagnostics
      do k = 1, nfsd
         d_afsd_weld(k) = c0
         do n = 1, ncat
            d_afsd_weld(k) = d_afsd_weld(k) + aicen(n)*d_afsdn_weld(k,n)
         end do ! n
      end do ! k

      end subroutine fsd_weld_thermo

!=======================================================================

      end module icepack_fsd

!=======================================================================

