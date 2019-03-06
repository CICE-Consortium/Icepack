!=========================================================================
!
! This module contains the subroutines required to define
! a floe size distribution tracer for sea ice
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
      use icepack_parameters, only: c0, c1, c2, c4, puny
      use icepack_parameters, only: pi, floeshape
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: icepack_init_fsd_bounds, icepack_init_fsd
!      public :: init_fsd, init_fsd_bounds,     &
!          frzmlt_bottom_lateral_fsd, &
!          renorm_mfstd, wave_dep_growth 

!      logical (kind=log_kind), public :: & 
!         write_diag_wave     ! if .true., write the lat/lons from find_wave to history file

!      real (kind=dbl_kind), public :: &
!        c_mrg            ! constant of proportionality for merging
                         ! see documentation for details

!      logical (kind=log_kind), public :: &
!         rdc_frzmlt      ! if true, only (1-oo) of frzmlt can melt (lat and bot)

!      integer(kind=int_kind), save, public ::  &
!         nfreq           ! number of frequencies in wave spectrum   

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
         alpha_mrg

!=======================================================================

      contains

!=======================================================================

!  Initialize ice fsd bounds (call whether or not restarting)
!  Define the bounds, midpoints and widths of floe size
!  categories in area and radius
!
!  Note that radius widths cannot be larger than twice previous
!  category width or floe merging will not have an effect
!
!  Note also that the bound of the lowest floe size category is used
!  to define the lead region width and the domain spacing for wave breaking
!
      subroutine icepack_init_fsd_bounds(ncat, nfsd, &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth, &  ! fsd size bin width in m (radius)
         c_fsd_range)      ! string for history output

      integer (kind=int_kind), intent(in) :: &
         ncat, & ! number of thickness categories
         nfsd    ! number of floe size categories

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

         call icepack_warnings_add(subname//' floe size categories not defined for nfsd')
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
         alpha_mrg     (nfsd, nfsd), & !
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

      ! scaled area for merging
      floe_area_binwidth = floe_area_h - floe_area_l
      area_lims(1:nfsd) = floe_area_l(1:nfsd)
      area_lims(nfsd+1) = floe_area_h(nfsd)
      area_lims_scaled = (area_lims - area_lims(1))/MAXVAL(floe_area_binwidth)

      area_scaled_h = area_lims_scaled(2:nfsd+1)
      area_scaled_l = area_lims_scaled(1:nfsd  )
      area_scaled_c = (area_scaled_h + area_scaled_l) / c2
      area_scaled_binwidth = area_scaled_h - area_scaled_l

      ! which floe sizes can combine during merging
      alpha_mrg(:,:) = -999
      do n  = 1, nfsd
         do m = 1, nfsd
            test = area_scaled_h(n) - area_scaled_c(m)
            do k = 1, nfsd
               if ((test >= area_scaled_l(k)) .and. (test < area_scaled_h(k))) then
                  alpha_mrg(n,m) = k + 1
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

            ! Write hin_max to character string
            write (c_fsd1, '(f6.3)') floe_rad(n-1)
            write (c_fsd2, '(f6.3)') floe_rad(n)

            ! Save character string to write to history file
            c_fsd_range(n)=c_fsd1//'m < fsd Cat '//c_nf//' < '//c_fsd2//'m'
         enddo

         write(warnstr,*) ' '
         call icepack_warnings_add(warnstr)

         ! write(*,*) &
         !    'floe size bin info: low, high, center, width, area_c'
         ! write(*,*) floe_rad_l
         ! write(*,*) floe_rad_h
         ! write(*,*) floe_rad_c
         ! write(*,*) floe_binwidth
         ! write(*,*) floe_area_l
         ! write(*,*) floe_area_h   
         ! write(*,*) floe_area_c

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

!         write(*,*)'init_fsd: initial number distribution of floes'
!         write(*,*) num_fsd(1:nfsd)
!         write(*,*)'init_fsd: initial fraction distribution of floes'
!         write(*,*) afsd(1:nfsd)
!         write(*,*) SUM(afsd(1:nfsd)) ! should be 1

      endif ! ice_ic

      end subroutine icepack_init_fsd

!=======================================================================

      end module icepack_fsd

!=======================================================================

