!  This module contains the subroutines required to fracture sea ice
!  by ocean surface waves
!
!  Theory based on:
!
!    Horvat, C., & Tziperman, E. (2015). A prognostic model of the sea-ice 
!    floe size and thickness distribution. The Cryosphere, 9(6), 2119–2134.
!    doi:10.5194/tc-9-2119-2015
!
!  and implementation described in:
!
!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018). An emergent
!    sea ice floe size distribution in a global coupled ocean--sea ice model. 
!    Journal of Geophysical Research: Oceans, 123(6), 4322–4337. 
!    doi:10.1029/2017JC013692
!
!  now with some modifications to allow direct input of ocean surface wave spectrum.
!
!  We calculate the fractures that would occur if waves enter a fully ice-covered 
!  region defined in one dimension in the direction of propagation, and then apply 
!  the outcome proportionally to the ice-covered fraction in each grid cell. Assuming
!  that sea ice flexes with the sea surface height field, strains are computed on this
!  sub-grid-scale 1D domain. If the strain between successive extrema exceeds a critical
!  value new floes are formed with diameters equal to the distance between the extrema.
!
!  authors: 2016-8 Lettie Roach, NIWA/VUW
!          
!
      module icepack_wavefracspec

      use icepack_kinds
      use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
      use icepack_parameters, only: bignum, puny, gravit, pi
      use icepack_tracers, only: nt_fsd
      use icepack_warnings, only: warnstr, icepack_warnings_add

      implicit none
      private
      public :: icepack_init_wave, icepack_step_wavefracture

      real (kind=dbl_kind), parameter  :: &
         swh_minval = 0.01_dbl_kind,  & ! minimum value of wave height (m)
         straincrit = 3.e-5_dbl_kind, & ! critical strain
         D          = 1.e4_dbl_kind,  & ! domain size
         dx         = c1,             & ! domain spacing
         threshold  = c10               ! peak-finding threshold -
                                        ! points are defined to be extrema if they
                                        ! are a local max or min over a distance 
                                        ! of 10m on both sides, based on the 
                                        ! observations of Toyota et al. (2011) who 
                                        ! find this to be the order of the smallest 
                                        ! floe size affected by wave fracture 
      integer (kind=int_kind) :: &
         nx = 10000                     ! number of points in domain

!=======================================================================

      contains

!=======================================================================
!
!  Initialize the wave spectrum and frequencies for the FSD
!
!  authors: 2018 Lettie Roach, NIWA/VUW

      subroutine icepack_init_wave(nfreq,                 &
                                   wave_spectrum_profile, &
                                   wavefreq, dwavefreq)

      integer(kind=int_kind), intent(in) :: &
         nfreq                    ! number of wave frequencies

      real(kind=dbl_kind), dimension(:), intent(out) :: &
         wave_spectrum_profile, & ! ocean surface wave spectrum as a function of frequency
	 		                ! power spectral density of surface elevation, E(f) (units m^2 s)
         wavefreq,              & ! wave frequencies (s^-1)
         dwavefreq                ! wave frequency bin widths (s^-1)

      ! local variables
      integer (kind=int_kind) :: k

      real(kind=dbl_kind), dimension(100) :: &
         wave_spectrum_data       ! default values for nfreq profile

      ! set for 25 frequencies

      wave_spectrum_data = c0

      ! FOR TESTING ONLY - do not use for actual runs!!
      wave_spectrum_data(1) = 0.00015429197810590267
      wave_spectrum_data(2) = 0.002913531381636858 
      wave_spectrum_data(3) = 0.02312942035496235
      wave_spectrum_data(4) = 0.07201970368623734
      wave_spectrum_data(5) = 0.06766948103904724 
      wave_spectrum_data(6) = 0.005527883302420378
      wave_spectrum_data(7) = 3.326293881400488e-05 
      wave_spectrum_data(8) = 6.815936703929992e-10 
      wave_spectrum_data(9) = 2.419401186610744e-20      

      do k = 1, nfreq
         wave_spectrum_profile(k) = wave_spectrum_data(k)
      enddo

      ! hardwired for wave coupling with NIWA version of Wavewatch
      ! From Wavewatch, f(n+1) = C*f(n) where C is a constant set by the user
      ! These freq are for C = 1.1
      wavefreq = (/0.04118,     0.045298,    0.0498278,   0.05481058,  0.06029164, &
                   0.06632081,  0.07295289,  0.08024818,  0.08827299,  0.09710029, &
                   0.10681032,  0.11749136,  0.1292405,   0.14216454,  0.15638101, &
                   0.17201911,  0.18922101,  0.20814312,  0.22895744,  0.25185317, &
                   0.27703848,  0.30474234,  0.33521661,  0.36873826,  0.40561208/)

      ! boundaries of bin n are at f(n)*sqrt(1/C) and f(n)*sqrt(C) 
      dwavefreq(:) = wavefreq(:)*(SQRT(1.1_dbl_kind) - SQRT(c1/1.1_dbl_kind))

      end subroutine icepack_init_wave

!=======================================================================
!
!  Calculate the change in the FSD arising from wave fracture
!
!  authors: 2017 Lettie Roach, NIWA/VUW
!
      function get_dafsd_wave(nfsd, afsd_init, fracture_hist, frac) &
                              result(d_afsd)

      integer (kind=int_kind), intent(in) :: &
         nfsd       ! number of floe size categories

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         afsd_init, fracture_hist

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         frac

      ! output
      real (kind=dbl_kind), dimension (nfsd) :: &
         d_afsd

      ! local variables
      real (kind=dbl_kind), dimension (nfsd) :: &
         loss, gain, omega

      integer (kind=int_kind) :: k

      character(len=*),parameter :: subname='(get_dafsd_wave)'

      do k = 1, nfsd
         ! fracture_hist is already normalized
         omega(k) = afsd_init(k)*SUM(fracture_hist(1:k-1)) 
      end do

      loss = omega

      do k =1,nfsd
         gain(k) = SUM(omega*frac(:,k)) 
      end do

      d_afsd(:) = gain(:) - loss(:)

      if (SUM(d_afsd(:)) > puny) then
         write(warnstr,*) subname, 'area not conserved, waves'
         call icepack_warnings_add(warnstr)
      endif

      WHERE (ABS(d_afsd).lt.puny) d_afsd = c0

      end  function get_dafsd_wave

!=======================================================================
!
!  Adaptive timestepping for wave fracture
!  See reference: Horvat & Tziperman (2017) JGR, Appendix A
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!
!
     function get_subdt_wave(nfsd, afsd_init, d_afsd) &
                              result(subdt)

      integer (kind=int_kind), intent(in) :: &
         nfsd       ! number of floe size categories

      real (kind=dbl_kind), dimension (nfsd), intent(in) :: &
         afsd_init, d_afsd ! floe size distribution tracer 

      ! output
      real (kind=dbl_kind) :: &
         subdt ! subcycle timestep (s)

      ! local variables
      real (kind=dbl_kind), dimension (nfsd) :: &
         check_dt ! to compute subcycle timestep (s)

      integer (kind=int_kind) :: k

      check_dt(:) = bignum 
      do k = 1, nfsd
          if (d_afsd(k) >  puny) check_dt(k) = (1-afsd_init(k))/d_afsd(k)
          if (d_afsd(k) < -puny) check_dt(k) = afsd_init(k)/ABS(d_afsd(k))
      end do 

      subdt = MINVAL(check_dt)

      end function get_subdt_wave

!=======================================================================
! 
!  Given fracture histogram computed from local wave spectrum, evolve 
!  the floe size distribution
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!
      subroutine icepack_step_wavefracture(wave_spec_type,   &
                  dt,            ncat,            nfsd,      &
                  nfreq,                                     &
                  aice,          vice,            aicen,     &
                  floe_rad_l,    floe_rad_c,                 &
                  wave_spectrum, wavefreq,        dwavefreq, &
                  trcrn,         d_afsd_wave)

      use icepack_fsd, only: icepack_cleanup_fsd

      character (len=char_len), intent(in) :: &
         wave_spec_type   ! type of wave spectrum forcing

      integer (kind=int_kind), intent(in) :: &
         nfreq,        & ! number of wave frequency categories
         ncat,         & ! number of thickness categories
         nfsd            ! number of floe size categories

      real (kind=dbl_kind), intent(in) :: &
         dt,           & ! time step
         aice,         & ! ice area fraction
         vice            ! ice volume per unit area

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen           ! ice area fraction (categories)

      real(kind=dbl_kind), dimension(:), intent(in) ::  &
         floe_rad_l,   & ! fsd size lower bound in m (radius)
         floe_rad_c      ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         wavefreq,     & ! wave frequencies (s^-1)
         dwavefreq       ! wave frequency bin widths (s^-1)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         wave_spectrum   ! ocean surface wave spectrum as a function of frequency
	 		 ! power spectral density of surface elevation, E(f) (units m^2 s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn           ! tracer array

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         d_afsd_wave     ! change in fsd due to waves

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         d_afsdn_wave    ! change in fsd due to waves, per category

      ! local variables
      integer (kind=int_kind) :: &  
         n, k, t, &
         nsubt ! number of subcycles 

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         afsdn           ! floe size and thickness distribution

      real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
         frac    

      real (kind=dbl_kind) :: &
         hbar         , & ! mean ice thickness
         elapsed_t    , & ! elapsed subcycling time
         subdt        , & ! subcycling time step
         cons_error       ! area conservation error

      real (kind=dbl_kind), dimension (nfsd) :: &
         fracture_hist, & ! fracture histogram
         afsd_init    , & ! tracer array
         afsd_tmp     , & ! tracer array
         d_afsd_tmp       ! change

      character(len=*),parameter :: subname='(icepack_step_wavefracture)'

      !------------------------------------

      ! initialize 
      d_afsd_wave    (:)   = c0
      d_afsdn_wave   (:,:) = c0
      fracture_hist  (:)   = c0

      ! if all ice is not in first floe size category
      if (.NOT. ALL(trcrn(nt_fsd,:).ge.c1-puny)) then
 
      ! do not try to fracture for minimal ice concentration or zero wave spectrum
      if ((aice > p01).and.(MAXVAL(wave_spectrum(:)) > puny)) then

         hbar = vice / aice

         ! calculate fracture histogram
         call wave_frac(nfsd, nfreq, wave_spec_type, &
                        floe_rad_l, floe_rad_c, &
                        wavefreq, dwavefreq, &
                        hbar, wave_spectrum, fracture_hist)

         ! if fracture occurs
         if (MAXVAL(fracture_hist) > puny) then
            ! protect against small numerical errors
            call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
            
            do n = 1, ncat
              
              afsd_init(:) = trcrn(nt_fsd:nt_fsd+nfsd-1,n)

              ! if there is ice, and a FSD, and not all ice is the smallest floe size 
              if ((aicen(n) > puny) .and. (SUM(afsd_init(:)) > puny) &
                                    .and.     (afsd_init(1) < c1)) then

                 afsd_tmp =  afsd_init

                  ! frac does not vary within subcycle
                  frac(:,:) = c0
                  do k = 2, nfsd
                     frac(k,1:k-1) = fracture_hist(1:k-1)
                  end do
                  do k = 1, nfsd
                     if (SUM(frac(k,:)) > c0) frac(k,:) = frac(k,:)/SUM(frac(k,:))
                  end do

                  ! adaptive sub-timestep
                  elapsed_t = c0
                  cons_error = c0
                  nsubt = 0
                  DO WHILE (elapsed_t < dt)
                     nsubt = nsubt + 1

                     ! if all floes in smallest category already, exit
                     if (afsd_tmp(1).ge.c1-puny) EXIT 

                     ! calculate d_afsd using current afstd
                     d_afsd_tmp = get_dafsd_wave(nfsd, afsd_tmp, fracture_hist, frac)
                     
                     ! check in case wave fracture struggles to converge
                     if (nsubt>100) then
                          print *, 'afsd_tmp ',afsd_tmp
                          print *, 'dafsd_tmp ',d_afsd_tmp
                          print *, 'subt ',nsubt
                          print *, &
                              'wave frac taking a while to converge....'
                     end if
 
                     ! required timestep
                     subdt = get_subdt_wave(nfsd, afsd_tmp, d_afsd_tmp)
                     subdt = MIN(subdt, dt)

                     ! update afsd
                     afsd_tmp = afsd_tmp + subdt * d_afsd_tmp(:) 

                     ! check conservation and negatives
                     if (MINVAL(afsd_tmp) < -puny) then
                        write(warnstr,*) subname, 'wb, <0 loop'
                        call icepack_warnings_add(warnstr)
                     endif
                     if (MAXVAL(afsd_tmp) > c1+puny) then
                         write(warnstr,*) subname, 'wb, >1 loop'
                        call icepack_warnings_add(warnstr)
                     endif

                     ! update time
                     elapsed_t = elapsed_t + subdt 

                  END DO ! elapsed_t < dt
 
                  ! In some cases---particularly for strong fracturing---the equation
                  ! for wave fracture does not quite conserve area.
                  ! With the dummy wave forcing, this happens < 2% of the time (in
                  ! 1997) and is always less than 10^-7.
                  ! Simply renormalizing may cause the first floe size 
                  ! category to reduce, which is not physically allowed
                  ! to happen. So we adjust here
                  cons_error = SUM(afsd_tmp) - c1

                  ! area loss: add to first category
                  if (cons_error.lt.c0) then
                      afsd_tmp(1) = afsd_tmp(1) - cons_error
                  else
                  ! area gain: take it from the largest possible category 
                  do k = nfsd, 1, -1
                     if (afsd_tmp(k).gt.cons_error) then
                        afsd_tmp(k) = afsd_tmp(k) - cons_error
                        EXIT
                     end if
                  end do
                  end if

                  ! update trcrn
                  trcrn(nt_fsd:nt_fsd+nfsd-1,n) = afsd_tmp/SUM(afsd_tmp)
                  call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
 

                  ! for diagnostics
                  d_afsdn_wave(:,n) = afsd_tmp(:) - afsd_init(:)  
                  d_afsd_wave (:)   = d_afsd_wave(:) + aicen(n)*d_afsdn_wave(:,n)

               endif ! aicen > puny
            enddo    ! n
         endif       ! fracture occurs

      endif          ! aice > p01
      end if         ! all small floes

     end subroutine icepack_step_wavefracture

!=======================================================================
!
!  Calculates functions to describe the change in the FSD when waves 
!  fracture ice, given a wave spectrum (1D frequency, nfreq (default 25)
!  frequency bins)
!
!  We calculate extrema and if these are successive maximum, 
!  minimum, maximum or vice versa, and have strain greater than a 
!  critical strain, break ice and create new floes with lengths equal
!  to these distances. Based on MatLab code written by Chris Horvat,
!  from Horvat & Tziperman (2015). 
!
!  Note that a realization of sea surface height requires a random phase.
!
!  authors: 2018 Lettie Roach, NIWA/VUW

      subroutine wave_frac(nfsd, nfreq, wave_spec_type, &
                           floe_rad_l, floe_rad_c, &
                           wavefreq, dwavefreq, &
                           hbar, spec_efreq, frac_local)

      integer (kind=int_kind), intent(in) :: &
         nfsd, &       ! number of floe size categories
         nfreq         ! number of wave frequency categories

      character (len=char_len), intent(in) :: &
        wave_spec_type ! type of wave spectrum forcing

      real (kind=dbl_kind),  intent(in) :: &
         hbar          ! mean ice thickness (m)

      real(kind=dbl_kind), dimension(:), intent(in) ::  &
         floe_rad_l, & ! fsd size lower bound in m (radius)
         floe_rad_c    ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         wavefreq,   & ! wave frequencies (s^-1) 
         dwavefreq,  & ! wave frequency bin widths (s^-1)
         spec_efreq    ! wave spectrum (m^2 s)

      real (kind=dbl_kind), dimension (nfsd), intent(out) :: &
         frac_local    ! fracturing histogram

      ! local variables

      integer (kind=int_kind) :: i, j, k

      integer, parameter :: &
         loopcts = 1  ! number of SSH realizations

      real (kind=dbl_kind), dimension(nfreq) :: &
         spec_elambda,             & ! spectrum as a function of wavelength (m^-1 s^-1)     
         reverse_spec_elambda,     & ! reversed
         reverse_lambda, lambda,   & ! wavelengths (m)
         reverse_dlambda, dlambda, & ! wavelength bin spacing (m)
         spec_coeff,               &
         phi, rand_array, summand

      real (kind=dbl_kind), dimension(2*nx) :: &
         fraclengths

      real (kind=dbl_kind), dimension(nx) :: &
         X,  &    ! spatial domain (m)
         eta      ! sea surface height field (m)

      real (kind=dbl_kind), dimension(nfsd) :: &
         frachistogram 

      logical (kind=log_kind) :: &
         e_stop        ! if true, stop and return zero omega and fsdformed

      e_stop = .false. ! if true, aborts fracture calc
  
      ! spatial domain
      do j = 1, nx
         X(j)= j*dx
      end do

      ! dispersion relation
      reverse_lambda (:) = gravit/(c2*pi*wavefreq (:)**2)
      reverse_dlambda(:) = gravit/(c2*pi*dwavefreq(:)**2)
      ! convert to lambda spectrum
      reverse_spec_elambda(:) = spec_efreq(:) &
                   *(p5 * (gravit/(c2*pi*reverse_lambda(:)**3) )**p5)
      ! reverse lambda
      lambda (:) = reverse_lambda (nfreq:1:-1)
      dlambda(:) = reverse_dlambda(nfreq:1:-1)
      spec_elambda(:) = reverse_spec_elambda(nfreq:1:-1) 
 
      ! spectral coefficients
      spec_coeff = sqrt(c2*spec_elambda*dlambda) 

      ! initialize fracture lengths
      fraclengths(:) = c0
     
      ! loop over n. realizations of SSH
      do i = 1, loopcts

         ! Phase for each Fourier component may be constant or
         ! a random phase that varies in each i loop
         ! See documentation for discussion
         if (trim(wave_spec_type)=='random') then
            call RANDOM_NUMBER(rand_array)
         else
            rand_array(:) = p5
         endif
         phi = c2*pi*rand_array
 
         do j = 1, nx
            ! SSH field in space (sum over wavelengths, no attenuation)
            summand = spec_coeff*COS(2*pi*X(j)/lambda+phi)
            eta(j)  = SUM(summand)
         end do
 
         if ((SUM(ABS(eta)) > puny).and.(hbar > puny)) then 
            call get_fraclengths(X, eta, fraclengths, hbar, e_stop)
         end if
      end do
 
      frachistogram(:) = c0

      if (.not. e_stop) then

         ! convert from diameter to radii
         fraclengths(:) = fraclengths(:)/c2

         ! bin into FS cats
         ! highest cat cannot be fractured into
         do j = 1, size(fraclengths)
            do k = 1, nfsd-1
               if ((fraclengths(j) >= floe_rad_l(k)) .and. &
                   (fraclengths(j) < floe_rad_l(k+1))) then
                  frachistogram(k) = frachistogram(k) + 1
               end if
            end do

         end do
      end if

      do k = 1, nfsd
         frac_local(k) = floe_rad_c(k)*frachistogram(k)
      end do

      ! normalize
      if (SUM(frac_local) /= c0) frac_local(:) = frac_local(:) / SUM(frac_local(:))

      end subroutine wave_frac

!===========================================================================
!
!  Given the (attenuated) sea surface height, find the strain across triplets
!  of max, min, max or min, max, min (local extrema within 10m).
!  If this strain is greater than the  critical strain, ice can fracture
!  and new floes are formed with sizes equal to the distances between
!  extrema. Based on MatLab code written by Chris Horvat,
!  from Horvat & Tziperman (2015). 
!
!  authors: 2016 Lettie Roach, NIWA/VUW
!
      subroutine get_fraclengths(X, eta, fraclengths, hbar, e_stop)

      real (kind=dbl_kind) :: &
         hbar             ! mean thickness (m)

      real (kind=dbl_kind), intent(in), dimension (nx) :: &
         X, &              ! spatial domain (m)
         eta               ! sea surface height field (m)

      real (kind=dbl_kind), intent(inout), dimension (2*nx) :: &
         fraclengths      ! the biggest number of fraclengths we could have is
                          ! two floe pieces created at each subgridpoint ie. 2*nx
                          ! This will never actually happen - most of the array
                          ! will be zeros

      logical (kind=log_kind), intent(inout) :: &
         e_stop           ! if true, stop and return zero omega and fsdformed

      ! local variables
      integer (kind=int_kind) :: &
         spcing,        & ! distance over which to search for extrema on each side of point
         j, k,          & ! indices to iterate over domain
         first, last,   & ! indices over which to search for extrema
         j_neg,         & ! nearest extrema backwards
         j_pos,         & ! nearest extrema forwards
         n_above          ! number of points where strain is above critical

      real (kind=dbl_kind), dimension(nx) :: &
         strain,        & ! the strain between triplets of extrema
         frac_size_one, & !
         frac_size_two

      logical (kind=log_kind), dimension(nx) :: &
         is_max, is_min,& ! arrays to hold whether each point is a local max or min
         is_extremum,   & ! or extremum
         is_triplet       ! or triplet of extrema

      real (kind=dbl_kind) :: &
         delta,         & ! difference in x between current and prev extrema
         delta_pos        ! difference in x between next and current extrema

      integer (kind=int_kind), dimension(1) :: &
         maxj, minj       ! indices of local max and min

      ! ------- equivalent of peakfinder2
      ! given eta and spcing, compute extremelocs in ascending order
      spcing = nint(threshold/dx)

      is_max = .false.
      is_min = .false.
      is_extremum = .false.
      is_triplet = .false.
      strain = c0
      frac_size_one = c0
      frac_size_two = c0
      j_neg = 0
      j_pos = 0      
      fraclengths(:) = c0

      ! search for local max and min within spacing
      ! on either side of each point

      do j = 1, nx

         first = MAX(1,j-spcing)
         last  = MIN(nx,j+spcing)

         maxj = MAXLOC(eta(first:last))
         minj = MINLOC(eta(first:last))

!         if (COUNT(eta(first:last) == MAXVAL(eta(first:last))) > 1) &
!                       stop 'more than one max'
!         if (COUNT(eta(first:last) == MINVAL(eta(first:last))) > 1) &
!                        stop 'more than one min'

         if (maxj(1)+first-1 == j) is_max(j) = .true.
         if (minj(1)+first-1 == j) is_min(j) = .true.

!         if (is_min(j).and.is_max(j)) then
!            print *, 'X ',X
!            print *, 'eta ',eta
!            print *, 'frst last' ,first, last
!            print *, 'maxj, minj ',maxj,minj
!            stop     'error in extrema'
!         end if
         if (is_min(j).or.is_max(j)) is_extremum(j) = .true.
      end do

      do j = 2, nx-1
         if (is_extremum(j)) then
            if (j == 2) then
               if (is_extremum(1)) j_neg = 1
            else
               do k = j-1, 1, -1
                  if (is_extremum(k)) then
                     j_neg = k
                     EXIT
                  end if
               end do
            end if
                       
            do k = j+1, nx
               if (is_extremum(k)) then
                  j_pos = k
                  EXIT
               end if
            end do
                       
            if ((j_neg > 0).and.(j_pos > 0)) then 
               if (is_max(j_neg).and.is_min(j).and.is_max(j_pos)) &
                  is_triplet(j) = .true.
               if (is_min(j_neg).and.is_max(j).and.is_min(j_pos)) &
                  is_triplet(j) = .true.
            end if

            ! calculate strain
            if (is_triplet(j)) then
               delta_pos = X(j_pos) - X(j    )
               delta     = X(j    ) - X(j_neg)

               strain(j) = p5*hbar*(eta(j_neg) - eta(j)) &
                                  / (delta*(delta+delta_pos))

               if (strain(j) > straincrit) then
                  frac_size_one(j) = X(j_pos) - X(j    )
                  frac_size_two(j) = X(j    ) - X(j_neg)
               end if
            end if
         end if

      end do

      n_above = COUNT(strain > straincrit)
      if (n_above > 0) then
         fraclengths(1:nx)      = frac_size_one(:)
         fraclengths(nx+1:2*nx) = frac_size_two(:)
         e_stop = .false.
      else
         e_stop = .true.
      end if

      end subroutine get_fraclengths

!=======================================================================
     
      end module icepack_wavefracspec

!=======================================================================


