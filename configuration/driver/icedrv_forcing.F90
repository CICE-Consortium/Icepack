!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
      module icedrv_forcing

      use icedrv_kinds
      use icedrv_domain_size, only: nx
      use icedrv_calendar, only: time, nyr, dayyr, mday, month, secday
      use icedrv_calendar, only: daymo, daycal, dt, yday, sec
      use icedrv_calendar, only: npt, use_leap_years, time0, year_init
      use icedrv_constants, only: nu_diag, nu_forcing, nu_open_clos
      use icedrv_constants, only: c0, c1, c2, c10, c100, p5, c4, c24
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_sea_freezing_temperature
      use icepack_intfc, only: icepack_init_wave
      use icedrv_system, only: icedrv_system_abort
      use icedrv_flux, only: zlvl, Tair, potT, rhoa, uatm, vatm, wind, &
         strax, stray, fsw, swvdr, swvdf, swidr, swidf, Qa, flw, frain, &
         fsnow, sst, sss, uocn, vocn, qdp, hmix, Tf, opening, closing, sstdat
#ifdef USE_NETCDF
      use netcdf
#endif

      implicit none
      private

      public :: init_forcing, get_forcing, interp_coeff, &
                interp_coeff_monthly, get_wave_spec

      integer (kind=int_kind) :: &
         ntime               ! number of data points in time

      integer (kind=int_kind), public :: &
         ycycle          , & ! number of years in forcing cycle
         fyear_init      , & ! first year of data in forcing cycle
         fyear           , & ! current year in forcing cycle
         fyear_final         ! last year in cycle

      real (kind=dbl_kind), allocatable :: &
            fsw_data(:), & ! field values at temporal data points
           cldf_data(:), &
          fsnow_data(:), &
           Tair_data(:), &
           uatm_data(:), &
           vatm_data(:), &
           wind_data(:), &
          strax_data(:), &
          stray_data(:), &
           rhum_data(:), &
             Qa_data(:), &
           rhoa_data(:), &
           potT_data(:), &
            flw_data(:), &
            qdp_data(:), &
            sst_data(:), &
            sss_data(:), &
           uocn_data(:), &
           vocn_data(:), &
          frain_data(:), &
          swvdr_data(:), &
          swvdf_data(:), &
          swidr_data(:), &
          swidf_data(:), &
           zlvl_data(:), &
           hmix_data(:), &
           open_data(:), &
           clos_data(:)

      real (kind=dbl_kind), dimension(nx) :: &
          sst_temp

      character(char_len), public :: &
         atm_data_format, & ! 'bin'=binary or 'nc'=netcdf
         ocn_data_format, & ! 'bin'=binary or 'nc'=netcdf
         bgc_data_format, & ! 'bin'=binary or 'nc'=netcdf
         atm_data_type,   & ! 'default', 'clim', 'CFS', 'MDF'
         ocn_data_type,   & ! 'default', 'SHEBA' 'MDF'
         bgc_data_type,   & ! 'default', 'ISPOL', 'NICE'
         lateral_flux_type,   & ! 'uniform_ice', 'open_water'
         atm_data_file,   & ! atmospheric forcing data file
         ocn_data_file,   & ! ocean forcing data file
         ice_data_file,   & ! ice forcing data file
         bgc_data_file,   & ! biogeochemistry forcing data file
         precip_units       ! 'mm_per_month', 'mm_per_sec', 'mks'

      character(char_len_long), public :: &
         data_dir           ! top directory for forcing data

      real (kind=dbl_kind), parameter, public :: &
         frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
         frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
         frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
         frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band

      logical (kind=log_kind), public :: &
         oceanmixed_ice        , & ! if true, use internal ocean mixed layer
         restore_ocn           , & ! restore sst if true
         precalc_forc              ! whether to precalculate forcing

      real (kind=dbl_kind), public :: &
         trest, &           ! restoring time scale (sec)
         trestore           ! restoring time scale (days)

      character (len=char_len_long), public :: &
         snw_ssp_table      ! snow table type 'test', 'snicar'

!=======================================================================

      contains

!=======================================================================

      subroutine init_forcing

! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the forcing data.

      integer (kind=int_kind) :: &
         i                ! index

      character(len=*), parameter :: subname='(init_forcing)'

      ! Initialize ntime and allocate data arrays
      if (precalc_forc) then
         if (trim(atm_data_type(1:3)) /= 'MDF') &
            call icedrv_system_abort(string=subname//&
            'precalc_forc should only be used with MDF atmosphere', &
            file=__FILE__,line=__LINE__)
         if (.not. ((trim(ocn_data_type(1:3)) == 'MDF') &
            .or. (trim(ocn_data_type(1:7)) == 'default'))) &
            call icedrv_system_abort(string=subname//&
            'precalc_forc should only be used with MDF ocean or'//&
            ' default ocean', file=__FILE__,line=__LINE__)
         ntime = npt
      else
         ntime = 8760
      endif
      allocate(fsw_data(ntime), cldf_data(ntime), fsnow_data(ntime), &
               Tair_data(ntime), uatm_data(ntime), vatm_data(ntime), &
               wind_data(ntime), strax_data(ntime), stray_data(ntime), &
               rhum_data(ntime), Qa_data(ntime), rhoa_data(ntime), &
               potT_data(ntime), flw_data(ntime), qdp_data(ntime), &
               sst_data(ntime), sss_data(ntime), uocn_data(ntime), &
               vocn_data(ntime), frain_data(ntime), swvdr_data(ntime), &
               swvdf_data(ntime), swidr_data(ntime), swidf_data(ntime), &
               zlvl_data(ntime), hmix_data(ntime), open_data(ntime), &
               clos_data(ntime))

      fyear       = fyear_init + mod(nyr-1,ycycle) ! current year
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      write (nu_diag,*) ' Initial forcing data year = ',fyear_init
      write (nu_diag,*) ' Final   forcing data year = ',fyear_final

     !-------------------------------------------------------------------
     ! Initialize forcing data to default values
     !-------------------------------------------------------------------

      ! many default forcing values are set in init_flux_atm
      i = 1 ! use first grid box value

          zlvl_data(:) = zlvl (i)    ! atmospheric data level (m)
          Tair_data(:) = Tair (i)    ! air temperature  (K)
          potT_data(:) = potT (i)    ! air potential temperature  (K)
          rhoa_data(:) = rhoa (i)    ! air density (kg/m^3)
          uatm_data(:) = uatm (i)    ! wind velocity components (m/s)
          vatm_data(:) = vatm (i)
          wind_data(:) = wind (i)    ! wind speed (m/s)
         strax_data(:) = strax(i)    ! wind stress components (N/m^2)
         stray_data(:) = stray(i)
           fsw_data(:) = fsw  (i)    ! incoming shortwave radiation (W/m^2)
         swvdr_data(:) = swvdr(i)    ! sw down, visible, direct  (W/m^2)
         swvdf_data(:) = swvdf(i)    ! sw down, visible, diffuse (W/m^2)
         swidr_data(:) = swidr(i)    ! sw down, near IR, direct  (W/m^2)
         swidf_data(:) = swidf(i)    ! sw down, near IR, diffuse (W/m^2)
            Qa_data(:) = Qa   (i)    ! specific humidity (kg/kg)
           flw_data(:) = flw  (i)    ! incoming longwave radiation (W/m^2)
         frain_data(:) = frain(i)    ! rainfall rate (kg/m^2 s)
         fsnow_data(:) = fsnow(i)    ! snowfall rate (kg/m^2 s)
           qdp_data(:) = qdp  (i)    ! deep ocean heat flux (W/m^2)
           sss_data(:) = sss  (i)    ! sea surface salinity
           hmix_data(:)= hmix (i)    ! ocean mixed layer depth (m)
          uocn_data(:) = uocn (i)    ! ocean current components (m/s)
          vocn_data(:) = vocn (i)
          cldf_data(:) = c0          ! cloud fraction

      if (trim(atm_data_type(1:4)) == 'CFS')   call atm_CFS
      if (trim(atm_data_type(1:4)) == 'clim')  call atm_climatological
      if (trim(atm_data_type(1:5)) == 'ISPOL') call atm_ISPOL
      if (trim(atm_data_type(1:4)) == 'NICE')  call atm_NICE
      if (trim(atm_data_type(1:3)) == 'MDF')   call atm_MDF
      if (trim(ocn_data_type(1:5)) == 'SHEBA') call ice_open_clos

      if (restore_ocn) then
        if (trestore == 0) then
          trest = dt        ! use data instantaneously
        else
          trest = real(trestore,kind=dbl_kind) * secday ! seconds
        end if
        sst_data(:) = sstdat(i)    ! default may be overwritten below
      else
        sst_data(:) = sst   (i)    ! default or restart value if not restoring
      endif

      if (trim(ocn_data_type(1:5)) == 'ISPOL') call ocn_ISPOL
      if (trim(ocn_data_type(1:4)) == 'NICE')  call ocn_NICE
      if (trim(ocn_data_type(1:3)) == 'MDF') call ocn_MDF

      call prepare_forcing (Tair_data,     fsw_data,      &
                            cldf_data,     &
                            frain_data,    fsnow_data,    &
                            Qa_data,       rhoa_data,     &
                            uatm_data,     vatm_data,     &
                            strax_data,    stray_data,    &
                            zlvl_data,     wind_data,     &
                            swvdr_data,    swvdf_data,    &
                            swidr_data,    swidf_data,    &
                            potT_data)

      end subroutine init_forcing

!=======================================================================

      subroutine get_forcing(timestep)

!ECH notes
! We will probably need to send in the time and work out what the data
! time slice is, instead of sending in the timestep.  This currently assumes
! the time step and the data both start Jan 1.

      integer (kind=int_kind), intent(in) :: &
         timestep         ! time step index

      integer (kind=int_kind) :: &
         i            , & ! data index
         recslot      , & ! spline slot for current record
         midmonth         ! middle day of month

      integer (kind=int_kind) :: &
         mlast, mnext     ! indices of bracketing time slices

      real (kind=dbl_kind) :: &
         c1intp, c2intp   ! interpolation coefficients

      integer (kind=int_kind) :: &
          recnum      , & ! record number for current data value
          maxrec      , & ! maximum number of data records
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          offndy          ! Julian day of first data record

      real (kind=dbl_kind) :: &
          sec6hr      , & ! number of seconds in 6 hours
          sec1hr      , & ! number of seconds in 1 hour
          offset          ! time to first data record since 1 Jan (s)

      character(len=*), parameter :: subname='(get_forcing)'

      if (precalc_forc) then
         ! Fill all grid boxes with same forcing data
         Tair (:) = Tair_data(timestep)
         Qa   (:) = Qa_data(timestep)
         uatm (:) = uatm_data(timestep)
         vatm (:) = vatm_data(timestep)
         fsnow(:) = fsnow_data(timestep)
         flw  (:) = flw_data(timestep)
         fsw  (:) = fsw_data(timestep)

         ! derived (or not otherwise set)
         potT (:) = potT_data(timestep)
         wind (:) = wind_data(timestep)
         strax(:) = strax_data(timestep)
         stray(:) = stray_data(timestep)
         rhoa (:) = rhoa_data(timestep)
         frain(:) = frain_data(timestep)
         swvdr(:) = swvdr_data(timestep)
         swvdf(:) = swvdf_data(timestep)
         swidr(:) = swidr_data(timestep)
         swidf(:) = swidf_data(timestep)

         ! Ocean forcing
         sst_temp(:) = sst_data(timestep)
         sss     (:) = sss_data(timestep)
         uocn    (:) = uocn_data(timestep)
         vocn    (:) = vocn_data(timestep)
         qdp     (:) = qdp_data(timestep)

      else
         if (trim(atm_data_type) == 'CFS') then
            ! calculate data index corresponding to current timestep
            i = mod(timestep-1,ntime)+1 ! repeat forcing cycle
            mlast = i
            mnext = mlast
            c1intp = c1
            c2intp = c0

            ! fill all grid boxes with the same forcing data
            Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
            Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
            uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
            vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
            fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)
            flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)
            fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)

            ! derived (or not otherwise set)
            potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
            wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
            strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
            stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
            rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
            frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)
            swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
            swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
            swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
            swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

         elseif (trim(atm_data_type) == 'clim') then
            midmonth = 15  ! assume data is given on 15th of every month
            recslot = 1                             ! latter half of month
            if (mday < midmonth) recslot = 2        ! first half of month
            if (recslot == 1) then
               mlast = month
               mnext = mod(month   ,12) + 1
            else ! recslot = 2
               mlast = mod(month+10,12) + 1
               mnext = month
            endif
            call interp_coeff_monthly(recslot, c1intp, c2intp)

            ! fill all grid boxes with the same forcing data
            Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
            Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
            uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
            vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
            fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)
            flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)
            fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)

            ! derived (or not otherwise set)
            potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
            wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
            strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
            stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
            rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
            frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)
            swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
            swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
            swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
            swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

         elseif (trim(atm_data_type) == 'ISPOL') then

         offndy = 0                              ! first data record (Julian day)
         offset = real(offndy,dbl_kind)*secday
         dataloc = 1                             ! data located at middle of interval
         maxrec = 365
         recslot = 2
         recnum = mod(int(yday)+maxrec-offndy-1,maxrec)+1
         mlast = mod(recnum+maxrec-2,maxrec) + 1
         mnext = mod(recnum-1,       maxrec) + 1
         call interp_coeff (recnum, recslot, secday, dataloc, &
                              c1intp, c2intp, offset)

         Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
         Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
         uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
         vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
         fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)

            ! derived (or not otherwise set)
            potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
            wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
            strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
            stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
            rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
            frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)

         sec6hr = secday/c4;                      ! seconds in 6 hours
         offndy = 0
         maxrec = 1460
         recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)
         recnum = mod(recnum+maxrec-4*offndy-1,maxrec)+1 ! data begins on 16 June 2004
         recslot = 2
         mlast = mod(recnum+maxrec-2,maxrec) + 1
         mnext = mod(recnum-1,       maxrec) + 1
         call interp_coeff (recnum, recslot, sec6hr, dataloc, &
                              c1intp, c2intp, offset)

         fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)
         flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)

            ! derived
            swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
            swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
            swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
            swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

         elseif (trim(atm_data_type) == 'NICE') then

         offndy = 0                              ! first data record (Julian day)
         offset = real(offndy,dbl_kind)*secday
         dataloc = 1                          ! data located in middle of interval
         maxrec = 365
         recslot = 2
         recnum = mod(int(yday)+maxrec-offndy-1,maxrec)+1
         mlast = mod(recnum+maxrec-2,maxrec) + 1
         mnext = mod(recnum-1,       maxrec) + 1
         call interp_coeff (recnum, recslot, secday, dataloc, &
                              c1intp, c2intp, offset)

         Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
         Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
         uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
         vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
         fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)

            ! derived (or not otherwise set)
            potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
            wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
            strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
            stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
            rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
            frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)

         sec6hr = secday/c4;                      ! seconds in 6 hours
         maxrec = 1460
         dataloc = 2                              ! data located at end of interval
         recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)
         recnum = mod(recnum+maxrec-4*offndy-1,maxrec)+1
         recslot = 2
         mlast = mod(recnum+maxrec-2,maxrec) + 1
         mnext = mod(recnum-1,       maxrec) + 1
         call interp_coeff (recnum, recslot, sec6hr, dataloc, &
                              c1intp, c2intp, offset)

         fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)
         flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)

            ! derived
            swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
            swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
            swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
            swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)

         endif

! possible bug:  is the ocean data also offset to the beginning of the field campaigns?

      if (trim(ocn_data_type) == 'ISPOL') then

         midmonth = 15  ! assume data is given on 15th of every month
         recslot = 1                             ! latter half of month
         if (mday < midmonth) recslot = 2        ! first half of month
         if (recslot == 1) then
            mlast = month
            mnext = mod(month   ,12) + 1
         else ! recslot = 2
            mlast = mod(month+10,12) + 1
            mnext = month
         endif
         call interp_coeff_monthly(recslot, c1intp, c2intp)

         sst_temp(:) = c1intp *  sst_data(mlast) + c2intp *  sst_data(mnext)
         sss     (:) = c1intp *  sss_data(mlast) + c2intp *  sss_data(mnext)
         uocn    (:) = c1intp * uocn_data(mlast) + c2intp * uocn_data(mnext)
         vocn    (:) = c1intp * vocn_data(mlast) + c2intp * vocn_data(mnext)
         qdp     (:) = c1intp *  qdp_data(mlast) + c2intp *  qdp_data(mnext)

      elseif (trim(ocn_data_type) == 'NICE') then

         dataloc = 2                          ! data located at end of interval
         maxrec = 365
         recslot = 2
         recnum = int(yday)
         mlast = mod(recnum+maxrec-2,maxrec) + 1
         mnext = mod(recnum-1,       maxrec) + 1
         call interp_coeff ( recnum, recslot, secday, dataloc, c1intp, c2intp)

         sst_temp(:) = c1intp *  sst_data(mlast) + c2intp *  sst_data(mnext)
         sss     (:) = c1intp *  sss_data(mlast) + c2intp *  sss_data(mnext)
         uocn    (:) = c1intp * uocn_data(mlast) + c2intp * uocn_data(mnext)
         vocn    (:) = c1intp * vocn_data(mlast) + c2intp * vocn_data(mnext)
         qdp     (:) = c1intp *  qdp_data(mlast) + c2intp *  qdp_data(mnext)
         hmix    (:) = c1intp * hmix_data(mlast) + c2intp * hmix_data(mnext)

      else

         ! use default values for all other data fields
         i = mod(timestep-1,ntime)+1 ! repeat forcing cycle
         mlast = i
         mnext = mlast
         c1intp = c1
         c2intp = c0

         sst_temp(:) = c1intp *  sst_data(mlast) + c2intp *  sst_data(mnext)
         sss     (:) = c1intp *  sss_data(mlast) + c2intp *  sss_data(mnext)
         uocn    (:) = c1intp * uocn_data(mlast) + c2intp * uocn_data(mnext)
         vocn    (:) = c1intp * vocn_data(mlast) + c2intp * vocn_data(mnext)
         qdp     (:) = c1intp *  qdp_data(mlast) + c2intp *  qdp_data(mnext)

      endif
   endif

      call finish_ocn_forcing(sst_temp)

      ! Lindsay SHEBA open/close dataset is hourly
      if (trim(ocn_data_type) == 'SHEBA') then

        sec1hr = secday/c24                      ! seconds in 1 hour
        maxrec = ntime
        recnum = 24*int(yday) - 23 + int(real(sec,kind=dbl_kind)/sec1hr)
        recslot = 2
        dataloc = 1                          ! data located at middle of interval
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff ( recnum, recslot, sec1hr, dataloc, c1intp, c2intp)

        opening(:) =   c1intp * open_data(mlast) + c2intp * open_data(mnext)
        closing(:) = -(c1intp * clos_data(mlast) + c2intp * clos_data(mnext))

      endif

      end subroutine get_forcing

!=======================================================================

      subroutine atm_climatological

      real (kind=dbl_kind), dimension(12) :: &
            fsw_clim, & ! field values at temporal data points
            flw_clim, &
           Tair_clim, &
           wind_clim, &
           rhum_clim, &
          fsnow_clim

      real (kind=dbl_kind) :: &
          Tffresh, qqqice, TTTice, rhos

      character(len=*), parameter :: subname='(atm_climatological)'

      ! Ice station meteorology from Lindsay (1998, J. Climate), Table 1, p. 325
      ! zlvl = c2 ! 2-m temperatures and wind speed

      data  fsw_clim /  0.0d0,   1.2d0,  31.5d0, 146.0d0, 263.3d0, 307.9d0, &
                      230.6d0, 134.7d0,  44.2d0,   2.6d0,   0.0d0,   0.0d0  /
      data  flw_clim /164.0d0, 160.5d0, 164.1d0, 188.1d0, 245.2d0, 291.2d0, &
                      303.9d0, 297.0d0, 263.8d0, 210.9d0, 177.0d0, 166.0d0  /
      data Tair_clim /-31.4d0, -32.8d0, -31.6d0, -24.1d0, -11.0d0,  -1.8d0, &
                       -0.1d0,  -1.4d0,  -8.0d0, -19.5d0, -27.6d0, -31.1d0  /
      data rhum_clim / 78.7d0,  78.4d0,  79.6d0,  82.1d0,  86.5d0,  91.7d0, &
                       95.1d0,  94.3d0,  90.7d0,  83.8d0,  80.1d0,  78.7d0  /
      data wind_clim /  4.4d0,   4.0d0,   4.0d0,   3.9d0,   3.9d0,   4.2d0, &
                        4.1d0,   4.2d0,   4.5d0,   4.2d0,   3.9d0,   4.0d0  /
!      data  shf_clim /  9.9d0,   8.4d0,   6.6d0,   0.1d0,  -5.8d0,  -1.6d0, &
!                        2.2d0,   1.2d0,   0.5d0,   2.0d0,   5.6d0,   7.0d0  /
!      data  lhf_clim /  1.3d0,   1.1d0,   1.1d0,   0.0d0,  -5.9d0, -10.3d0, &
!                       -6.5d0,  -6.7d0,  -3.9d0,  -0.1d0,   1.0d0,   1.1d0  /

      ! Semtner (1976, JPO) snowfall spec., p. 383 in m/s snow volume (.4 m/yr)
      data fsnow_clim/ 3.17d-9, 3.17d-9, 3.17d-9, 3.17d-9, 1.90d-8,    0.0d0, &
                           0.0d0, 1.63d-8, 4.89d-8, 4.89d-8, 3.17d-9, 3.17d-9 /

      !-----------------------------------------------------------------
      ! query icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(Tffresh_out=Tffresh, qqqice_out=qqqice, &
           TTTice_out=TTTice, rhos_out=rhos)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------

       fsw_data (1:12) =  fsw_clim (1:12)
       flw_data (1:12) =  flw_clim (1:12)
      rhum_data (1:12) = rhum_clim (1:12)
      wind_data (1:12) = wind_clim (1:12)

      rhoa_data (1:12) = 1.275_dbl_kind ! air density (kg/m^3)
      Tair_data (1:12) = Tair_clim (1:12) + Tffresh
      uatm_data (1:12) = wind_clim (1:12)
      vatm_data (1:12) = c0

      ! Qa = rhum * saturation humidity (1.275 kg/m^3 = air density)
        Qa_data (1:12) = (rhum_clim(1:12)/c100)*qqqice &
                       * exp(-TTTice/Tair_data(1:12))/rhoa_data(1:12)

      fsnow_data(1:12) = rhos*fsnow_clim(1:12) ! convert vol -> mass flux
      frain_data(1:12) = c0

      ! 6 W/m2 warming of mixed layer from deep ocean
        qdp_data(:) = -6.0 ! 2 W/m2 from deep + 4 W/m2 counteracting larger
                              ! SH+LH with bulk transfer than in MU 71
      ! Warn that this overwrites default and namelist value
      write(nu_diag,*) subname
      write(nu_diag,*) 'WARNING: atm_data_type = clim overwrites '//&
         'oceanic heat flux convergence from default or namelist'

      end subroutine atm_climatological

!=======================================================================

      subroutine atm_CFS

      integer (kind=int_kind) :: &
         nt             ! loop index

      real (kind=dbl_kind) :: &
         dlwsfc,  &     ! downwelling longwave (W/m2)
         dswsfc,  &     ! downwelling shortwave (W/m2)
         windu10, &     ! wind components (m/s)
         windv10, &     !
         temp2m,  &     ! 2m air temperature (K)
         spechum ,&     ! specific humidity (kg/kg)
         precip         ! precipitation (kg/m2/s)

      character (char_len_long) string1
      character (char_len_long) filename
      character(len=*), parameter :: subname='(atm_CFS)'

!      atm_data_file = 'cfsv2_2015_220_70_01hr.txt'
      filename = trim(data_dir)//'/CFS/'//trim(atm_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')
      read (nu_forcing, *) string1 ! headers
      read (nu_forcing, *) string1 ! units

      do nt = 1, ntime
         read (nu_forcing, '(6(f10.5,1x),2(f10.8,1x))') &
         dswsfc, dlwsfc, windu10, windv10, temp2m, spechum, precip

           flw_data(nt) = dlwsfc
           fsw_data(nt) = dswsfc
          uatm_data(nt) = windu10
          vatm_data(nt) = windv10
          Tair_data(nt) = temp2m
          potT_data(nt) = temp2m
            Qa_data(nt) = spechum
         fsnow_data(nt) = precip
      enddo

      close (nu_forcing)

      end subroutine atm_CFS

!=======================================================================

      subroutine prepare_forcing (Tair,     fsw,      &
                                  cldf,     &
                                  frain,    fsnow,    &
                                  Qa,       rhoa,     &
                                  uatm,     vatm,     &
                                  strax,    stray,    &
                                  zlvl,     wind,     &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  potT)

      ! this routine acts on the data fields prior to interpolation

      real (kind=dbl_kind), dimension(ntime), &
         intent(inout) :: &
         Tair    , & ! air temperature  (K)
         fsw     , & ! incoming shortwave radiation (W/m^2)
         cldf    , & ! cloud fraction
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow   , & ! snowfall rate (kg/m^2 s)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         strax   , & ! wind stress components (N/m^2)
         stray   , &
         zlvl    , & ! atm level height (m)
         wind    , & ! wind speed (m/s)
!        flw     , & ! incoming longwave radiation (W/m^2)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         potT        ! air potential temperature  (K)

      ! local variables

      integer (kind=int_kind) :: &
         nt

      logical (kind=log_kind) :: &
         calc_strair

      real (kind=dbl_kind), parameter :: &
         lapse_rate = 0.0065_dbl_kind      ! (K/m) lapse rate over sea level

      real (kind=dbl_kind) :: &
         precip_factor, zlvl0, &
         Tffresh

      character(len=*), parameter :: subname='(prepare_forcing)'

      zlvl0 = c10 ! default

      !-----------------------------------------------------------------
      ! query icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_query_parameters(Tffresh_out=Tffresh)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! convert precipitation units to kg/m^2 s
      !-----------------------------------------------------------------
      if (trim(precip_units) == 'mm_per_month') then
         precip_factor = 12._dbl_kind/(secday*dayyr)
      elseif (trim(precip_units) == 'mm_per_day') then
         precip_factor = c1/secday
      elseif (trim(precip_units) == 'mm_per_sec' .or. &
              trim(precip_units) == 'mks') then
         precip_factor = c1    ! mm/sec = kg/m^2 s
      endif

      do nt = 1, ntime

      !-----------------------------------------------------------------
      ! make sure interpolated values are physically realistic
      !-----------------------------------------------------------------
         cldf (nt) = max(min(cldf(nt),c1),c0)
         fsw  (nt) = max(fsw(nt),c0)
         fsnow(nt) = max(fsnow(nt),c0)
         rhoa (nt) = max(rhoa(nt),c0)
         Qa   (nt) = max(Qa(nt),c0)

      !-----------------------------------------------------------------
      ! calculations specific to datasets
      !-----------------------------------------------------------------

         if (trim(atm_data_type) == 'CFS') then
            ! precip is in kg/m^2/s
            zlvl0 = c10
            ! downward longwave as in Parkinson and Washington (1979)
!            call longwave_parkinson_washington(Tair(nt), cldf(nt), flw(nt))

         elseif (trim(atm_data_type) == 'clim') then
            ! precip is in kg/m^2/s
            zlvl0 = c2

         elseif (trim(atm_data_type) == 'ISPOL' .or. &
                 trim(atm_data_type) == 'NICE') then
            Tair(nt) = Tair(nt) -  lapse_rate*8.0_dbl_kind

         endif                     ! atm_data_type

! this longwave depends on the current ice aice and sst and so can not be
! computed ahead of time
!            ! longwave based on Rosati and Miyakoda, JPO 18, p. 1607 (1988)
!            call longwave_rosati_miyakoda(cldf(i,j), Tsfc(i,j), &
!                                          aice(i,j), sst(i,j),  &
!                                          Qa(i,j),   Tair(i,j), &
!                                          hm(i,j),   flw(i,j))

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

         zlvl(nt) = zlvl0
         potT(nt) = Tair(nt)

         ! divide shortwave into spectral bands
         swvdr(nt) = fsw(nt)*frcvdr        ! visible direct
         swvdf(nt) = fsw(nt)*frcvdf        ! visible diffuse
         swidr(nt) = fsw(nt)*frcidr        ! near IR direct
         swidf(nt) = fsw(nt)*frcidf        ! near IR diffuse

         ! precipitation
         fsnow(nt) = fsnow(nt) * precip_factor

         ! determine whether precip is rain or snow
         ! HadGEM forcing provides separate snowfall and rainfall rather
         ! than total precipitation
!         if (trim(atm_data_type) /= 'hadgem') then
            frain(nt) = c0
            if (Tair(nt) >= Tffresh) then
                frain(nt) = fsnow(nt)
                fsnow(nt) = c0
            endif
!         endif

         if (calc_strair) then
            wind (nt) = sqrt(uatm(nt)**2 + vatm(nt)**2)
            strax(nt) = c0
            stray(nt) = c0
         ! else  ! strax, stray, wind are read from files
         endif                   ! calc_strair

      enddo ! ntime

      end subroutine prepare_forcing

!=======================================================================

      subroutine interp_coeff_monthly (recslot, c1intp, c2intp)

! Compute coefficients for interpolating monthly data to current time step.

      integer (kind=int_kind), intent(in) :: &
          recslot         ! slot (1 or 2) for current record

      real (kind=dbl_kind), intent(inout) :: &
         c1intp, c2intp   ! interpolation coefficients

      ! local variables

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2           ! seconds elapsed at month midpoint

      real (kind=dbl_kind) :: &
          daymid(0:13)     ! month mid-points

      character(len=*), parameter :: subname='(interp_coeff_monthly)'

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = 14._dbl_kind - daymo(12)  ! Dec 15, 0 sec

      ! make time cyclic
      tt = mod(time/secday,dayyr)

      ! Find neighboring times

      if (recslot==2) then      ! first half of month
        t2 = daycal(month) + daymid(month)   ! midpoint, current month
        if (month == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
        else
          t1 = daycal(month-1) + daymid(month-1) ! midpoint, previous month
        endif
      else                      ! second half of month
        t1 = daycal(month) + daymid(month)    ! midpoint, current month
        t2 = daycal(month+1) + daymid(month+1)! day 15 of next month (0 sec)
      endif

      ! Compute coefficients
      c1intp = (t2 - tt) / (t2 - t1)
      c2intp =  c1 - c1intp

      end subroutine interp_coeff_monthly

!=======================================================================

      subroutine interp_coeff (recnum, recslot, secint, dataloc, &
                               c1intp, c2intp, offset)

! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.

      integer (kind=int_kind), intent(in) :: &
          recnum      , & ! record number for current data value
          recslot     , & ! spline slot for current record
          dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

      real (kind=dbl_kind), intent(in) :: &
          secint          ! seconds in data interval

      real (kind=dbl_kind), intent(inout) :: &
         c1intp, c2intp   ! interpolation coefficients

      real (kind=dbl_kind), intent(in), optional :: &
          offset          ! amount of time data is offset (s)

      ! local variables

      real (kind=dbl_kind) :: &
          secyr            ! seconds in a year

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2       , & ! seconds elapsed at data points
          rcnum            ! recnum => dbl_kind

      character(len=*), parameter :: subname='(interp_coeff)'

      secyr = dayyr * secday         ! seconds in a year
      if (present(offset)) then
         tt = mod(time-offset,secyr)
      else
         tt = mod(time,secyr)
      endif

      ! Find neighboring times
      rcnum = real(recnum,kind=dbl_kind)
      if (recslot==2) then           ! current record goes in slot 2
         if (dataloc==1) then        ! data located at middle of interval
            t2 = (rcnum-p5)*secint
         else                        !  data located at end of interval
            t2 = rcnum*secint
         endif
         t1 = t2 - secint            !  - 1 interval
      else                           ! recslot = 1
         if (dataloc==1) then        ! data located at middle of interval
            t1 = (rcnum-p5)*secint
         else
            t1 = rcnum*secint        ! data located at end of interval
         endif
         t2 = t1 + secint            !  + 1 interval
      endif

      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      end subroutine interp_coeff

!=======================================================================

      subroutine atm_ISPOL
      ! there is data for 366 days, but we only use 365

      integer (kind=int_kind) :: &
         i        ! index

      real (kind=dbl_kind), dimension(366) :: &
         tair , & ! air temperature
         qa   , & ! specific humidity
         uatm , & ! wind velocity components
         vatm , &
         fsnow    ! snowfall rate
         ! aday   ! not used

      real (kind=dbl_kind), dimension(1464) :: &
         fsw  , & ! shortwave
         flw      ! longwave
         ! atime  ! not used

      character (char_len_long) filename

      character(len=*), parameter :: subname='(atm_ISPOL)'

!      atm_data_file = 'ISPOL_atm_forcing.txt'
      filename = trim(data_dir)//'/ISPOL_2004/'//trim(atm_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) tair
      read(nu_forcing,*) qa
      read(nu_forcing,*) fsw
      read(nu_forcing,*) flw
      read(nu_forcing,*) uatm
      read(nu_forcing,*) vatm
      read(nu_forcing,*) fsnow
!      read(nu_forcing,*) aday
!      read(nu_forcing,*) atime

      do i = 1, 366 ! daily
         Tair_data (i) = tair (i)
         Qa_data   (i) = qa   (i)
         uatm_data (i) = uatm (i)
         vatm_data (i) = vatm (i)
         fsnow_data(i) = fsnow(i)
      end do
      do i = 1, 1464 ! 6hr, 1464/4=366 days
         fsw_data  (i) = fsw  (i)
         flw_data  (i) = flw  (i)
      end do

      close(nu_forcing)

      end subroutine atm_ISPOL

!=======================================================================

      subroutine atm_NICE
      ! there is data for 366 days, but we only use 365

      integer (kind=int_kind) :: &
         i        ! index

      real (kind=dbl_kind), dimension(366) :: &
         tair , & ! air temperature
         qa   , & ! specific humidity
         uatm , & ! wind velocity components
         vatm , &
         fsnow    ! snowfall rate
         ! aday   ! not used

      real (kind=dbl_kind), dimension(1464) :: &
         fsw  , & ! shortwave
         flw      ! longwave
         ! atime  ! not used

      character (char_len_long) filename

      character(len=*), parameter :: subname='(atm_NICE)'

!      atm_data_file = 'NICE_atm_forcing.txt'
      filename = trim(data_dir)//'/NICE_2015/'//trim(atm_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) tair
      read(nu_forcing,*) qa
      read(nu_forcing,*) fsw
      read(nu_forcing,*) flw
      read(nu_forcing,*) uatm
      read(nu_forcing,*) vatm
      read(nu_forcing,*) fsnow
!      read(nu_forcing,*) aday
!      read(nu_forcing,*) atime

      do i = 1, 366
         Tair_data (i) = tair (i)
         Qa_data   (i) = qa   (i)
         uatm_data (i) = uatm (i)
         vatm_data (i) = vatm (i)
         fsnow_data(i) = fsnow(i)
      end do
      do i = 1, 1464
         fsw_data  (i) = fsw  (i)
         flw_data  (i) = flw  (i)
      end do

      close(nu_forcing)

      end subroutine atm_NICE

!=======================================================================

      subroutine atm_MDF

      integer (kind=int_kind) :: &
         nt,      &  ! timestep index for Icepack arrays
         i,       &  ! index for forcing data arrays
         bound,   &  ! bound for subsetting data
         dimlen,  &  ! length of the data arrays
         ncid,    &  ! NetCDF file id
         dimid,   &  ! NetCDF dimension id
         status,  &  ! NetCDF status flag
         varid       ! NetCDF variable id

      integer (kind=8), allocatable :: &
         data_time(:)   ! array for time array in forcing data

      integer (kind=8), dimension(ntime) :: &
         model_time ! array for Icepack minutely time

      real (kind=dbl_kind) :: &
         model_time0    ! start time for model in seconds since 1970

      real (kind=dbl_kind), allocatable :: &
         data(:)  ! data array from file

      character (char_len) :: &
         calendar_type, &  ! data calendar type
         varname

      character (char_len_long) :: &
         filename, &
         time_basis     ! time basis for data

      integer (kind=int_kind), dimension(ntime, 2) :: &
         data_sections  ! 2D array for indices corresponding
                        ! to which data values should be averaged to
                        ! create the model forcing values

      real (kind=dbl_kind), parameter :: &
         Gregorian_year = 365.2425, &  ! days in Gregorian year per cf standard
         model_miss_val = -9999.00     ! missing value for internal use

      character(len=*), parameter :: subname='(atm_MDF)'

      filename = trim(data_dir)//'/MDF/'//trim(atm_data_file)

      if (atm_data_format /= 'nc') then
         call icedrv_system_abort(string=subname//&
         ' ERROR: only NetCDF input implemented for atm_MDF', &
         file=__FILE__,line=__LINE__)
      else
#ifdef USE_NETCDF
         ! Open forcing file
         status = nf90_open(trim(filename), nf90_nowrite, ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(&
            string=subname//'Couldnt open netcdf file', &
                           file=__FILE__,line=__LINE__)
         
         ! Create array for the time step values in seconds since 1970
         ! CF standard calendar is Gregorian
         ! May have strange behavior if dt is not an integer
         model_time0 = (year_init - 1970) * Gregorian_year * 24 * 3600 + time0
         do nt = 1, ntime
            model_time(nt) = int(model_time0 + dt * nt, kind=8)
         enddo

         ! Read, average, and interpolate forcing data from each variable
         ! Moving average forcing values into model arrays
         call load_var_MDF("tas", Tair_data, ncid, model_time)
         call load_var_MDF("hus", Qa_data, ncid, model_time)
         call load_var_MDF("uas", uatm_data, ncid, model_time)
         call load_var_MDF("vas", vatm_data, ncid, model_time)
         call load_var_MDF("rlds", flw_data, ncid, model_time)
         call load_var_MDF("rsds", fsw_data, ncid, model_time)

         ! Precipitation data is optional
         ! Check whether snowfall rate is in the dataset
         status = nf90_inq_varid(ncid, "prsn", varid)
         if (status /= nf90_noerr) then
            write(nu_diag,*) subname
            write(nu_diag,*) 'WARNING: no snowfall rate in forcing.'// &
               ' Icepack will assume snowfall rate is zero.'
            fsnow_data(:) = c0
         else
            call load_var_MDF("prsn", fsnow_data, ncid, model_time)
         endif
         ! Check whether total precipitation rate is in the dataset
         status = nf90_inq_varid(ncid, "pr", varid)
         if (status /= nf90_noerr) then
            write(nu_diag,*) subname
            write(nu_diag,*) 'WARNING: no precipitation rate in '// &
               'forcing. Icepack will assume rainfall rate is zero.'
            frain_data(:) = c0
         else
            call load_var_MDF("pr", frain_data, ncid, model_time)
         endif
         ! rainfall is total precipitation minus snowfall
         ! subject to the constraint that rainfall is non-negative
         do nt = 1, ntime
            frain_data(nt) = max(frain_data(nt) - fsnow_data(nt), c0)
         enddo

#else
         call icedrv_system_abort(string=subname//&
         ' ERROR: atm_data_format = "nc" requires USE_NETCDF', &
         file=__FILE__,line=__LINE__)
#endif
      endif

      end subroutine atm_MDF

!=======================================================================

#ifdef USE_NETCDF
      subroutine MDF_average(data_var_name, model_var_arr, &
         data_var_len, ncid, data_sections, model_miss_val)

      character(len=*), intent(in) :: &
         data_var_name  ! Name of the variable in the MDF forcing file

      real (kind=dbl_kind), dimension(ntime), intent(out) :: &
         model_var_arr  ! array to place averaged forcing data in

      integer (kind=int_kind), intent(in) :: &
         data_var_len, &   ! Size of data array in MDF forcing file
         ncid              ! NetCDF file id

      integer (kind=int_kind), dimension(ntime, 2) :: &
         data_sections     ! indices for which data values to average

      real (kind=dbl_kind), intent(in) :: &
         model_miss_val ! for when there is no data in a time step

      ! Local variables
      real (kind=dbl_kind), dimension(data_var_len) :: &
         data_var_arr      ! array for data from forcing file

      real (kind=dbl_kind) :: &
         work, &           ! variable for averaging
         data_miss_val, &  ! value of missing data
         count             ! counter for data to average

      integer (kind=int_kind) :: &
         status, &         ! NetCDF status flag
         nt,     &         ! timestep index for Icepack arrays
         i,      &         ! index for forcing data arrays
         varid             ! NetCDF variable id

      character(len=*), parameter :: subname='(MDF_average)'

      ! Allocate get data and missing value from file
      status = nf90_inq_varid(ncid, trim(data_var_name), varid)
      if (status /= nf90_noerr) call icedrv_system_abort(&
            string=subname//'Couldnt get '//data_var_name//' var id', &
                           file=__FILE__,line=__LINE__)
      status = nf90_get_att(ncid, varid, "missing_value", data_miss_val)
      if (status /= nf90_noerr) call icedrv_system_abort(&
            string=subname//'Couldnt get '//data_var_name//' missing value', &
                           file=__FILE__,line=__LINE__)
      status = nf90_get_var(ncid, varid, data_var_arr)
      if (status /= nf90_noerr) call icedrv_system_abort(&
            string=subname//'Couldnt get '//data_var_name//' values', &
                           file=__FILE__,line=__LINE__)

      ! For each model time point average non-missing data values
      do nt = 1, ntime
         count = 0
         work = c0
         do i = data_sections(nt, 1), data_sections(nt, 2)
            if (data_var_arr(i) /= data_miss_val) then
               work = work + data_var_arr(i)
               count = count + 1
            endif
         end do
         if (count > 0) then
            model_var_arr(nt) = work / count
         else
            model_var_arr(nt) = model_miss_val
         endif
      end do

      end subroutine MDF_average
#endif

!=======================================================================

      subroutine MDF_interpolate(model_var_arr, model_miss_val)

      real (kind=dbl_kind), dimension(ntime), intent(inout) :: &
         model_var_arr  ! array to place averaged forcing data in

      real (kind=dbl_kind), intent(in) :: &
         model_miss_val ! for when there is no data in a time step

      integer (kind=int_kind) :: &
         mlast,         &  ! index of last present data
         nt, m,         &  ! model timestep indices
         count             ! counter for missing values

      character(len=*), parameter :: subname='(MDF_interpolate)'

      ! Interpolate, extrapolate for first and last values
      if (model_var_arr(1) == model_miss_val) then
         mlast = 0
      else
         mlast = 1
      endif
      do nt = 2, ntime
         if (model_var_arr(nt) == model_miss_val) then
            ! Do nothing (i.e., allow nt to increment) unless we're at end
            if (nt == ntime) then
               do m = mlast + 1, nt
                  model_var_arr(m) = model_var_arr(mlast)
               end do
            endif
         else if ((nt - mlast) == 1) then
            ! No missing data, increment mlast
            mlast = nt
         else
            ! If we're at the start extrapolate to fill
            if (mlast==0) then
               do m = mlast + 1, nt - 1
                  model_var_arr(m) = model_var_arr(nt)
               end do
            else
               ! Interpolate missing data
               do m = mlast + 1, nt - 1
                  model_var_arr(m) = model_var_arr(mlast) &
                     + (model_var_arr(nt) - model_var_arr(mlast)) &
                     * (m - mlast) / (nt - mlast)
               end do
            endif
            mlast = nt
         endif
      end do

      end subroutine MDF_interpolate
!=======================================================================

      subroutine load_var_MDF(data_var_name, model_var_arr, ncid, &
                                 model_time)

      character(len=*), intent(in) :: &
         data_var_name  ! Name of the variable in the MDF forcing file

      real (kind=dbl_kind), dimension(ntime), intent(out) :: &
         model_var_arr  ! array to place forcing data in

      integer (kind=int_kind), intent(in) :: &
         ncid           ! NetCDF file id

      integer (kind=8), dimension(ntime), intent(in) :: &
         model_time     ! model time array

      ! Local variables

      integer (kind=int_kind) :: &
         nt,      &  ! timestep index for Icepack arrays
         i,       &  ! index for forcing data arrays
         bound,   &  ! bound for subsetting data
         dimlen,  &  ! length of the data arrays
         dimid,   &  ! NetCDF dimension id
         status,  &  ! NetCDF status flag
         nvardims,&  ! number of dimensions for variable
         varid       ! NetCDF variable id

      integer (kind=8), allocatable :: &
         data_time(:)   ! array for time array in forcing data

      integer, dimension(1) :: &
         vardimids      ! dimension id for variable

      character (char_len) :: &
         calendar_type, &  ! data calendar type
         dimname           ! name for variables dimension

      character (char_len_long) :: &
         time_basis     ! time basis for data

      integer (kind=int_kind), dimension(ntime, 2) :: &
         data_sections  ! 2D array for indices corresponding
                        ! to which data values should be averaged to
                        ! create the model forcing values

      real (kind=dbl_kind), parameter :: &
         Gregorian_year = 365.2425, &  ! days in Gregorian year per cf standard
         model_miss_val = -9999.00     ! missing value for internal use      

      character(len=*), parameter :: subname='(load_var_MDF)'

#ifdef USE_NETCDF

      ! Get varid and missing value from file
      status = nf90_inq_varid(ncid, trim(data_var_name), varid)
      if (status /= nf90_noerr) call icedrv_system_abort(&
            string=subname//'Couldnt get '//data_var_name//' var id', &
                           file=__FILE__,line=__LINE__)
      ! Get information about the dimension for this variable
      status = nf90_inquire_variable(ncid, varid, ndims=nvardims)
      if (nvardims /= 1) call icedrv_system_abort(&
            string=subname//data_var_name//' has more than 1 dimension', &
                           file=__FILE__,line=__LINE__)
      status = nf90_inquire_variable(ncid, varid, dimids=vardimids)
      status = nf90_inquire_dimension(ncid, vardimids(1), &
            name=dimname, len=dimlen)
      ! Check that dimname matches pattern
      if (dimname(1:4) /= 'time') call icedrv_system_abort(&
            string=subname//data_var_name//' dimension name is not timeXXXX', &
                           file=__FILE__,line=__LINE__)
      if (verify(trim(dimname(5:)), "0123456789") /= 0) call icedrv_system_abort(&
            string=subname//data_var_name//' dimension name is not timeXXXX', &
                           file=__FILE__,line=__LINE__)
      allocate (data_time(dimlen))
      ! Check that cadence variable exists and calendars match
      status = nf90_inq_varid(ncid, trim(dimname), varid)
      if (status /= nf90_noerr) call icedrv_system_abort(&
         string=subname//'Couldnt get '//trim(dimname)//' var id', &
                        file=__FILE__,line=__LINE__)
      status = nf90_get_att(ncid, varid, "calendar", calendar_type)
      if (status /= nf90_noerr) call icedrv_system_abort(&
         string=subname//'Couldnt get calendar attribute', &
                        file=__FILE__,line=__LINE__)
      ! In future this check could be replaced with calendar matching
      if (calendar_type /= "standard" .or. .not. use_leap_years) then
         call icedrv_system_abort(&
         string=subname//'Forcing calendar not standard or not using leap years',&
         file=__FILE__,line=__LINE__)
      endif
      ! Get the time array
      !! Note, in the file the value is actually unsigned, need to make sure this
      ! doesn't cause issues since Fortran 90 doesn't support unsigned ints.
      status = nf90_get_var(ncid, varid, data_time)
      if (status /= nf90_noerr) call icedrv_system_abort(&
         string=subname//'Couldnt read '//trim(dimname)//' values', &
                        file=__FILE__,line=__LINE__)
      ! Convert the data time from minutes into seconds for compatability w/ icepack
      data_time = data_time * 60

      ! Check that the time basis in forcing file is correct
      status = nf90_get_att(ncid, varid, "units", time_basis)
      if (status /= nf90_noerr) call icedrv_system_abort(&
         string=subname//'Couldnt get '//trim(dimname)//' units', &
                        file=__FILE__,line=__LINE__)
      if (time_basis /= "minutes since 1970-01-01 00:00:00") then
         call icedrv_system_abort(&
         string=subname//'Time basis is not minutes since 1970',&
         file=__FILE__,line=__LINE__)
      endif
      
      ! Check that we are not extrapolating forcing outside of time bounds
      if (model_time(1) < data_time(1)) call icedrv_system_abort(&
      string=subname//'Simulation starts before forcing',&
      file=__FILE__,line=__LINE__)
      if (model_time(ntime) > data_time(dimlen)) call icedrv_system_abort(&
      string=subname//'Simulation ends after forcing',&
      file=__FILE__,line=__LINE__)

      ! data_sections is a 2D array where the first dimension
      ! is the same length as model_time. The 1D array at each index
      ! contains the start and stop indices of the data to be averaged
      ! into each model timestep
      ! Get the first start index
      bound = model_time(1) - (model_time(2) - model_time(1))/2
      i = 1
      do while (data_time(i) < bound)
         i = i + 1
      end do
      data_sections(1, 1) = i
      do nt = 1, ntime - 1
         ! Bound is halfway between this time step and the next
         bound = (model_time(nt + 1) + model_time(nt))/2
         do while (data_time(i) < bound)
            i = i + 1
         end do ! i - 1 is now the last element in timestep nt
         data_sections(nt, 2) = i - 1
         data_sections(nt + 1, 1) = i
      end do
      ! Get the last index
      bound = model_time(ntime) + (model_time(ntime) - model_time(ntime - 1))/2
      i = dimlen
      do while (data_time(i) > bound)
         i = i - 1
      end do
      data_sections(ntime, 2) = i

      ! Moving average forcing values into model arrays
      call MDF_average(data_var_name, model_var_arr, dimlen, ncid, &
         data_sections, model_miss_val)
      ! Linearly interpolate missing values
      call MDF_interpolate(model_var_arr, model_miss_val)

#else
      call icedrv_system_abort(string=subname//&
      ' load_var_MDF requires USE_NETCDF', &
      file=__FILE__,line=__LINE__)
#endif

      end subroutine load_var_MDF

!=======================================================================

      subroutine ocn_NICE

      integer (kind=int_kind) :: &
         i

      real (kind=dbl_kind), dimension(365) :: &
         t   , &  ! sea surface temperature
         s   , &  ! sea surface salinity
         hblt, &  ! mixed layer depth
         u   , &  ! ocean current, x
         v   , &  ! ocean current, y
         dhdx, &  ! sea surface slope
         dhdy, &  ! sea surface slope
         qdp      ! deep ocean heat flux

      character (char_len_long) filename

      character(len=*), parameter :: subname='(ocn_NICE)'

!      ocn_data_file = 'oceanmixed_daily_3.txt'
      filename = trim(data_dir)//'/NICE_2015/'//trim(ocn_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) t
      read(nu_forcing,*) s
      read(nu_forcing,*) hblt
      read(nu_forcing,*) u
      read(nu_forcing,*) v
      read(nu_forcing,*) dhdx  ! not used for Icepack
      read(nu_forcing,*) dhdy  ! not used for Icepack
      read(nu_forcing,*) qdp

      close(nu_forcing)

      do i = 1, 365 ! daily
         sst_data (i) = t   (i)
         sss_data (i) = s   (i)
         hmix_data(i) = hblt(i)
         uocn_data(i) = u   (i)
         vocn_data(i) = v   (i)
         qdp_data (i) = qdp (i)
      end do

      end subroutine ocn_NICE

!=======================================================================

      subroutine ocn_ISPOL

      integer (kind=int_kind) :: &
         i

      real (kind=dbl_kind), dimension(12) :: &
         t   , &  ! sea surface temperature
         s   , &  ! sea surface salinity
         hblt, &  ! mixed layer depth
         u   , &  ! ocean current, x
         v   , &  ! ocean current, y
         dhdx, &  ! sea surface slope
         dhdy, &  ! sea surface slope
         qdp      ! deep ocean heat flux

      character (char_len_long) filename

      character(len=*), parameter :: subname='(ocn_ISPOL)'

!      ocn_data_file = 'pop_frc.gx1v3.051202_but_hblt_from_010815_ispol.txt'
      filename = trim(data_dir)//'/ISPOL_2004/'//trim(ocn_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_forcing, file=filename, form='formatted')

      read(nu_forcing,*) t
      read(nu_forcing,*) s
      read(nu_forcing,*) hblt
      read(nu_forcing,*) u
      read(nu_forcing,*) v
      read(nu_forcing,*) dhdx  ! not used for Icepack
      read(nu_forcing,*) dhdy  ! not used for Icepack
      read(nu_forcing,*) qdp

      close(nu_forcing)

      do i = 1, 12 ! monthly
         sst_data (i) = t   (i)
         sss_data (i) = s   (i)
         hmix_data(i) = hblt(i)
         uocn_data(i) = u   (i)
         vocn_data(i) = v   (i)
         qdp_data (i) = qdp (i)
      end do

    end subroutine ocn_ISPOL

!=======================================================================

      subroutine ocn_MDF

      integer (kind=int_kind) :: &
         nt,      &  ! timestep index for Icepack arrays
         i,       &  ! index for forcing data arrays
         bound,   &  ! bound for subsetting data
         dimlen,  &  ! length of the data arrays
         ncid,    &  ! NetCDF file id
         dimid,   &  ! NetCDF dimension id
         status,  &  ! NetCDF status flag
         varid       ! NetCDF variable id

      integer (kind=8), allocatable :: &
         data_time(:)   ! array for time array in forcing data

      integer (kind=8), dimension(ntime) :: &
         model_time ! array for Icepack minutely time

      real (kind=dbl_kind) :: &
         model_time0    ! start time for model in seconds since 1970

      real (kind=dbl_kind), allocatable :: &
         data(:)  ! data array from file

      character (char_len) :: &
         calendar_type, &  ! data calendar type
         varname

      character (char_len_long) :: &
         filename, &
         time_basis     ! time basis for data

      integer (kind=int_kind), dimension(ntime, 2) :: &
         data_sections  ! 2D array for indices corresponding
                        ! to which data values should be averaged to
                        ! create the model forcing values

      real (kind=dbl_kind), parameter :: &
         Gregorian_year = 365.2425, &  ! days in Gregorian year per cf standard
         model_miss_val = -9999.00, &  ! missing value for internal use
         leg4_end_time  = 1596034800, &! end of leg 4 in seconds since 1970
         leg5_start_time= 1598451600   ! start of leg 5 in seconds since 1970

      character(len=*), parameter :: subname='(ocn_MDF)'

      filename = trim(data_dir)//'/MDF/'//trim(ocn_data_file)

      if (ocn_data_format /= 'nc') then
         call icedrv_system_abort(string=subname//&
         ' ERROR: only NetCDF input implemented for ocn_MDF', &
         file=__FILE__,line=__LINE__)
      else
#ifdef USE_NETCDF
         ! Open forcing file
         status = nf90_open(trim(filename), nf90_nowrite, ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(&
            string=subname//'Couldnt open netcdf file', &
                           file=__FILE__,line=__LINE__)
         
         ! Create array for the time step values in seconds since 1970
         ! CF standard calendar is Gregorian
         ! May have strange behavior if dt is not an integer
         model_time0 = (year_init - 1970) * Gregorian_year * 24 * 3600 + time0
         do nt = 1, ntime
            model_time(nt) = int(model_time0 + dt * nt, kind=8)
         enddo

         ! Warn if simulation includes leg 4-5 transition
         if ((model_time(1) < leg5_start_time) .and. &
             (model_time(ntime) > leg4_end_time) .and. &
             (index(ocn_data_file, 'MOSAiC') > 0)) then
            write(nu_diag,*) subname
            write(nu_diag,*) 'WARNING: Forcing may be from MOSAiC '// &
               'and time includes MOSAIC leg 4-5 repositioning. ' // &
               'If so, forcing interpolation is not valid.'
         endif

         ! Read, average, and interpolate forcing data from each variable
         ! Moving average forcing values into model arrays
         call load_var_MDF("so", sss_data, ncid, model_time)
         call load_var_MDF("mlotst", hmix_data, ncid, model_time)
         call load_var_MDF("hfsot", qdp_data, ncid, model_time)
         call load_var_MDF("tos", sst_data, ncid, model_time)

#else
         call icedrv_system_abort(string=subname//&
         ' ERROR: ocn_data_format = "nc" requires USE_NETCDF', &
         file=__FILE__,line=__LINE__)
#endif
      endif

      end subroutine ocn_MDF

!=======================================================================

      subroutine finish_ocn_forcing(sst_temp)

! Compute ocean freezing temperature Tf based on tfrz_option
! 'minus1p8'         Tf = -1.8 C
! 'constant'         Tf = Tocnfrz
! 'linear_salt'      Tf = -depressT * sss
! 'mushy'            Tf conforms with mushy layer thermo (ktherm=2)

      real (kind=dbl_kind), dimension(nx), intent(in)  :: &
          sst_temp

      ! local variables

      integer (kind=int_kind) :: &
         i           ! horizontal indices

      character(len=*), parameter :: subname='(finish_ocn_forcing)'

      do i = 1, nx
         sss (i) = max (sss(i), c0)
         hmix(i) = max(hmix(i), c0)
         Tf  (i) = icepack_sea_freezing_temperature(sss(i))
         if (restore_ocn) sst(i) = sst(i) + (sst_temp(i)-sst(i))*dt/trest
      enddo
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      end subroutine finish_ocn_forcing

!=======================================================================

     subroutine ice_open_clos


      integer (kind=int_kind) :: i

      real (kind=dbl_kind) :: xtime

      character (char_len_long) filename

!      ice_data_file = 'open_clos_lindsay.dat'
      filename = trim(data_dir)//'/SHEBA/'//trim(ice_data_file)

      write (nu_diag,*) 'Reading ',filename

      open (nu_open_clos, file=filename, form='formatted')

      ! hourly data
      do i=1,ntime
         read(nu_open_clos,*) xtime, open_data(i), clos_data(i)
      enddo

      close (nu_open_clos)

     end subroutine ice_open_clos

!=======================================================================

      subroutine get_wave_spec

      use icedrv_arrays_column, only: wave_spectrum, wave_sig_ht, &
                                   dwavefreq, wavefreq
      use icedrv_domain_size, only: nfreq

      ! local variables
      integer (kind=int_kind) :: &
         k

      real(kind=dbl_kind), dimension(nfreq) :: &
         wave_spectrum_profile  ! wave spectrum

       wave_spectrum(:,:) = c0

      ! wave spectrum and frequencies
      ! get hardwired frequency bin info and a dummy wave spectrum profile

      call icepack_init_wave(nfreq=nfreq,                 &
                             wave_spectrum_profile=wave_spectrum_profile, &
                             wavefreq=wavefreq, dwavefreq=dwavefreq)

      do k = 1, nfreq
          wave_spectrum(:,k) = wave_spectrum_profile(k)
      enddo

      end subroutine get_wave_spec

!=======================================================================

      end module icedrv_forcing

!=======================================================================
