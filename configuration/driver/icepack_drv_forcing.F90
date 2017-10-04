!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
      module icepack_drv_forcing

      use icepack_kinds_mod
      use icepack_drv_domain_size, only: ncat, nx
      use icepack_drv_calendar, only: nyr, days_per_year, dayyr, month, &
                              daymo, daycal, dt
!      use ice_calendar, only: istep, istep1, time, time_forc, year_init, &
!                              sec, mday, nyr, yday
      use icepack_drv_constants, only: nu_diag, nu_forcing, secday
      use icepack_intfc_shared, only: calc_strair

      implicit none
      private
      public :: init_forcing, get_forcing, read_clim_data, &
          interpolate_data, interp_coeff, interp_coeff_monthly, &
          read_data_point
!                , read_clim_data_nc, &
!                read_data_nc_point
      save

      integer (kind=int_kind), public :: &
         ntime           , & ! number of data points in time
         ycycle          , & ! number of years in forcing cycle
         fyear_init      , & ! first year of data in forcing cycle
         fyear           , & ! current year in forcing cycle
         fyear_final         ! last year in cycle

!      real (kind=dbl_kind), public  :: &
!           c1intp, c2intp , & ! interpolation coefficients
!           ftime              ! forcing time (for restart)

!      integer (kind=int_kind) :: &
!           oldrecnum = 0  , & ! old record number (save between steps)
!           oldrecnum4X = 0    !

!      real (kind=dbl_kind), dimension(ntime) ::
      real (kind=dbl_kind), dimension(8760) :: & ! hardwired for now
            fsw_data, & ! field values at temporal data points
           cldf_data, &
          fsnow_data, &
           Tair_data, &
           uatm_data, &
           vatm_data, &
           wind_data, &
          strax_data, &
          stray_data, &
             Qa_data, &
           rhoa_data, &
           potT_data, &
            flw_data, &
            sst_data, &
            sss_data, & 
           uocn_data, &
           vocn_data, &
         sublim_data, &
          frain_data, &
           zlvl_data, &
          swvdr_data, &
          swvdf_data, &
          swidr_data, &
          swidf_data

      real (kind=dbl_kind), public  :: &
           c1intp, c2intp , & ! interpolation coefficients
           ftime              ! forcing time (for restart)

      character(char_len), public :: & 
         atm_data_format, & ! 'bin'=binary or 'nc'=netcdf
         ocn_data_format, & ! 'bin'=binary or 'nc'=netcdf
         bgc_data_format, & ! 'bin'=binary or 'nc'=netcdf
         atm_data_type,   & ! 'default', 'monthly', 'ncar', 
                            ! 'LYq' or 'hadgem' or 'oned'
         ocn_data_type,   & ! 'default', 'clim', 'ncar', 'oned'
         bgc_data_type,   & ! 'default', 'clim', 'ncar', 'oned',
                            ! 'hadgem_sst' or 'hadgem_sst_uvocn'
         precip_units       ! 'mm_per_month', 'mm_per_sec', 'mks'
 
      character(char_len_long), public :: & 
         data_dir           ! top directory for forcing data

      ! as in the dummy atm (latm)
      real (kind=dbl_kind), parameter, public :: &
         frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
         frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
         frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
         frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band

      logical (kind=log_kind), public :: &
         restore_sst                 ! restore sst if true

      integer (kind=int_kind), public :: &
         trestore                    ! restoring time scale (days)

      real (kind=dbl_kind), public :: & 
         trest                       ! restoring time scale (sec)
      logical (kind=log_kind), public :: &
         dbug             ! prints debugging output if true

!=======================================================================

      contains

!=======================================================================

      subroutine init_forcing

! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the forcing data.

      use icepack_constants, only: c0
      use icepack_drv_flux, only: zlvl, Tair, potT, rhoa, uatm, vatm, wind, &
         strax, stray, fsw, swvdr, swvdf, swidr, swidf, Qa, flw, frain, &
         fsnow, sst, sss, uocn, vocn

      integer (kind=int_kind) :: &
         i                ! index

      fyear       = fyear_init + mod(nyr-1,ycycle) ! current year
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      write (nu_diag,*) ' Initial forcing data year = ',fyear_init
      write (nu_diag,*) ' Final   forcing data year = ',fyear_final

    !-------------------------------------------------------------------
    ! Initialize forcing data to default values
    !-------------------------------------------------------------------
! maybe these should all be zero, and set defaults in atm_GOFS to be clear
! which ones are defaults and which are being read in

      ! many default forcing values are set in init_flux_atm
      i = 1 ! use first grid box value
          zlvl_data(:) = zlvl (i)    ! atmospheric level height (m)
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
           sst_data(:) = sst  (i)    ! sea surface temperature
           sss_data(:) = sst  (i)    ! sea surface salinity
          uocn_data(:) = uocn (i)    ! wind velocity components (m/s)
          vocn_data(:) = vocn (i)

          cldf_data(:) = c0     ! cloud fraction

      if (trim(atm_data_type) == 'GOFS') call atm_GOFS
      if (trim(atm_data_type) == 'ISPOL') call atm_ISPOL


      call prepare_forcing (Tair_data,     fsw_data,      &    
                            cldf_data,     flw_data,      &
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
! We will probably need to send in the time and working out what the data
! time slice is, instead of sending in the timestep.  This currently assumes
! the time step and the data both start Jan 1.
! Interpolate if necessary - for now, this assumes the data and timesteps match.

!      use icepack_constants, only: c0
      use icepack_drv_flux, only: zlvl, Tair, potT, rhoa, uatm, vatm, wind, &
         strax, stray, fsw, swvdr, swvdf, swidr, swidf, Qa, flw, frain, &
         fsnow, sst, sss, uocn, vocn
      use icepack_intfc_shared, only: restore_bgc

      integer (kind=int_kind), intent(in) :: &
         timestep         ! time step index

      integer (kind=int_kind) :: &
         i                ! data index

      if (trim(atm_data_type) == 'default') return

      ! calculate data index corresponding to current timestep
      i = mod(timestep-1,ntime)+1 ! repeat forcing cycle

      ! fill all grid boxes with the same forcing data
      flw  (:) =   flw_data(i)
      fsw  (:) =   fsw_data(i)
      Tair (:) =  Tair_data(i)
      Qa   (:) =    Qa_data(i)
      fsnow(:) = fsnow_data(i)

      zlvl (:) = zlvl_data (i)    ! atmospheric level height (m)
      Tair (:) = Tair_data (i)    ! air temperature  (K)
      potT (:) = potT_data (i)    ! air potential temperature  (K)
      rhoa (:) = rhoa_data (i)    ! air density (kg/m^3)
      uatm (:) = uatm_data (i)    ! wind velocity components (m/s)
      vatm (:) = vatm_data (i)    
      wind (:) = wind_data (i)    ! wind speed (m/s)
      strax(:) = strax_data(i)    ! wind stress components (N/m^2)
      stray(:) = stray_data(i)   
      fsw  (:) = fsw_data  (i)    ! incoming shortwave radiation (W/m^2)
      swvdr(:) = swvdr_data(i)    ! sw down, visible, direct  (W/m^2)
      swvdf(:) = swvdf_data(i)    ! sw down, visible, diffuse (W/m^2)
      swidf(:) = swidr_data(i)    ! sw down, near IR, direct  (W/m^2)
      swidf(:) = swidf_data(i)    ! sw down, near IR, diffuse (W/m^2)
      Qa   (:) = Qa_data   (i)    ! specific humidity (kg/kg)
      flw  (:) = flw_data  (i)    ! incoming longwave radiation (W/m^2)
      frain(:) = frain_data(i)    ! rainfall rate (kg/m^2 s)
      fsnow(:) = fsnow_data(i)    ! snowfall rate (kg/m^2 s)

      if (trim(ocn_data_type) == 'default') return

      sst  (:) = sst_data  (i)    ! sea surface temperature
      sss  (:) = sss_data  (i)    ! sea surface salinity
      uocn (:) = uocn_data (i)    ! wind velocity components (m/s)
      vocn (:) = vocn_data (i) 


!cn we need trest here...
!cn it is normally read from forcing_nml
      if (restore_sst .or. restore_bgc) then
         if (trestore == 0) then
            trest = dt        ! use data instantaneously
         else
            trest = real(trestore,kind=dbl_kind) * secday ! seconds
         endif
      endif

!for debugging, for now
if (i==8760) then
write (nu_diag,*) flw
write (nu_diag,*) fsw
write (nu_diag,*) Tair
write (nu_diag,*) Qa
write (nu_diag,*) fsnow
write (nu_diag,*) frain
write (nu_diag,*) zlvl
write (nu_diag,*) potT
write (nu_diag,*) rhoa
write (nu_diag,*) uatm
write (nu_diag,*) vatm
write (nu_diag,*) wind
write (nu_diag,*) strax
write (nu_diag,*) stray
write (nu_diag,*) swvdr
write (nu_diag,*) swvdf
write (nu_diag,*) swidr
write (nu_diag,*) swidf
write (nu_diag,*) sst
write (nu_diag,*) uocn
write (nu_diag,*) vocn
endif

      end subroutine get_forcing

!=======================================================================

      subroutine atm_GOFS

      integer (kind=int_kind) :: &
         nu_navy, &     ! unit number
         nt             ! loop index

      real (kind=dbl_kind) :: &
         dlwsfc,  &     ! downwelling longwave (W/m2)
         dswsfc,  &     ! downwelling shortwave (W/m2)
         ltntht,  &     ! latent heat (W/m2)
         sensht,  &     ! sensible heat (W/m2)
         temp2m,  &     ! 2m air temperature (K)
         pottmp,  &     ! potential temperature (K) (=2m air temperature)
         spechum ,&     ! specific humidity (kg/kg)
         precip         ! precipitation (kg/m2/s)

      character (char_len_long) string1
      character (char_len_long) filename

      nu_navy = 12
      filename = trim(data_dir)//'/cfsv2_2015_220_70_01hr.ascii'

      write (nu_diag,*) 'Reading ',filename

      open (nu_navy, file=filename, form='formatted')
      read (nu_navy, *) string1 ! headers
      read (nu_navy, *) string1 ! units

      ntime = 8760 ! one year
      do nt = 1, ntime
         read (nu_navy, '(6(f10.5,1x),2(f10.8,1x))') &
         dlwsfc, dswsfc, ltntht, sensht, temp2m, pottmp, spechum, precip

           flw_data(nt) = dlwsfc
           fsw_data(nt) = dswsfc
          Tair_data(nt) = temp2m
          potT_data(nt) = pottmp
            Qa_data(nt) = spechum
         fsnow_data(nt) = precip
      enddo

      close (nu_navy)

!      write(nu_diag,*), 'GOFS data', &
!         dlwsfc, dswsfc, ltntht, sensht, temp2m, pottmp, spechum, precip

      end subroutine atm_GOFS

!=======================================================================

    subroutine atm_ISPOL             ! current forcing year

      character (char_len_long) filename
      
      filename = &
          trim(data_dir)//'ISPOL_atm_forcing.ascii'
      
      !write (nu_diag,*) ' '
      !write (nu_diag,*) 'Atmospheric data file:'
      !write (nu_diag,*) trim(filename)

#if 0
!this is from nicoles subroutine ISPOL_data
! Defines atmospheric data fields for Antarctic Weddell sea location 

! authors: Nicole Jeffery, LANL
!
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_constants, only: c0, c1, c2, c4, p001, p01, p1, secday, &
          field_loc_center, field_type_scalar, p5
      use ice_flux, only: uatm, vatm, Tair, fsw,  Qa, qdp, rhoa, &
          frain, fsnow, flw
      use ice_grid, only:  tmask
      use ice_diagnostics, only: latpnt, lonpnt 
#ifdef ncdf
      use netcdf
#endif

!local parameters

      character (char_len_long) :: & 
         met_file,   &    ! netcdf filename
         fieldname        ! field name in netcdf file

      integer (kind=int_kind) :: &
         fid              ! file id for netCDF file 

      real (kind=dbl_kind):: &
         work             ! temporary variable
 
      real (kind=dbl_kind) :: &
          vmin, vmax

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: &
         status           ! status flag

      integer (kind=int_kind) :: &
         iblk             ! block index

      real (kind=dbl_kind) :: & ! used to determine specific humidity
         Temp               , & ! air temperature (K)
         rh                 , & ! relative humidity (%)
         Psat               , & ! saturation vapour pressure (hPa)
         ws                     ! saturation mixing ratio

      real (kind=dbl_kind), dimension(2), save :: &
         Tair_data_p      , &      ! air temperature (K) for interpolation
         Qa_data_p,  fsnow_data_p, &
         fsw_data_p, flw_data_p, &
         uatm_data_p, vatm_data_p
         
      real (kind=dbl_kind), parameter :: & ! coefficients for Hyland-Wexler Qa 
         ps1 = 0.58002206e4_dbl_kind,    & ! (K) 
         ps2 = 1.3914993_dbl_kind,       & !
         ps3 = 0.48640239e-1_dbl_kind,   & ! (K^-1) 
         ps4 = 0.41764768e-4_dbl_kind,   & ! (K^-2)
         ps5 = 0.14452093e-7_dbl_kind,   & ! (K^-3)
         ps6 = 6.5459673_dbl_kind,       & !
         ws1 = 621.97_dbl_kind,          & ! for saturation mixing ratio 
         Pair = 1020._dbl_kind,          & ! Sea level pressure (hPa) 
         lapse_rate = 0.0065_dbl_kind      ! (K/m) lapse rate over sea level
    
      ! for interpolation of hourly data                
      integer (kind=int_kind) :: &
          i, j, k     , &
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          recnum4X    , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          sec_day        !  fix time to noon

       real (kind=dbl_kind) :: &
         hour_angle, &
         solar_time, &
         declin    , &
         cosZ      , &
         year_day  , &
         e, d      , &
         sw0       , &
         deg2rad   , &
         fsw_pnt   , &
         sumsw0    , &
         Qa_pnt                

      real (kind=dbl_kind) :: &
          sec1hr              ! number of seconds in 1 hour

      logical (kind=log_kind) :: readm, read1
                  
      diag = .false.   ! write diagnostic information 
   
#ifdef ncdf 
      if (trim(atm_data_format) == 'nc') then     ! read nc file
      
     !-------------------------------------------------------------------
     ! data from NCEP_DOE Reanalysis 2 and Bareiss et al 2008
     ! daily data located at the end of the 24-hour period. 
     !-------------------------------------------------------------------

      dataloc = 2                          ! data located at end of interval
      sec1hr = secday                      ! seconds in day
      maxrec = 366                         ! 

      ! current record number
      recnum = int(yday)   

      ! Compute record numbers for surrounding data (2 on each side)
      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!     ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ixp = -99
      call interp_coeff (recnum, recslot, sec1hr, dataloc)

      read1 = .false.
      if (istep==1 .or. oldrecnum .ne. recnum) read1 = .true.
      
      ! Daily 2m Air temperature 1991
                                                
        met_file = atm_file 
        fieldname='Tair' 

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, Tair_data_p, &
                    field_loc_center, field_type_scalar)

        Tair(:,:,:) =  c1intp * Tair_data_p(1) &
                       + c2intp * Tair_data_p(2) &
                     - lapse_rate*8.0_dbl_kind

        fieldname='Qa'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, Qa_data_p, &
                    field_loc_center, field_type_scalar)

        Qa_pnt= c1intp * Qa_data_p(1) &
                          + c2intp * Qa_data_p(2) 
        Qa(:,:,:) = Qa_pnt

        fieldname='uatm'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, uatm_data_p, &
                    field_loc_center, field_type_scalar)

        uatm(:,:,:) =  c1intp * uatm_data_p(1) &
                          + c2intp * uatm_data_p(2) 

        fieldname='vatm'
 
        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, vatm_data_p, &
                    field_loc_center, field_type_scalar)

        vatm(:,:,:) =  c1intp * vatm_data_p(1) &
                          + c2intp * vatm_data_p(2)
 
        fieldname='fsnow'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, fsnow_data_p, &
                    field_loc_center, field_type_scalar)

        fsnow(:,:,:) =  (c1intp * fsnow_data_p(1) + &
                          + c2intp * fsnow_data_p(2)) 
      
        !-----------------------------
        !fsw and flw are every 6 hours
        !------------------------------
        dataloc = 2                          ! data located at end of interval
        sec1hr = secday/c4                   ! seconds in 6 hours
        maxrec = 1464                        ! 366*4

      ! current record number
        recnum4X = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec1hr)   

      ! Compute record numbers for surrounding data (2 on each side)
      ixm = mod(recnum4X+maxrec-2,maxrec) + 1
      ixx = mod(recnum4X-1,       maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ixp = -99
      call interp_coeff (recnum4X, recslot, sec1hr, dataloc)

      read1 = .false.
      if (istep==1 .or. oldrecnum4X .ne. recnum4X) read1 = .true.

        fieldname='fsw'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, fsw_data_p, &
                    field_loc_center, field_type_scalar)

        fsw(:,:,:) =  c1intp * fsw_data_p(1) &
                          + c2intp * fsw_data_p(2)

        fieldname='flw' 

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, flw_data_p, &
                    field_loc_center, field_type_scalar)

        flw(:,:,:) =  c1intp * flw_data_p(1) &
                          + c2intp * flw_data_p(2) 
     endif  !nc
#else      
    
      uatm(:,:,:) = c0              !wind velocity (m/s)
      vatm(:,:,:) = c0
      fsw(:,:,:)  = c0 
      fsnow (:,:,:) = c0          

#endif

      !flw   given cldf and Tair  calculated in prepare_forcing

      !-----------------------------
      ! fixed data
      ! May not be needed
      !-----------------------------
        rhoa (:,:,:) = 1.3_dbl_kind ! air density (kg/m^3)
        cldf(:,:,:) =  c1  !0.25_dbl_kind ! cloud fraction
        frain(:,:,:) = c0            ! this is available in hourlymet_rh file
  
      ! Save record number for next time step
      oldrecnum = recnum
      oldrecnum4X = recnum4X


#endif

      
    end subroutine atm_ISPOL
    
!=======================================================================


      subroutine prepare_forcing (Tair,     fsw,      &    
                                  cldf,     flw,      &
                                  frain,    fsnow,    &
                                  Qa,       rhoa,     &
                                  uatm,     vatm,     &
                                  strax,    stray,    &
                                  zlvl,     wind,     &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  potT)

      use icepack_constants, only: c0, c1, c10, secday, Tffresh
 
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
         flw     , & ! incoming longwave radiation (W/m^2)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         potT        ! air potential temperature  (K)

      ! local variables

      integer (kind=int_kind) :: &
         nt

      real (kind=dbl_kind) :: workx, worky, &
         precip_factor, zlvl0

      zlvl0 = c10 ! default

      !-----------------------------------------------------------------
      ! convert precipitation units to kg/m^2 s
      !-----------------------------------------------------------------
      if (trim(precip_units) == 'mm_per_month') then
         precip_factor = 12._dbl_kind/(secday*days_per_year) 
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

         if (trim(atm_data_type) == 'GOFS') then
            ! precip is in kg/m^2/s
            zlvl0 = c10
            ! downward longwave as in Parkinson and Washington (1979)
!            call longwave_parkinson_washington(Tair(nt), cldf(nt), flw(nt))
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
               wind(nt) = sqrt(uatm(nt)**2 + vatm(nt)**2)
         ! else  ! strax, stray, wind are read from files
         endif                   ! calc_strair

      enddo ! ntime

      end subroutine prepare_forcing

!=======================================================================
      subroutine read_clim_data (readflag, recd, ixm, ixx, ixp, &
                                 data_file, field_data)


! Read data needed for interpolation, as in read_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.

#if 0
      use ice_diagnostics, only: check_step
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite
#endif
      logical (kind=log_kind),intent(in) :: readflag

      integer (kind=int_kind), intent(in) :: &
        recd            , & ! baseline record number
        ixm,ixx,ixp         ! record numbers of 3 data values
                            ! relative to recd

      character (char_len_long), intent(in) ::  data_file

!cn       integer (kind=int_kind), intent(in) :: &
!cn            field_loc, &      ! location of field on staggered grid
!cn            field_type        ! type of field (scalar, vector, angle)

      !cn real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
      real (kind=dbl_kind), dimension(nx,2), &
        intent(out) :: &
        field_data         ! 2 values needed for interpolation

      ! local variables
      integer (kind=int_kind) :: &
        nbits          , & ! = 32 for single precision, 64 for double
        nrec           , & ! record number to read
        arg            , & ! value of time argument in field_data
        i

#if 0
      call ice_timer_start(timer_readwrite)  ! reading/writing

      nbits = 64                ! double precision data
      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug)) &
        write(nu_diag,*) '  ', trim(data_file)
#endif
      if (readflag) then

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         !cn call ice_open (nu_forcing, data_file, nbits)
         open(nu_forcing,file=data_file,form='unformatted')
!cn this is probably not going to work? We need to get rec straight...
!         elseif (atype == 'rda4') then
!            allocate(work_gr(nx_global,ny_global))
!            read(nu,rec=nrec) work_gr
!            work_g1 = work_gr
!            deallocate(work_gr)
!cn in restart we do:
!            read(nu_forcing) (field_data(i,arg),i=1,nx) !cn





         arg = 0
         if (ixm /= -99) then
            arg = 1
            nrec = recd + ixm
            !cn call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
            !cn               'rda8', dbug, field_loc, field_type)
            read(nu_forcing,rec=nrec) field_data(:,arg) !cn
         endif

         arg = arg + 1
         nrec = recd + ixx
         !cn call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
         !cn               'rda8', dbug, field_loc, field_type)
         read(nu_forcing,rec=nrec) field_data(:,arg) !cn
            
         if (ixp /= -99) then
            arg = arg + 1
            nrec = recd + ixp
            !cn call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
            !cn               'rda8', dbug, field_loc, field_type)
            read(nu_forcing,rec=nrec) field_data(:,arg) !cn
         endif

#if 0
         if (my_task == master_task) close (nu_forcing)
#endif
      endif                     ! readflag

#if 0

      call ice_timer_stop(timer_readwrite)  ! reading/writing

#endif
      end subroutine read_clim_data

!=======================================================================

      subroutine interpolate_data (field_data, field)

! Linear interpolation

! author: Elizabeth C. Hunke, LANL

      !cn use ice_domain, only: nblocks

      !cn real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
      real (kind=dbl_kind), dimension(nx,2), &
        intent(in) :: &
        field_data    ! 2 values used for interpolation

      !cn real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
      real (kind=dbl_kind), dimension(nx), &
        intent(out) :: &
        field         ! interpolated field
      ! local variables

      integer (kind=int_kind) :: i
      do i = 1, nx
        field(i) = c1intp * field_data(i,1) &
            + c2intp * field_data(i,2)
      enddo
#if 0
      ! local variables

      integer (kind=int_kind) :: i,j, iblk

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            field(i,j,iblk) = c1intp * field_data(i,j,1,iblk) &
                            + c2intp * field_data(i,j,2,iblk)
         enddo
         enddo
      enddo
#endif
      end subroutine interpolate_data

!=======================================================================
      subroutine interp_coeff (recnum, recslot, secint, dataloc)

! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.

      use icepack_drv_constants, only: c1, p5, secday

      integer (kind=int_kind), intent(in) :: &
          recnum      , & ! record number for current data value
          recslot     , & ! spline slot for current record
          dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

      real (kind=dbl_kind), intent(in) :: &
          secint                    ! seconds in data interval

      ! local variables

      real (kind=dbl_kind) :: &
          secyr            ! seconds in a year

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2       , & ! seconds elapsed at data points
          rcnum            ! recnum => dbl_kind

      secyr = dayyr * secday         ! seconds in a year
      tt = mod(ftime,secyr)

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

      subroutine interp_coeff_monthly (recslot)
! Compute coefficients for interpolating monthly data to current time step.

      use icepack_drv_constants, only: c1, secday

      integer (kind=int_kind), intent(in) :: &
          recslot         ! slot (1 or 2) for current record

      ! local variables

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2           ! seconds elapsed at month midpoint

      real (kind=dbl_kind) :: &
          daymid(0:13)     ! month mid-points

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = 14._dbl_kind - daymo(12)  ! Dec 15, 0 sec

      ! make time cyclic
      tt = mod(ftime/secday,dayyr)

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

      subroutine read_data_point (flag, recd, yr, ixm, ixx, ixp, &
                            maxrec, data_file, fieldname, field_data, &
                            field_loc, field_type)
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ixp value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.
!
      use icepack_drv_constants, only: c0
      !cn use icepack_drv_diagnostics, only: check_step
      !cn use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite

      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                , & ! baseline record number
         yr                  , & ! year of forcing data
         ixm, ixx, ixp       , & ! record numbers of 3 data values
                                 ! relative to recd
         maxrec                  ! maximum record value

      character (char_len_long), intent(in) :: &
         data_file               ! data file to be read

      character (char_len), intent(in) :: &
         fieldname               ! field name in netCDF file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(2), &
         intent(out) :: &
         field_data              ! 2 values needed for interpolation
#if 0

#ifdef ncdf 
      integer (kind=int_kind) :: &
         nrec             , & ! record number to read
         n2, n4           , & ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
         arg              , & ! value of time argument in field_data
         fid                  ! file id for netCDF routines

      call ice_timer_start(timer_readwrite)  ! reading/writing

      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= -99) then
         ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  !call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                    ! call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1

      ! write(nu_diag,*) '!! read_data_nc !!!', trim(data_file)
      ! write(nu_diag,*) 'istep  ', istep
      ! write(nu_diag,*) 'fyear_final  ', fyear_final
      ! write(nu_diag,*) 'fyear_init  ', fyear_init
      ! write(nu_diag,*) 'ixm, ixx, ixp  ', ixm, ixx, ixp
      ! write(nu_diag,*) 'maxrec ', maxrec
      ! write(nu_diag,*) 'fieldname  ', fieldname
 
            call ice_open_nc (data_file, fid)

            arg = 1
            nrec = recd + n2

            call ice_read_nc & 
                 (fid, nrec, fieldname, field_data(arg), dbug, &
                  field_loc, field_type)

            !if (ixx==1) call ice_close_nc(fid)
            call ice_close_nc(fid)
         endif                  ! ixm ne -99

        ! always read ixx data from data file for current year
        ! call file_year (data_file, yr)
         call ice_open_nc (data_file, fid)

         arg = arg + 1
         nrec = recd + ixx

         call ice_read_nc & 
              (fid, nrec, fieldname, field_data(arg), dbug, &
               field_loc, field_type)

         if (ixp /= -99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  call ice_close_nc(fid)
                  !call file_year (data_file, yr+1)
                  call ice_open_nc (data_file, fid)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     call ice_close_nc(fid)
                    ! call file_year (data_file, fyear_init)
                     call ice_open_nc (data_file, fid)

                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4

            call ice_read_nc & 
                 (fid, nrec, fieldname, field_data(arg), dbug, &
                  field_loc, field_type)
         endif                  ! ixp /= -99

         call ice_close_nc(fid)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

#else
      field_data = c0 ! to satisfy intent(out) attribute
#endif

#endif
      end subroutine read_data_point

!=======================================================================

      end module icepack_drv_forcing

!=======================================================================
