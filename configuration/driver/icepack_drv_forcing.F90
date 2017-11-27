!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
      module icepack_drv_forcing

      use icepack_drv_kinds
      use icepack_drv_domain_size, only: ncat, nx
      use icepack_drv_calendar, only: time, nyr, dayyr, mday, month, &
         daymo, daycal, &
         dt, yday , days_per_year
      use icepack_drv_constants, only: nu_diag, nu_forcing, secday
      use icepack_drv_parameters, only: calc_strair

      implicit none
      private
      public :: init_forcing, get_forcing, interp_coeff, interp_coeff_monthly

      save

      integer (kind=int_kind), parameter :: &
         ntime = 8760        ! number of data points in time

      integer (kind=int_kind), public :: &
         ycycle          , & ! number of years in forcing cycle
         fyear_init      , & ! first year of data in forcing cycle
         fyear           , & ! current year in forcing cycle
         fyear_final         ! last year in cycle

      real (kind=dbl_kind), dimension(ntime) :: &
            fsw_data, & ! field values at temporal data points
           cldf_data, &
          fsnow_data, &
           Tair_data, &
           uatm_data, &
           vatm_data, &
           wind_data, &
          strax_data, &
          stray_data, &
           rhum_data, &
             Qa_data, &
           rhoa_data, &
           potT_data, &
            flw_data, &
            qdp_data, &
            sst_data, &
            sss_data, & 
           uocn_data, &
           vocn_data, &
         sublim_data, &
          frain_data, &
          swvdr_data, &
          swvdf_data, &
          swidr_data, &
          swidf_data, &
           zlvl_data, &
          hmix_data

      real (kind=dbl_kind), dimension(nx) :: &
          sst_temp

      character(char_len), public :: & 
         atm_data_format, & ! 'bin'=binary or 'nc'=netcdf
         ocn_data_format, & ! 'bin'=binary or 'nc'=netcdf
         bgc_data_format, & ! 'bin'=binary or 'nc'=netcdf
         atm_data_type,   & ! 'default', 'clim', 'CFS'
         ocn_data_type,   & ! 'default'
         bgc_data_type,   & ! 'default'
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
         dbug             ! prints debugging output if true

      logical (kind=log_kind), public :: &
         restore_ocn               ! restore sst if true

      real (kind=dbl_kind), public :: & 
         trest, &                    ! restoring time scale (sec)
         trestore                    ! restoring time scale (days)

!=======================================================================

      contains

!=======================================================================

      subroutine init_forcing

! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the forcing data.

      use icepack_constants, only: c0
      use icepack_drv_flux, only: zlvl, Tair, potT, rhoa, uatm, vatm, wind, &
         strax, stray, fsw, swvdr, swvdf, swidr, swidf, Qa, flw, frain, &
         fsnow, sst, sss, uocn, vocn, qdp

      integer (kind=int_kind) :: &
         i                ! index

      fyear       = fyear_init + mod(nyr-1,ycycle) ! current year
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      write (nu_diag,*) ' Initial forcing data year = ',fyear_init
      write (nu_diag,*) ' Final   forcing data year = ',fyear_final

    !-------------------------------------------------------------------
    ! Initialize forcing data to default values
    !-------------------------------------------------------------------
! maybe these should all be zero, and set defaults in particular source 
! routines to be clear which ones are defaults and which are being read in

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
           sst_data(:) = sst  (i)    ! sea surface temperature
           sss_data(:) = sss  (i)    ! sea surface salinity
          uocn_data(:) = uocn (i)    ! wind velocity components (m/s)
          vocn_data(:) = vocn (i)

          cldf_data(:) = c0     ! cloud fraction

      if (trim(atm_data_type(1:4)) == 'CFS') call atm_CFS
      if (trim(atm_data_type(1:4)) == 'clim') call atm_climatological

      if (trim(atm_data_type) == 'ISPOL') call atm_ISPOL
      if (trim(atm_data_type) == 'NICE') call atm_NICE

      !cn if (restore_ocn .or. restore_bgc) then
      if (restore_ocn) then
        if (trestore == 0) then
          trest = dt        ! use data instantaneously
        else
          trest = real(trestore,kind=dbl_kind) * secday ! seconds
        end if
      endif

      if (trim(ocn_data_type) == 'ISPOL') then
        call ocn_ISPOL
      endif


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
! We will probably need to send in the time and work out what the data
! time slice is, instead of sending in the timestep.  This currently assumes
! the time step and the data both start Jan 1.

      use icepack_constants, only: c0, c1
      use icepack_drv_flux, only: zlvl, Tair, potT, rhoa, uatm, vatm, wind, &
         strax, stray, fsw, swvdr, swvdf, swidr, swidf, Qa, flw, frain, &
         fsnow, sst, sss, uocn, vocn, qdp, hmix

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

      integer (kind=int_kind) :: &  !cn ispol stuff
          recnum, dataloc, maxrec

      real (kind=dbl_kind) :: &
          sec6hr
  
      !values for default case
      i = mod(timestep-1,ntime)+1 ! repeat forcing cycle
      mlast = i
      mnext = mlast
      c1intp = c1
      c2intp = c0

      if (trim(atm_data_type) == 'CFS') then
         ! calculate data index corresponding to current timestep
         i = mod(timestep-1,ntime)+1 ! repeat forcing cycle
         mlast = i
         mnext = mlast
         c1intp = c1
         c2intp = c0
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
        
      elseif (trim(atm_data_type) == 'ISPOL' .or. &
          trim(atm_data_type) == 'NICE') then

        !we have a years worth of data, some in days, some in 6hr or quarter day
        !it also looks like the data Nicole starts on June 17, not Jan 1 ?????
        !also need to repeat as above

        dataloc = 2                          ! data located at end of interval
        maxrec = 366  
        recslot = 2
        recnum = int(yday)   
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff ( recnum, recslot, secday, dataloc, c1intp, c2intp)

        Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
        Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
        uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
        vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
        fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)
        
        sec6hr = secday/4;                      ! seconds in 6 hours
        maxrec = 1464  
        recnum = int(yday*4)
        mlast = mod(recnum+maxrec-2,maxrec) + 1
        mnext = mod(recnum-1,       maxrec) + 1
        call interp_coeff (recnum, recslot, sec6hr, dataloc, c1intp, c2intp)

        fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)
        flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)

      endif


!cn this is called from get_forcing_ocn in cice...
      if(trim(ocn_data_type) == 'ISPOL') then
        
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
        sst_temp(:) = c1intp *   sst_data(mlast) + c2intp *   sst_data(mnext)
        sss(:) = c1intp *   sss_data(mlast) + c2intp *   sss_data(mnext)
        uocn(:) = c1intp *  uocn_data(mlast) + c2intp *  uocn_data(mnext)
        vocn(:) = c1intp *  vocn_data(mlast) + c2intp *  vocn_data(mnext)
        do i = 1, nx
          sss (i) = max (sss(i), c0) 
          hmix(i) = max(hmix(i), c0)           
        end do

        call finish_ocn_forcing(sst_temp)

      elseif (trim(ocn_data_type) == 'NICE') then
        
!cn the nice stuff seems to be more complicated than ispol....        
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
        sst_temp(:) = c1intp *   sst_data(mlast) + c2intp *   sst_data(mnext)
        sss(:) = c1intp *   sss_data(mlast) + c2intp *   sss_data(mnext)
        uocn(:) = c1intp *  uocn_data(mlast) + c2intp *  uocn_data(mnext)
        vocn(:) = c1intp *  vocn_data(mlast) + c2intp *  vocn_data(mnext)
        do i = 1, nx
          sss (i) = max (sss(i), c0) 
          hmix(i) = max(hmix(i), c0)           
        end do

        call finish_ocn_forcing(sst_temp)

      else
        

      endif

      ! fill all grid boxes with the same forcing data
      flw  (:) = c1intp *   flw_data(mlast) + c2intp *   flw_data(mnext)
      Tair (:) = c1intp *  Tair_data(mlast) + c2intp *  Tair_data(mnext)
      potT (:) = c1intp *  potT_data(mlast) + c2intp *  potT_data(mnext)
      rhoa (:) = c1intp *  rhoa_data(mlast) + c2intp *  rhoa_data(mnext)
      uatm (:) = c1intp *  uatm_data(mlast) + c2intp *  uatm_data(mnext)
      vatm (:) = c1intp *  vatm_data(mlast) + c2intp *  vatm_data(mnext)
      wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
      strax(:) = c1intp * strax_data(mlast) + c2intp * strax_data(mnext)
      stray(:) = c1intp * stray_data(mlast) + c2intp * stray_data(mnext)
      wind (:) = c1intp *  wind_data(mlast) + c2intp *  wind_data(mnext)
      fsw  (:) = c1intp *   fsw_data(mlast) + c2intp *   fsw_data(mnext)
      swvdr(:) = c1intp * swvdr_data(mlast) + c2intp * swvdr_data(mnext)
      swvdf(:) = c1intp * swvdf_data(mlast) + c2intp * swvdf_data(mnext)
      swidr(:) = c1intp * swidr_data(mlast) + c2intp * swidr_data(mnext)
      swidf(:) = c1intp * swidf_data(mlast) + c2intp * swidf_data(mnext)
      Qa   (:) = c1intp *    Qa_data(mlast) + c2intp *    Qa_data(mnext)
      frain(:) = c1intp * frain_data(mlast) + c2intp * frain_data(mnext)
      fsnow(:) = c1intp * fsnow_data(mlast) + c2intp * fsnow_data(mnext)

        qdp(:) = c1intp *   qdp_data(mlast) + c2intp *   qdp_data(mnext)

! for debugging
!if (timestep==8760.or.timestep==8761) then
if (0==1) then ! off
write (nu_diag,*) 'timestep',timestep, mlast, mnext
write (nu_diag,*) 'index',mlast,mnext
write (nu_diag,*) 'flw',flw
write (nu_diag,*) 'fsw',fsw
write (nu_diag,*) 'Tair',Tair
write (nu_diag,*) 'Qa',Qa
write (nu_diag,*) 'fsnow',fsnow
write (nu_diag,*) 'frain',frain
write (nu_diag,*) 'potT',potT
write (nu_diag,*) 'rhoa',rhoa
write (nu_diag,*) 'uatm',uatm
write (nu_diag,*) 'vatm',vatm
write (nu_diag,*) 'wind',wind
write (nu_diag,*) 'strax',strax
write (nu_diag,*) 'stray',stray
write (nu_diag,*) 'swvdr',swvdr
write (nu_diag,*) 'swvdf',swvdf
write (nu_diag,*) 'swidr',swidr
write (nu_diag,*) 'swidf',swidf
write (nu_diag,*) 'sst',sst
write (nu_diag,*) 'sss',sss
write (nu_diag,*) 'uocn',uocn
write (nu_diag,*) 'vocn',vocn
write (nu_diag,*) 'qdp',qdp
endif

      end subroutine get_forcing

!=======================================================================

      subroutine atm_climatological

      use icepack_constants, only: c0, c1, c2, c100, qqqice, TTTice, & 
         rhos, Tffresh 

      real (kind=dbl_kind), dimension(12) :: &
            fsw_clim, & ! field values at temporal data points
            flw_clim, &
           Tair_clim, &
           wind_clim, &
           rhum_clim, &
          fsnow_clim

      ! Ice station meteorology from Lindsay (1998, J. Climate), Table 1, p. 325
      ! zlvl = c2 ! 2-m temperatures and wind speed

      data  fsw_clim /  0.0,   1.2,  31.5, 146.0, 263.3, 307.9, &
                      230.6, 134.7,  44.2,   2.6,   0.0,   0.0  /
      data  flw_clim /164.0, 160.5, 164.1, 188.1, 245.2, 291.2, &
                      303.9, 297.0, 263.8, 210.9, 177.0, 166.0  /
      data Tair_clim /-31.4, -32.8, -31.6, -24.1, -11.0,  -1.8, &
                       -0.1,  -1.4,  -8.0, -19.5, -27.6, -31.1  /
      data rhum_clim / 78.7,  78.4,  79.6,  82.1,  86.5,  91.7, &
                       95.1,  94.3,  90.7,  83.8,  80.1,  78.7  /
      data wind_clim /  4.4,   4.0,   4.0,   3.9,   3.9,   4.2, &
                        4.1,   4.2,   4.5,   4.2,   3.9,   4.0  /
!      data  shf_clim /  9.9,   8.4,   6.6,   0.1,  -5.8,  -1.6, &
!                        2.2,   1.2,   0.5,   2.0,   5.6,   7.0  /
!      data  lhf_clim /  1.3,   1.1,   1.1,   0.0,  -5.9, -10.3, &
!                       -6.5,  -6.7,  -3.9,  -0.1,   1.0,   1.1  /

      ! Semtner (1976, JPO) snowfall spec., p. 383 in m/s snow volume (.4 m/yr)
      data fsnow_clim/ 3.17e-9, 3.17e-9, 3.17e-9, 3.17e-9, 1.90e-8,    0.0, &
                           0.0, 1.63e-8, 4.89e-8, 4.89e-8, 3.17e-9, 3.17e-9 /

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
        qdp_data(1:12) = -6.0 ! 2 W/m2 from deep + 4 W/m2 counteracting larger
                              ! SH+LH with bulk transfer than in MU 71

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
         pottmp,  &     ! potential temperature (K) (=2m air temperature)
         spechum ,&     ! specific humidity (kg/kg)
         precip         ! precipitation (kg/m2/s)

      character (char_len_long) string1
      character (char_len_long) filename

      filename = trim(data_dir)//'/CFS/cfsv2_2015_220_70_01hr.txt'

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

!      write(nu_diag,*), 'CFS data', &
!         dswsfc, dlwsfc, windu10, windv10, temp2m, spechum, precip

      end subroutine atm_CFS

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

      ! this routine acts on the data fields prior to interpolation

      use icepack_constants, only: c0, c1, c2, c10, secday, Tffresh

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

      real (kind=dbl_kind), parameter :: &    
         lapse_rate = 0.0065_dbl_kind      ! (K/m) lapse rate over sea level

      real (kind=dbl_kind) :: workx, worky, &
         precip_factor, zlvl0

      zlvl0 = c10 ! default

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

      use icepack_constants, only: c1, secday

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
                               c1intp, c2intp)

! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.

      use icepack_constants, only: c1, p5, secday

      integer (kind=int_kind), intent(in) :: &
          recnum      , & ! record number for current data value
          recslot     , & ! spline slot for current record
          dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

      real (kind=dbl_kind), intent(in) :: &
          secint                    ! seconds in data interval

      real (kind=dbl_kind), intent(inout) :: &
         c1intp, c2intp   ! interpolation coefficients

      ! local variables

      real (kind=dbl_kind) :: &
          secyr            ! seconds in a year

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2       , & ! seconds elapsed at data points
          rcnum            ! recnum => dbl_kind

      secyr = dayyr * secday         ! seconds in a year
      tt = mod(time,secyr)

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

      integer (kind=int_kind) :: &
         nu_ispol,&     ! unit number
         i

      real (kind=dbl_kind), dimension(366) :: &
          tair, &
          qa, &
          uatm, &
          vatm, &
          fsnow, &
          aday
      real (kind=dbl_kind), dimension(1464) :: &
          fsw, &
          flw, &
          atime

      character (char_len_long) filename
      
      filename = &
          trim(data_dir)//'ISPOL_atm_forcing.txt'

      write (nu_diag,*) 'Reading ',filename

      open (nu_ispol, file=filename, form='formatted')

      read(nu_ispol,*) tair
      read(nu_ispol,*) qa
      read(nu_ispol,*) fsw
      read(nu_ispol,*) flw
      read(nu_ispol,*) uatm
      read(nu_ispol,*) vatm
      read(nu_ispol,*) fsnow
      read(nu_ispol,*) aday
      read(nu_ispol,*) atime

      do i = 1, 366 !daily
        Tair_data(i) = tair(i)
        Qa_data(i) = qa(i)
        uatm_data(i) = uatm(i)
        vatm_data(i) = vatm(i)
        fsnow_data(i) = fsnow(i)
      end do
      do i = 1, 1464 ! 6hr, 1464/4=366 days
        fsw_data(i) = fsw(i)
        flw_data(i) = flw(i)
      end do

      !write(*,*) tair
      !write(*,*) qa
      !write(*,*) fsw
      !write(*,*) flw
      !write(*,*) uatm
      !write(*,*) vatm
      !write(*,*) fsnow
      !write(*,*) aday
      !write(*,*) atime
      
      close(nu_ispol)

      !write (nu_diag,*) ' '
      !write (nu_diag,*) 'Atmospheric data file:'
      !write (nu_diag,*) trim(filename)
      
    end subroutine atm_ISPOL

!=======================================================================

    subroutine atm_NICE

      integer (kind=int_kind) :: &
         nu_nice,&     ! unit number
         i

      real (kind=dbl_kind), dimension(366) :: &
          tair, &
          qa, &
          uatm, &
          vatm, &
          fsnow, &
          aday
      real (kind=dbl_kind), dimension(1464) :: &
          fsw, &
          flw, &
          atime

      character (char_len_long) filename
      
      filename = &
          trim(data_dir)//'NICE_atm_forcing.txt'

      write (nu_diag,*) 'Reading ',filename

      open (nu_nice, file=filename, form='formatted')

      read(nu_nice,*) tair
      read(nu_nice,*) qa
      read(nu_nice,*) fsw
      read(nu_nice,*) flw
      read(nu_nice,*) uatm
      read(nu_nice,*) vatm
      read(nu_nice,*) fsnow
      read(nu_nice,*) aday
      read(nu_nice,*) atime
      do i = 1, 366
        Tair_data(i) = tair(i)
        Qa_data(i) = qa(i)
        uatm_data(i) = uatm(i)
        vatm_data(i) = vatm(i)
        fsnow_data(i) = fsnow(i)
      end do
      do i = 1, 1464
        fsw_data(i) = fsw(i)
        flw_data(i) = flw(i)
      end do

      !write(*,*) tair
      !write(*,*) qa
      !write(*,*) fsw
      !write(*,*) flw
      !write(*,*) uatm
      !write(*,*) vatm
      !write(*,*) fsnow
      !write(*,*) aday
      !write(*,*) atime
      
      close(nu_nice)

!cn there is probably more to do here, see below...


#if 0
      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year

      atm_file = &
          trim(atm_data_dir)//'NICE_atm_forcing.nc'


      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Atmospheric data file:'
         write (nu_diag,*) trim(atm_file)
      endif                     ! master_task


!from nicoles      subroutine NICE_data(dt)

! Defines ocean data fields for NICE-2015 Arctic location (Nansen Basin)

! authors: Nicole Jeffery, LANL

! Restore sst if desired. sst is updated with surface fluxes in ice_ocean.F.

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, p5, secday, &
          field_loc_center, field_type_scalar, p5
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_domain_size, only: max_blocks
      use ice_flux, only: sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty, &
            qdp, hmix
      use ice_restart_shared, only: restart
      use ice_grid, only: hm, tmask, umask
      use ice_colpkg, only: colpkg_liquidus_temperature
      use ice_diagnostics, only: latpnt, lonpnt
#ifdef ncdf
      use netcdf
#endif

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

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

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: &
         status           ! status flag

      integer (kind=int_kind) :: &
         iblk             ! block index

      ! for interpolation of hourly data
      integer (kind=int_kind) :: &
          i, j, k, n  , &
          ixm,ixx,ixp , & ! record numbers
          recnum      , & ! record number
          recnum4X    , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          sec_day        !  fix time to noon

      real (kind=dbl_kind) :: &
          sec1hr              ! number of seconds in 1 hour

      logical (kind=log_kind) :: readm, read1

      diag = .false.   ! write diagnostic information
     !-------------------------------------------------------------------
     ! NICE_2015/oceanmixed_daily.nc
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

      do n = nfld, 1, -1
        !$OMP PARALLEL DO PRIVATE(iblk,i,j)
        do iblk = 1, nblocks
        ! use sst_data arrays as temporary work space until n=1
        if (ixm /= -99) then
          sst_data(:,:,1,iblk) = ocn_frc_d(:,:,iblk,n,ixm)
          sst_data(:,:,2,iblk) = ocn_frc_d(:,:,iblk,n,recnum)
        else
          sst_data(:,:,1,iblk) = ocn_frc_d(:,:,iblk,n,recnum)
          sst_data(:,:,2,iblk) = ocn_frc_d(:,:,iblk,n,ixp)
        endif
        enddo
        !$OMP END PARALLEL DO

        call interpolate_data (sst_data,work1)
        ! masking by hm is necessary due to NaNs in the data file
        do j = 1, ny_block 
          do i = 1, nx_block 
            if (n == 2) sss    (i,j,:) = c0
            if (n == 3) hmix   (i,j,:) = c0
            if (n == 4) uocn   (i,j,:) = c0
            if (n == 5) vocn   (i,j,:) = c0
            if (n == 6) ss_tltx(i,j,:) = c0
            if (n == 7) ss_tlty(i,j,:) = c0
            if (n == 8) qdp    (i,j,:) = c0
            do iblk = 1, nblocks
              if (hm(i,j,iblk) == c1) then
                if (n == 2) sss    (i,j,iblk) = work1(i,j,iblk)
                if (n == 3) hmix   (i,j,iblk) = work1(i,j,iblk)
                if (n == 4) uocn   (i,j,iblk) = work1(i,j,iblk)
                if (n == 5) vocn   (i,j,iblk) = work1(i,j,iblk)
                if (n == 6) ss_tltx(i,j,iblk) = work1(i,j,iblk)
                if (n == 7) ss_tlty(i,j,iblk) = work1(i,j,iblk)
                if (n == 8) qdp    (i,j,iblk) = work1(i,j,iblk)
              endif
            enddo
          enddo
        enddo
      enddo

      do j = 1, ny_block 
         do i = 1, nx_block 
            sss (i,j,:) = max (sss(i,j,:), c0) 
            hmix(i,j,:) = max(hmix(i,j,:), c0) 
         enddo 
      enddo 

      call ocn_freezing_temperature

      if (restore_ocn) then
        do j = 1, ny_block 
         do i = 1, nx_block 
           sst(i,j,:) = sst(i,j,:) + (work1(i,j,:)-sst(i,j,:))*dt/trest 
         enddo 
        enddo 
!     else sst is only updated in ice_ocean.F
      endif

      ! initialize sst properly on first step
      if (istep1 <= 1 .and. .not. (restart)) then
        call interpolate_data (sst_data,sst)
        !$OMP PARALLEL DO PRIVATE(iblk,i,j)
        do iblk = 1, nblocks
         do j = 1, ny_block 
          do i = 1, nx_block 
            if (hm(i,j,iblk) == c1) then
              sst(i,j,iblk) =  max (sst(i,j,iblk), Tf(i,j,iblk)) 
            else
              sst(i,j,iblk) = c0
            endif
          enddo 
         enddo 
        enddo 
        !$OMP END PARALLEL DO
      endif

      if (dbug) then
         if (my_task == master_task)  &
               write (nu_diag,*) 'ocn_data_ncar'
           vmin = global_minval(Tf,distrb_info,tmask)
           vmax = global_maxval(Tf,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'Tf',vmin,vmax
           vmin = global_minval(sst,distrb_info,tmask)
           vmax = global_maxval(sst,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'sst',vmin,vmax
           vmin = global_minval(sss,distrb_info,tmask)
           vmax = global_maxval(sss,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'sss',vmin,vmax
           vmin = global_minval(hmix,distrb_info,tmask)
           vmax = global_maxval(hmix,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'hmix',vmin,vmax
           vmin = global_minval(uocn,distrb_info,umask)
           vmax = global_maxval(uocn,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'uocn',vmin,vmax
           vmin = global_minval(vocn,distrb_info,umask)
           vmax = global_maxval(vocn,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'vocn',vmin,vmax
           vmin = global_minval(ss_tltx,distrb_info,umask)
           vmax = global_maxval(ss_tltx,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'ss_tltx',vmin,vmax
           vmin = global_minval(ss_tlty,distrb_info,umask)
           vmax = global_maxval(ss_tlty,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'ss_tlty',vmin,vmax
           vmin = global_minval(qdp,distrb_info,tmask)
           vmax = global_maxval(qdp,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'qdp',vmin,vmax
      endif


#endif           

    end subroutine atm_NICE

!=======================================================================

    subroutine ocn_ISPOL
!cn cice fills ocn_frc_m here, I am not sure what that does...


      integer (kind=int_kind) :: &
         nu_ispol,&     ! unit number
         i

      real (kind=dbl_kind), dimension(12) :: &
          t, &  !probably temperature, Tf?
          s, &  !probably sss_data
          hblt, &  !probably hmix
          u, &  !probably uocn_data seems to be zeroed out anyway??
          v, &  !probably vocn_data seems to be zeroed out anyway??
          dhdx, &  !probably ss_tltx
          dhdy, &  !probably ss_tlty 
          qdp  !probably heat flux

      character (char_len_long) filename
      
      filename = &
          trim(data_dir)//'pop_frc.gx1v3.051202_but_hblt_from_010815_ispol.txt'

      write (nu_diag,*) 'Reading ',filename

      open (nu_ispol, file=filename, form='formatted')

      read(nu_ispol,*) t
      read(nu_ispol,*) s
      read(nu_ispol,*) hblt
      read(nu_ispol,*) u
      read(nu_ispol,*) v
      read(nu_ispol,*) dhdx
      read(nu_ispol,*) dhdy
      read(nu_ispol,*) qdp
      do i = 1, 12 ! monthly
        !t(i)
        sss_data(i) = s(i)
        hmix_data(i) = hblt(i)
        uocn_data(i) = u(i)
        vocn_data(i) = v(i)
        !dhdx(i)
        !dhdy(1)
        qdp_data(i) = qdp(i)
      end do

      !write(*,*) t
      !write(*,*) s
      !write(*,*) hblt
      !write(*,*) u
      !write(*,*) v
      !write(*,*) dhdx
      !write(*,*) dhdy
      !write(*,*) qdp
      
      !stop

    end subroutine ocn_ISPOL

!=======================================================================

      subroutine finish_ocn_forcing(sst_temp)

 ! Compute ocean freezing temperature Tf based on tfrz_option
 ! 'minus1p8'         Tf = -1.8 C (default)
 ! 'linear_salt'      Tf = -depressT * sss
 ! 'mushy'            Tf conforms with mushy layer thermo (ktherm=2)

      use icepack_therm_shared, only: icepack_sea_freezing_temperature
      use icepack_drv_flux, only: sss, Tf, sst

      real (kind=dbl_kind), dimension(nx), intent(in)  :: &
          sst_temp

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk           ! horizontal indices

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do i = 1, nx
        Tf(i) = icepack_sea_freezing_temperature(sss(i))
        if (restore_ocn) then
          sst(i) = sst(i) + (sst_temp(i)-sst(i))*dt/trest
        endif
      enddo
      !$OMP END PARALLEL DO

    end subroutine finish_ocn_forcing


!=======================================================================

      end module icepack_drv_forcing

!=======================================================================
