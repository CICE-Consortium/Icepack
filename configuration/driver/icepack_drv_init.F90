!=======================================================================

! parameter and variable initializations
!
! authors Elizabeth C. Hunke, LANL

      module icepack_drv_init

      use icepack_drv_kinds
      use icepack_drv_domain_size, only: nx
      use icepack_intfc, only: icepack_init_constants

      implicit none
      private
      public :: input_data, init_grid2, init_state

      character(len=char_len_long), public, save :: &
         ice_ic      ! method of ice cover initialization
                     ! 'default' or 'none' => conditions specified in code
                     ! restart = .true. overwrites default initial
                     !    condition using filename given by ice_ic

      real (kind=dbl_kind), dimension (nx), public, save :: &
         TLON   , & ! longitude of temp pts (radians)
         TLAT       ! latitude of temp pts (radians)

      logical (kind=log_kind), &
         dimension (nx), public, save :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

!=======================================================================

      contains

!=======================================================================

! Namelist variables, set to default values; may be altered
! at run time
!
! author Elizabeth C. Hunke, LANL

      subroutine input_data

      use icepack_drv_parameters, only: ustar_min, albicev, albicei, albsnowv, albsnowi
      use icepack_drv_parameters, only: ahmax, shortwave, albedo_type, R_ice, R_pnd
      use icepack_drv_parameters, only: R_snw, dT_mlt, rsnw_mlt
      use icepack_drv_parameters, only: kstrength, krdg_partic, krdg_redist, mu_rdg
      use icepack_drv_parameters, only: atmbndy, calc_strair, formdrag, highfreq, natmiter
      use icepack_drv_parameters, only: kitd, kcatbound, hs0, dpscale, frzpnd
      use icepack_drv_parameters, only: rfracmin, rfracmax, pndaspect, hs1, hp1
      use icepack_drv_parameters, only: ktherm, calc_Tsfc, conduct
      use icepack_drv_arrays_column, only: oceanmixed_ice
      use icepack_drv_constants, only: c0, c1, puny, ice_stdout, nu_diag, nu_diag_out, nu_nml
      use icepack_drv_diagnostics, only: diag_file, nx_names
      use icepack_drv_domain_size, only: nilyr, nslyr, max_ntrcr, ncat, n_aero
      use icepack_drv_calendar, only: year_init, istep0
      use icepack_drv_calendar, only: dumpfreq, diagfreq
      use icepack_drv_calendar, only: npt, dt, ndtd, days_per_year, use_leap_years
      use icepack_drv_restart_shared, only: restart, restart_dir, restart_file
      use icepack_drv_flux, only: update_ocn_f, l_mpond_fresh, cpl_bgc
      use icepack_drv_flux, only: default_season
      use icepack_drv_forcing, only: precip_units,    fyear_init,      ycycle
      use icepack_drv_forcing, only: atm_data_type,   ocn_data_type,   bgc_data_type
      use icepack_drv_forcing, only: atm_data_format, ocn_data_format, bgc_data_format
      use icepack_drv_forcing, only: data_dir,        dbug
      use icepack_drv_forcing, only: restore_ocn, trestore

      use icepack_drv_tracers, only: tr_iage, tr_FY, tr_lvl, tr_pond
      use icepack_drv_tracers, only: tr_pond_cesm, tr_pond_lvl, tr_pond_topo
      use icepack_drv_tracers, only: tr_aero, nt_Tsfc, nt_qice, nt_qsno, nt_sice
      use icepack_drv_tracers, only: nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd
      use icepack_drv_tracers, only: nt_aero, ntrcr
      use icepack_drv_parameters, only: a_rapid_mode, Rac_rapid_mode
      use icepack_drv_parameters, only: aspect_rapid_mode, dSdt_slow_mode
      use icepack_drv_parameters, only: phi_c_slow_mode, phi_i_mushy
      use icepack_drv_parameters, only: tfrz_option, kalg, fbot_xfer_type

      ! local variables

      character (32) :: &
         nml_filename = 'icepack_in' ! namelist input file name

      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        n,         & ! loop index
        diag_len     ! length of diag file

      character (len=char_len) :: diag_file_names

      character (len=6) :: chartmp
      character (len=32) :: str
      character (len=20) :: format_str

      logical :: exists

      real (kind=real_kind) :: rpcesm, rplvl, rptopo 
      real (kind=dbl_kind) :: Cf

      !-----------------------------------------------------------------
      ! Namelist variables.
      !-----------------------------------------------------------------

      namelist /setup_nml/ &
        days_per_year,  use_leap_years, year_init,       istep0,        &
        dt,             npt,            ndtd,                           &
        ice_ic,         restart,        restart_dir,     restart_file,  &
        dumpfreq,    &
        diagfreq,       diag_file,                      &
        cpl_bgc

      namelist /grid_nml/ &
        kcatbound

      namelist /thermo_nml/ &
        kitd,           ktherm,          conduct,                       &
        a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,             &
        dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy

      namelist /dynamics_nml/ &
        kstrength,      krdg_partic,    krdg_redist,    mu_rdg,         &
        Cf

      namelist /shortwave_nml/ &
        shortwave,      albedo_type,                                    &
        albicev,        albicei,         albsnowv,      albsnowi,       &
        ahmax,          R_ice,           R_pnd,         R_snw,          &
        dT_mlt,         rsnw_mlt,        kalg

      namelist /ponds_nml/ &
        hs0,            dpscale,         frzpnd,                        &
        rfracmin,       rfracmax,        pndaspect,     hs1,            &
        hp1

      namelist /forcing_nml/ &
        atmbndy,         calc_strair,     calc_Tsfc,       &
        update_ocn_f,    l_mpond_fresh,   ustar_min,       &
        fbot_xfer_type,  oceanmixed_ice,                   &
        formdrag,        highfreq,        natmiter,        &
        tfrz_option,     default_season,                   &
        precip_units,    fyear_init,      ycycle,          &
        atm_data_type,   ocn_data_type,   bgc_data_type,   &
        atm_data_format, ocn_data_format, bgc_data_format, &
        data_dir,        trestore,        restore_ocn

      namelist /tracer_nml/   &
        tr_iage,      &
        tr_FY,        &
        tr_lvl,       &
        tr_pond_cesm, &
        tr_pond_lvl,  &
        tr_pond_topo, &
        tr_aero

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      days_per_year = 365    ! number of days in a year
      use_leap_years= .false.! if true, use leap years (Feb 29)
      year_init = 0          ! initial year
      istep0 = 0             ! no. of steps taken in previous integrations,
                             ! real (dumped) or imagined (to set calendar)
      dt = 3600.0_dbl_kind   ! time step, s      
      npt = 99999            ! total number of time steps (dt) 
      diagfreq = 24          ! how often diag output is written
      diag_file = 'ice_diag' ! history file name prefix
      cpl_bgc = .false.      ! 
      dumpfreq='y'           ! restart frequency option
      restart = .false.      ! if true, read restart files for initialization
      restart_dir  = './'    ! write to executable dir for default
      restart_file = 'iced'  ! restart file name prefix
      ice_ic       = 'default'      ! specified in code

      kitd = 1           ! type of itd conversions (0 = delta, 1 = linear)
      kcatbound = 1      ! category boundary formula (0 = old, 1 = new, etc)
      ndtd = 1           ! dynamic time steps per thermodynamic time step
      kstrength = 1          ! 1 = Rothrock 75 strength, 0 = Hibler 79
      krdg_partic = 1        ! 1 = new participation, 0 = Thorndike et al 75
      krdg_redist = 1        ! 1 = new redistribution, 0 = Hibler 80
      mu_rdg = 3             ! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
      Cf = 17.0_dbl_kind     ! ratio of ridging work to PE change in ridging 
      shortwave = 'default'  ! 'default' or 'dEdd' (delta-Eddington)
      albedo_type = 'default'! or 'constant'
      ktherm = 1             ! 0 = 0-layer, 1 = BL99, 2 = mushy thermo
      conduct = 'bubbly'     ! 'MU71' or 'bubbly' (Pringle et al 2007)
      calc_Tsfc = .true.     ! calculate surface temperature
      update_ocn_f = .false. ! include fresh water and salt fluxes for frazil
      ustar_min = 0.005      ! minimum friction velocity for ocean heat flux (m/s)
      l_mpond_fresh = .false.     ! logical switch for including meltpond freshwater
                                  ! flux feedback to ocean model
      fbot_xfer_type = 'constant' ! transfer coefficient type for ocn heat flux
      R_ice     = 0.00_dbl_kind   ! tuning parameter for sea ice
      R_pnd     = 0.00_dbl_kind   ! tuning parameter for ponded sea ice
      R_snw     = 1.50_dbl_kind   ! tuning parameter for snow over sea ice
      dT_mlt    = 1.5_dbl_kind    ! change in temp to give non-melt to melt change
                                  ! in snow grain radius
      rsnw_mlt  = 1500._dbl_kind  ! maximum melting snow grain radius
      kalg      = 0.60_dbl_kind   ! algae absorption coefficient for 0.5 m thick layer
                                  ! 0.5 m path of 75 mg Chl a / m2
      hp1       = 0.01_dbl_kind   ! critical pond lid thickness for topo ponds
      hs0       = 0.03_dbl_kind   ! snow depth for transition to bare sea ice (m)
      hs1       = 0.03_dbl_kind   ! snow depth for transition to bare pond ice (m)
      dpscale   = c1              ! alter e-folding time scale for flushing 
      frzpnd    = 'cesm'          ! melt pond refreezing parameterization
      rfracmin  = 0.15_dbl_kind   ! minimum retained fraction of meltwater
      rfracmax  = 0.85_dbl_kind   ! maximum retained fraction of meltwater
      pndaspect = 0.8_dbl_kind    ! ratio of pond depth to area fraction
      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR
      ahmax     = 0.3_dbl_kind    ! thickness above which ice albedo is constant (m)
      atmbndy   = 'default'       ! or 'constant'

      default_season = 'winter' ! default forcing data
      fyear_init = 1998
      ycycle = 1

      atm_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      atm_data_type   = 'default'
      calc_strair     = .true.    ! calculate wind stress
      formdrag        = .false.   ! calculate form drag
      highfreq        = .false.   ! calculate high frequency RASM coupling
      natmiter        = 5         ! number of iterations for atm boundary layer calcs
      precip_units    = 'mks'     ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
      tfrz_option     = 'mushy'   ! freezing temp formulation
      oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
      ocn_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      ocn_data_type   = 'default'
      bgc_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      bgc_data_type   = 'default'
      data_dir    = ' '
!      dbug      = .false.         ! true writes diagnostics for input forcing
      restore_ocn     = .false.   ! restore sst if true
      trestore        = 90        ! restoring timescale, days (0 instantaneous)

      ! extra tracers
      tr_iage      = .false. ! ice age
      tr_FY        = .false. ! ice age
      tr_lvl       = .false. ! level ice 
      tr_pond_cesm = .false. ! CESM melt ponds
      tr_pond_lvl  = .false. ! level-ice melt ponds
      tr_pond_topo = .false. ! explicit melt ponds (topographic)
      tr_aero      = .false. ! aerosols

      ! mushy layer gravity drainage physics
      a_rapid_mode      =  0.5e-3_dbl_kind ! channel radius for rapid drainage mode (m)
      Rac_rapid_mode    =    10.0_dbl_kind ! critical Rayleigh number
      aspect_rapid_mode =     1.0_dbl_kind ! aspect ratio (larger is wider)
      dSdt_slow_mode    = -1.5e-7_dbl_kind ! slow mode drainage strength (m s-1 K-1)
      phi_c_slow_mode   =    0.05_dbl_kind ! critical liquid fraction porosity cutoff
      phi_i_mushy       =    0.85_dbl_kind ! liquid fraction of congelation ice

      !-----------------------------------------------------------------
      ! read from input file name from command line if it exists,
      ! otherwise the default is icepack_in
      !-----------------------------------------------------------------

      if ( command_argument_count() == 1 ) then
        call get_command_argument(1,nml_filename)
      endif

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
        nml_error = -1
      else
        nml_error =  1
      endif
      
      do while (nml_error > 0)
        print*,'Reading namelist file   ',nml_filename
        print*,'Reading setup_nml'
        read(nu_nml, nml=setup_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading tracer_nml'
        read(nu_nml, nml=tracer_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading thermo_nml'
        read(nu_nml, nml=thermo_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading shortwave_nml'
        read(nu_nml, nml=shortwave_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading ponds_nml'
        read(nu_nml, nml=ponds_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading forcing_nml'
        read(nu_nml, nml=forcing_nml,iostat=nml_error)
        if (nml_error /= 0) exit
      end do
      if (nml_error == 0) close(nu_nml)
      if (nml_error /= 0) then
        write(ice_stdout,*) 'error reading namelist'
        stop
      endif
      close(nu_nml)
      
      !-----------------------------------------------------------------
      ! set up diagnostics output and resolve conflicts
      !-----------------------------------------------------------------
      
      write(ice_stdout,*) 'Diagnostic output will be in files '
      write(ice_stdout,*)'    ','icepack.runlog.timestamp'
! tcraig, see below, no longer opened, using icepack.runlog.timestamp for "6"
!      write(ice_stdout,*)'    ',trim(diag_file)

      do n = 1,nx
        write(nx_names(n),'(a,i2.2)') 'point_',n
      enddo      
      nx_names(1) = 'icefree'
      nx_names(2) = 'slab'
      nx_names(3) = 'full_ITD'
      nx_names(4) = 'land'

      diag_len = len(trim(diag_file))
      do n = 1,nx
        diag_file_names=' '
        write(diag_file_names,'(a,a,a)') trim(diag_file),'.',trim(nx_names(n))
        write(ice_stdout,*)'    ',trim(diag_file_names)
        open(nu_diag_out+n-1, file=diag_file_names, status='unknown')
      end do

! tcraig, want nu_diag == ice_stdout == 6 to go to the icepack.runlog file
!      open (nu_diag, file=diag_file, status='unknown')
      
      write(nu_diag,*) '-----------------------------------'
      write(nu_diag,*) '  ICEPACK model diagnostic output  '
      write(nu_diag,*) '-----------------------------------'
      write(nu_diag,*) ' '



#ifndef ncdf
      ! netcdf is unavailable
!      atm_data_format = 'bin'
!      ocn_data_format = 'bin' 
#endif

      if (ncat == 1 .and. kitd == 1) then
         write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1.'
         write (nu_diag,*) 'Use kitd = 0 (delta function ITD) with kcatbound = 0'
         write (nu_diag,*) 'or for column configurations use kcatbound = -1'
         stop
      endif

      if (ncat /= 1 .and. kcatbound == -1) then
            write (nu_diag,*) &
               'WARNING: ITD required for ncat > 1'
            write (nu_diag,*) &
               'WARNING: Setting kitd and kcatbound to default values'
         kitd = 1
         kcatbound = 0
      endif

      rpcesm = c0
      rplvl  = c0
      rptopo = c0
      if (tr_pond_cesm) rpcesm = c1
      if (tr_pond_lvl ) rplvl  = c1
      if (tr_pond_topo) rptopo = c1

      tr_pond = .false. ! explicit melt ponds
      if (rpcesm + rplvl + rptopo > puny) tr_pond = .true.

      if (rpcesm + rplvl + rptopo > c1 + puny) then
            write (nu_diag,*) 'WARNING: Must use only one melt pond scheme'
            stop
      endif

      if (tr_pond_lvl .and. .not. tr_lvl) then
            write (nu_diag,*) 'WARNING: tr_pond_lvl=T but tr_lvl=F'
            write (nu_diag,*) 'WARNING: Setting tr_lvl=T'
         tr_lvl = .true.
      endif

      if (tr_pond_lvl .and. abs(hs0) > puny) then
            write (nu_diag,*) 'WARNING: tr_pond_lvl=T and hs0/=0'
            write (nu_diag,*) 'WARNING: Setting hs0=0'
         hs0 = c0
      endif

      if (tr_pond_cesm .and. trim(frzpnd) /= 'cesm') then
            write (nu_diag,*) 'WARNING: tr_pond_cesm=T'
            write (nu_diag,*) 'WARNING: frzpnd, dpscale not used'
         frzpnd = 'cesm'
      endif

      if (trim(shortwave) /= 'dEdd' .and. tr_pond .and. calc_tsfc) then
            write (nu_diag,*) 'WARNING: Must use dEdd shortwave'
            write (nu_diag,*) 'WARNING: with tr_pond and calc_tsfc=T.'
            write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
         shortwave = 'dEdd'
      endif

      if (tr_aero .and. n_aero==0) then
            write (nu_diag,*) 'WARNING: aerosols activated but'
            write (nu_diag,*) 'WARNING: not allocated in tracer array.'
            write (nu_diag,*) 'WARNING: Activate in compilation script.'
         stop
      endif

      if (tr_aero .and. trim(shortwave) /= 'dEdd') then
            write (nu_diag,*) 'WARNING: aerosols activated but dEdd'
            write (nu_diag,*) 'WARNING: shortwave is not.'
            write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
         shortwave = 'dEdd'
      endif

      rfracmin = min(max(rfracmin,c0),c1)
      rfracmax = min(max(rfracmax,c0),c1)

      if (ktherm == 2 .and. .not. calc_Tsfc) then
            write (nu_diag,*) 'WARNING: ktherm = 2 and calc_Tsfc = F'
            write (nu_diag,*) 'WARNING: Setting calc_Tsfc = T'
         calc_Tsfc = .true.
      endif

      if (ktherm == 1 .and. trim(tfrz_option) /= 'linear_salt') then
         write (nu_diag,*) &
         'WARNING: ktherm = 1 and tfrz_option = ',trim(tfrz_option)
         write (nu_diag,*) &
         'WARNING: For consistency, set tfrz_option = linear_salt'
      endif
      if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
         write (nu_diag,*) &
         'WARNING: ktherm = 2 and tfrz_option = ',trim(tfrz_option)
         write (nu_diag,*) &
         'WARNING: For consistency, set tfrz_option = mushy'
      endif

      if (formdrag) then
      if (trim(atmbndy) == 'constant') then
            write (nu_diag,*) 'WARNING: atmbndy = constant not allowed with formdrag'
            write (nu_diag,*) 'WARNING: Setting atmbndy = default'
         atmbndy = 'default'
      endif

      if (.not. calc_strair) then
            write (nu_diag,*) 'WARNING: formdrag=T but calc_strair=F'
            write (nu_diag,*) 'WARNING: Setting calc_strair=T'
         calc_strair = .true.
      endif

      if (tr_pond_cesm) then
            write (nu_diag,*) 'ERROR: formdrag=T but frzpnd=cesm' 
         stop
      endif

      if (.not. tr_lvl) then
            write (nu_diag,*) 'WARNING: formdrag=T but tr_lvl=F'
            write (nu_diag,*) 'WARNING: Setting tr_lvl=T'
         tr_lvl = .true.
      endif
      endif

      if (trim(fbot_xfer_type) == 'Cdn_ocn' .and. .not. formdrag)  then
            write (nu_diag,*) 'WARNING: formdrag=F but fbot_xfer_type=Cdn_ocn'
            write (nu_diag,*) 'WARNING: Setting fbot_xfer_type = constant'
         fbot_xfer_type = 'constant'
      endif

      call icepack_init_constants(Cf_in=Cf)

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------

         write(nu_diag,*) ' Document ice_in namelist parameters:'
         write(nu_diag,*) ' ==================================== '
         write(nu_diag,*) ' '
         write(nu_diag,1020) ' days_per_year             = ', days_per_year
         write(nu_diag,1010) ' use_leap_years            = ', use_leap_years
         write(nu_diag,1020) ' year_init                 = ', year_init
         write(nu_diag,1020) ' istep0                    = ', istep0
         write(nu_diag,1000) ' dt                        = ', dt
         write(nu_diag,1020) ' npt                       = ', npt
         write(nu_diag,1020) ' diagfreq                  = ', diagfreq
         write(nu_diag,1030) ' dumpfreq                  = ', &
                               trim(dumpfreq)
         write(nu_diag,1010) ' restart                   = ', restart
         write(nu_diag,*)    ' restart_dir               = ', &
                               trim(restart_dir)
         write(nu_diag,*)    ' restart_file              = ', &
                               trim(restart_file)
         write(nu_diag,*)    ' ice_ic                    = ', &
                               trim(ice_ic)
         write(nu_diag,1020) ' kitd                      = ', kitd
         write(nu_diag,1020) ' kcatbound                 = ', &
                               kcatbound
         write(nu_diag,1020) ' ndtd                      = ', ndtd
         write(nu_diag,1020) ' kstrength                 = ', kstrength
         write(nu_diag,1020) ' krdg_partic               = ', &
                               krdg_partic
         write(nu_diag,1020) ' krdg_redist               = ', &
                               krdg_redist
         if (krdg_redist == 1) &
         write(nu_diag,1000) ' mu_rdg                    = ', mu_rdg
         if (kstrength == 1) &
         write(nu_diag,1000) ' Cf                        = ', Cf
         write(nu_diag,1030) ' shortwave                 = ', &
                               trim(shortwave)
         if (cpl_bgc) then
             write(nu_diag,1000) ' BGC coupling is switched ON'
         else
             write(nu_diag,1000) ' BGC coupling is switched OFF'
          endif

         if (trim(shortwave) == 'dEdd') then
         write(nu_diag,1000) ' R_ice                     = ', R_ice
         write(nu_diag,1000) ' R_pnd                     = ', R_pnd
         write(nu_diag,1000) ' R_snw                     = ', R_snw
         write(nu_diag,1000) ' dT_mlt                    = ', dT_mlt
         write(nu_diag,1000) ' rsnw_mlt                  = ', rsnw_mlt
         write(nu_diag,1000) ' kalg                      = ', kalg
         write(nu_diag,1000) ' hp1                       = ', hp1
         write(nu_diag,1000) ' hs0                       = ', hs0
         else
         write(nu_diag,1030) ' albedo_type               = ', &
                               trim(albedo_type)
         write(nu_diag,1000) ' albicev                   = ', albicev
         write(nu_diag,1000) ' albicei                   = ', albicei
         write(nu_diag,1000) ' albsnowv                  = ', albsnowv
         write(nu_diag,1000) ' albsnowi                  = ', albsnowi
         write(nu_diag,1000) ' ahmax                     = ', ahmax
         endif

         write(nu_diag,1000) ' rfracmin                  = ', rfracmin
         write(nu_diag,1000) ' rfracmax                  = ', rfracmax
         if (tr_pond_lvl) then
         write(nu_diag,1000) ' hs1                       = ', hs1
         write(nu_diag,1000) ' dpscale                   = ', dpscale
         write(nu_diag,1030) ' frzpnd                    = ', trim(frzpnd)
         endif
         if (tr_pond .and. .not. tr_pond_lvl) &
         write(nu_diag,1000) ' pndaspect                 = ', pndaspect

         write(nu_diag,1020) ' ktherm                    = ', ktherm
         if (ktherm == 1) &
         write(nu_diag,1030) ' conduct                   = ', conduct
         if (ktherm == 2) then
         write(nu_diag,1005) ' a_rapid_mode              = ', a_rapid_mode
         write(nu_diag,1005) ' Rac_rapid_mode            = ', Rac_rapid_mode
         write(nu_diag,1005) ' aspect_rapid_mode         = ', aspect_rapid_mode
         write(nu_diag,1005) ' dSdt_slow_mode            = ', dSdt_slow_mode
         write(nu_diag,1005) ' phi_c_slow_mode           = ', phi_c_slow_mode
         write(nu_diag,1005) ' phi_i_mushy               = ', phi_i_mushy
         endif

         write(nu_diag,1030) ' atmbndy                   = ', &
                               trim(atmbndy)
         write(nu_diag,1010) ' formdrag                  = ', formdrag
         write(nu_diag,1010) ' highfreq                  = ', highfreq
         write(nu_diag,1020) ' natmiter                  = ', natmiter
         write(nu_diag,1010) ' calc_strair               = ', calc_strair
         write(nu_diag,1010) ' calc_Tsfc                 = ', calc_Tsfc

         write(nu_diag,*)    ' atm_data_type             = ', &
                               trim(atm_data_type)
         write(nu_diag,*)    ' ocn_data_type             = ', &
                               trim(ocn_data_type)
         write(nu_diag,*)    ' bgc_data_type             = ', &
                               trim(bgc_data_type)

         write(nu_diag,1010) ' update_ocn_f              = ', update_ocn_f
         write(nu_diag,1010) ' l_mpond_fresh             = ', l_mpond_fresh
         write(nu_diag,1005) ' ustar_min                 = ', ustar_min
         write(nu_diag, *)   ' fbot_xfer_type            = ', &
                               trim(fbot_xfer_type)
         write(nu_diag,1010) ' oceanmixed_ice            = ', &
                               oceanmixed_ice
         write(nu_diag,*)    ' tfrz_option               = ', &
                               trim(tfrz_option)

         write (nu_diag,*) ' '

         write(nu_diag,1010) ' restore_ocn               = ', &
             restore_ocn
         !if (restore_ice .or. restore_ocn) &
         if ( restore_ocn) &
         write(nu_diag,1005) ' trestore                  = ', trestore

         ! tracers
         write(nu_diag,1010) ' tr_iage                   = ', tr_iage
         write(nu_diag,1010) ' tr_FY                     = ', tr_FY
         write(nu_diag,1010) ' tr_lvl                    = ', tr_lvl
         write(nu_diag,1010) ' tr_pond_cesm              = ', tr_pond_cesm
         write(nu_diag,1010) ' tr_pond_lvl               = ', tr_pond_lvl
         write(nu_diag,1010) ' tr_pond_topo              = ', tr_pond_topo
         write(nu_diag,1010) ' tr_aero                   = ', tr_aero

         nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
         ntrcr = 1             ! count tracers, starting with Tsfc = 1

         nt_qice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! qice in nilyr layers

         nt_qsno = ntrcr + 1
         ntrcr = ntrcr + nslyr ! qsno in nslyr layers

         nt_sice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! sice in nilyr layers

         nt_iage = max_ntrcr
         if (tr_iage) then
             ntrcr = ntrcr + 1
             nt_iage = ntrcr   ! chronological ice age
         endif

         nt_FY = max_ntrcr
         if (tr_FY) then
             ntrcr = ntrcr + 1
             nt_FY = ntrcr     ! area of first year ice
         endif

         nt_alvl = max_ntrcr
         nt_vlvl = max_ntrcr
         if (tr_lvl) then
             ntrcr = ntrcr + 1
             nt_alvl = ntrcr
             ntrcr = ntrcr + 1
             nt_vlvl = ntrcr
         endif

         nt_apnd = max_ntrcr
         nt_hpnd = max_ntrcr
         nt_ipnd = max_ntrcr
         if (tr_pond) then            ! all explicit melt pond schemes
             ntrcr = ntrcr + 1
             nt_apnd = ntrcr
             ntrcr = ntrcr + 1
             nt_hpnd = ntrcr
             if (tr_pond_lvl) then
                 ntrcr = ntrcr + 1    ! refrozen pond ice lid thickness
                 nt_ipnd = ntrcr      ! on level-ice ponds (if frzpnd='hlid')
             endif
             if (tr_pond_topo) then
                 ntrcr = ntrcr + 1    ! 
                 nt_ipnd = ntrcr      ! refrozen pond ice lid thickness
             endif
         endif

         nt_aero = max_ntrcr - 4*n_aero
         if (tr_aero) then
             nt_aero = ntrcr + 1
             ntrcr = ntrcr + 4*n_aero ! 4 dEdd layers, n_aero species
         endif
              
         if (ntrcr > max_ntrcr-1) then
            write(nu_diag,*) 'max_ntrcr-1 < number of namelist tracers'
            write(nu_diag,*) 'max_ntrcr-1 = ',max_ntrcr-1,' ntrcr = ',ntrcr
            stop
         endif                               

         write(nu_diag,*) ' '
         write(nu_diag,1020) 'max_ntrcr = ', max_ntrcr
         write(nu_diag,1020) 'ntrcr = ', ntrcr
         write(nu_diag,*) ' '
         write(nu_diag,1020)'nt_sice = ', nt_sice
         write(nu_diag,1020)'nt_qice = ', nt_qice
         write(nu_diag,1020)'nt_qsno = ', nt_qsno
         write(nu_diag,*)' '
         write(nu_diag,1020)'ncat', ncat
         write(nu_diag,1020)'nilyr', nilyr
         write(nu_diag,1020)'nslyr', nslyr
         write(nu_diag,*)' '
         write(nu_diag,1020)'nx', nx
         write(nu_diag,*)' '

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character
 1040    format (a30,2x,6i6)   ! integer
 1050    format (a30,2x,6a6)   ! character

      if (formdrag) then
         if (nt_apnd==0) then
            write(nu_diag,*)'ERROR: nt_apnd:',nt_apnd
            stop
         elseif (nt_hpnd==0) then
            write(nu_diag,*)'ERROR: nt_hpnd:',nt_hpnd
            stop
         elseif (nt_ipnd==0) then
            write(nu_diag,*)'ERROR: nt_ipnd:',nt_ipnd
            stop
         elseif (nt_alvl==0) then
            write(nu_diag,*)'ERROR: nt_alvl:',nt_alvl
            stop
         elseif (nt_vlvl==0) then
            write(nu_diag,*)'ERROR: nt_vlvl:',nt_vlvl
            stop
         endif
      endif

      end subroutine input_data

!=======================================================================

! Horizontal grid initialization:
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_grid2

      use icepack_drv_constants, only: c0, puny
      use icepack_drv_constants, only: pi, p5, c1

      integer :: i

      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      TLAT(:) = p5*pi  ! pi/2, North pole
      TLON(:) = c0

      do i = 2, nx
         TLAT(i) = TLAT(i-1) - p5*pi/180._dbl_kind ! half-deg increments
      enddo

      tmask(:) = .true.
      tmask(nx) = .false.  ! land in last grid cell

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

         lmask_n(:) = .false.
         lmask_s(:) = .false.

         do i = 1, nx
            if (TLAT(i) >= -puny) lmask_n(i) = .true. ! N. Hem.
            if (TLAT(i) <  -puny) lmask_s(i) = .true. ! S. Hem.
         enddo

      end subroutine init_grid2

!=======================================================================

! Initialize state for the itd model
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_state

      use icepack_intfc, only: icepack_aggregate
      use icepack_drv_parameters, only: heat_capacity
      use icepack_drv_constants, only: c0, c1, nu_diag
      use icepack_drv_domain_size, only: ncat, nilyr, nslyr, max_ntrcr, n_aero
      use icepack_drv_flux, only: sst, Tf, Tair, salinz, Tmltz
      use icepack_drv_state, only: trcr_depend, aicen, trcrn, vicen, vsnon
      use icepack_drv_state, only: aice0, aice, vice, vsno, trcr, aice_init
      use icepack_drv_state, only: n_trcr_strata, nt_strata, trcr_base
      use icepack_drv_tracers, only: ntrcr
      use icepack_drv_tracers, only: tr_iage, tr_FY, tr_lvl, tr_aero
      use icepack_drv_tracers, only: tr_pond_cesm, tr_pond_lvl, tr_pond_topo
      use icepack_drv_tracers, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_fy
      use icepack_drv_tracers, only: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero

      integer (kind=int_kind) :: &
         i           , & ! horizontal indes
         k           , & ! vertical index
         it              ! tracer index

      !-----------------------------------------------------------------
      ! Check number of layers in ice and snow.
      !-----------------------------------------------------------------
         if (nilyr < 1) then
            write (nu_diag,*) 'nilyr =', nilyr
            write (nu_diag,*) 'Must have at least one ice layer'
            stop
         endif

         if (nslyr < 1) then
            write (nu_diag,*) 'nslyr =', nslyr
            write (nu_diag,*) 'Must have at least one snow layer'
            stop
         endif

         if (.not.heat_capacity) then

            write (nu_diag,*) 'WARNING - Zero-layer thermodynamics'

            if (nilyr > 1) then
               write (nu_diag,*) 'nilyr =', nilyr
               write (nu_diag,*)        &
                    'Must have nilyr = 1 if ktherm = 0'
               stop
            endif

            if (nslyr > 1) then
               write (nu_diag,*) 'nslyr =', nslyr
               write (nu_diag,*)        &
                    'Must have nslyr = 1 if heat_capacity = F'
               stop
            endif

         endif   ! heat_capacity = F

      !-----------------------------------------------------------------
      ! Set tracer types
      !-----------------------------------------------------------------

      trcr_depend(nt_Tsfc) = 0 ! ice/snow surface temperature
      do k = 1, nilyr
         trcr_depend(nt_sice + k - 1) = 1 ! volume-weighted ice salinity
         trcr_depend(nt_qice + k - 1) = 1 ! volume-weighted ice enthalpy
      enddo
      do k = 1, nslyr
         trcr_depend(nt_qsno + k - 1) = 2 ! volume-weighted snow enthalpy
      enddo
      if (tr_iage) trcr_depend(nt_iage)  = 1   ! volume-weighted ice age
      if (tr_FY)   trcr_depend(nt_FY)    = 0   ! area-weighted first-year ice area
      if (tr_lvl)  trcr_depend(nt_alvl)  = 0   ! level ice area
      if (tr_lvl)  trcr_depend(nt_vlvl)  = 1   ! level ice volume
      if (tr_pond_cesm) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
      endif
      if (tr_pond_lvl) then
                   trcr_depend(nt_apnd)  = 2+nt_alvl   ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_pond_topo) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_aero) then ! volume-weighted aerosols
         do it = 1, n_aero
            trcr_depend(nt_aero+(it-1)*4  ) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+1) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+2) = 1 ! ice
            trcr_depend(nt_aero+(it-1)*4+3) = 1 ! ice
         enddo
      endif

      do it = 1, ntrcr
         ! mask for base quantity on which tracers are carried
         if (trcr_depend(it) == 0) then      ! area
            trcr_base(it,1) = c1 
         elseif (trcr_depend(it) == 1) then  ! ice volume
            trcr_base(it,2) = c1 
         elseif (trcr_depend(it) == 2) then  ! snow volume
            trcr_base(it,3) = c1 
         else
            trcr_base(it,1) = c1    ! default: ice area
            trcr_base(it,2) = c0
            trcr_base(it,3) = c0
         endif

         ! initialize number of underlying tracer layers
         n_trcr_strata(it) = 0
         ! default indices of underlying tracer layers
         nt_strata   (it,1) = 0
         nt_strata   (it,2) = 0
      enddo

      if (tr_pond_cesm) then
         n_trcr_strata(nt_hpnd)   = 1       ! melt pond depth
         nt_strata    (nt_hpnd,1) = nt_apnd ! on melt pond area
      endif
      if (tr_pond_lvl) then
         n_trcr_strata(nt_apnd)   = 1       ! melt pond area
         nt_strata    (nt_apnd,1) = nt_alvl ! on level ice area
         n_trcr_strata(nt_hpnd)   = 2       ! melt pond depth
         nt_strata    (nt_hpnd,2) = nt_apnd ! on melt pond area
         nt_strata    (nt_hpnd,1) = nt_alvl ! on level ice area
         n_trcr_strata(nt_ipnd)   = 2       ! refrozen pond lid
         nt_strata    (nt_ipnd,2) = nt_apnd ! on melt pond area
         nt_strata    (nt_ipnd,1) = nt_alvl ! on level ice area
      endif
      if (tr_pond_topo) then
         n_trcr_strata(nt_hpnd)   = 1       ! melt pond depth
         nt_strata    (nt_hpnd,1) = nt_apnd ! on melt pond area
         n_trcr_strata(nt_ipnd)   = 1       ! refrozen pond lid
         nt_strata    (nt_ipnd,1) = nt_apnd ! on melt pond area
      endif

      !-----------------------------------------------------------------
      ! Set state variables
      !-----------------------------------------------------------------

      call set_state_var (nx,          ice_ic,       &
          TLON  (:),   TLAT (:),     &
          Tair  (:),   sst  (:),     &
          Tf    (:),                 &
          salinz(:,:), Tmltz(:,:),   &
          aicen (:,:), trcrn(:,:,:), &
          vicen (:,:), vsnon(:,:))

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, nx
         aice(i) = c0
         vice(i) = c0
         vsno(i) = c0
         do it = 1, max_ntrcr
            trcr(i,it) = c0
         enddo

         if (tmask(i)) &
         call icepack_aggregate (ncat,               &
                                aicen(i,:),  &
                                trcrn(i,1:ntrcr,:), &
                                vicen(i,:),   &
                                vsnon(i,:),   &
                                aice (i),   &
                                trcr (i,1:ntrcr),   &
                                vice (i),   &
                                vsno (i),   &
                                aice0(i),   &
                                ntrcr,               &
                                trcr_depend  (1:ntrcr),&
                                trcr_base    (1:ntrcr,:),&
                                n_trcr_strata(1:ntrcr),&
                                nt_strata    (1:ntrcr,:))

         aice_init(i) = aice(i)

      enddo
      !$OMP END PARALLEL DO

      end subroutine init_state

!=======================================================================

! Initialize state in each ice thickness category
!
! authors: Elizabeth Hunke, LANL

      subroutine set_state_var (nx,       ice_ic,   &
                                TLON,     TLAT, &
                                Tair,     sst,  &
                                Tf,       &
                                salinz,   Tmltz, &
                                aicen,    trcrn, &
                                vicen,    vsnon)

      use icepack_drv_arrays_column, only: hin_max
      use icepack_intfc, only: icepack_init_trcr
      use icepack_drv_constants, only: c0, c1, c2, c3, p2, p5, rhoi, rhos, Lfresh
      use icepack_drv_constants, only: cp_ice, cp_ocn, Tsmelt, Tffresh, puny
      use icepack_drv_domain_size, only: nilyr, nslyr, max_ntrcr, ncat
      use icepack_drv_tracers, only: tr_brine, tr_lvl
      use icepack_drv_tracers, only: nt_Tsfc, nt_qice, nt_qsno, nt_sice
      use icepack_drv_tracers, only: nt_fbri, nt_alvl, nt_vlvl
!      use icepack_drv_forcing, only: atm_data_type

      integer (kind=int_kind), intent(in) :: &
         nx          ! number of grid cells

      character(len=char_len_long), intent(in) :: & 
         ice_ic      ! method of ice cover initialization  ! not used for now

      real (kind=dbl_kind), dimension (nx), &
         intent(in) :: &
         TLON   , & ! longitude of temperature pts (radians)  ! not used for now
         TLAT       ! latitude of temperature pts (radians)

      real (kind=dbl_kind), dimension (nx), intent(in) :: &
         Tair        ! air temperature  (K)

      ! ocean values may be redefined here, unlike in CICE
      real (kind=dbl_kind), dimension (nx), intent(inout) :: &
         Tf      , & ! freezing temperature (C) 
         sst         ! sea surface temperature (C) 

      real (kind=dbl_kind), dimension (nx,nilyr), &
         intent(in) :: &
         salinz  , & ! initial salinity profile
         Tmltz       ! initial melting temperature profile

      real (kind=dbl_kind), dimension (nx,ncat), &
         intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx,max_ntrcr,ncat), &
         intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i           , & ! horizontal indices
         k           , & ! ice layer index
         n           , & ! thickness category index
         it              ! tracer index

      real (kind=dbl_kind) :: &
         Tsfc, sum, hbar

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      real (kind=dbl_kind), parameter :: &
         hsno_init = 0.25_dbl_kind   ! initial snow thickness (m)

      ! Initialize state variables.
      ! If restarting, these values are overwritten.

      do n = 1, ncat
         do i = 1, nx
            aicen(i,n) = c0
            vicen(i,n) = c0
            vsnon(i,n) = c0
            trcrn(i,nt_Tsfc,n) = Tf(i)  ! surface temperature 
            if (max_ntrcr >= 2) then
               do it = 2, max_ntrcr
                  trcrn(i,it,n) = c0
               enddo
            endif
            if (tr_lvl)   trcrn(i,nt_alvl,n) = c1
            if (tr_lvl)   trcrn(i,nt_vlvl,n) = c1
            if (tr_brine) trcrn(i,nt_fbri,n) = c1
            do k = 1, nilyr
               trcrn(i,nt_sice+k-1,n) = salinz(i,k)
            enddo
            do k = 1, nslyr
               trcrn(i,nt_qsno+k-1,n) = -rhos * Lfresh
            enddo
         enddo
         ainit(n) = c0
         hinit(n) = c0
      enddo

      !-----------------------------------------------------------------
      ! For Icepack testing, the grid vector is populated with several
      ! different ice distributions, including ice-free, a single-
      ! thickness slab, a full thickness distribution (as in CICE),
      ! and land
      !-----------------------------------------------------------------

      i = 1  ! ice-free
             ! already initialized above 

      !-----------------------------------------------------------------

      i = 2  ! 2-m slab, no snow
      ainit(3) = c1  ! assumes we are using the default ITD boundaries
      hinit(3) = c2
      do n = 1, ncat
        ! ice volume, snow volume
        aicen(i,n) = ainit(n)
        vicen(i,n) = hinit(n) * ainit(n) ! m
        vsnon(i,n) = c0
        ! tracers
        call icepack_init_trcr(Tair(i),     Tf(i),      &
            salinz(i,:), Tmltz(i,:), &
            Tsfc,                        &
            nilyr,         nslyr,        &
            qin(:),        qsn(:))
        
        ! surface temperature
        trcrn(i,nt_Tsfc,n) = Tsfc ! deg C
        ! ice enthalpy, salinity 
        do k = 1, nilyr
          trcrn(i,nt_qice+k-1,n) = qin(k)
          trcrn(i,nt_sice+k-1,n) = salinz(i,k)
        enddo
        ! snow enthalpy
        do k = 1, nslyr
          trcrn(i,nt_qsno+k-1,n) = qsn(k)
        enddo               ! nslyr
        ! brine fraction
        if (tr_brine) trcrn(i,nt_fbri,n) = c1
      enddo                  ! ncat
      
      !-----------------------------------------------------------------

      i = 3  ! full thickness distribution
      ! initial category areas in cells with ice
      hbar = c3  ! initial ice thickness with greatest area
      ! Note: the resulting average ice thickness 
      ! tends to be less than hbar due to the
      ! nonlinear distribution of ice thicknesses 

      sum = c0
      do n = 1, ncat
        if (n < ncat) then
          hinit(n) = p5*(hin_max(n-1) + hin_max(n)) ! m
        else                ! n=ncat
          hinit(n) = (hin_max(n-1) + c1) ! m
        endif
        ! parabola, max at h=hbar, zero at h=0, 2*hbar
        ainit(n) = max(c0, (c2*hbar*hinit(n) - hinit(n)**2))
        sum = sum + ainit(n)
      enddo
      do n = 1, ncat
        ainit(n) = ainit(n) / (sum + puny/ncat) ! normalize
      enddo
      
      do n = 1, ncat
        ! ice volume, snow volume
        aicen(i,n) = ainit(n)
        vicen(i,n) = hinit(n) * ainit(n) ! m
        vsnon(i,n) = min(aicen(i,n)*hsno_init,p2*vicen(i,n))
        ! tracers
        call icepack_init_trcr(Tair(i),     Tf(i),      &
            salinz(i,:), Tmltz(i,:), &
            Tsfc,                        &
            nilyr,         nslyr,        &
            qin(:),        qsn(:))
        
        ! surface temperature
        trcrn(i,nt_Tsfc,n) = Tsfc ! deg C
        ! ice enthalpy, salinity 
        do k = 1, nilyr
          trcrn(i,nt_qice+k-1,n) = qin(k)
          trcrn(i,nt_sice+k-1,n) = salinz(i,k)
        enddo
        ! snow enthalpy
        do k = 1, nslyr
          trcrn(i,nt_qsno+k-1,n) = qsn(k)
        enddo               ! nslyr
        ! brine fraction
        if (tr_brine) trcrn(i,nt_fbri,n) = c1
      enddo                  ! ncat
      
      !-----------------------------------------------------------------
      
      ! land
      ! already initialized above (tmask = 0)
      i = 4
      sst(i) = c0
      Tf(i) = c0

    end subroutine set_state_var

!=======================================================================

  end module icepack_drv_init

!=======================================================================
