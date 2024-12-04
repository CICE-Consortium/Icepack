!=======================================================================

! parameter and variable initializations
!
! authors Elizabeth C. Hunke, LANL

      module icedrv_init

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, ice_stdout, nu_diag_out, nu_nml
      use icedrv_constants, only: c0, c1, c2, c3, p2, p5
      use icedrv_domain_size, only: nx
      use icedrv_flux, only: sst_init
      use icepack_intfc, only: icepack_init_parameters
      use icepack_intfc, only: icepack_init_fsd
      use icepack_intfc, only: icepack_init_tracer_flags
      use icepack_intfc, only: icepack_init_tracer_sizes
      use icepack_intfc, only: icepack_init_tracer_indices
      use icepack_intfc, only: icepack_init_enthalpy
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_indices
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icedrv_system, only: icedrv_system_abort, icedrv_system_flush

      implicit none
      private
      public :: input_data, init_grid2, init_state, init_fsd

      character(len=char_len_long), public :: &
         ice_ic      ! method of ice cover initialization
                     ! 'default' or 'none' => conditions specified in code
                     ! restart = .true. overwrites default initial
                     !    condition using filename given by ice_ic

      real (kind=dbl_kind), dimension (nx), public :: &
         TLON   , & ! longitude of temp pts (radians)
         TLAT       ! latitude of temp pts (radians)

      logical (kind=log_kind), &
         dimension (nx), public :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      real (kind=dbl_kind), public :: &
         hi_init_slab,   & ! initial ice thickness for slab cell (nx=2)
         hsno_init_slab, & ! initial snow thickness for slab cell (nx=2)
         hbar_init_itd,  & ! hbar for ice thickness for itd cell (nx=3)
         hsno_init_itd     ! initial snow thickness for itd cell (nx=3)

!=======================================================================

      contains

!=======================================================================

! Namelist variables, set to default values; may be altered
! at run time
!
! author Elizabeth C. Hunke, LANL

      subroutine input_data

      use icedrv_diagnostics, only: diag_file, nx_names
      use icedrv_domain_size, only: nilyr, nslyr, nblyr, max_ntrcr, ncat
      use icedrv_domain_size, only: n_iso, n_aero, nfsd
      use icedrv_calendar, only: year_init, istep0
      use icedrv_calendar, only: dumpfreq, diagfreq, dump_last
      use icedrv_calendar, only: npt, dt, ndtd, days_per_year, use_leap_years
      use icedrv_history, only: history_format
      use icedrv_restart_shared, only: restart, restart_dir, restart_file, restart_format
      use icedrv_flux, only: l_mpond_fresh, cpl_bgc
      use icedrv_flux, only: default_season
      use icedrv_flux, only: sss_fixed, qdp_fixed, hmix_fixed
      use icedrv_forcing, only: precip_units,    fyear_init,      ycycle
      use icedrv_forcing, only: atm_data_type,   ocn_data_type,   bgc_data_type
      use icedrv_forcing, only: atm_data_file,   ocn_data_file,   bgc_data_file
      use icedrv_forcing, only: ice_data_file
      use icedrv_forcing, only: atm_data_format, ocn_data_format, bgc_data_format
      use icedrv_forcing, only: data_dir
      use icedrv_forcing, only: oceanmixed_ice, restore_ocn, trestore
      use icedrv_forcing, only: snw_ssp_table, lateral_flux_type
      use icedrv_forcing, only: precalc_forc

      ! local variables

      character (32) :: &
         nml_filename = 'icepack_in' ! namelist input file name

      integer (kind=int_kind) :: &
         nml_error, & ! namelist i/o error flag
         n            ! loop index

      character (len=char_len) :: diag_file_names
      character (len=char_len), dimension(4) :: nx_names_default

      real (kind=dbl_kind) :: ustar_min, albicev, albicei, albsnowv, albsnowi, &
         ahmax, R_ice, R_pnd, R_snw, dT_mlt, rsnw_mlt, ksno, hi_min, Tliquidus_max, &
         mu_rdg, hs0, dpscale, rfracmin, rfracmax, pndaspect, hs1, hp1, &
         a_rapid_mode, Rac_rapid_mode, aspect_rapid_mode, dSdt_slow_mode, &
         phi_c_slow_mode, phi_i_mushy, kalg, emissivity, floediam, hfrazilmin, &
         rsnw_fall, rsnw_tmax, rhosnew, rhosmin, rhosmax, &
         windmin, drhosdwind, snwlvlfac

      integer (kind=int_kind) :: ktherm, kstrength, krdg_partic, krdg_redist, &
         natmiter, kitd, kcatbound

      character (len=char_len) :: shortwave, albedo_type, conduct, fbot_xfer_type, &
         cpl_frazil, congel_freeze, tfrz_option, saltflux_option, &
         frzpnd, atmbndy, wave_spec_type, snwredist, snw_aging_table

      logical (kind=log_kind) :: sw_redist, use_smliq_pnd, snwgrain, update_ocn_f
      real (kind=dbl_kind)    :: sw_frac, sw_dtemp

      ! Flux convergence tolerance
      real (kind=dbl_kind) :: atmiter_conv

      ! Ice reference salinity for fluxes
      real (kind=dbl_kind) :: ice_ref_salinity

      logical (kind=log_kind) :: calc_Tsfc, formdrag, highfreq, calc_strair, calc_dragio
      logical (kind=log_kind) :: conserv_check

      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_pond, tr_snow
      logical (kind=log_kind) :: tr_iso, tr_aero, tr_fsd
      logical (kind=log_kind) :: tr_pond_lvl, tr_pond_topo, wave_spec
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_FY
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
                                 nt_smice, nt_smliq, nt_rhos, nt_rsnw, &
                                 nt_aero, nt_fsd, nt_isosno, nt_isoice

      real (kind=real_kind) :: rplvl, rptopo
      real (kind=dbl_kind) :: Cf, puny

      character(len=*), parameter :: subname='(input_data)'

      !-----------------------------------------------------------------
      ! Namelist variables
      !-----------------------------------------------------------------

      namelist /setup_nml/ &
        days_per_year,  use_leap_years, year_init,       istep0,        &
        dt,             npt,            ndtd,            dump_last,     &
        ice_ic,         restart,        restart_dir,     restart_file,  &
        restart_format, &
        dumpfreq,       diagfreq,       diag_file,       cpl_bgc,       &
        conserv_check,  history_format,                                 &
        hi_init_slab,   hsno_init_slab, hbar_init_itd,   hsno_init_itd, &
        sst_init

      namelist /grid_nml/ &
        kcatbound

      namelist /thermo_nml/ &
        kitd,           ktherm,          ksno,     conduct,             &
        a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,             &
        dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy,                   &
        floediam,       hfrazilmin,      Tliquidus_max,    hi_min

      namelist /dynamics_nml/ &
        kstrength,      krdg_partic,    krdg_redist,    mu_rdg,         &
        Cf

      namelist /shortwave_nml/ &
        shortwave,      albedo_type,                                    &
        albicev,        albicei,         albsnowv,      albsnowi,       &
        ahmax,          R_ice,           R_pnd,         R_snw,          &
        sw_redist,      sw_frac,         sw_dtemp,                      &
        dT_mlt,         rsnw_mlt,        kalg,          snw_ssp_table

      namelist /ponds_nml/ &
        hs0,            dpscale,         frzpnd,                        &
        rfracmin,       rfracmax,        pndaspect,     hs1,            &
        hp1
      namelist /snow_nml/ &
        snwredist,      snwgrain,       rsnw_fall,     rsnw_tmax,      &
        rhosnew,        rhosmin,        rhosmax,       snwlvlfac,      &
        windmin,        drhosdwind,     use_smliq_pnd, snw_aging_table

      namelist /forcing_nml/ &
        atmbndy,         calc_strair,     calc_Tsfc,       &
        update_ocn_f,    l_mpond_fresh,   ustar_min,       &
        fbot_xfer_type,  oceanmixed_ice,  emissivity,      &
        formdrag,        highfreq,        natmiter,        &
        atmiter_conv,    calc_dragio,     congel_freeze,   &
        tfrz_option,     saltflux_option, ice_ref_salinity, &
        default_season,  wave_spec_type,  cpl_frazil,      &
        precip_units,    fyear_init,      ycycle,          &
        atm_data_type,   ocn_data_type,   bgc_data_type,   &
        lateral_flux_type,                                 &
        atm_data_file,   ocn_data_file,   bgc_data_file,   &
        ice_data_file,                                     &
        atm_data_format, ocn_data_format, bgc_data_format, &
        data_dir,        trestore,        restore_ocn,     &
        sss_fixed,       qdp_fixed,       hmix_fixed,     &
        precalc_forc

      namelist /tracer_nml/   &
        tr_iage,      &
        tr_FY,        &
        tr_lvl,       &
        tr_pond_lvl,  &
        tr_pond_topo, &
        tr_snow,      &
        tr_aero,      &
        tr_fsd,       &
        tr_iso

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(ustar_min_out=ustar_min, Cf_out=Cf, &
           albicev_out=albicev, albicei_out=albicei, ksno_out = ksno,   &
           albsnowv_out=albsnowv, albsnowi_out=albsnowi, hi_min_out=hi_min, &
           natmiter_out=natmiter, ahmax_out=ahmax, shortwave_out=shortwave, &
           atmiter_conv_out = atmiter_conv, calc_dragio_out=calc_dragio, &
           albedo_type_out=albedo_type, R_ice_out=R_ice, R_pnd_out=R_pnd, &
           R_snw_out=R_snw, dT_mlt_out=dT_mlt, rsnw_mlt_out=rsnw_mlt, &
           kstrength_out=kstrength, krdg_partic_out=krdg_partic, &
           krdg_redist_out=krdg_redist, mu_rdg_out=mu_rdg, &
           atmbndy_out=atmbndy, calc_strair_out=calc_strair, &
           formdrag_out=formdrag, highfreq_out=highfreq, &
           emissivity_out=emissivity, &
           kitd_out=kitd, kcatbound_out=kcatbound, hs0_out=hs0, &
           dpscale_out=dpscale, frzpnd_out=frzpnd, &
           rfracmin_out=rfracmin, rfracmax_out=rfracmax, &
           pndaspect_out=pndaspect, hs1_out=hs1, hp1_out=hp1, &
           ktherm_out=ktherm, calc_Tsfc_out=calc_Tsfc, &
           floediam_out=floediam, hfrazilmin_out=hfrazilmin, &
           update_ocn_f_out = update_ocn_f, cpl_frazil_out = cpl_frazil, &
           conduct_out=conduct, a_rapid_mode_out=a_rapid_mode, &
           Rac_rapid_mode_out=Rac_rapid_mode, &
           aspect_rapid_mode_out=aspect_rapid_mode, &
           dSdt_slow_mode_out=dSdt_slow_mode, &
           phi_c_slow_mode_out=phi_c_slow_mode, Tliquidus_max_out=Tliquidus_max, &
           phi_i_mushy_out=phi_i_mushy, conserv_check_out=conserv_check, &
           congel_freeze_out=congel_freeze, &
           tfrz_option_out=tfrz_option, saltflux_option_out=saltflux_option, &
           ice_ref_salinity_out=ice_ref_salinity, kalg_out=kalg, &
           fbot_xfer_type_out=fbot_xfer_type, puny_out=puny, &
           wave_spec_type_out=wave_spec_type, &
           sw_redist_out=sw_redist, sw_frac_out=sw_frac, sw_dtemp_out=sw_dtemp, &
           snwredist_out=snwredist, use_smliq_pnd_out=use_smliq_pnd, &
           snwgrain_out=snwgrain, rsnw_fall_out=rsnw_fall, rsnw_tmax_out=rsnw_tmax, &
           rhosnew_out=rhosnew, rhosmin_out = rhosmin, rhosmax_out=rhosmax, &
           windmin_out=windmin, drhosdwind_out=drhosdwind, snwlvlfac_out=snwlvlfac, &
           snw_aging_table_out=snw_aging_table)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

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
      dump_last=.false.      ! restart at end of run
      restart = .false.      ! if true, read restart files for initialization
      restart_dir  = './'    ! write to executable dir for default
      restart_file = 'iced'  ! restart file name prefix
      restart_format = 'bin' ! default restart format is binary, other option 'nc'
      history_format = 'none'     ! if 'nc', write history files. Otherwise do nothing
      ice_ic       = 'default'    ! initial conditions are specified in the code
                                  ! otherwise, the filename for reading restarts
      hi_init_slab   = c2            ! initial ice thickness for slab cell (nx=2)
      hsno_init_slab = c0            ! initial snow thickness for slab cell (nx=2)
      hbar_init_itd  = c3            ! hbar for ice thickness for itd cell (nx=3)
      hsno_init_itd  = 0.25_dbl_kind ! initial snow thickness for itd cell (nx=3)
      sst_init       = -1.8_dbl_kind ! initial mixed layer temperature (all cells)
      ndtd = 1               ! dynamic time steps per thermodynamic time step
      l_mpond_fresh = .false.     ! logical switch for including meltpond freshwater
                                  ! flux feedback to ocean model
      default_season  = 'winter'  ! default forcing data, if data is not read in
      fyear_init      = 1998      ! initial forcing year
      ycycle          = 1         ! number of years in forcing cycle
      atm_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      atm_data_type   = 'default' ! source of atmospheric forcing data
      atm_data_file   = ' '       ! atmospheric forcing data file
      precip_units    = 'mks'     ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
      oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
      wave_spec_type  = 'none'    ! type of wave spectrum forcing
      ocn_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      ocn_data_type   = 'default' ! source of ocean forcing data
      ocn_data_file   = ' '       ! ocean forcing data file
      ice_data_file   = ' '       ! ice forcing data file (opening, closing)
      lateral_flux_type   = 'uniform_ice' ! if 'uniform_ice' assume closing
                                           ! fluxes in uniform ice
      bgc_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      bgc_data_type   = 'default' ! source of BGC forcing data
      bgc_data_file   = ' '       ! biogeochemistry forcing data file
      data_dir    = ' '           ! root location of data files
      restore_ocn     = .false.   ! restore sst if true
      trestore        = 90        ! restoring timescale, days (0 instantaneous)
      snw_ssp_table   = 'test'    ! snow table type, test or snicar
      precalc_forc    = .false.   ! whether to precalculate forcing

      ! extra tracers
      tr_iage      = .false. ! ice age
      tr_FY        = .false. ! ice age
      tr_lvl       = .false. ! level ice
      tr_pond_lvl  = .false. ! level-ice melt ponds
      tr_pond_topo = .false. ! topographic melt ponds
      tr_snow      = .false. ! snow tracers (wind redistribution, metamorphosis)
      tr_aero      = .false. ! aerosols
      tr_fsd       = .false. ! floe size distribution
      tr_iso       = .false. ! isotopes

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: open file '// &
            trim(nml_filename), &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading setup_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: setup_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=setup_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: setup_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading grid_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: grid_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=grid_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: grid_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading thermo_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: thermo_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=thermo_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: thermo_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading tracer_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: tracer_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=tracer_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: tracer_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading shortwave_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: shortwave_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=shortwave_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: shortwave_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading ponds_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: ponds_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=ponds_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: ponds_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading snow_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: snow_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=snow_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: snow_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      write(nu_diag,*) subname,' Reading forcing_nml'
      rewind(unit=nu_nml, iostat=nml_error)
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: forcing_nml rewind ', &
            file=__FILE__, line=__LINE__)
      endif
      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=forcing_nml,iostat=nml_error)
      end do
      if (nml_error /= 0) then
         call icedrv_system_abort(string=subname//'ERROR: forcing_nml reading ', &
            file=__FILE__, line=__LINE__)
      endif

      close(nu_nml)

      !-----------------------------------------------------------------
      ! set up diagnostics output and resolve conflicts
      !-----------------------------------------------------------------

      write(ice_stdout,*) 'Diagnostic output will be in files '
      write(ice_stdout,*)'    ','icepack.runlog.timestamp'

      do n = 1,nx
         write(nx_names(n),'(a,i2.2)') 'point_',n
      enddo
      nx_names_default(1) = 'icefree'
      nx_names_default(2) = 'slab'
      nx_names_default(3) = 'full_ITD'
      nx_names_default(4) = 'land'
      do n = 1,nx
         nx_names(n) = nx_names_default(n)
      enddo

      do n = 1,nx
         diag_file_names=' '
         write(diag_file_names,'(a,a,a)') trim(diag_file),'.',trim(nx_names(n))
         write(ice_stdout,*)'    ',trim(diag_file_names)
         open(nu_diag_out+n-1, file=diag_file_names, status='unknown')
      end do

      write(nu_diag,*) '-----------------------------------'
      write(nu_diag,*) '  ICEPACK model diagnostic output  '
      write(nu_diag,*) '-----------------------------------'
      write(nu_diag,*) ' '

      if (ncat == 1 .and. kitd == 1) then
         write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1.'
         write (nu_diag,*) 'Use kitd = 0 (delta function ITD) with kcatbound = 0'
         write (nu_diag,*) 'or for column configurations use kcatbound = -1'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (ncat /= 1 .and. kcatbound == -1) then
         write (nu_diag,*) 'WARNING: ITD required for ncat > 1'
         write (nu_diag,*) 'WARNING: Setting kitd and kcatbound to default values'
         kitd = 1
         kcatbound = 0
      endif

      rplvl  = c0
      rptopo = c0
      if (tr_pond_lvl ) rplvl  = c1
      if (tr_pond_topo) rptopo = c1

      tr_pond = .false. ! explicit melt ponds
      if (rplvl + rptopo > puny) tr_pond = .true.

      if (rplvl + rptopo > c1 + puny) then
         write (nu_diag,*) 'WARNING: Must use only one melt pond scheme'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
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

      if (trim(shortwave(1:4)) /= 'dEdd' .and. tr_pond .and. calc_tsfc) then
         write (nu_diag,*) 'WARNING: Must use dEdd shortwave'
         write (nu_diag,*) 'WARNING: with tr_pond and calc_tsfc=T.'
         write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
         shortwave = 'dEdd'
      endif

      if (snwredist(1:3) == 'ITD' .and. .not. tr_snow) then
         write (nu_diag,*) 'WARNING: snwredist on but tr_snow=F'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (snwredist(1:4) == 'bulk' .and. .not. tr_lvl) then
         write (nu_diag,*) 'WARNING: snwredist=bulk but tr_lvl=F'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (snwredist(1:6) == 'ITDrdg' .and. .not. tr_lvl) then
         write (nu_diag,*) 'WARNING: snwredist=ITDrdg but tr_lvl=F'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (use_smliq_pnd .and. .not. snwgrain) then
         write (nu_diag,*) 'WARNING: use_smliq_pnd = T but'
         write (nu_diag,*) 'WARNING: snow metamorphosis not used'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (use_smliq_pnd .and. .not. tr_snow) then
         write (nu_diag,*) 'WARNING: use_smliq_pnd = T but'
         write (nu_diag,*) 'WARNING: snow tracers are not active'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (snwgrain .and. .not. tr_snow) then
         write (nu_diag,*) 'WARNING: snwgrain=T but tr_snow=F'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (trim(snw_aging_table) /= 'test') then
         write (nu_diag,*) 'WARNING: snw_aging_table /= test'
         write (nu_diag,*) 'WARNING: netcdf not available'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (tr_iso .and. n_iso==0) then
         write (nu_diag,*) 'WARNING: isotopes activated but'
         write (nu_diag,*) 'WARNING: not allocated in tracer array.'
         write (nu_diag,*) 'WARNING: Activate in compilation script.'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (tr_aero .and. n_aero==0) then
         write (nu_diag,*) 'WARNING: aerosols activated but'
         write (nu_diag,*) 'WARNING: not allocated in tracer array.'
         write (nu_diag,*) 'WARNING: Activate in compilation script.'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (restart_format /= 'bin' .and. restart_format /= 'nc') then
         write (nu_diag,*) 'WARNING: restart_format value unknown '//trim(restart_format)
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (history_format /= 'none' .and. history_format /= 'nc') then
         write (nu_diag,*) 'WARNING: history_format value unknown '//trim(history_format)
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (tr_aero .and. trim(shortwave(1:4)) /= 'dEdd') then
         write (nu_diag,*) 'WARNING: aerosols activated but dEdd'
         write (nu_diag,*) 'WARNING: shortwave is not.'
         write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
         shortwave = 'dEdd'
      endif

      if (snwgrain .and. trim(shortwave(1:4)) /= 'dEdd') then
         write (nu_diag,*) 'WARNING: snow grain radius activated but'
         write (nu_diag,*) 'WARNING: dEdd shortwave is not.'
      endif

      if (snwredist(1:4) /= 'none' .and. trim(shortwave(1:4)) /= 'dEdd') then
         write (nu_diag,*) 'WARNING: snow redistribution activated but'
         write (nu_diag,*) 'WARNING: dEdd shortwave is not.'
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

      if (ktherm == 1 .and. trim(saltflux_option) /= 'constant') then
         write (nu_diag,*) &
         'WARNING: ktherm = 1 and saltflux_option = ',trim(saltflux_option)
         write (nu_diag,*) &
         'WARNING: For consistency, set saltflux_option = constant'
      endif

      if (ktherm == 0) then
         write (nu_diag,*) 'WARNING: ktherm = 0 zero-layer thermodynamics'
         write (nu_diag,*) 'WARNING: has been deprecated'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (formdrag) then
      if (trim(atmbndy) == 'constant') then
         write (nu_diag,*) 'WARNING: atmbndy = constant not allowed with formdrag'
         write (nu_diag,*) 'WARNING: Setting atmbndy = similarity'
         atmbndy = 'similarity'
      endif

      if (.not. calc_strair) then
         write (nu_diag,*) 'WARNING: formdrag=T but calc_strair=F'
         write (nu_diag,*) 'WARNING: Setting calc_strair=T'
         calc_strair = .true.
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

      wave_spec = .false.
      if (tr_fsd .and. (trim(wave_spec_type) /= 'none')) wave_spec = .true.
      if (tr_fsd .and. (trim(wave_spec_type) == 'none')) then
         write (nu_diag,*) 'WARNING: tr_fsd=T but wave_spec=F - not recommended'
      end if

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
         write(nu_diag,1030) ' dumpfreq                  = ', trim(dumpfreq)
         write(nu_diag,1010) ' dump_last                 = ', dump_last
         write(nu_diag,1010) ' restart                   = ', restart
         write(nu_diag,1030) ' restart_dir               = ', trim(restart_dir)
         write(nu_diag,1030) ' restart_file              = ', trim(restart_file)
         write(nu_diag,1030) ' restart_format            = ', trim(restart_format)
         write(nu_diag,1030) ' history_format            = ', trim(history_format)
         write(nu_diag,1030) ' ice_ic                    = ', trim(ice_ic)
         write(nu_diag,1005) ' hi_init_slab              = ', hi_init_slab
         write(nu_diag,1005) ' hsno_init_slab            = ', hsno_init_slab
         write(nu_diag,1005) ' hbar_init_itd             = ', hbar_init_itd
         write(nu_diag,1005) ' hsno_init_itd             = ', hsno_init_itd
         write(nu_diag,1005) ' sst_init                  = ', sst_init
         write(nu_diag,1010) ' conserv_check             = ', conserv_check
         write(nu_diag,1020) ' kitd                      = ', kitd
         write(nu_diag,1020) ' kcatbound                 = ', kcatbound
         write(nu_diag,1020) ' ndtd                      = ', ndtd
         write(nu_diag,1020) ' kstrength                 = ', kstrength
         write(nu_diag,1020) ' krdg_partic               = ', krdg_partic
         write(nu_diag,1020) ' krdg_redist               = ', krdg_redist
         if (krdg_redist == 1) &
         write(nu_diag,1000) ' mu_rdg                    = ', mu_rdg
         if (kstrength == 1) &
         write(nu_diag,1000) ' Cf                        = ', Cf
         write(nu_diag,1000) ' ksno                      = ', ksno
         write(nu_diag,1030) ' shortwave                 = ', trim(shortwave)
         if (cpl_bgc) then
             write(nu_diag,1000) ' BGC coupling is switched ON'
         else
             write(nu_diag,1000) ' BGC coupling is switched OFF'
         endif

         if (trim(shortwave(1:4)) == 'dEdd') then
         write(nu_diag,1000) ' R_ice                     = ', R_ice
         write(nu_diag,1000) ' R_pnd                     = ', R_pnd
         write(nu_diag,1000) ' R_snw                     = ', R_snw
         write(nu_diag,1000) ' dT_mlt                    = ', dT_mlt
         write(nu_diag,1000) ' rsnw_mlt                  = ', rsnw_mlt
         write(nu_diag,1000) ' kalg                      = ', kalg
         write(nu_diag,1000) ' hp1                       = ', hp1
         write(nu_diag,1000) ' hs0                       = ', hs0
         else
         write(nu_diag,1030) ' albedo_type               = ', trim(albedo_type)
         write(nu_diag,1000) ' albicev                   = ', albicev
         write(nu_diag,1000) ' albicei                   = ', albicei
         write(nu_diag,1000) ' albsnowv                  = ', albsnowv
         write(nu_diag,1000) ' albsnowi                  = ', albsnowi
         write(nu_diag,1000) ' ahmax                     = ', ahmax
         endif

         if (trim(shortwave) == 'dEdd_snicar_ad') then
         write(nu_diag,1030) ' snw_ssp_table             = ', trim(snw_ssp_table)
         endif

         write(nu_diag,1010) ' sw_redist                 = ', sw_redist
         write(nu_diag,1005) ' sw_frac                   = ', sw_frac
         write(nu_diag,1005) ' sw_dtemp                  = ', sw_dtemp

         write(nu_diag,1000) ' rfracmin                  = ', rfracmin
         write(nu_diag,1000) ' rfracmax                  = ', rfracmax
         if (tr_pond_lvl) then
         write(nu_diag,1000) ' hs1                       = ', hs1
         write(nu_diag,1000) ' dpscale                   = ', dpscale
         write(nu_diag,1030) ' frzpnd                    = ', trim(frzpnd)
         endif
         if (tr_pond .and. .not. tr_pond_lvl) &
         write(nu_diag,1000) ' pndaspect                 = ', pndaspect

         if (tr_snow) then
         write(nu_diag,1030) ' snwredist                 = ', trim(snwredist)
         write(nu_diag,1010) ' snwgrain                  = ', snwgrain
         write(nu_diag,1010) ' use_smliq_pnd             = ', use_smliq_pnd
         write(nu_diag,1030) ' snw_aging_table           = ', trim(snw_aging_table)
         write(nu_diag,1000) ' rsnw_fall                 = ', rsnw_fall
         write(nu_diag,1000) ' rsnw_tmax                 = ', rsnw_tmax
         write(nu_diag,1000) ' rhosnew                   = ', rhosnew
         write(nu_diag,1000) ' rhosmin                   = ', rhosmin
         write(nu_diag,1000) ' rhosmax                   = ', rhosmax
         write(nu_diag,1000) ' windmin                   = ', windmin
         write(nu_diag,1000) ' drhosdwind                = ', drhosdwind
         write(nu_diag,1000) ' snwlvlfac                 = ', snwlvlfac
         endif

         write(nu_diag,1020) ' ktherm                    = ', ktherm
         if (ktherm == 1) &
         write(nu_diag,1030) ' conduct                   = ', trim(conduct)
         write(nu_diag,1005) ' emissivity                = ', emissivity
         if (ktherm == 2) then
         write(nu_diag,1005) ' a_rapid_mode              = ', a_rapid_mode
         write(nu_diag,1005) ' Rac_rapid_mode            = ', Rac_rapid_mode
         write(nu_diag,1005) ' aspect_rapid_mode         = ', aspect_rapid_mode
         write(nu_diag,1005) ' dSdt_slow_mode            = ', dSdt_slow_mode
         write(nu_diag,1005) ' phi_c_slow_mode           = ', phi_c_slow_mode
         write(nu_diag,1005) ' phi_i_mushy               = ', phi_i_mushy
         write(nu_diag,1005) ' Tliquidus_max             = ', Tliquidus_max
         endif

         write(nu_diag,1030) ' atmbndy                   = ', trim(atmbndy)
         write(nu_diag,1010) ' formdrag                  = ', formdrag
         write(nu_diag,1010) ' highfreq                  = ', highfreq
         write(nu_diag,1020) ' natmiter                  = ', natmiter
         write(nu_diag,1005) ' atmiter_conv              = ', atmiter_conv
         write(nu_diag,1010) ' calc_strair               = ', calc_strair
         write(nu_diag,1010) ' calc_Tsfc                 = ', calc_Tsfc
         write(nu_diag,1010) ' calc_dragio               = ', calc_dragio
         write(nu_diag,1005) ' floediam                  = ', floediam
         write(nu_diag,1005) ' hfrazilmin                = ', hfrazilmin

         write(nu_diag,1030) ' atm_data_type             = ', trim(atm_data_type)
         write(nu_diag,1030) ' ocn_data_type             = ', trim(ocn_data_type)
         write(nu_diag,1030) ' bgc_data_type             = ', trim(bgc_data_type)

         write(nu_diag,*)    '  lateral_flux_type        = ', trim(lateral_flux_type)

         write(nu_diag,1030) ' atm_data_file             = ', trim(atm_data_file)
         write(nu_diag,1030) ' ocn_data_file             = ', trim(ocn_data_file)
         write(nu_diag,1030) ' bgc_data_file             = ', trim(bgc_data_file)
         write(nu_diag,1030) ' ice_data_file             = ', trim(ice_data_file)
         write(nu_diag,1010) ' precalc_forc              = ', precalc_forc

         if (trim(atm_data_type)=='default') &
         write(nu_diag,1030) ' default_season            = ', trim(default_season)

         write(nu_diag,1030) ' cpl_frazil                = ', trim(cpl_frazil)
         write(nu_diag,1010) ' update_ocn_f              = ', update_ocn_f
         write(nu_diag,1010) ' wave_spec                 = ', wave_spec
         if (wave_spec) &
         write(nu_diag,1030) ' wave_spec_type            = ', trim(wave_spec_type)
         write(nu_diag,1010) ' l_mpond_fresh             = ', l_mpond_fresh
         write(nu_diag,1005) ' ustar_min                 = ', ustar_min
         write(nu_diag,1005) ' hi_min                    = ', hi_min
         write(nu_diag,1030) ' fbot_xfer_type            = ', trim(fbot_xfer_type)
         write(nu_diag,1010) ' oceanmixed_ice            = ', oceanmixed_ice
         write(nu_diag,1030) ' congel_freeze             = ', trim(congel_freeze)
         write(nu_diag,1030) ' tfrz_option               = ', trim(tfrz_option)
         write(nu_diag,*)    ' saltflux_option           = ', &
                               trim(saltflux_option)
         if (trim(saltflux_option) == 'constant') then
            write(nu_diag,1005)    ' ice_ref_salinity          = ', ice_ref_salinity
         endif
         write(nu_diag,1010) ' restore_ocn               = ', restore_ocn
         if (restore_ocn) &
         write(nu_diag,1005) ' trestore                  = ', trestore

         ! tracers
         write(nu_diag,1010) ' tr_iage                   = ', tr_iage
         write(nu_diag,1010) ' tr_FY                     = ', tr_FY
         write(nu_diag,1010) ' tr_lvl                    = ', tr_lvl
         write(nu_diag,1010) ' tr_pond_lvl               = ', tr_pond_lvl
         write(nu_diag,1010) ' tr_pond_topo              = ', tr_pond_topo
         write(nu_diag,1010) ' tr_snow                   = ', tr_snow
         write(nu_diag,1010) ' tr_aero                   = ', tr_aero
         write(nu_diag,1010) ' tr_fsd                    = ', tr_fsd

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
             nt_alvl = ntrcr   ! area of level ice
             ntrcr = ntrcr + 1
             nt_vlvl = ntrcr   ! volume of level ice
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

         nt_smice = max_ntrcr
         nt_smliq = max_ntrcr
         nt_rhos  = max_ntrcr
         nt_rsnw  = max_ntrcr
         if (tr_snow) then
            nt_smice = ntrcr + 1
            ntrcr = ntrcr + nslyr     ! mass of ice in nslyr snow layers
            nt_smliq = ntrcr + 1
            ntrcr = ntrcr + nslyr     ! mass of liquid in nslyr snow layers
            nt_rhos = ntrcr + 1
            ntrcr = ntrcr + nslyr     ! snow density in nslyr layers
            nt_rsnw = ntrcr + 1
            ntrcr = ntrcr + nslyr     ! snow grain radius in nslyr layers
         endif

         nt_fsd = max_ntrcr
         if (tr_fsd) then
             nt_fsd = ntrcr + 1       ! floe size distribution
             ntrcr = ntrcr + nfsd
         end if

         nt_isosno = max_ntrcr
         nt_isoice = max_ntrcr
         if (tr_iso) then             ! isotopes
             nt_isosno = ntrcr + 1
             ntrcr = ntrcr + n_iso    ! n_iso species in snow
             nt_isoice = ntrcr + 1
             ntrcr = ntrcr + n_iso    ! n_iso species in ice
         end if

         nt_aero = max_ntrcr - 4*n_aero
         if (tr_aero) then
             nt_aero = ntrcr + 1
             ntrcr = ntrcr + 4*n_aero ! 4 dEdd layers, n_aero species
         endif

         if (ntrcr > max_ntrcr-1) then
            write(nu_diag,*) 'max_ntrcr-1 < number of namelist tracers'
            write(nu_diag,*) 'max_ntrcr-1 = ',max_ntrcr-1,' ntrcr = ',ntrcr
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         endif

         write(nu_diag,*) ' '
         write(nu_diag,1020) 'max_ntrcr = ', max_ntrcr
         write(nu_diag,1020) 'ntrcr = ', ntrcr
         write(nu_diag,*) ' '
         write(nu_diag,1020) 'nt_sice = ', nt_sice
         write(nu_diag,1020) 'nt_qice = ', nt_qice
         write(nu_diag,1020) 'nt_qsno = ', nt_qsno
         write(nu_diag,*)' '
         write(nu_diag,1020) 'ncat    = ', ncat
         write(nu_diag,1020) 'nilyr   = ', nilyr
         write(nu_diag,1020) 'nslyr   = ', nslyr
         write(nu_diag,1020) 'nblyr   = ', nblyr
         write(nu_diag,1020) 'nfsd    = ', nfsd
         write(nu_diag,1020) 'n_iso   = ', n_iso
         write(nu_diag,1020) 'n_aero  = ', n_aero
         write(nu_diag,*)' '
         write(nu_diag,1020) 'nx      = ', nx
         write(nu_diag,*)' '

         call icedrv_system_flush(nu_diag)

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f10.6) ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,    a)    ! character
 1040    format (a30,2x,6i6)   ! integer
 1050    format (a30,2x,6a6)   ! character

      if (formdrag) then
         if (nt_apnd==0) then
            write(nu_diag,*)'ERROR: nt_apnd:',nt_apnd
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         elseif (nt_hpnd==0) then
            write(nu_diag,*)'ERROR: nt_hpnd:',nt_hpnd
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         elseif (nt_ipnd==0) then
            write(nu_diag,*)'ERROR: nt_ipnd:',nt_ipnd
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         elseif (nt_alvl==0) then
            write(nu_diag,*)'ERROR: nt_alvl:',nt_alvl
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         elseif (nt_vlvl==0) then
            write(nu_diag,*)'ERROR: nt_vlvl:',nt_vlvl
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         endif
      endif

      !-----------------------------------------------------------------
      ! set Icepack values
      !-----------------------------------------------------------------

      call icepack_init_parameters(ustar_min_in=ustar_min, Cf_in=Cf, &
           albicev_in=albicev, albicei_in=albicei, ksno_in=ksno, &
           albsnowv_in=albsnowv, albsnowi_in=albsnowi, hi_min_in=hi_min, &
           natmiter_in=natmiter, ahmax_in=ahmax, shortwave_in=shortwave, &
           atmiter_conv_in = atmiter_conv, calc_dragio_in=calc_dragio, &
           albedo_type_in=albedo_type, R_ice_in=R_ice, R_pnd_in=R_pnd, &
           R_snw_in=R_snw, dT_mlt_in=dT_mlt, rsnw_mlt_in=rsnw_mlt, &
           kstrength_in=kstrength, krdg_partic_in=krdg_partic, &
           krdg_redist_in=krdg_redist, mu_rdg_in=mu_rdg, &
           atmbndy_in=atmbndy, calc_strair_in=calc_strair, &
           formdrag_in=formdrag, highfreq_in=highfreq, &
           emissivity_in=emissivity, &
           kitd_in=kitd, kcatbound_in=kcatbound, hs0_in=hs0, &
           dpscale_in=dpscale, frzpnd_in=frzpnd, &
           rfracmin_in=rfracmin, rfracmax_in=rfracmax, &
           pndaspect_in=pndaspect, hs1_in=hs1, hp1_in=hp1, &
           floediam_in=floediam, hfrazilmin_in=hfrazilmin, &
           ktherm_in=ktherm, calc_Tsfc_in=calc_Tsfc, &
           conduct_in=conduct, a_rapid_mode_in=a_rapid_mode, &
           update_ocn_f_in=update_ocn_f, cpl_frazil_in=cpl_frazil, &
           Rac_rapid_mode_in=Rac_rapid_mode, &
           aspect_rapid_mode_in=aspect_rapid_mode, &
           dSdt_slow_mode_in=dSdt_slow_mode, &
           phi_c_slow_mode_in=phi_c_slow_mode, Tliquidus_max_in=Tliquidus_max, &
           phi_i_mushy_in=phi_i_mushy, conserv_check_in=conserv_check, &
           congel_freeze_in=congel_freeze, &
           tfrz_option_in=tfrz_option, saltflux_option_in=saltflux_option, &
           ice_ref_salinity_in=ice_ref_salinity, kalg_in=kalg, &
           fbot_xfer_type_in=fbot_xfer_type, &
           wave_spec_type_in=wave_spec_type, wave_spec_in=wave_spec, &
           sw_redist_in=sw_redist, sw_frac_in=sw_frac, sw_dtemp_in=sw_dtemp, &
           snwredist_in=snwredist, use_smliq_pnd_in=use_smliq_pnd, &
           snw_aging_table_in=snw_aging_table, &
           snwgrain_in=snwgrain, rsnw_fall_in=rsnw_fall, rsnw_tmax_in=rsnw_tmax, &
           rhosnew_in=rhosnew, rhosmin_in=rhosmin, rhosmax_in=rhosmax, &
           windmin_in=windmin, drhosdwind_in=drhosdwind, snwlvlfac_in=snwlvlfac)
      call icepack_init_tracer_sizes(ntrcr_in=ntrcr, &
           ncat_in=ncat, nilyr_in=nilyr, nslyr_in=nslyr, nblyr_in=nblyr, &
           nfsd_in=nfsd, n_iso_in=n_iso, n_aero_in=n_aero)
      call icepack_init_tracer_flags(tr_iage_in=tr_iage, &
           tr_FY_in=tr_FY, tr_lvl_in=tr_lvl, tr_aero_in=tr_aero, &
           tr_iso_in=tr_iso, tr_snow_in=tr_snow, &
           tr_pond_in=tr_pond, &
           tr_pond_lvl_in=tr_pond_lvl, &
           tr_pond_topo_in=tr_pond_topo, tr_fsd_in=tr_fsd)
      call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, &
           nt_sice_in=nt_sice, nt_qice_in=nt_qice, &
           nt_qsno_in=nt_qsno, nt_iage_in=nt_iage, &
           nt_fy_in=nt_fy, nt_alvl_in=nt_alvl, nt_vlvl_in=nt_vlvl, &
           nt_apnd_in=nt_apnd, nt_hpnd_in=nt_hpnd, nt_ipnd_in=nt_ipnd, &
           nt_smice_in=nt_smice, nt_smliq_in=nt_smliq, &
           nt_rhos_in=nt_rhos, nt_rsnw_in=nt_rsnw, &
           nt_aero_in=nt_aero, nt_fsd_in=nt_fsd, &
           nt_isosno_in=nt_isosno, nt_isoice_in=nt_isoice)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      end subroutine input_data

!=======================================================================

! Horizontal grid initialization:
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_grid2

      integer :: i
      real (kind=dbl_kind) :: pi, puny
      character(len=*), parameter :: subname='(init_grid2)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(pi_out=pi,puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      TLAT(:) = p5*pi  ! pi/2, North pole
      TLON(:) = c0

      do i = 2, nx
         TLAT(i) = TLAT(i-1) - p5*pi/180._dbl_kind ! half-deg increments
      enddo

      tmask(:) = .true.

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
      use icedrv_domain_size, only: ncat, nilyr, nslyr, nblyr, max_ntrcr
      use icedrv_domain_size, only: n_iso, n_aero, nfsd
      use icedrv_flux, only: sst, Tf, Tair, salinz, Tmltz
      use icedrv_state, only: trcr_depend, aicen, trcrn, vicen, vsnon
      use icedrv_state, only: aice0, aice, vice, vsno, trcr, aice_init
      use icedrv_state, only: n_trcr_strata, nt_strata, trcr_base

      integer (kind=int_kind) :: &
         i           , & ! horizontal indes
         k           , & ! vertical index
         it              ! tracer index

      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_aero, tr_fsd, tr_iso
      logical (kind=log_kind) :: tr_pond_lvl, tr_pond_topo, tr_snow
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_fy
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
                                 nt_smice, nt_smliq, nt_rhos, nt_rsnw, &
                                 nt_aero, nt_fsd, nt_isosno, nt_isoice

      character(len=*), parameter :: subname='(init_state)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
         call icepack_query_tracer_flags(tr_iage_out=tr_iage, &
              tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, &
              tr_iso_out=tr_iso, tr_snow_out=tr_snow, &
              tr_pond_lvl_out=tr_pond_lvl, &
              tr_pond_topo_out=tr_pond_topo, tr_fsd_out=tr_fsd)
         call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, &
              nt_sice_out=nt_sice, nt_qice_out=nt_qice, &
              nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_fy_out=nt_fy, &
              nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, &
              nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
              nt_ipnd_out=nt_ipnd, &
              nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, &
              nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw, &
              nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice, &
              nt_aero_out=nt_aero, nt_fsd_out=nt_fsd)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Check number of layers in ice and snow.
      !-----------------------------------------------------------------
         if (nilyr < 1) then
            write (nu_diag,*) 'nilyr =', nilyr
            write (nu_diag,*) 'Must have at least one ice layer'
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         endif

         if (nslyr < 1) then
            write (nu_diag,*) 'nslyr =', nslyr
            write (nu_diag,*) 'Must have at least one snow layer'
            call icedrv_system_abort(file=__FILE__,line=__LINE__)
         endif

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
      if (tr_snow) then
         do k = 1, nslyr
            trcr_depend(nt_smice + k - 1) = 2          ! ice mass in snow
            trcr_depend(nt_smliq + k - 1) = 2          ! liquid mass in snow
            trcr_depend(nt_rhos  + k - 1) = 2          ! effective snow density
            trcr_depend(nt_rsnw  + k - 1) = 2          ! snow radius
         enddo
      endif
      if (tr_fsd) then
         do it = 1, nfsd
            trcr_depend(nt_fsd + it - 1) = 0    ! area-weighted floe size distribution
         enddo
      endif
      if (tr_iso) then  ! isotopes
         do it = 1, n_iso
            trcr_depend(nt_isosno) = 2          ! snow
            trcr_depend(nt_isoice) = 1          ! ice
         enddo
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

      call set_state_var (nx, &
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
         call icepack_aggregate(trcrn=trcrn(i,1:ntrcr,:),     &
                                aicen=aicen(i,:),             &
                                vicen=vicen(i,:),             &
                                vsnon=vsnon(i,:),             &
                                trcr=trcr (i,1:ntrcr),        &
                                aice=aice (i),                &
                                vice=vice (i),                &
                                vsno=vsno (i),                &
                                aice0=aice0(i),               &
                                trcr_depend=trcr_depend(1:ntrcr),     &
                                trcr_base=trcr_base    (1:ntrcr,:),   &
                                n_trcr_strata=n_trcr_strata(1:ntrcr), &
                                nt_strata=nt_strata    (1:ntrcr,:), &
                                Tf = Tf(i))

         aice_init(i) = aice(i)

      enddo
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      end subroutine init_state

!=======================================================================

! Initialize state in each ice thickness category
!
! authors: Elizabeth Hunke, LANL

      subroutine set_state_var (nx, &
                                Tair,     sst,    &
                                Tf,               &
                                salinz,   Tmltz,  &
                                aicen,    trcrn,  &
                                vicen,    vsnon)

      use icedrv_arrays_column, only: hin_max
      use icedrv_domain_size, only: nilyr, nslyr, max_ntrcr, ncat, nfsd

      integer (kind=int_kind), intent(in) :: &
         nx          ! number of grid cells

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         Tair       ! air temperature  (K)

      ! ocean values may be redefined here, unlike in CICE
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         Tf     , & ! freezing temperature (C)
         sst        ! sea surface temperature (C)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         salinz , & ! initial salinity profile
         Tmltz      ! initial melting temperature profile

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i     , & ! horizontal indices
         k     , & ! ice layer index
         n     , & ! thickness category index
         it        ! tracer index

      real (kind=dbl_kind) :: &
         Tsfc, sum, hbar, &
         rhos, Lfresh, puny, rsnw_fall

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      logical (kind=log_kind) :: tr_brine, tr_lvl, tr_fsd, tr_snow
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice, nt_fsd
      integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl
      integer (kind=int_kind) :: nt_rhos, nt_rsnw, nt_smice, nt_smliq

      character(len=*), parameter :: subname='(set_state_var)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl, &
           tr_fsd_out=tr_fsd, tr_snow_out=tr_snow)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fsd_out=nt_fsd, &
           nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, &
           nt_rsnw_out=nt_rsnw, nt_rhos_out=nt_rhos, &
           nt_smice_out=nt_smice, nt_smliq_out=nt_smliq)
      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh, puny_out=puny, &
           rsnw_fall_out=rsnw_fall)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Initialize state variables.
      ! If restarting, these values are overwritten.
      !-----------------------------------------------------------------

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
      if (i <= nx) then ! need to set ocean parameters
         sst(i) = sst_init
      endif

      !-----------------------------------------------------------------

      i = 2  ! 100% ice concentration slab, thickness and snow from namelist
      if (i <= nx) then
         sst(i) = sst_init ! initial ocean temperature
         do n = 1, ncat
            if (hi_init_slab <= hin_max(n)) then
               ainit(n) = c1
               hinit(n) = hi_init_slab
               exit
            endif
         enddo
         if (hi_init_slab > hin_max(ncat)) then
            ainit(ncat) = c1
            hinit(ncat) = hi_init_slab
         endif
      do n = 1, ncat
         ! ice volume, snow volume
         aicen(i,n) = ainit(n)
         vicen(i,n) = hinit(n) * ainit(n) ! m
         vsnon(i,n) = hsno_init_slab * ainit(n)
         ! tracers
         call icepack_init_enthalpy(Tair = Tair(i),     &
                                Tf       = Tf(i),       &
                                Sprofile = salinz(i,:), &
                                Tprofile = Tmltz(i,:),  &
                                Tsfc     = Tsfc,        &
                                qin=qin(:), qsn=qsn(:))

         ! floe size distribution
         if (tr_fsd) call icepack_init_fsd(ice_ic=ice_ic, &
                                  afsd=trcrn(i,nt_fsd:nt_fsd+nfsd-1,n))
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
         ! snow radius, effective density, ice and liquid mass content
         if (tr_snow) then
            do k = 1, nslyr
               trcrn(i,nt_rsnw +k-1,n) = rsnw_fall
               trcrn(i,nt_rhos +k-1,n) = rhos
               trcrn(i,nt_smice+k-1,n) = rhos
               trcrn(i,nt_smliq+k-1,n) = c0
            enddo               ! nslyr
         endif
      enddo                  ! ncat
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      endif  ! (i <= nx)
      !-----------------------------------------------------------------

      i = 3  ! full thickness distribution
      if (i <= nx) then
      sst(i) = sst_init
      ! initial category areas in cells with ice
      hbar = hbar_init_itd  ! initial ice thickness with greatest area
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
         vsnon(i,n) = min(aicen(i,n)*hsno_init_itd,p2*vicen(i,n))
         ! tracers
         call icepack_init_enthalpy(Tair = Tair(i),     &
                                Tf       = Tf(i),       &
                                Sprofile = salinz(i,:), &
                                Tprofile = Tmltz(i,:),  &
                                Tsfc     = Tsfc,        &
                                qin=qin(:), qsn=qsn(:))
         ! floe size distribution
         if (tr_fsd) call icepack_init_fsd(ice_ic=ice_ic, &
                                  afsd=trcrn(i,nt_fsd:nt_fsd+nfsd-1,n))

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
         ! snow radius, effective density, ice and liquid mass content
         if (tr_snow) then
            do k = 1, nslyr
               trcrn(i,nt_rsnw +k-1,n) = rsnw_fall
               trcrn(i,nt_rhos +k-1,n) = rhos
               trcrn(i,nt_smice+k-1,n) = rhos
               trcrn(i,nt_smliq+k-1,n) = c0
            enddo               ! nslyr
         endif
      enddo                  ! ncat
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      endif  ! (i <= nx)

      !-----------------------------------------------------------------

      ! land
      ! already initialized above (tmask = 0)
      i = 4
      if (i <= nx) then
         tmask(i) = .false.
         sst(i) = c0
         Tf(i) = c0
      endif

      end subroutine set_state_var

!=======================================================================

!  Initialize floe size distribution tracer (call prior to reading restart data)

      subroutine init_fsd

      use icedrv_arrays_column, only: wavefreq, dwavefreq, wave_sig_ht, &
         wave_spectrum, d_afsd_newi, d_afsd_latg, d_afsd_latm, &
         d_afsd_wave, d_afsd_weld

      wavefreq       (:)   = c0
      dwavefreq      (:)   = c0
      wave_sig_ht    (:)   = c0
      wave_spectrum  (:,:) = c0
      d_afsd_newi    (:,:) = c0
      d_afsd_latg    (:,:) = c0
      d_afsd_latm    (:,:) = c0
      d_afsd_wave    (:,:) = c0
      d_afsd_weld    (:,:) = c0

      end subroutine init_fsd

!=======================================================================

  end module icedrv_init

!=======================================================================
