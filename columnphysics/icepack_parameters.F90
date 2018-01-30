!  SVN:$Id: icepack_parameters.F90 1226 2017-05-22 22:45:03Z tcraig $
!=========================================================================
!
! flags for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_parameters

      use icepack_kinds
      use icepack_constants, only: c3, c0, c1, p5, p1
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private

      public :: icepack_init_parameters
      public :: icepack_query_parameters
      public :: icepack_write_parameters

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         ktherm          ! type of thermodynamics
                         ! 0 = 0-layer approximation
                         ! 1 = Bitz and Lipscomb 1999
                         ! 2 = mushy layer theory

      character (char_len), public :: &
         conduct, &      ! 'MU71' or 'bubbly'
         fbot_xfer_type  ! transfer coefficient type for ice-ocean heat flux

      logical (kind=log_kind), public :: &
         heat_capacity, &! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics
         calc_Tsfc   ,  &! if true, calculate surface temperature
                         ! if false, Tsfc is computed elsewhere and
                         ! atmos-ice fluxes are provided to CICE
         solve_zsal  ,  &! if true, update salinity profile from solve_S_dt
         modal_aero      ! if true, use modal aerosal optical properties
                         ! only for use with tr_aero or tr_zaero

      real (kind=dbl_kind), public :: &
         dts_b,   &      ! zsalinity timestep
         ustar_min       ! minimum friction velocity for ice-ocean heat flux

      ! mushy thermo
      real(kind=dbl_kind), public :: &
         a_rapid_mode      , & ! channel radius for rapid drainage mode (m)
         Rac_rapid_mode    , & ! critical Rayleigh number for rapid drainage mode
         aspect_rapid_mode , & ! aspect ratio for rapid drainage mode (larger=wider)
         dSdt_slow_mode    , & ! slow mode drainage strength (m s-1 K-1)
         phi_c_slow_mode   , & ! liquid fraction porosity cutoff for slow mode
         phi_i_mushy           ! liquid fraction of congelation ice


!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

      character (len=char_len), public :: &
         shortwave, & ! shortwave method, 'ccsm3' or 'dEdd'
         albedo_type  ! albedo parameterization, 'ccsm3' or 'constant'
                      ! shortwave='dEdd' overrides this parameter

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind), public :: &
         albicev  , & ! visible ice albedo for h > ahmax
         albicei  , & ! near-ir ice albedo for h > ahmax
         albsnowv , & ! cold snow albedo, visible
         albsnowi , & ! cold snow albedo, near IR
         ahmax        ! thickness above which ice albedo is constant (m)

      ! dEdd tuning parameters, set in namelist
      real (kind=dbl_kind), public :: &
         R_ice    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt   , & ! change in temp for non-melt to melt snow grain 
                      ! radius change (C)
         rsnw_mlt , & ! maximum melting snow grain radius (10^-6 m)
         kalg         ! algae absorption coefficient for 0.5 m thick layer

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: & ! defined in namelist 
         kstrength  , & ! 0 for simple Hibler (1979) formulation 
                        ! 1 for Rothrock (1975) pressure formulation 
         krdg_partic, & ! 0 for Thorndike et al. (1975) formulation 
                        ! 1 for exponential participation function 
         krdg_redist    ! 0 for Hibler (1980) formulation 
                        ! 1 for exponential redistribution function 

      real (kind=dbl_kind), public :: &  
         mu_rdg         ! gives e-folding scale of ridged ice (m^.5) 
                        ! (krdg_redist = 1) 

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

      character (len=char_len), public :: &
         atmbndy ! atmo boundary method, 'default' ('ccsm3') or 'constant'

      logical (kind=log_kind), public :: &
         calc_strair, &  ! if true, calculate wind stress components
         formdrag,    &  ! if true, calculate form drag
         highfreq        ! if true, use high frequency coupling

      integer (kind=int_kind), public :: &
         natmiter        ! number of iterations for boundary layer calculations

      character(len=char_len), public :: &
         tfrz_option              ! form of ocean freezing temperature
                                  ! 'minus1p8' = -1.8 C
                                  ! 'linear_salt' = -depressT * sss
                                  ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         kitd        , & ! type of itd conversions
                         !   0 = delta function
                         !   1 = linear remap
         kcatbound       !   0 = old category boundary formula
                         !   1 = new formula giving round numbers
                         !   2 = WMO standard
                         !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         hs0             ! snow depth for transition to bare sea ice (m)

      ! level-ice ponds
      character (len=char_len), public :: &
         frzpnd          ! pond refreezing parameterization

      real (kind=dbl_kind), public :: &
         dpscale, &      ! alter e-folding time scale for flushing 
         rfracmin, &     ! minimum retained fraction of meltwater
         rfracmax, &     ! maximum retained fraction of meltwater
         pndaspect, &    ! ratio of pond depth to pond fraction
         hs1             ! tapering parameter for snow on pond ice

      ! topo ponds
      real (kind=dbl_kind), public :: &
         hp1             ! critical parameter for pond ice thickness

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

      character(char_len), public :: &          
         bgc_flux_type      ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006' 

      logical (kind=log_kind), public :: &
         z_tracers,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae        ! if .true., algal absorption of Shortwave is computed in the

      logical (kind=log_kind), public :: & 
         skl_bgc         ! if true, solve skeletal biochemistry

      real (kind=dbl_kind), public :: & 
         grid_o      , & ! for bottom flux        
         l_sk        , & ! characteristic diffusive scale (zsalinity) (m)
         phi_snow    , & ! porosity of snow
         initbio_frac    ! fraction of ocean tracer concentration used to initialize tracer 

      real (kind=dbl_kind), public :: & 
         grid_oS     , & ! for bottom flux (zsalinity)
         l_skS           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)

      real (kind=dbl_kind), public :: &
         fr_resp           , &   ! fraction of algal growth lost due to respiration        
         algal_vel         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol        , &   ! solubility fraction
         frazil_scav             ! fraction or multiple of bgc concentrated in frazil ice

      !-----------------------------------------------------------------
      ! From algal_dyn in icepack_algae.F90
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         T_max            , & ! maximum temperature (C)
         fsal             , & ! Salinity limitation (ppt)
         op_dep_min       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s       , & ! fraction of grazing spilled or slopped
         fr_graze_e       , & ! fraction of assimilation excreted 
         fr_mort2min      , & ! fractionation of mortality to Am
         fr_dFe           , & ! fraction of remineralized nitrogen (in units of algal iron)
         k_nitrif         , & ! nitrification rate (1/day)           
         t_iron_conv      , & ! desorption loss pFe to dFe (day)
         max_loss         , & ! restrict uptake to % of remaining value 
         max_dfe_doc1     , & ! max ratio of dFe to saccharides in the ice (nM Fe/muM C)    
         fr_resp_s        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS         , & ! fraction conversion given high yield
         t_sk_conv        , & ! Stefels conversion time (d)
         t_sk_ox              ! DMS oxidation time (d)

!=======================================================================

      contains

!=======================================================================

! subroutine to set the column package internal parameters

      subroutine icepack_init_parameters(   &
           ktherm_in, conduct_in, fbot_xfer_type_in, calc_Tsfc_in, dts_b_in, ustar_min_in, a_rapid_mode_in, &
           Rac_rapid_mode_in, aspect_rapid_mode_in, dSdt_slow_mode_in, phi_c_slow_mode_in, &
           phi_i_mushy_in, shortwave_in, albedo_type_in, albicev_in, albicei_in, albsnowv_in, &
           albsnowi_in, ahmax_in, R_ice_in, R_pnd_in, R_snw_in, dT_mlt_in, rsnw_mlt_in, &
           kalg_in, kstrength_in, krdg_partic_in, krdg_redist_in, mu_rdg_in, &
           atmbndy_in, calc_strair_in, formdrag_in, highfreq_in, natmiter_in, &
           tfrz_option_in, kitd_in, kcatbound_in, hs0_in, frzpnd_in, &
           dpscale_in, rfracmin_in, rfracmax_in, pndaspect_in, hs1_in, hp1_in, &
           bgc_flux_type_in, z_tracers_in, scale_bgc_in, solve_zbgc_in, dEdd_algae_in, &
           modal_aero_in, skl_bgc_in, solve_zsal_in, grid_o_in, l_sk_in, &
           initbio_frac_in, grid_oS_in, l_skS_in, &
           phi_snow_in, heat_capacity_in, &
           fr_resp_in, algal_vel_in, R_dFe2dust_in, dustFe_sol_in, &
           T_max_in, fsal_in, op_dep_min_in, fr_graze_s_in, fr_graze_e_in, fr_mort2min_in, &
           fr_dFe_in, k_nitrif_in, t_iron_conv_in, max_loss_in, max_dfe_doc1_in, fr_resp_s_in, &
           y_sk_DMS_in, t_sk_conv_in, t_sk_ox_in, frazil_scav_in)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in), optional :: &
             ktherm_in          ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(in), optional :: &
             conduct_in, &      ! 'MU71' or 'bubbly'
             fbot_xfer_type_in  ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(in), optional :: &
             heat_capacity_in, &! if true, ice has nonzero heat capacity
                                ! if false, use zero-layer thermodynamics
             calc_Tsfc_in       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(in), optional :: &
             dts_b_in,   &      ! zsalinity timestep
             ustar_min_in       ! minimum friction velocity for ice-ocean heat flux
 
        ! mushy thermo
        real(kind=dbl_kind), intent(in), optional :: &
             a_rapid_mode_in      , & ! channel radius for rapid drainage mode (m)
             Rac_rapid_mode_in    , & ! critical Rayleigh number for rapid drainage mode
             aspect_rapid_mode_in , & ! aspect ratio for rapid drainage mode (larger=wider)
             dSdt_slow_mode_in    , & ! slow mode drainage strength (m s-1 K-1)
             phi_c_slow_mode_in   , & ! liquid fraction porosity cutoff for slow mode
             phi_i_mushy_in           ! liquid fraction of congelation ice
        
!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

        character (len=char_len), intent(in), optional :: &
             shortwave_in, & ! shortwave method, 'ccsm3' or 'dEdd'
             albedo_type_in  ! albedo parameterization, 'ccsm3' or 'constant'
                             ! shortwave='dEdd' overrides this parameter

        ! baseline albedos for ccsm3 shortwave, set in namelist
        real (kind=dbl_kind), intent(in), optional :: &
             albicev_in  , & ! visible ice albedo for h > ahmax
             albicei_in  , & ! near-ir ice albedo for h > ahmax
             albsnowv_in , & ! cold snow albedo, visible
             albsnowi_in , & ! cold snow albedo, near IR
             ahmax_in        ! thickness above which ice albedo is constant (m)
        
        ! dEdd tuning parameters, set in namelist
        real (kind=dbl_kind), intent(in), optional :: &
             R_ice_in    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
             R_pnd_in    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
             R_snw_in    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
             dT_mlt_in   , & ! change in temp for non-melt to melt snow grain 
                             ! radius change (C)
             rsnw_mlt_in , & ! maximum melting snow grain radius (10^-6 m)
             kalg_in         ! algae absorption coefficient for 0.5 m thick layer

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in), optional :: & ! defined in namelist 
             kstrength_in  , & ! 0 for simple Hibler (1979) formulation 
                               ! 1 for Rothrock (1975) pressure formulation 
             krdg_partic_in, & ! 0 for Thorndike et al. (1975) formulation 
                               ! 1 for exponential participation function 
             krdg_redist_in    ! 0 for Hibler (1980) formulation 
                               ! 1 for exponential redistribution function 
 
        real (kind=dbl_kind), intent(in), optional :: &  
             mu_rdg_in         ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

        character (len=char_len), intent(in), optional :: &
             atmbndy_in ! atmo boundary method, 'default' ('ccsm3') or 'constant'
        
        logical (kind=log_kind), intent(in), optional :: &
             calc_strair_in, &  ! if true, calculate wind stress components
             formdrag_in,    &  ! if true, calculate form drag
             highfreq_in        ! if true, use high frequency coupling
        
        integer (kind=int_kind), intent(in), optional :: &
             natmiter_in        ! number of iterations for boundary layer calculations
        
!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

        character(len=char_len), intent(in), optional :: &
             tfrz_option_in              ! form of ocean freezing temperature
                                         ! 'minus1p8' = -1.8 C
                                         ! 'linear_salt' = -depressT * sss
                                         ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in), optional :: &
             kitd_in        , & ! type of itd conversions
                                !   0 = delta function
                                !   1 = linear remap
             kcatbound_in       !   0 = old category boundary formula
                                !   1 = new formula giving round numbers
                                !   2 = WMO standard
                                !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

     character(char_len), intent(in), optional :: &     
        bgc_flux_type_in    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      

      logical (kind=log_kind), intent(in), optional :: &
         z_tracers_in,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_in,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_in,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_in,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_in        ! if .true., use modal aerosol formulation in shortwave
        
      logical (kind=log_kind), intent(in), optional :: & 
         skl_bgc_in,        &   ! if true, solve skeletal biochemistry
         solve_zsal_in          ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), intent(in), optional :: & 
         grid_o_in      , & ! for bottom flux        
         l_sk_in        , & ! characteristic diffusive scale (zsalinity) (m)
         initbio_frac_in, & ! fraction of ocean tracer concentration used to initialize tracer 
         phi_snow_in        ! snow porosity at the ice/snow interface 

      real (kind=dbl_kind), intent(in), optional :: & 
         grid_oS_in     , & ! for bottom flux (zsalinity)
         l_skS_in           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(in), optional :: &
         fr_resp_in           , &   ! fraction of algal growth lost due to respiration
         algal_vel_in         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_in        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_in        , &   ! solubility fraction
         T_max_in            , & ! maximum temperature (C)
         fsal_in             , & ! Salinity limitation (ppt)
         op_dep_min_in       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s_in       , & ! fraction of grazing spilled or slopped
         fr_graze_e_in       , & ! fraction of assimilation excreted 
         fr_mort2min_in      , & ! fractionation of mortality to Am
         fr_dFe_in           , & ! fraction of remineralized nitrogen 
                                    ! (in units of algal iron)
         k_nitrif_in         , & ! nitrification rate (1/day)            
         t_iron_conv_in      , & ! desorption loss pFe to dFe (day)
         max_loss_in         , & ! restrict uptake to % of remaining value 
         max_dfe_doc1_in     , & ! max ratio of dFe to saccharides in the ice 
                                    ! (nM Fe/muM C)    
         fr_resp_s_in        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS_in         , & ! fraction conversion given high yield
         t_sk_conv_in        , & ! Stefels conversion time (d)
         t_sk_ox_in          , & ! DMS oxidation time (d)
         frazil_scav_in          ! scavenging fraction or multiple in frazil ice

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

        real (kind=dbl_kind), intent(in), optional :: &
             hs0_in             ! snow depth for transition to bare sea ice (m)
        
        ! level-ice ponds
        character (len=char_len), intent(in), optional :: &
             frzpnd_in          ! pond refreezing parameterization
        
        real (kind=dbl_kind), intent(in), optional :: &
             dpscale_in, &      ! alter e-folding time scale for flushing 
             rfracmin_in, &     ! minimum retained fraction of meltwater
             rfracmax_in, &     ! maximum retained fraction of meltwater
             pndaspect_in, &    ! ratio of pond depth to pond fraction
             hs1_in             ! tapering parameter for snow on pond ice
        
        ! topo ponds
        real (kind=dbl_kind), intent(in), optional :: &
             hp1_in             ! critical parameter for pond ice thickness

        character(len=*),parameter :: subname='(icepack_init_parameters)'

        if (present(ktherm_in)               ) ktherm        = ktherm_in
        if (present(conduct_in)              ) conduct       = conduct_in
        if (present(fbot_xfer_type_in)       ) fbot_xfer_type    = fbot_xfer_type_in
        if (present(heat_capacity_in)        ) heat_capacity     = heat_capacity_in
        if (present(calc_Tsfc_in)            ) calc_Tsfc         = calc_Tsfc_in
        if (present(dts_b_in)                ) dts_b             = dts_b_in
        if (present(ustar_min_in)            ) ustar_min         = ustar_min_in
        if (present(a_rapid_mode_in)         ) a_rapid_mode      = a_rapid_mode_in
        if (present(Rac_rapid_mode_in)       ) Rac_rapid_mode    = Rac_rapid_mode_in
        if (present(aspect_rapid_mode_in)    ) aspect_rapid_mode = aspect_rapid_mode_in
        if (present(dSdt_slow_mode_in)       ) dSdt_slow_mode    = dSdt_slow_mode_in
        if (present(phi_c_slow_mode_in)      ) phi_c_slow_mode   = phi_c_slow_mode_in
        if (present(phi_i_mushy_in)          ) phi_i_mushy       = phi_i_mushy_in
        if (present(shortwave_in)            ) shortwave     = shortwave_in
        if (present(albedo_type_in)          ) albedo_type   = albedo_type_in
        if (present(albicev_in)              ) albicev       = albicev_in
        if (present(albicei_in)              ) albicei       = albicei_in
        if (present(albsnowv_in)             ) albsnowv      = albsnowv_in
        if (present(albsnowi_in)             ) albsnowi      = albsnowi_in
        if (present(ahmax_in)                ) ahmax         = ahmax_in
        if (present(R_ice_in)                ) R_ice         = R_ice_in
        if (present(R_pnd_in)                ) R_pnd         = R_pnd_in
        if (present(R_snw_in)                ) R_snw         = R_snw_in
        if (present(dT_mlt_in)               ) dT_mlt        = dT_mlt_in
        if (present(rsnw_mlt_in)             ) rsnw_mlt      = rsnw_mlt_in
        if (present(kalg_in)                 ) kalg          = kalg_in
        if (present(kstrength_in)            ) kstrength     = kstrength_in
        if (present(krdg_partic_in)          ) krdg_partic   = krdg_partic_in
        if (present(krdg_redist_in)          ) krdg_redist   = krdg_redist_in
        if (present(mu_rdg_in)               ) mu_rdg        = mu_rdg_in
        if (present(atmbndy_in)              ) atmbndy       = atmbndy_in
        if (present(calc_strair_in)          ) calc_strair   = calc_strair_in
        if (present(formdrag_in)             ) formdrag      = formdrag_in
        if (present(highfreq_in)             ) highfreq      = highfreq_in
        if (present(natmiter_in)             ) natmiter      = natmiter_in
        if (present(tfrz_option_in)          ) tfrz_option   = tfrz_option_in
        if (present(kitd_in)                 ) kitd          = kitd_in
        if (present(kcatbound_in)            ) kcatbound     = kcatbound_in
        if (present(hs0_in)                  ) hs0           = hs0_in
        if (present(frzpnd_in)               ) frzpnd        = frzpnd_in
        if (present(dpscale_in)              ) dpscale       = dpscale_in
        if (present(rfracmin_in)             ) rfracmin      = rfracmin_in
        if (present(rfracmax_in)             ) rfracmax      = rfracmax_in
        if (present(pndaspect_in)            ) pndaspect     = pndaspect_in
        if (present(hs1_in)                  ) hs1           = hs1_in
        if (present(hp1_in)                  ) hp1           = hp1_in
        if (present(bgc_flux_type_in)        ) bgc_flux_type = bgc_flux_type_in
        if (present(z_tracers_in)            ) z_tracers     = z_tracers_in
        if (present(scale_bgc_in)            ) scale_bgc     = scale_bgc_in
        if (present(solve_zbgc_in)           ) solve_zbgc    = solve_zbgc_in
        if (present(dEdd_algae_in)           ) dEdd_algae    = dEdd_algae_in
        if (present(modal_aero_in)           ) modal_aero    = modal_aero_in
        if (present(skl_bgc_in)              ) skl_bgc       = skl_bgc_in
        if (present(solve_zsal_in)           ) solve_zsal    = solve_zsal_in
        if (present(grid_o_in)               ) grid_o        = grid_o_in
        if (present(l_sk_in)                 ) l_sk          = l_sk_in
        if (present(initbio_frac_in)         ) initbio_frac  = initbio_frac_in
        if (present(grid_oS_in)              ) grid_oS       = grid_oS_in
        if (present(l_skS_in)                ) l_skS         = l_skS_in
        if (present(phi_snow_in)             ) phi_snow      = phi_snow_in
        if (present(fr_resp_in)              ) fr_resp       = fr_resp_in
        if (present(algal_vel_in)            ) algal_vel     = algal_vel_in
        if (present(R_dFe2dust_in)           ) R_dFe2dust    = R_dFe2dust_in
        if (present(dustFe_sol_in)           ) dustFe_sol    = dustFe_sol_in
        if (present(T_max_in)                ) T_max         = T_max_in
        if (present(fsal_in)                 ) fsal          = fsal_in
        if (present(op_dep_min_in)           ) op_dep_min    = op_dep_min_in
        if (present(fr_graze_s_in)           ) fr_graze_s    = fr_graze_s_in
        if (present(fr_graze_e_in)           ) fr_graze_e    = fr_graze_e_in
        if (present(fr_mort2min_in)          ) fr_mort2min   = fr_mort2min_in
        if (present(fr_dFe_in)               ) fr_dFe        = fr_dFe_in
        if (present(k_nitrif_in)             ) k_nitrif      = k_nitrif_in
        if (present(t_iron_conv_in)          ) t_iron_conv   = t_iron_conv_in
        if (present(max_loss_in)             ) max_loss      = max_loss_in
        if (present(max_dfe_doc1_in)         ) max_dfe_doc1  = max_dfe_doc1_in
        if (present(fr_resp_s_in)            ) fr_resp_s     = fr_resp_s_in
        if (present(y_sk_DMS_in)             ) y_sk_DMS      = y_sk_DMS_in
        if (present(t_sk_conv_in)            ) t_sk_conv     = t_sk_conv_in
        if (present(t_sk_ox_in)              ) t_sk_ox       = t_sk_ox_in
        if (present(frazil_scav_in)          ) frazil_scav   = frazil_scav_in

      end subroutine icepack_init_parameters

!=======================================================================

! subroutine to query the column package internal parameters

      subroutine icepack_query_parameters(   &
           ktherm_out, conduct_out, fbot_xfer_type_out, calc_Tsfc_out, dts_b_out, ustar_min_out, a_rapid_mode_out, &
           Rac_rapid_mode_out, aspect_rapid_mode_out, dSdt_slow_mode_out, phi_c_slow_mode_out, &
           phi_i_mushy_out, shortwave_out, albedo_type_out, albicev_out, albicei_out, albsnowv_out, &
           albsnowi_out, ahmax_out, R_ice_out, R_pnd_out, R_snw_out, dT_mlt_out, rsnw_mlt_out, &
           kalg_out, kstrength_out, krdg_partic_out, krdg_redist_out, mu_rdg_out, &
           atmbndy_out, calc_strair_out, formdrag_out, highfreq_out, natmiter_out, &
           tfrz_option_out, kitd_out, kcatbound_out, hs0_out, frzpnd_out, &
           dpscale_out, rfracmin_out, rfracmax_out, pndaspect_out, hs1_out, hp1_out, &
           bgc_flux_type_out, z_tracers_out, scale_bgc_out, solve_zbgc_out, dEdd_algae_out, &
           modal_aero_out, skl_bgc_out, solve_zsal_out, grid_o_out, l_sk_out, &
           initbio_frac_out, grid_oS_out, l_skS_out, &
           phi_snow_out, heat_capacity_out, &
           fr_resp_out, algal_vel_out, R_dFe2dust_out, dustFe_sol_out, &         
           T_max_out, fsal_out, op_dep_min_out, fr_graze_s_out, fr_graze_e_out, fr_mort2min_out, &        
           fr_dFe_out, k_nitrif_out, t_iron_conv_out, max_loss_out, max_dfe_doc1_out, fr_resp_s_out, &          
           y_sk_DMS_out, t_sk_conv_out, t_sk_ox_out, frazil_scav_out)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(out), optional :: &
             ktherm_out         ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(out), optional :: &
             conduct_out, &     ! 'MU71' or 'bubbly'
             fbot_xfer_type_out ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(out), optional :: &
             heat_capacity_out,&! if true, ice has nonzero heat capacity
                                ! if false, use zero-layer thermodynamics
             calc_Tsfc_out      ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(out), optional :: &
             dts_b_out,   &      ! zsalinity timestep
             ustar_min_out       ! minimum friction velocity for ice-ocean heat flux
 
        ! mushy thermo
        real(kind=dbl_kind), intent(out), optional :: &
             a_rapid_mode_out      , & ! channel radius for rapid drainage mode (m)
             Rac_rapid_mode_out    , & ! critical Rayleigh number for rapid drainage mode
             aspect_rapid_mode_out , & ! aspect ratio for rapid drainage mode (larger=wider)
             dSdt_slow_mode_out    , & ! slow mode drainage strength (m s-1 K-1)
             phi_c_slow_mode_out   , & ! liquid fraction porosity cutoff for slow mode
             phi_i_mushy_out           ! liquid fraction of congelation ice
        
!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

        character (len=char_len), intent(out), optional :: &
             shortwave_out, & ! shortwave method, 'ccsm3' or 'dEdd'
             albedo_type_out  ! albedo parameterization, 'ccsm3' or 'constant'
                             ! shortwave='dEdd' overrides this parameter

        ! baseline albedos for ccsm3 shortwave, set in namelist
        real (kind=dbl_kind), intent(out), optional :: &
             albicev_out  , & ! visible ice albedo for h > ahmax
             albicei_out  , & ! near-ir ice albedo for h > ahmax
             albsnowv_out , & ! cold snow albedo, visible
             albsnowi_out , & ! cold snow albedo, near IR
             ahmax_out        ! thickness above which ice albedo is constant (m)
        
        ! dEdd tuning parameters, set in namelist
        real (kind=dbl_kind), intent(out), optional :: &
             R_ice_out    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
             R_pnd_out    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
             R_snw_out    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
             dT_mlt_out   , & ! change in temp for non-melt to melt snow grain 
                             ! radius change (C)
             rsnw_mlt_out , & ! maximum melting snow grain radius (10^-6 m)
             kalg_out         ! algae absorption coefficient for 0.5 m thick layer

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(out), optional :: & ! defined in namelist 
             kstrength_out  , & ! 0 for simple Hibler (1979) formulation 
                                ! 1 for Rothrock (1975) pressure formulation 
             krdg_partic_out, & ! 0 for Thorndike et al. (1975) formulation 
                                ! 1 for exponential participation function 
             krdg_redist_out    ! 0 for Hibler (1980) formulation 
                                ! 1 for exponential redistribution function 
 
        real (kind=dbl_kind), intent(out), optional :: &  
             mu_rdg_out         ! gives e-folding scale of ridged ice (m^.5) 
                                ! (krdg_redist = 1) 

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

        character (len=char_len), intent(out), optional :: &
             atmbndy_out ! atmo boundary method, 'default' ('ccsm3') or 'constant'
        
        logical (kind=log_kind), intent(out), optional :: &
             calc_strair_out, &  ! if true, calculate wind stress components
             formdrag_out,    &  ! if true, calculate form drag
             highfreq_out        ! if true, use high frequency coupling
        
        integer (kind=int_kind), intent(out), optional :: &
             natmiter_out        ! number of iterations for boundary layer calculations
        
!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

        character(len=char_len), intent(out), optional :: &
             tfrz_option_out              ! form of ocean freezing temperature
                                         ! 'minus1p8' = -1.8 C
                                         ! 'linear_salt' = -depressT * sss
                                         ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(out), optional :: &
             kitd_out        , & ! type of itd conversions
                                !   0 = delta function
                                !   1 = linear remap
             kcatbound_out       !   0 = old category boundary formula
                                !   1 = new formula giving round numbers
                                !   2 = WMO standard
                                !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

     character(char_len), intent(out), optional :: &     
        bgc_flux_type_out    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      

      logical (kind=log_kind), intent(out), optional :: &
         z_tracers_out,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_out,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_out,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_out,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_out        ! if .true., use modal aerosol formulation in shortwave
        
      logical (kind=log_kind), intent(out), optional :: & 
         skl_bgc_out,        &   ! if true, solve skeletal biochemistry
         solve_zsal_out          ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), intent(out), optional :: & 
         grid_o_out      , & ! for bottom flux        
         l_sk_out        , & ! characteristic diffusive scale (zsalinity) (m)
         initbio_frac_out, & ! fraction of ocean tracer concentration used to initialize tracer 
         phi_snow_out        ! snow porosity at the ice/snow interface 

      real (kind=dbl_kind), intent(out), optional :: & 
         grid_oS_out     , & ! for bottom flux (zsalinity)
         l_skS_out           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(out), optional :: &
         fr_resp_out           , &   ! fraction of algal growth lost due to respiration
         algal_vel_out         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_out        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_out        , &   ! solubility fraction
         T_max_out            , & ! maximum temperature (C)
         fsal_out             , & ! Salinity limitation (ppt)
         op_dep_min_out       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s_out       , & ! fraction of grazing spilled or slopped
         fr_graze_e_out       , & ! fraction of assimilation excreted 
         fr_mort2min_out      , & ! fractionation of mortality to Am
         fr_dFe_out           , & ! fraction of remineralized nitrogen 
                                    ! (in units of algal iron)
         k_nitrif_out         , & ! nitrification rate (1/day)            
         t_iron_conv_out      , & ! desorption loss pFe to dFe (day)
         max_loss_out         , & ! restrict uptake to % of remaining value 
         max_dfe_doc1_out     , & ! max ratio of dFe to saccharides in the ice 
                                    ! (nM Fe/muM C)    
         fr_resp_s_out        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS_out         , & ! fraction conversion given high yield
         t_sk_conv_out        , & ! Stefels conversion time (d)
         t_sk_ox_out          , & ! DMS oxidation time (d)
         frazil_scav_out          ! scavenging fraction or multiple in frazil ice

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

        real (kind=dbl_kind), intent(out), optional :: &
             hs0_out             ! snow depth for transition to bare sea ice (m)
        
        ! level-ice ponds
        character (len=char_len), intent(out), optional :: &
             frzpnd_out          ! pond refreezing parameterization
        
        real (kind=dbl_kind), intent(out), optional :: &
             dpscale_out, &      ! alter e-folding time scale for flushing 
             rfracmin_out, &     ! minimum retained fraction of meltwater
             rfracmax_out, &     ! maximum retained fraction of meltwater
             pndaspect_out, &    ! ratio of pond depth to pond fraction
             hs1_out             ! tapering parameter for snow on pond ice
        
        ! topo ponds
        real (kind=dbl_kind), intent(out), optional :: &
             hp1_out             ! critical parameter for pond ice thickness

        character(len=*),parameter :: subname='(icepack_query_parameters)'

        if (present(ktherm_out)               ) ktherm_out        = ktherm
        if (present(conduct_out)              ) conduct_out       = conduct
        if (present(fbot_xfer_type_out)       ) fbot_xfer_type_out    = fbot_xfer_type
        if (present(heat_capacity_out)        ) heat_capacity_out     = heat_capacity
        if (present(calc_Tsfc_out)            ) calc_Tsfc_out         = calc_Tsfc
        if (present(dts_b_out)                ) dts_b_out             = dts_b
        if (present(ustar_min_out)            ) ustar_min_out         = ustar_min
        if (present(a_rapid_mode_out)         ) a_rapid_mode_out      = a_rapid_mode
        if (present(Rac_rapid_mode_out)       ) Rac_rapid_mode_out    = Rac_rapid_mode
        if (present(aspect_rapid_mode_out)    ) aspect_rapid_mode_out = aspect_rapid_mode
        if (present(dSdt_slow_mode_out)       ) dSdt_slow_mode_out    = dSdt_slow_mode
        if (present(phi_c_slow_mode_out)      ) phi_c_slow_mode_out   = phi_c_slow_mode
        if (present(phi_i_mushy_out)          ) phi_i_mushy_out       = phi_i_mushy
        if (present(shortwave_out)            ) shortwave_out     = shortwave
        if (present(albedo_type_out)          ) albedo_type_out   = albedo_type
        if (present(albicev_out)              ) albicev_out       = albicev
        if (present(albicei_out)              ) albicei_out       = albicei
        if (present(albsnowv_out)             ) albsnowv_out      = albsnowv
        if (present(albsnowi_out)             ) albsnowi_out      = albsnowi
        if (present(ahmax_out)                ) ahmax_out         = ahmax
        if (present(R_ice_out)                ) R_ice_out         = R_ice
        if (present(R_pnd_out)                ) R_pnd_out         = R_pnd
        if (present(R_snw_out)                ) R_snw_out         = R_snw
        if (present(dT_mlt_out)               ) dT_mlt_out        = dT_mlt
        if (present(rsnw_mlt_out)             ) rsnw_mlt_out      = rsnw_mlt
        if (present(kalg_out)                 ) kalg_out          = kalg
        if (present(kstrength_out)            ) kstrength_out     = kstrength
        if (present(krdg_partic_out)          ) krdg_partic_out   = krdg_partic
        if (present(krdg_redist_out)          ) krdg_redist_out   = krdg_redist
        if (present(mu_rdg_out)               ) mu_rdg_out        = mu_rdg
        if (present(atmbndy_out)              ) atmbndy_out       = atmbndy
        if (present(calc_strair_out)          ) calc_strair_out   = calc_strair
        if (present(formdrag_out)             ) formdrag_out      = formdrag
        if (present(highfreq_out)             ) highfreq_out      = highfreq
        if (present(natmiter_out)             ) natmiter_out      = natmiter
        if (present(tfrz_option_out)          ) tfrz_option_out   = tfrz_option
        if (present(kitd_out)                 ) kitd_out          = kitd
        if (present(kcatbound_out)            ) kcatbound_out     = kcatbound
        if (present(hs0_out)                  ) hs0_out           = hs0
        if (present(frzpnd_out)               ) frzpnd_out        = frzpnd
        if (present(dpscale_out)              ) dpscale_out       = dpscale
        if (present(rfracmin_out)             ) rfracmin_out      = rfracmin
        if (present(rfracmax_out)             ) rfracmax_out      = rfracmax
        if (present(pndaspect_out)            ) pndaspect_out     = pndaspect
        if (present(hs1_out)                  ) hs1_out           = hs1
        if (present(hp1_out)                  ) hp1_out           = hp1
        if (present(bgc_flux_type_out)        ) bgc_flux_type_out = bgc_flux_type
        if (present(z_tracers_out)            ) z_tracers_out     = z_tracers
        if (present(scale_bgc_out)            ) scale_bgc_out     = scale_bgc
        if (present(solve_zbgc_out)           ) solve_zbgc_out    = solve_zbgc
        if (present(dEdd_algae_out)           ) dEdd_algae_out    = dEdd_algae
        if (present(modal_aero_out)           ) modal_aero_out    = modal_aero
        if (present(skl_bgc_out)              ) skl_bgc_out       = skl_bgc
        if (present(solve_zsal_out)           ) solve_zsal_out    = solve_zsal
        if (present(grid_o_out)               ) grid_o_out        = grid_o
        if (present(l_sk_out)                 ) l_sk_out          = l_sk
        if (present(initbio_frac_out)         ) initbio_frac_out  = initbio_frac
        if (present(grid_oS_out)              ) grid_oS_out       = grid_oS
        if (present(l_skS_out)                ) l_skS_out         = l_skS
        if (present(phi_snow_out)             ) phi_snow_out      = phi_snow
        if (present(fr_resp_out)              ) fr_resp_out       = fr_resp
        if (present(algal_vel_out)            ) algal_vel_out     = algal_vel
        if (present(R_dFe2dust_out)           ) R_dFe2dust_out    = R_dFe2dust
        if (present(dustFe_sol_out)           ) dustFe_sol_out    = dustFe_sol
        if (present(T_max_out)                ) T_max_out         = T_max
        if (present(fsal_out)                 ) fsal_out          = fsal
        if (present(op_dep_min_out)           ) op_dep_min_out    = op_dep_min
        if (present(fr_graze_s_out)           ) fr_graze_s_out    = fr_graze_s
        if (present(fr_graze_e_out)           ) fr_graze_e_out    = fr_graze_e
        if (present(fr_mort2min_out)          ) fr_mort2min_out   = fr_mort2min
        if (present(fr_dFe_out)               ) fr_dFe_out        = fr_dFe
        if (present(k_nitrif_out)             ) k_nitrif_out      = k_nitrif
        if (present(t_iron_conv_out)          ) t_iron_conv_out   = t_iron_conv
        if (present(max_loss_out)             ) max_loss_out      = max_loss
        if (present(max_dfe_doc1_out)         ) max_dfe_doc1_out  = max_dfe_doc1
        if (present(fr_resp_s_out)            ) fr_resp_s_out     = fr_resp_s
        if (present(y_sk_DMS_out)             ) y_sk_DMS_out      = y_sk_DMS
        if (present(t_sk_conv_out)            ) t_sk_conv_out     = t_sk_conv
        if (present(t_sk_ox_out)              ) t_sk_ox_out       = t_sk_ox
        if (present(frazil_scav_out)          ) frazil_scav_out   = frazil_scav

      end subroutine icepack_query_parameters

!=======================================================================

! subroutine to write the column package internal parameters

      subroutine icepack_write_parameters(iounit)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             iounit   ! unit number for output

        character(len=*),parameter :: subname='(icepack_write_parameters)'

        write(iounit,*) subname
        write(iounit,*) "  ktherm        = ", ktherm
        write(iounit,*) "  conduct       = ", conduct
        write(iounit,*) "  fbot_xfer_type    = ", fbot_xfer_type
        write(iounit,*) "  heat_capacity     = ", heat_capacity
        write(iounit,*) "  calc_Tsfc         = ", calc_Tsfc
        write(iounit,*) "  dts_b             = ", dts_b
        write(iounit,*) "  ustar_min         = ", ustar_min
        write(iounit,*) "  a_rapid_mode      = ", a_rapid_mode
        write(iounit,*) "  Rac_rapid_mode    = ", Rac_rapid_mode
        write(iounit,*) "  aspect_rapid_mode = ", aspect_rapid_mode
        write(iounit,*) "  dSdt_slow_mode    = ", dSdt_slow_mode
        write(iounit,*) "  phi_c_slow_mode   = ", phi_c_slow_mode
        write(iounit,*) "  phi_i_mushy       = ", phi_i_mushy
        write(iounit,*) "  shortwave     = ", shortwave
        write(iounit,*) "  albedo_type   = ", albedo_type
        write(iounit,*) "  albicev       = ", albicev
        write(iounit,*) "  albicei       = ", albicei
        write(iounit,*) "  albsnowv      = ", albsnowv
        write(iounit,*) "  albsnowi      = ", albsnowi
        write(iounit,*) "  ahmax         = ", ahmax
        write(iounit,*) "  R_ice         = ", R_ice
        write(iounit,*) "  R_pnd         = ", R_pnd
        write(iounit,*) "  R_snw         = ", R_snw
        write(iounit,*) "  dT_mlt        = ", dT_mlt
        write(iounit,*) "  rsnw_mlt      = ", rsnw_mlt
        write(iounit,*) "  kalg          = ", kalg
        write(iounit,*) "  kstrength     = ", kstrength
        write(iounit,*) "  krdg_partic   = ", krdg_partic
        write(iounit,*) "  krdg_redist   = ", krdg_redist
        write(iounit,*) "  mu_rdg        = ", mu_rdg
        write(iounit,*) "  atmbndy       = ", atmbndy
        write(iounit,*) "  calc_strair   = ", calc_strair
        write(iounit,*) "  formdrag      = ", formdrag
        write(iounit,*) "  highfreq      = ", highfreq
        write(iounit,*) "  natmiter      = ", natmiter
        write(iounit,*) "  tfrz_option   = ", tfrz_option
        write(iounit,*) "  kitd          = ", kitd
        write(iounit,*) "  kcatbound     = ", kcatbound
        write(iounit,*) "  hs0           = ", hs0
        write(iounit,*) "  frzpnd        = ", frzpnd
        write(iounit,*) "  dpscale       = ", dpscale
        write(iounit,*) "  rfracmin      = ", rfracmin
        write(iounit,*) "  rfracmax      = ", rfracmax
        write(iounit,*) "  pndaspect     = ", pndaspect
        write(iounit,*) "  hs1           = ", hs1
        write(iounit,*) "  hp1           = ", hp1
        write(iounit,*) "  bgc_flux_type = ", bgc_flux_type
        write(iounit,*) "  z_tracers     = ", z_tracers
        write(iounit,*) "  scale_bgc     = ", scale_bgc
        write(iounit,*) "  solve_zbgc    = ", solve_zbgc
        write(iounit,*) "  dEdd_algae    = ", dEdd_algae
        write(iounit,*) "  modal_aero    = ", modal_aero
        write(iounit,*) "  skl_bgc       = ", skl_bgc
        write(iounit,*) "  solve_zsal    = ", solve_zsal
        write(iounit,*) "  grid_o        = ", grid_o
        write(iounit,*) "  l_sk          = ", l_sk
        write(iounit,*) "  initbio_frac  = ", initbio_frac
        write(iounit,*) "  grid_oS       = ", grid_oS
        write(iounit,*) "  l_skS         = ", l_skS
        write(iounit,*) "  phi_snow      = ", phi_snow
        write(iounit,*) "  fr_resp       = ", fr_resp
        write(iounit,*) "  algal_vel     = ", algal_vel
        write(iounit,*) "  R_dFe2dust    = ", R_dFe2dust
        write(iounit,*) "  dustFe_sol    = ", dustFe_sol
        write(iounit,*) "  T_max         = ", T_max
        write(iounit,*) "  fsal          = ", fsal
        write(iounit,*) "  op_dep_min    = ", op_dep_min
        write(iounit,*) "  fr_graze_s    = ", fr_graze_s
        write(iounit,*) "  fr_graze_e    = ", fr_graze_e
        write(iounit,*) "  fr_mort2min   = ", fr_mort2min
        write(iounit,*) "  fr_dFe        = ", fr_dFe
        write(iounit,*) "  k_nitrif      = ", k_nitrif
        write(iounit,*) "  t_iron_conv   = ", t_iron_conv
        write(iounit,*) "  max_loss      = ", max_loss
        write(iounit,*) "  max_dfe_doc1  = ", max_dfe_doc1
        write(iounit,*) "  fr_resp_s     = ", fr_resp_s
        write(iounit,*) "  y_sk_DMS      = ", y_sk_DMS
        write(iounit,*) "  t_sk_conv     = ", t_sk_conv
        write(iounit,*) "  t_sk_ox       = ", t_sk_ox
        write(iounit,*) "  frazil_scav   = ", frazil_scav

      end subroutine icepack_write_parameters

!=======================================================================

    end module icepack_parameters

!=======================================================================
