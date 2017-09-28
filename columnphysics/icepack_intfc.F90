!  SVN:$Id: icepack_intfc.F90 1227 2017-05-22 22:49:10Z tcraig $
!=========================================================================
!
! flags and interface routines for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_intfc

      use icepack_kinds_mod    ! kinds
      use icepack_intfc_shared ! namelist and other parameters

      use icepack_itd, only: icepack_init_itd
      use icepack_itd, only: icepack_init_itd_hist
      use icepack_itd, only: icepack_aggregate

      use icepack_mechred, only: icepack_step_ridge
      use icepack_mechred, only: icepack_ice_strength

      use icepack_shortwave, only: icepack_prep_radiation
      use icepack_shortwave, only: icepack_step_radiation

      use icepack_brine, only: icepack_init_hbrine
      use icepack_brine, only: icepack_init_zsalinity

      use icepack_zbgc , only: icepack_init_bgc
      use icepack_zbgc , only: icepack_init_zbgc
      use icepack_zbgc , only: icepack_init_bgc_trcr
      use icepack_zbgc , only: icepack_biogeochemistry
      use icepack_zbgc , only: icepack_init_OceanConcArray
      use icepack_zbgc , only: icepack_init_ocean_conc

      use icepack_atmo , only: icepack_atm_boundary
      use icepack_ocean, only: icepack_ocn_mixed_layer

      use icepack_therm_vertical, only: icepack_step_therm1
      use icepack_therm_itd     , only: icepack_step_therm2
      use icepack_therm_shared  , only: icepack_ice_temperature
      use icepack_therm_shared  , only: icepack_snow_temperature
      use icepack_therm_shared  , only: icepack_liquidus_temperature
      use icepack_therm_shared  , only: icepack_sea_freezing_temperature
      use icepack_therm_shared  , only: icepack_enthalpy_snow
      use icepack_therm_shared  , only: icepack_init_thermo
      use icepack_therm_shared  , only: icepack_init_trcr

      use icepack_orbital , only: icepack_init_orbit

      use icepack_warnings, only: icepack_clear_warnings
      use icepack_warnings, only: icepack_get_warnings
      use icepack_warnings, only: icepack_print_warnings

      implicit none

      public

      ! initialization
      public :: &
           icepack_init_parameters, &
           icepack_init_tracer_flags, &
           icepack_init_tracer_indices, &
           icepack_init_tracer_numbers

!=======================================================================

      contains

!=======================================================================
!     Initialization routines
!=======================================================================
! subroutine to set the column package internal parameters

      subroutine icepack_init_parameters(&
           ktherm_in, &
           conduct_in, &
           fbot_xfer_type_in, &
           calc_Tsfc_in, &
           ustar_min_in, &
           a_rapid_mode_in, &
           Rac_rapid_mode_in, &
           aspect_rapid_mode_in, &
           dSdt_slow_mode_in, &
           phi_c_slow_mode_in, &
           phi_i_mushy_in, &
           shortwave_in, &
           albedo_type_in, &
           albicev_in, &
           albicei_in, &
           albsnowv_in, &
           albsnowi_in, &
           ahmax_in, &
           R_ice_in, &
           R_pnd_in, &
           R_snw_in, &
           dT_mlt_in, &
           rsnw_mlt_in, &
           kalg_in, &
           kstrength_in, &
           krdg_partic_in, &
           krdg_redist_in, &
           mu_rdg_in, &
           Cf_in, &
           atmbndy_in, &
           calc_strair_in, &
           formdrag_in, &
           highfreq_in, &
           natmiter_in, &
           oceanmixed_ice_in, &
           tfrz_option_in, &
           kitd_in, &
           kcatbound_in, &
           hs0_in, &
           frzpnd_in, &
           dpscale_in, &
           rfracmin_in, &
           rfracmax_in, &
           pndaspect_in, &
           hs1_in, &
           hp1_in, &
         ! bgc_data_dir_in, &
         ! sil_data_type_in, &
         ! nit_data_type_in, &
         ! fe_data_type_in, &
           bgc_flux_type_in, &
           z_tracers_in, &
           scale_bgc_in, &
           solve_zbgc_in, &
           dEdd_algae_in, &
           modal_aero_in, &
           skl_bgc_in, &
           solve_zsal_in, &
           grid_o_in, &
           l_sk_in, &
           grid_o_t_in, &
           initbio_frac_in, &
           frazil_scav_in, &
           grid_oS_in, &
           l_skS_in, &
           phi_snow_in, &
           ratio_Si2N_diatoms_in, &
           ratio_Si2N_sp_in, &
           ratio_Si2N_phaeo_in, &
           ratio_S2N_diatoms_in, &
           ratio_S2N_sp_in, &      
           ratio_S2N_phaeo_in, &   
           ratio_Fe2C_diatoms_in, & 
           ratio_Fe2C_sp_in, &     
           ratio_Fe2C_phaeo_in, &  
           ratio_Fe2N_diatoms_in, & 
           ratio_Fe2N_sp_in, &     
           ratio_Fe2N_phaeo_in, &  
           ratio_Fe2DON_in, &       
           ratio_Fe2DOC_s_in, &     
           ratio_Fe2DOC_l_in, &     
           fr_resp_in, &            
           tau_min_in, &            
           tau_max_in, &            
           algal_vel_in, &          
           R_dFe2dust_in, &         
           dustFe_sol_in, &         
           chlabs_diatoms_in, &    
           chlabs_sp_in, &         
           chlabs_phaeo_in, &      
           alpha2max_low_diatoms_in, &  
           alpha2max_low_sp_in, &       
           alpha2max_low_phaeo_in, &    
           beta2max_diatoms_in, & 
           beta2max_sp_in, &       
           beta2max_phaeo_in, &    
           mu_max_diatoms_in, &   
           mu_max_sp_in, &         
           mu_max_phaeo_in, &      
           grow_Tdep_diatoms_in, &
           grow_Tdep_sp_in, &      
           grow_Tdep_phaeo_in, &   
           fr_graze_diatoms_in, & 
           fr_graze_sp_in, &       
           fr_graze_phaeo_in, &    
           mort_pre_diatoms_in, & 
           mort_pre_sp_in, &       
           mort_pre_phaeo_in, &    
           mort_Tdep_diatoms_in, &
           mort_Tdep_sp_in, &       
           mort_Tdep_phaeo_in, &    
           k_exude_diatoms_in, &  
           k_exude_sp_in, &         
           k_exude_phaeo_in, &      
           K_Nit_diatoms_in, &    
           K_Nit_sp_in, &           
           K_Nit_phaeo_in, &        
           K_Am_diatoms_in, &     
           K_Am_sp_in, &             
           K_Am_phaeo_in, &          
           K_Sil_diatoms_in, &    
           K_Sil_sp_in, &            
           K_Sil_phaeo_in, &         
           K_Fe_diatoms_in, &     
           K_Fe_sp_in, &             
           K_Fe_phaeo_in, &           
           f_don_protein_in, &    
           kn_bac_protein_in, &   
           f_don_Am_protein_in, & 
           f_doc_s_in, &            
           f_doc_l_in, &               
           f_exude_s_in, &          
           f_exude_l_in, &           
           k_bac_s_in, &            
           k_bac_l_in, &             
           T_max_in, &              
           fsal_in, &               
           op_dep_min_in, &         
           fr_graze_s_in, &         
           fr_graze_e_in, &         
           fr_mort2min_in, &        
           fr_dFe_in, &             
           k_nitrif_in, &           
           t_iron_conv_in, &        
           max_loss_in, &           
           max_dfe_doc1_in, &       
           fr_resp_s_in, &          
           y_sk_DMS_in, &           
           t_sk_conv_in, &          
           t_sk_ox_in, &             
           algaltype_diatoms_in, &   
           algaltype_sp_in, &       
           algaltype_phaeo_in, &    
           nitratetype_in, &        
           ammoniumtype_in, &       
           silicatetype_in, &       
           dmspptype_in, &          
           dmspdtype_in, &          
           humtype_in, &            
           doctype_s_in, &          
           doctype_l_in, &          
           dontype_protein_in, &     
           fedtype_1_in, &           
           feptype_1_in, &           
           zaerotype_bc1_in, &       
           zaerotype_bc2_in, &       
           zaerotype_dust1_in, &     
           zaerotype_dust2_in, &     
           zaerotype_dust3_in, &     
           zaerotype_dust4_in, &     
           ratio_C2N_diatoms_in, &   
           ratio_C2N_sp_in, &        
           ratio_C2N_phaeo_in, &     
           ratio_chl2N_diatoms_in, & 
           ratio_chl2N_sp_in, &      
           ratio_chl2N_phaeo_in, &   
           F_abs_chl_diatoms_in, &   
           F_abs_chl_sp_in, &        
           F_abs_chl_phaeo_in, &
           ratio_C2N_proteins_in)
           !restore_bgc_in)

        use icepack_intfc_shared, only: &
             ktherm, &
             conduct, &
             fbot_xfer_type, &
             calc_Tsfc, &
             ustar_min, &
             a_rapid_mode, &
             Rac_rapid_mode, &
             aspect_rapid_mode, &
             dSdt_slow_mode, &
             phi_c_slow_mode, &
             phi_i_mushy, &
             shortwave, &
             albedo_type, &
             albicev, &
             albicei, &
             albsnowv, &
             albsnowi, &
             ahmax, &
             R_ice, &
             R_pnd, &
             R_snw, &
             dT_mlt, &
             rsnw_mlt, &
             kalg, &
             kstrength, &
             krdg_partic, &
             krdg_redist, &
             mu_rdg, &
             Cf, &
             atmbndy, &
             calc_strair, &
             formdrag, &
             highfreq, &
             natmiter, &
             oceanmixed_ice, &
             tfrz_option, &
             kitd, &
             kcatbound, &
             hs0, &
             frzpnd, &
             dpscale, &
             rfracmin, &
             rfracmax, &
             pndaspect, &
             hs1, &
             hp1, &
           ! bgc_data_dir, &
           ! sil_data_type, &
           ! nit_data_type, &
           ! fe_data_type, &
             bgc_flux_type, &
             z_tracers, &
             scale_bgc, &
             solve_zbgc, &
             dEdd_algae, &
             modal_aero, &
             skl_bgc, &
             solve_zsal, &
             grid_o, &
             l_sk, &
             grid_o_t, &
             initbio_frac, &
             frazil_scav, &
             grid_oS, &
             l_skS, &
             phi_snow, &
             ratio_Si2N_diatoms, & 
             ratio_Si2N_sp     , &
             ratio_Si2N_phaeo  , &
             ratio_S2N_diatoms , & 
             ratio_S2N_sp      , &
             ratio_S2N_phaeo   , &
             ratio_Fe2C_diatoms, & 
             ratio_Fe2C_sp     , &
             ratio_Fe2C_phaeo  , &
             ratio_Fe2N_diatoms, & 
             ratio_Fe2N_sp     , &
             ratio_Fe2N_phaeo  , &
             ratio_Fe2DON      , & 
             ratio_Fe2DOC_s    , & 
             ratio_Fe2DOC_l    , & 
             fr_resp           , & 
             tau_min           , & 
             tau_max           , & 
             algal_vel         , & 
             R_dFe2dust        , & 
             dustFe_sol        , & 
             chlabs_diatoms    , &
             chlabs_sp         , &
             chlabs_phaeo      , &
             alpha2max_low_diatoms , & 
             alpha2max_low_sp      , & 
             alpha2max_low_phaeo   , & 
             beta2max_diatoms , &
             beta2max_sp      , & 
             beta2max_phaeo   , & 
             mu_max_diatoms   , &
             mu_max_sp        , & 
             mu_max_phaeo     , & 
             grow_Tdep_diatoms, &
             grow_Tdep_sp     , & 
             grow_Tdep_phaeo  , & 
             fr_graze_diatoms , &
             fr_graze_sp      , & 
             fr_graze_phaeo   , & 
             mort_pre_diatoms , &
             mort_pre_sp      , & 
             mort_pre_phaeo   , & 
             mort_Tdep_diatoms, &
             mort_Tdep_sp     , &  
             mort_Tdep_phaeo  , &  
             k_exude_diatoms  , &
             k_exude_sp       , &  
             k_exude_phaeo    , &  
             K_Nit_diatoms    , &
             K_Nit_sp         , &  
             K_Nit_phaeo      , &  
             K_Am_diatoms     , &
             K_Am_sp          , &   
             K_Am_phaeo       , &   
             K_Sil_diatoms    , &
             K_Sil_sp         , &   
             K_Sil_phaeo      , &   
             K_Fe_diatoms     , &
             K_Fe_sp          , &   
             K_Fe_phaeo       , &    
             f_don_protein    , &
             kn_bac_protein   , &
             f_don_Am_protein , &
             f_doc_s            , &
             f_doc_l            , &   
             f_exude_s          , &
             f_exude_l          , & 
             k_bac_s            , &
             k_bac_l            , & 
             T_max              , &
             fsal               , &
             op_dep_min         , &
             fr_graze_s         , &
             fr_graze_e         , &
             fr_mort2min        , &
             fr_dFe             , &
             k_nitrif           , &
             t_iron_conv        , &
             max_loss           , &
             max_dfe_doc1       , &
             fr_resp_s          , &
             y_sk_DMS           , &
             t_sk_conv          , &
             t_sk_ox            , & 
             algaltype_diatoms  , & 
             algaltype_sp       , &
             algaltype_phaeo    , &
             nitratetype        , &
             ammoniumtype       , &
             silicatetype       , &
             dmspptype          , &
             dmspdtype          , &
             humtype            , &
             doctype_s          , &
             doctype_l          , &
             dontype_protein    , & 
             fedtype_1          , & 
             feptype_1          , & 
             zaerotype_bc1      , & 
             zaerotype_bc2      , & 
             zaerotype_dust1    , & 
             zaerotype_dust2    , & 
             zaerotype_dust3    , & 
             zaerotype_dust4    , & 
             ratio_C2N_diatoms  , & 
             ratio_C2N_sp       , & 
             ratio_C2N_phaeo    , & 
             ratio_chl2N_diatoms, & 
             ratio_chl2N_sp     , & 
             ratio_chl2N_phaeo  , & 
             F_abs_chl_diatoms  , & 
             F_abs_chl_sp       , & 
             F_abs_chl_phaeo    , & 
             ratio_C2N_proteins
            !restore_bgc

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             ktherm_in          ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(in) :: &
             conduct_in, &      ! 'MU71' or 'bubbly'
             fbot_xfer_type_in  ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(in) :: &
             calc_Tsfc_in       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(in) :: &
             ustar_min_in       ! minimum friction velocity for ice-ocean heat flux
 
        ! mushy thermo
        real(kind=dbl_kind), intent(in) :: &
             a_rapid_mode_in      , & ! channel radius for rapid drainage mode (m)
             Rac_rapid_mode_in    , & ! critical Rayleigh number for rapid drainage mode
             aspect_rapid_mode_in , & ! aspect ratio for rapid drainage mode (larger=wider)
             dSdt_slow_mode_in    , & ! slow mode drainage strength (m s-1 K-1)
             phi_c_slow_mode_in   , & ! liquid fraction porosity cutoff for slow mode
             phi_i_mushy_in           ! liquid fraction of congelation ice
        
!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             shortwave_in, & ! shortwave method, 'default' ('ccsm3') or 'dEdd'
             albedo_type_in  ! albedo parameterization, 'default' ('ccsm3') or 'constant'
                             ! shortwave='dEdd' overrides this parameter

        ! baseline albedos for ccsm3 shortwave, set in namelist
        real (kind=dbl_kind), intent(in) :: &
             albicev_in  , & ! visible ice albedo for h > ahmax
             albicei_in  , & ! near-ir ice albedo for h > ahmax
             albsnowv_in , & ! cold snow albedo, visible
             albsnowi_in , & ! cold snow albedo, near IR
             ahmax_in        ! thickness above which ice albedo is constant (m)
        
        ! dEdd tuning parameters, set in namelist
        real (kind=dbl_kind), intent(in) :: &
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

        integer (kind=int_kind), intent(in) :: & ! defined in namelist 
             kstrength_in  , & ! 0 for simple Hibler (1979) formulation 
                               ! 1 for Rothrock (1975) pressure formulation 
             krdg_partic_in, & ! 0 for Thorndike et al. (1975) formulation 
                               ! 1 for exponential participation function 
             krdg_redist_in    ! 0 for Hibler (1980) formulation 
                               ! 1 for exponential redistribution function 
 
        real (kind=dbl_kind), intent(in) :: &  
             mu_rdg_in, &      ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 
             Cf_in             ! ratio of ridging work to PE change in ridging (kstrength = 1)

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             atmbndy_in ! atmo boundary method, 'default' ('ccsm3') or 'constant'
        
        logical (kind=log_kind), intent(in) :: &
             calc_strair_in, &  ! if true, calculate wind stress components
             formdrag_in,    &  ! if true, calculate form drag
             highfreq_in        ! if true, use high frequency coupling
        
        integer (kind=int_kind), intent(in) :: &
             natmiter_in        ! number of iterations for boundary layer calculations
        
!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

        logical (kind=log_kind), intent(in) :: &
             oceanmixed_ice_in           ! if true, use ocean mixed layer
        
        character(len=char_len), intent(in) :: &
             tfrz_option_in              ! form of ocean freezing temperature
                                         ! 'minus1p8' = -1.8 C
                                         ! 'linear_salt' = -depressT * sss
                                         ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
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

   !  character(char_len_long), intent(in) :: & 
   !     bgc_data_dir_in   ! directory for biogeochemistry data

     character(char_len), intent(in) :: &     
        bgc_flux_type_in    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      
    !    sil_data_type_in  , & ! 'default', 'clim'
    !    nit_data_type_in  , & ! 'default', 'clim'   
    !    fe_data_type_in   , & ! 'default', 'clim'      

      logical (kind=log_kind), intent(in) :: &
         z_tracers_in,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_in,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_in,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_in,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_in        ! if .true., use modal aerosol formulation in shortwave
        
      logical (kind=log_kind), intent(in) :: & 
         skl_bgc_in,        &   ! if true, solve skeletal biochemistry
         solve_zsal_in          ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), intent(in) :: & 
         grid_o_in      , & ! for bottom flux        
         l_sk_in        , & ! characteristic diffusive scale (zsalinity) (m)
         grid_o_t_in    , & ! top grid point length scale 
         initbio_frac_in, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav_in , & ! multiple of ocean tracer concentration due to frazil scavenging
         phi_snow_in        ! snow porosity at the ice/snow interface 

      real (kind=dbl_kind), intent(in) :: & 
         grid_oS_in     , & ! for bottom flux (zsalinity)
         l_skS_in           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(in) :: &
         ratio_Si2N_diatoms_in, &   ! algal Si to N (mol/mol)
         ratio_Si2N_sp_in     , &
         ratio_Si2N_phaeo_in  , &
         ratio_S2N_diatoms_in , &   ! algal S  to N (mol/mol)
         ratio_S2N_sp_in      , &
         ratio_S2N_phaeo_in   , &
         ratio_Fe2C_diatoms_in, &   ! algal Fe to C  (umol/mol)
         ratio_Fe2C_sp_in     , &
         ratio_Fe2C_phaeo_in  , &
         ratio_Fe2N_diatoms_in, &   ! algal Fe to N  (umol/mol)
         ratio_Fe2N_sp_in     , &
         ratio_Fe2N_phaeo_in  , &
         ratio_Fe2DON_in      , &   ! Fe to N of DON (nmol/umol)
         ratio_Fe2DOC_s_in    , &   ! Fe to C of DOC (nmol/umol) saccharids
         ratio_Fe2DOC_l_in    , &   ! Fe to C of DOC (nmol/umol) lipids
         fr_resp_in           , &   ! fraction of algal growth lost due to respiration
         tau_min_in           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
         tau_max_in           , &   ! long time mobile to stationary exchanges (s) = 2 days
         algal_vel_in         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_in        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_in        , &   ! solubility fraction
         chlabs_diatoms_in   , & ! chl absorption (1/m/(mg/m^3))
         chlabs_sp_in        , & !
         chlabs_phaeo_in     , & !
         alpha2max_low_diatoms_in , & ! light limitation (1/(W/m^2))  
         alpha2max_low_sp_in      , & 
         alpha2max_low_phaeo_in   , & 
         beta2max_diatoms_in , & ! light inhibition (1/(W/m^2))  
         beta2max_sp_in      , & 
         beta2max_phaeo_in   , & 
         mu_max_diatoms_in   , & ! maximum growth rate (1/day)       
         mu_max_sp_in        , & 
         mu_max_phaeo_in     , & 
         grow_Tdep_diatoms_in, & ! Temperature dependence of growth (1/C)
         grow_Tdep_sp_in     , & 
         grow_Tdep_phaeo_in  , & 
         fr_graze_diatoms_in , & ! Fraction grazed
         fr_graze_sp_in      , & 
         fr_graze_phaeo_in   , & 
         mort_pre_diatoms_in , & ! Mortality (1/day)
         mort_pre_sp_in      , & 
         mort_pre_phaeo_in   , & 
         mort_Tdep_diatoms_in, & ! T dependence of mortality (1/C) 
         mort_Tdep_sp_in     , &  
         mort_Tdep_phaeo_in  , &  
         k_exude_diatoms_in  , & ! algal exudation (1/d)
         k_exude_sp_in       , &  
         k_exude_phaeo_in    , &  
         K_Nit_diatoms_in    , & ! nitrate half saturation (mmol/m^3)
         K_Nit_sp_in         , &  
         K_Nit_phaeo_in      , &  
         K_Am_diatoms_in     , & ! ammonium half saturation (mmol/m^3)
         K_Am_sp_in          , &   
         K_Am_phaeo_in       , &   
         K_Sil_diatoms_in    , & ! silicate half saturation (mmol/m^3)
         K_Sil_sp_in         , &   
         K_Sil_phaeo_in      , &   
         K_Fe_diatoms_in     , & ! iron half saturation (nM)
         K_Fe_sp_in          , &   
         K_Fe_phaeo_in       , &    
         f_don_protein_in    , & ! fraction of spilled grazing to proteins            
         kn_bac_protein_in   , & ! Bacterial degredation of DON (1/d)                  
         f_don_Am_protein_in , & ! fraction of remineralized DON to ammonium          
         f_doc_s_in          , & ! fraction of mortality to DOC 
         f_doc_l_in          , &   
         f_exude_s_in        , & ! fraction of exudation to DOC
         f_exude_l_in        , & 
         k_bac_s_in          , & ! Bacterial degredation of DOC (1/d)
         k_bac_l_in          , & 
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
         t_sk_ox_in          , &   ! DMS oxidation time (d)
         algaltype_diatoms_in  , & ! mobility type
         algaltype_sp_in       , & !
         algaltype_phaeo_in    , & !
         nitratetype_in        , & !
         ammoniumtype_in       , & !
         silicatetype_in       , & !
         dmspptype_in          , & !
         dmspdtype_in          , & !
         humtype_in            , & !
         doctype_s_in          , & !
         doctype_l_in          , & !
         dontype_protein_in    , & !
         fedtype_1_in          , & !
         feptype_1_in          , & !
         zaerotype_bc1_in      , & !
         zaerotype_bc2_in      , & !
         zaerotype_dust1_in    , & !
         zaerotype_dust2_in    , & !
         zaerotype_dust3_in    , & !
         zaerotype_dust4_in    , & !
         ratio_C2N_diatoms_in  , & ! algal C to N ratio (mol/mol)
         ratio_C2N_sp_in       , & !
         ratio_C2N_phaeo_in    , & !
         ratio_chl2N_diatoms_in, & ! algal chlorophyll to N ratio (mg/mmol)
         ratio_chl2N_sp_in     , & !
         ratio_chl2N_phaeo_in  , & !
         F_abs_chl_diatoms_in  , & ! scales absorbed radiation for dEdd
         F_abs_chl_sp_in       , & !
         F_abs_chl_phaeo_in    , & !
         ratio_C2N_proteins_in     ! ratio of C to N in proteins (mol/mol)       

     !logical (kind=log_kind), intent(in) :: & 
     !   restore_bgc_in      ! if true, restore nitrate

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

        real (kind=dbl_kind), intent(in) :: &
             hs0_in             ! snow depth for transition to bare sea ice (m)
        
        ! level-ice ponds
        character (len=char_len), intent(in) :: &
             frzpnd_in          ! pond refreezing parameterization
        
        real (kind=dbl_kind), intent(in) :: &
             dpscale_in, &      ! alter e-folding time scale for flushing 
             rfracmin_in, &     ! minimum retained fraction of meltwater
             rfracmax_in, &     ! maximum retained fraction of meltwater
             pndaspect_in, &    ! ratio of pond depth to pond fraction
             hs1_in             ! tapering parameter for snow on pond ice
        
        ! topo ponds
        real (kind=dbl_kind), intent(in) :: &
             hp1_in             ! critical parameter for pond ice thickness
        
        ktherm = ktherm_in
        conduct = conduct_in
        fbot_xfer_type = fbot_xfer_type_in
        calc_Tsfc = calc_Tsfc_in
        ustar_min = ustar_min_in
        a_rapid_mode = a_rapid_mode_in
        Rac_rapid_mode = Rac_rapid_mode_in
        aspect_rapid_mode = aspect_rapid_mode_in
        dSdt_slow_mode = dSdt_slow_mode_in
        phi_c_slow_mode = phi_c_slow_mode_in
        phi_i_mushy = phi_i_mushy_in
        shortwave = shortwave_in
        albedo_type = albedo_type_in
        albicev = albicev_in
        albicei = albicei_in
        albsnowv = albsnowv_in
        albsnowi = albsnowi_in
        ahmax = ahmax_in
        R_ice = R_ice_in
        R_pnd = R_pnd_in
        R_snw = R_snw_in
        dT_mlt = dT_mlt_in
        rsnw_mlt = rsnw_mlt_in
        kalg = kalg_in
        kstrength = kstrength_in
        krdg_partic = krdg_partic_in
        krdg_redist = krdg_redist_in
        mu_rdg = mu_rdg_in
        Cf = Cf_in
        atmbndy = atmbndy_in
        calc_strair = calc_strair_in
        formdrag = formdrag_in
        highfreq = highfreq_in
        natmiter = natmiter_in
        oceanmixed_ice = oceanmixed_ice_in
        tfrz_option = tfrz_option_in
        kitd = kitd_in
        kcatbound = kcatbound_in
        hs0 = hs0_in
        frzpnd = frzpnd_in
        dpscale = dpscale_in
        rfracmin = rfracmin_in
        rfracmax = rfracmax_in
        pndaspect = pndaspect_in
        hs1 = hs1_in
        hp1 = hp1_in
     !  bgc_data_dir = bgc_data_dir_in
     !  sil_data_type= sil_data_type_in
     !  nit_data_type = nit_data_type_in
     !  fe_data_type = fe_data_type_in
        bgc_flux_type = bgc_flux_type_in
        z_tracers = z_tracers_in
        scale_bgc = scale_bgc_in
        solve_zbgc = solve_zbgc_in
        dEdd_algae = dEdd_algae_in
        skl_bgc = skl_bgc_in
        grid_o = grid_o_in
        l_sk = l_sk_in
        grid_o_t = grid_o_t_in
        initbio_frac = initbio_frac_in
        frazil_scav = frazil_scav_in
        grid_oS = grid_oS_in
        l_skS = l_skS_in
        phi_snow = phi_snow_in
     !  restore_bgc = restore_bgc_in
        ratio_Si2N_diatoms= ratio_Si2N_diatoms_in 
        ratio_Si2N_sp     = ratio_Si2N_sp_in
        ratio_Si2N_phaeo  = ratio_Si2N_phaeo_in
        ratio_S2N_diatoms = ratio_S2N_diatoms_in
        ratio_S2N_sp      = ratio_S2N_sp_in
        ratio_S2N_phaeo   = ratio_S2N_phaeo_in
        ratio_Fe2C_diatoms= ratio_Fe2C_diatoms_in 
        ratio_Fe2C_sp     = ratio_Fe2C_sp_in
        ratio_Fe2C_phaeo  = ratio_Fe2C_phaeo_in
        ratio_Fe2N_diatoms= ratio_Fe2N_diatoms_in 
        ratio_Fe2N_sp     = ratio_Fe2N_sp_in
        ratio_Fe2N_phaeo  = ratio_Fe2N_phaeo_in
        ratio_Fe2DON      = ratio_Fe2DON_in
        ratio_Fe2DOC_s    = ratio_Fe2DOC_s_in
        ratio_Fe2DOC_l    = ratio_Fe2DOC_l_in
        fr_resp           = fr_resp_in
        tau_min           = tau_min_in
        tau_max           = tau_max_in
        algal_vel         = algal_vel_in
        R_dFe2dust        = R_dFe2dust_in
        dustFe_sol        = dustFe_sol_in
        chlabs_diatoms    = chlabs_diatoms_in
        chlabs_sp         = chlabs_sp_in
        chlabs_phaeo      = chlabs_phaeo_in
        alpha2max_low_diatoms = alpha2max_low_diatoms_in
        alpha2max_low_sp      = alpha2max_low_sp_in
        alpha2max_low_phaeo   = alpha2max_low_phaeo_in
        beta2max_diatoms = beta2max_diatoms_in
        beta2max_sp      = beta2max_sp_in
        beta2max_phaeo   = beta2max_phaeo_in
        mu_max_diatoms   = mu_max_diatoms_in
        mu_max_sp        = mu_max_sp_in
        mu_max_phaeo     = mu_max_phaeo_in
        grow_Tdep_diatoms= grow_Tdep_diatoms_in
        grow_Tdep_sp     = grow_Tdep_sp_in
        grow_Tdep_phaeo  = grow_Tdep_phaeo_in
        fr_graze_diatoms = fr_graze_diatoms_in
        fr_graze_sp      = fr_graze_sp_in
        fr_graze_phaeo   = fr_graze_phaeo_in
        mort_pre_diatoms = mort_pre_diatoms_in
        mort_pre_sp      = mort_pre_sp_in
        mort_pre_phaeo   = mort_pre_phaeo_in
        mort_Tdep_diatoms= mort_Tdep_diatoms_in
        mort_Tdep_sp     = mort_Tdep_sp_in
        mort_Tdep_phaeo  = mort_Tdep_phaeo_in
        k_exude_diatoms  = k_exude_diatoms_in
        k_exude_sp       = k_exude_sp_in
        k_exude_phaeo    = k_exude_phaeo_in
        K_Nit_diatoms    = K_Nit_diatoms_in
        K_Nit_sp         = K_Nit_sp_in
        K_Nit_phaeo      = K_Nit_phaeo_in
        K_Am_diatoms     = K_Am_diatoms_in
        K_Am_sp          = K_Am_sp_in
        K_Am_phaeo       = K_Am_phaeo_in
        K_Sil_diatoms    = K_Sil_diatoms_in
        K_Sil_sp         = K_Sil_sp_in
        K_Sil_phaeo      = K_Sil_phaeo_in
        K_Fe_diatoms     = K_Fe_diatoms_in
        K_Fe_sp          = K_Fe_sp_in
        K_Fe_phaeo       = K_Fe_phaeo_in
        f_don_protein    = f_don_protein_in
        kn_bac_protein   = kn_bac_protein_in
        f_don_Am_protein = f_don_Am_protein_in
        f_doc_s          = f_doc_s_in
        f_doc_l          = f_doc_l_in
        f_exude_s        = f_exude_s_in
        f_exude_l        = f_exude_l_in
        k_bac_s          = k_bac_s_in
        k_bac_l          = k_bac_l_in
        T_max            = T_max_in
        fsal             = fsal_in
        op_dep_min       = op_dep_min_in
        fr_graze_s       = fr_graze_s_in
        fr_graze_e       = fr_graze_e_in
        fr_mort2min      = fr_mort2min_in
        fr_dFe           = fr_dFe_in
        k_nitrif         = k_nitrif_in
        t_iron_conv      = t_iron_conv_in
        max_loss         = max_loss_in
        max_dfe_doc1     = max_dfe_doc1_in
        fr_resp_s        = fr_resp_s_in
        y_sk_DMS         = y_sk_DMS_in
        t_sk_conv        = t_sk_conv_in
        t_sk_ox          = t_sk_ox_in
        algaltype_diatoms  = algaltype_diatoms_in
        algaltype_sp       = algaltype_sp_in
        algaltype_phaeo    = algaltype_phaeo_in
        nitratetype        = nitratetype_in
        ammoniumtype       = ammoniumtype_in
        silicatetype       = silicatetype_in
        dmspptype          = dmspptype_in
        dmspdtype          = dmspdtype_in
        humtype            = humtype_in
        doctype_s          = doctype_s_in
        doctype_l          = doctype_l_in
        dontype_protein    = dontype_protein_in
        fedtype_1          = fedtype_1_in
        feptype_1          = feptype_1_in
        zaerotype_bc1      = zaerotype_bc1_in
        zaerotype_bc2      = zaerotype_bc2_in
        zaerotype_dust1    = zaerotype_dust1_in
        zaerotype_dust2    = zaerotype_dust2_in
        zaerotype_dust3    = zaerotype_dust3_in
        zaerotype_dust4    = zaerotype_dust4_in
        ratio_C2N_diatoms  = ratio_C2N_diatoms_in
        ratio_C2N_sp       = ratio_C2N_sp_in
        ratio_C2N_phaeo    = ratio_C2N_phaeo_in
        ratio_chl2N_diatoms= ratio_chl2N_diatoms_in
        ratio_chl2N_sp     = ratio_chl2N_sp_in
        ratio_chl2N_phaeo  = ratio_chl2N_phaeo_in
        F_abs_chl_diatoms  = F_abs_chl_diatoms_in
        F_abs_chl_sp       = F_abs_chl_sp_in
        F_abs_chl_phaeo    = F_abs_chl_phaeo_in
        ratio_C2N_proteins = ratio_C2N_proteins_in

      end subroutine icepack_init_parameters

!=======================================================================
! set tracer active flags

      subroutine icepack_init_tracer_flags(&
           tr_iage_in      , & ! if .true., use age tracer
           tr_FY_in        , & ! if .true., use first-year area tracer
           tr_lvl_in       , & ! if .true., use level ice tracer
           tr_pond_in      , & ! if .true., use melt pond tracer
           tr_pond_cesm_in , & ! if .true., use cesm pond tracer
           tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
           tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
           tr_aero_in      , & ! if .true., use aerosol tracers
           tr_brine_in     , & ! if .true., brine height differs from ice thickness
           tr_bgc_S_in     , & ! if .true., use zsalinity
           tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
           tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
           tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
           tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
           tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
           tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
           tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
           tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
           tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
           tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
           tr_bgc_hum_in   , & ! if .true., hum as tracer 
           tr_bgc_PON_in)      ! if .true., PON as product tracer 


        use icepack_intfc_tracers, only: &
             tr_iage      , & ! if .true., use age tracer
             tr_FY        , & ! if .true., use first-year area tracer
             tr_lvl       , & ! if .true., use level ice tracer
             tr_pond      , & ! if .true., use melt pond tracer
             tr_pond_cesm , & ! if .true., use cesm pond tracer
             tr_pond_lvl  , & ! if .true., use level-ice pond tracer
             tr_pond_topo , & ! if .true., use explicit topography-based ponds
             tr_aero      , & ! if .true., use aerosol tracers
             tr_brine     , & ! if .true., brine height differs from ice thickness
             tr_bgc_S     , & ! if .true., use zsalinity
             tr_zaero     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe    , & ! if .true., Fe as product tracer 
             tr_bgc_hum   , & ! if .true., hum as product tracer 
             tr_bgc_PON       ! if .true., PON as product tracer 


        logical, intent(in) :: &
             tr_iage_in      , & ! if .true., use age tracer
             tr_FY_in        , & ! if .true., use first-year area tracer
             tr_lvl_in       , & ! if .true., use level ice tracer
             tr_pond_in      , & ! if .true., use melt pond tracer
             tr_pond_cesm_in , & ! if .true., use cesm pond tracer
             tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
             tr_aero_in      , & ! if .true., use aerosol tracers
             tr_brine_in     , & ! if .true., brine height differs from ice thickness
             tr_bgc_S_in     , & ! if .true., use zsalinity
             tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
             tr_bgc_hum_in   , & ! if .true., hum as product tracer 
             tr_bgc_PON_in       ! if .true., PON as product tracer 

        tr_iage      = tr_iage_in
        tr_FY        = tr_FY_in
        tr_lvl       = tr_lvl_in
        tr_pond      = tr_pond_in
        tr_pond_cesm = tr_pond_cesm_in
        tr_pond_lvl  = tr_pond_lvl_in
        tr_pond_topo = tr_pond_topo_in
        tr_aero      = tr_aero_in
        tr_brine     = tr_brine_in
        tr_bgc_S     = tr_bgc_S_in
        tr_zaero     = tr_zaero_in 
        tr_bgc_Nit   = tr_bgc_Nit_in
        tr_bgc_N     = tr_bgc_N_in 
        tr_bgc_DON   = tr_bgc_DON_in
        tr_bgc_C     = tr_bgc_C_in 
        tr_bgc_chl   = tr_bgc_chl_in
        tr_bgc_Am    = tr_bgc_Am_in
        tr_bgc_Sil   = tr_bgc_Sil_in
        tr_bgc_DMS   = tr_bgc_DMS_in
        tr_bgc_Fe    = tr_bgc_Fe_in 
        tr_bgc_hum   = tr_bgc_hum_in
        tr_bgc_PON   = tr_bgc_PON_in 

      end subroutine icepack_init_tracer_flags

!=======================================================================

      subroutine icepack_init_tracer_indices(&
           nt_Tsfc_in, & ! ice/snow temperature
           nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
           nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
           nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
           nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
           nt_iage_in, & ! volume-weighted ice age
           nt_FY_in, & ! area-weighted first-year ice area
           nt_alvl_in, & ! level ice area fraction
           nt_vlvl_in, & ! level ice volume fraction
           nt_apnd_in, & ! melt pond area fraction
           nt_hpnd_in, & ! melt pond depth
           nt_ipnd_in, & ! melt pond refrozen lid thickness
           nt_aero_in, & ! starting index for aerosols in ice 
           nt_zaero_in,   & !  black carbon and other aerosols
           nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
           nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
           nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
           nt_bgc_DOC_in, & !  dissolved organic carbon
           nt_bgc_DON_in, & !  dissolved organic nitrogen
           nt_bgc_DIC_in, & !  dissolved inorganic carbon
           nt_bgc_Fed_in, & !  dissolved iron
           nt_bgc_Fep_in, & !  particulate iron
           nt_bgc_Nit_in, & ! nutrients  
           nt_bgc_Am_in,  & ! 
           nt_bgc_Sil_in, & !
           nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
           nt_bgc_DMSPd_in,&! 
           nt_bgc_DMS_in, & ! 
           nt_bgc_hum_in, & ! 
           nt_bgc_PON_in, & ! zooplankton and detritus  
           nlt_zaero_in,  & !  black carbon and other aerosols
           nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
           nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
           nlt_bgc_chl_in,& ! diatoms, phaeocystis, pico/small 
           nlt_bgc_DOC_in,& !  dissolved organic carbon
           nlt_bgc_DON_in,& !  dissolved organic nitrogen
           nlt_bgc_DIC_in,& !  dissolved inorganic carbon
           nlt_bgc_Fed_in,& !  dissolved iron
           nlt_bgc_Fep_in,& !  particulate iron
           nlt_bgc_Nit_in,& ! nutrients  
           nlt_bgc_Am_in, & ! 
           nlt_bgc_Sil_in,& !
           nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
           nlt_bgc_DMSPd_in,&! 
           nlt_bgc_DMS_in,& ! 
           nlt_bgc_hum_in,& ! 
           nlt_bgc_PON_in,& ! zooplankton and detritus  
           nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
           nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
           nlt_chl_sw_in, & ! points to total chla in trcrn_sw
           nlt_zaero_sw_in,&! black carbon and dust in trcrn_sw
                            ! Index Dimensions: 
           n_algae, n_algalC, & !
           n_algalchl, n_DOC, & !
           n_DON,n_DIC,n_dFe, & !
           n_pFe, n_aerosols, & !
           bio_index_o_in,    & ! nlt index to fixed data array
           bio_index_in,      & ! nlt index to nt index
           nbtrcr)

        use icepack_intfc_tracers, only: &
             nt_Tsfc, & ! ice/snow temperature
             nt_qice, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno, & ! volume-weighted snow enthalpy (in layers)
             nt_sice, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage, & ! volume-weighted ice age
             nt_FY, & ! area-weighted first-year ice area
             nt_alvl, & ! level ice area fraction
             nt_vlvl, & ! level ice volume fraction
             nt_apnd, & ! melt pond area fraction
             nt_hpnd, & ! melt pond depth
             nt_ipnd, & ! melt pond refrozen lid thickness
             nt_aero, & ! starting index for aerosols in ice
             nt_zaero,   & !  black carbon and other aerosols
             nt_bgc_N ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl, & ! diatoms, phaeocystis, pico/small 
             nt_bgc_DOC, & !  dissolved organic carbon
             nt_bgc_DON, & !  dissolved organic nitrogen
             nt_bgc_DIC, & !  dissolved inorganic carbon
             nt_bgc_Fed, & !  dissolved iron
             nt_bgc_Fep, & !  particulate iron
             nt_bgc_Nit, & ! nutrients  
             nt_bgc_Am,  & ! 
             nt_bgc_Sil, & !
             nt_bgc_DMSPp,&! trace gases (skeletal layer)
             nt_bgc_DMSPd,&! 
             nt_bgc_DMS, & ! 
             nt_bgc_hum, & ! 
             nt_bgc_PON, & ! zooplankton and detritus 
             nlt_zaero,  & !  black carbon and other aerosols
             nlt_bgc_N , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl,& ! diatoms, phaeocystis, pico/small 
             nlt_bgc_DOC,& !  dissolved organic carbon
             nlt_bgc_DON,& !  dissolved organic nitrogen
             nlt_bgc_DIC,& !  dissolved inorganic carbon
             nlt_bgc_Fed,& !  dissolved iron
             nlt_bgc_Fep,& !  particulate iron
             nlt_bgc_Nit,& ! nutrients  
             nlt_bgc_Am, & ! 
             nlt_bgc_Sil,& !
             nlt_bgc_DMSPp,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd,&! 
             nlt_bgc_DMS,& ! 
             nlt_bgc_hum,& ! 
             nlt_bgc_PON,& ! zooplankton and detritus   
             nt_zbgc_frac,&! fraction of tracer in the mobile phase
             nt_bgc_S,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw, & ! points to total chla in trcrn_sw
             nlt_zaero_sw,&! black carbon and dust in trcrn_sw
             bio_index_o,& !
             bio_index 
        
        integer, intent(in) :: &
             nt_Tsfc_in, & ! ice/snow temperature
             nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_in, & ! volume-weighted ice age
             nt_FY_in, & ! area-weighted first-year ice area
             nt_alvl_in, & ! level ice area fraction
             nt_vlvl_in, & ! level ice volume fraction
             nt_apnd_in, & ! melt pond area fraction
             nt_hpnd_in, & ! melt pond depth
             nt_ipnd_in, & ! melt pond refrozen lid thickness
             nt_aero_in, & ! starting index for aerosols in ice
             nt_bgc_Nit_in, & ! nutrients  
             nt_bgc_Am_in,  & ! 
             nt_bgc_Sil_in, & !
             nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_in,&! 
             nt_bgc_DMS_in, & ! 
             nt_bgc_hum_in, & ! 
             nt_bgc_PON_in, & ! zooplankton and detritus   
             nlt_bgc_Nit_in,& ! nutrients  
             nlt_bgc_Am_in, & ! 
             nlt_bgc_Sil_in,& !
             nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_in,&! 
             nlt_bgc_DMS_in,& ! 
             nlt_bgc_hum_in,& ! 
             nlt_bgc_PON_in,& ! zooplankton and detritus  
             nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
             nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_in    ! points to total chla in trcrn_sw

       integer, intent(in) :: &
             n_algae,    & !  Dimensions
             n_algalC,   & !
             n_algalchl, & !
             n_DOC,      & !
             n_DON,      & !
             n_DIC,      & !
             n_dFe,      & !
             n_pFe,      & ! 
             n_aerosols, & !
             nbtrcr

        integer (kind=int_kind), dimension(:), intent(in) :: &
             bio_index_o_in, & 
             bio_index_in  

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_in   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DOC_in, & !  dissolved organic carbon
             nlt_bgc_DOC_in   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DON_in, & !  dissolved organic nitrogen
             nlt_bgc_DON_in   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DIC_in, & ! dissolved inorganic carbon
             nlt_bgc_DIC_in   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_Fed_in, & !  dissolved iron
             nt_bgc_Fep_in, & !  particulate iron
             nlt_bgc_Fed_in,& !  dissolved iron
             nlt_bgc_Fep_in   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_zaero_in,   & !  black carbon and other aerosols
             nlt_zaero_in,  & !  black carbon and other aerosols
             nlt_zaero_sw_in  ! black carbon and dust in trcrn_sw

        ! local
        integer (kind=int_kind) :: k

        nt_Tsfc = nt_Tsfc_in
        nt_qice = nt_qice_in
        nt_qsno = nt_qsno_in
        nt_sice = nt_sice_in
        nt_fbri = nt_fbri_in
        nt_iage = nt_iage_in
        nt_FY = nt_FY_in
        nt_alvl = nt_alvl_in
        nt_vlvl = nt_vlvl_in
        nt_apnd = nt_apnd_in
        nt_hpnd = nt_hpnd_in
        nt_ipnd = nt_ipnd_in
        nt_aero = nt_aero_in
        nt_bgc_Nit = nt_bgc_Nit_in
        nt_bgc_Am  = nt_bgc_Am_in
        nt_bgc_Sil = nt_bgc_Sil_in
        nt_bgc_DMSPp=nt_bgc_DMSPp_in
        nt_bgc_DMSPd=nt_bgc_DMSPd_in
        nt_bgc_DMS = nt_bgc_DMS_in
        nt_bgc_hum = nt_bgc_hum_in
        nt_bgc_PON = nt_bgc_PON_in
        nlt_bgc_Nit = nlt_bgc_Nit_in
        nlt_bgc_Am  = nlt_bgc_Am_in
        nlt_bgc_Sil = nlt_bgc_Sil_in
        nlt_bgc_DMSPp=nlt_bgc_DMSPp_in
        nlt_bgc_DMSPd=nlt_bgc_DMSPd_in
        nlt_bgc_DMS = nlt_bgc_DMS_in
        nlt_bgc_hum = nlt_bgc_hum_in
        nlt_bgc_PON = nlt_bgc_PON_in
        nlt_chl_sw  = nlt_chl_sw_in
        nt_zbgc_frac=nt_zbgc_frac_in
        nt_bgc_S   = nt_bgc_S_in

        nt_bgc_N(:)    = 0
        nt_bgc_C(:)    = 0
        nt_bgc_chl(:)  = 0
        nlt_bgc_N(:)   = 0
        nlt_bgc_C(:)   = 0
        nlt_bgc_chl(:) = 0
        nt_bgc_DOC(:)  = 0
        nlt_bgc_DOC(:) = 0
        nt_bgc_DIC(:)  = 0
        nlt_bgc_DIC(:) = 0
        nt_bgc_DON(:)  = 0
        nlt_bgc_DON(:) = 0
        nt_bgc_Fed(:)  = 0
        nt_bgc_Fep(:)  = 0
        nlt_bgc_Fed(:) = 0
        nlt_bgc_Fep(:) = 0
        nt_zaero(:)    = 0
        nlt_zaero(:)   = 0
        nlt_zaero_sw(:)= 0
        bio_index(:)   = 0
        bio_index_o(:) = 0

        do k = 1, nbtrcr
           bio_index_o(k)= bio_index_o_in(k)
           bio_index(k)  = bio_index_in(k)
        enddo
        do k = 1, n_algae
           nt_bgc_N(k) = nt_bgc_N_in(k) 
           nlt_bgc_N(k)= nlt_bgc_N_in(k) 
        enddo
        do k = 1, n_algalC
           nt_bgc_C(k) = nt_bgc_C_in(k) 
           nlt_bgc_C(k)= nlt_bgc_C_in(k) 
        enddo
        do k = 1, n_algalchl
           nt_bgc_chl(k) = nt_bgc_chl_in(k) 
           nlt_bgc_chl(k)= nlt_bgc_chl_in(k) 
        enddo
        do k = 1, n_DOC
           nt_bgc_DOC(k) = nt_bgc_DOC_in(k) 
           nlt_bgc_DOC(k)= nlt_bgc_DOC_in(k) 
        enddo
        do k = 1, n_DON
           nt_bgc_DON(k) = nt_bgc_DON_in(k) 
           nlt_bgc_DON(k)= nlt_bgc_DON_in(k) 
        enddo
        do k = 1, n_DIC
           nt_bgc_DIC(k) = nt_bgc_DIC_in(k) 
           nlt_bgc_DIC(k)= nlt_bgc_DIC_in(k) 
        enddo
        do k = 1, n_dFe  
           nt_bgc_Fed(k) = nt_bgc_Fed_in(k) 
           nlt_bgc_Fed(k)= nlt_bgc_Fed_in(k) 
        enddo
        do k = 1, n_pFe  
           nt_bgc_Fep(k) = nt_bgc_Fep_in(k) 
           nlt_bgc_Fep(k)= nlt_bgc_Fep_in(k) 
        enddo
        do k = 1, n_aerosols
           nt_zaero(k)    = nt_zaero_in(k)   
           nlt_zaero(k)   = nlt_zaero_in(k)   
           nlt_zaero_sw(k)= nlt_zaero_sw_in(k)   
        enddo

      end subroutine icepack_init_tracer_indices

!=======================================================================
! set the number of column tracers

      subroutine icepack_init_tracer_numbers(&
         ntrcr_in, nbtrcr_in, nbtrcr_sw_in)

      use icepack_intfc_tracers, only: &
         ntrcr, nbtrcr, nbtrcr_sw

      integer (kind=int_kind), intent(in) :: &
         ntrcr_in  , &! number of tracers in use
         nbtrcr_in , &! number of bio tracers in use
         nbtrcr_sw_in ! number of shortwave bio tracers in use
        
         ntrcr     = ntrcr_in
         nbtrcr    = nbtrcr_in
         nbtrcr_sw = nbtrcr_sw_in

      end subroutine icepack_init_tracer_numbers

!=======================================================================

      end module icepack_intfc

!=======================================================================
