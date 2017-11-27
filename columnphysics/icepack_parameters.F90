!  SVN:$Id: icepack_parameters.F90 1226 2017-05-22 22:45:03Z tcraig $
!=========================================================================
!
! flags for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_parameters

      use icepack_kinds
      use icepack_constants, only: c3, c0, c1, p5, p1

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

      real (kind=dbl_kind), parameter, public :: &
         saltmax = 3.2_dbl_kind,   & ! max salinity at ice base for BL99 (ppt)
         ! phi_init and dSin0_frazil are used for mushy thermo, ktherm=2
         phi_init = 0.75_dbl_kind, & ! initial liquid fraction of frazil
         min_salin = p1          , & ! threshold for brine pocket treatment 
         salt_loss =0.4_dbl_kind, &  ! fraction of salt retained in zsalinity 
         min_bgc        = 0.01_dbl_kind, & ! fraction of ocean bgc concentration in surface melt 
         dSin0_frazil = c3 ! bulk salinity reduction of newly formed frazil

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

      real (kind=dbl_kind), parameter, public :: &
         hi_ssl = 0.050_dbl_kind, & ! ice surface scattering layer thickness (m)
         hs_ssl = 0.040_dbl_kind    ! snow surface scattering layer thickness (m)

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
         mu_rdg, &      ! gives e-folding scale of ridged ice (m^.5) 
                        ! (krdg_redist = 1) 
         Cf             ! ratio of ridging work to PE change in ridging (kstrength = 1)

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

!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

      logical (kind=log_kind), public :: &
         oceanmixed_ice           ! if true, use ocean mixed layer

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

      !-----------------------------------------------------------------
      ! dimensions
      !-----------------------------------------------------------------
      integer (kind=int_kind), parameter, public :: &
         max_algae  =   3       , & ! maximum number of algal types 
         max_dic    =   1       , & ! maximum number of dissolved inorganic carbon types 
         max_doc    =   3       , & ! maximum number of dissolved organic carbon types
         max_don    =   1       , & ! maximum number of dissolved organic nitrogen types
         max_fe     =   2       , & ! maximum number of iron types
         nmodal1    =   10      , & ! dimension for modal aerosol radiation parameters
         nmodal2    =   8       , & ! dimension for modal aerosol radiation parameters
         max_aero   =   6       , & ! maximum number of aerosols 
         max_nbtrcr = max_algae*2 & ! algal nitrogen and chlorophyll
                    + max_dic     & ! dissolved inorganic carbon
                    + max_doc     & ! dissolved organic carbon
                    + max_don     & ! dissolved organic nitrogen
                    + 5           & ! nitrate, ammonium, silicate, PON, and humics
                    + 3           & ! DMSPp, DMSPd, DMS
                    + max_fe*2    & ! dissolved Fe and  particulate Fe
                    + max_aero      ! aerosols

      !-----------------------------------------------------------------
      ! namelist
      !-----------------------------------------------------------------
      character(char_len_long), public :: & 
         bgc_data_dir   ! directory for biogeochemistry data

      character(char_len), public :: &          
         sil_data_type  , & ! 'default', 'clim'
         nit_data_type  , & ! 'default', 'clim'   
         fe_data_type   , & ! 'default', 'clim'      
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
         grid_o_t    , & ! top grid point length scale 
         phi_snow    , & ! porosity of snow
         initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

      real (kind=dbl_kind), public :: & 
         grid_oS     , & ! for bottom flux (zsalinity)
         l_skS           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)

      logical (kind=log_kind), public :: & 
         restore_bgc      ! if true, restore nitrate

      !-----------------------------------------------------------------
      ! From icepack_zbgc_shared.F90
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         ratio_Si2N_diatoms, &   ! algal Si to N (mol/mol)
         ratio_Si2N_sp     , &
         ratio_Si2N_phaeo  , &
         ratio_S2N_diatoms , &   ! algal S  to N (mol/mol)
         ratio_S2N_sp      , &
         ratio_S2N_phaeo   , &
         ratio_Fe2C_diatoms, &   ! algal Fe to C  (umol/mol)
         ratio_Fe2C_sp     , &
         ratio_Fe2C_phaeo  , &
         ratio_Fe2N_diatoms, &   ! algal Fe to N  (umol/mol)
         ratio_Fe2N_sp     , &
         ratio_Fe2N_phaeo  , &
         ratio_Fe2DON      , &   ! Fe to N of DON (nmol/umol)
         ratio_Fe2DOC_s    , &   ! Fe to C of DOC (nmol/umol) saccharids
         ratio_Fe2DOC_l    , &   ! Fe to C of DOC (nmol/umol) lipids
         fr_resp           , &   ! fraction of algal growth lost due to respiration        
         tau_min           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
         tau_max           , &   ! long time mobile to stationary exchanges (s) = 2 days
         algal_vel         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol              ! solubility fraction

      !-----------------------------------------------------------------
      ! From algal_dyn in icepack_algae.F90
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         chlabs_diatoms   , & ! chl absorption (1/m/(mg/m^3))
         chlabs_sp        , & !
         chlabs_phaeo     , & !
         alpha2max_low_diatoms , & ! light limitation (1/(W/m^2))  
         alpha2max_low_sp      , & 
         alpha2max_low_phaeo   , & 
         beta2max_diatoms , & ! light inhibition (1/(W/m^2))  
         beta2max_sp      , & 
         beta2max_phaeo   , & 
         mu_max_diatoms   , & ! maximum growth rate (1/day)       
         mu_max_sp        , & 
         mu_max_phaeo     , & 
         grow_Tdep_diatoms, & ! Temperature dependence of growth (1/C)
         grow_Tdep_sp     , & 
         grow_Tdep_phaeo  , & 
         fr_graze_diatoms , & ! Fraction grazed
         fr_graze_sp      , & 
         fr_graze_phaeo   , & 
         mort_pre_diatoms , & ! Mortality (1/day)
         mort_pre_sp      , & 
         mort_pre_phaeo   , & 
         mort_Tdep_diatoms, & ! T dependence of mortality (1/C)
         mort_Tdep_sp     , &  
         mort_Tdep_phaeo  , &  
         k_exude_diatoms  , & ! algal exudation (1/d)
         k_exude_sp       , &  
         k_exude_phaeo    , &  
         K_Nit_diatoms    , & ! nitrate half saturation (mmol/m^3)
         K_Nit_sp         , &  
         K_Nit_phaeo      , &  
         K_Am_diatoms     , & ! ammonium half saturation (mmol/m^3)
         K_Am_sp          , &   
         K_Am_phaeo       , &   
         K_Sil_diatoms    , & ! silicate half saturation (mmol/m^3)
         K_Sil_sp         , &   
         K_Sil_phaeo      , &   
         K_Fe_diatoms     , & ! iron half saturation (nM)
         K_Fe_sp          , &   
         K_Fe_phaeo       , &    
         f_don_protein    , & ! fraction of spilled grazing to proteins           
         kn_bac_protein   , & ! Bacterial degredation of DON (1/d)                
         f_don_Am_protein , & ! fraction of remineralized DON to ammonium         
         f_doc_s          , & ! fraction of mortality to DOC 
         f_doc_l          , &   
         f_exude_s        , & ! fraction of exudation to DOC
         f_exude_l        , & 
         k_bac_s          , & ! Bacterial degredation of DOC (1/d)
         k_bac_l          , & 
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

      !-----------------------------------------------------------------
      ! former parameters now in namelist
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         algaltype_diatoms  , & ! mobility type
         algaltype_sp       , & !
         algaltype_phaeo    , & !
         nitratetype        , & !
         ammoniumtype       , & !
         silicatetype       , & !
         dmspptype          , & !
         dmspdtype          , & !
         humtype            , & !
         doctype_s          , & !
         doctype_l          , & !
         dontype_protein    , & !
         fedtype_1          , & !
         feptype_1          , & !
         zaerotype_bc1      , & !
         zaerotype_bc2      , & !
         zaerotype_dust1    , & !
         zaerotype_dust2    , & !
         zaerotype_dust3    , & !
         zaerotype_dust4    , & !
         ratio_C2N_diatoms  , & ! algal C to N ratio (mol/mol)
         ratio_C2N_sp       , & !
         ratio_C2N_phaeo    , & !
         ratio_chl2N_diatoms, & ! algal chlorophyll to N ratio (mg/mmol)
         ratio_chl2N_sp     , & !
         ratio_chl2N_phaeo  , & !
         F_abs_chl_diatoms  , & ! scales absorbed radiation for dEdd
         F_abs_chl_sp       , & !
         F_abs_chl_phaeo    , & !
         ratio_C2N_proteins     ! ratio of C to N in proteins (mol/mol)       

      !-----------------------------------------------------------------
      ! Transport type 
      !-----------------------------------------------------------------
      ! In delta Eddington, algal particles are assumed to cause no
      ! significant scattering (Brieglib and Light), only absorption
      ! in the visible spectral band (200-700 nm)
      ! Algal types: Diatoms, flagellates, Phaeocycstis
      ! DOC        : Proteins, EPS, Lipids
      !-----------------------------------------------------------------
      real (kind=dbl_kind), parameter, dimension(max_dic), public :: &
         dictype   = (/-c1/)  ! not in namelist

      real (kind=dbl_kind), dimension(max_algae), public :: &
         algaltype   ! tau_min for both retention and release

      real (kind=dbl_kind), dimension(max_doc), public :: &
         doctype 

      real (kind=dbl_kind), dimension(max_don), public :: &
         dontype  

      real (kind=dbl_kind), dimension(max_fe), public :: &
         fedtype 

      real (kind=dbl_kind), dimension(max_fe), public :: &
         feptype  

      !------------------------------------------------------------
      ! Aerosol order and type should be consistent with order/type 
      ! specified in delta Eddington:  1) hydrophobic black carbon;
      ! 2) hydrophilic black carbon; 3) dust (0.05-0.5 micron);
      ! 4) dust (0.5-1.25 micron); 5) dust (1.25-2.5 micron);
      ! 6) dust (2.5-5 micron) 
      !-------------------------------------------------------------
      real (kind=dbl_kind), dimension(max_aero), public :: &
         zaerotype  

      !-----------------------------------------------------------------
      ! Forcing input, history and diagnostic output
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         rhosi     = 940.0_dbl_kind, & ! average sea ice density
                                       ! Cox and Weeks, 1982: 919-974 kg/m^2
         sk_l      = 0.03_dbl_kind     ! skeletal layer thickness (m)

      real (kind=dbl_kind), dimension(max_algae), public :: &
         R_C2N     ,      & ! algal C to N (mole/mole) 
         R_chl2N   ,      & ! 3 algal chlorophyll to N (mg/mmol)
         F_abs_chl          ! to scale absorption in Dedd

      real (kind=dbl_kind), dimension(max_don), public :: &  ! increase compare to algal R_Fe2C
         R_C2N_DON 

!=======================================================================
!=======================================================================

      contains

!=======================================================================

! subroutine to set the column package internal parameters

      subroutine icepack_init_parameters(   &
           ktherm_in, conduct_in, fbot_xfer_type_in, calc_Tsfc_in, ustar_min_in, a_rapid_mode_in, &
           Rac_rapid_mode_in, aspect_rapid_mode_in, dSdt_slow_mode_in, phi_c_slow_mode_in, &
           phi_i_mushy_in, shortwave_in, albedo_type_in, albicev_in, albicei_in, albsnowv_in, &
           albsnowi_in, ahmax_in, R_ice_in, R_pnd_in, R_snw_in, dT_mlt_in, rsnw_mlt_in, &
           kalg_in, kstrength_in, krdg_partic_in, krdg_redist_in, mu_rdg_in, Cf_in, &
           atmbndy_in, calc_strair_in, formdrag_in, highfreq_in, natmiter_in, &
           oceanmixed_ice_in, tfrz_option_in, kitd_in, kcatbound_in, hs0_in, frzpnd_in, &
           dpscale_in, rfracmin_in, rfracmax_in, pndaspect_in, hs1_in, hp1_in, &
         ! bgc_data_dir_in, sil_data_type_in, nit_data_type_in, fe_data_type_in, &
           bgc_flux_type_in, z_tracers_in, scale_bgc_in, solve_zbgc_in, dEdd_algae_in, &
           modal_aero_in, skl_bgc_in, solve_zsal_in, grid_o_in, l_sk_in, &
           grid_o_t_in, initbio_frac_in, frazil_scav_in, grid_oS_in, l_skS_in, &
           phi_snow_in, ratio_Si2N_diatoms_in, ratio_Si2N_sp_in, ratio_Si2N_phaeo_in, &
           ratio_S2N_diatoms_in, ratio_S2N_sp_in, ratio_S2N_phaeo_in, ratio_Fe2C_diatoms_in, & 
           ratio_Fe2C_sp_in, ratio_Fe2C_phaeo_in, ratio_Fe2N_diatoms_in, ratio_Fe2N_sp_in, &     
           ratio_Fe2N_phaeo_in, ratio_Fe2DON_in, ratio_Fe2DOC_s_in, ratio_Fe2DOC_l_in, &     
           fr_resp_in, tau_min_in, tau_max_in, algal_vel_in, R_dFe2dust_in, dustFe_sol_in, &         
           chlabs_diatoms_in, chlabs_sp_in, chlabs_phaeo_in, alpha2max_low_diatoms_in, &
           alpha2max_low_sp_in, alpha2max_low_phaeo_in, beta2max_diatoms_in, beta2max_sp_in, &       
           beta2max_phaeo_in, mu_max_diatoms_in, mu_max_sp_in, mu_max_phaeo_in, &      
           grow_Tdep_diatoms_in, grow_Tdep_sp_in, grow_Tdep_phaeo_in, &   
           fr_graze_diatoms_in, fr_graze_sp_in, fr_graze_phaeo_in, &    
           mort_pre_diatoms_in, mort_pre_sp_in, mort_pre_phaeo_in, &    
           mort_Tdep_diatoms_in, mort_Tdep_sp_in, mort_Tdep_phaeo_in, &    
           k_exude_diatoms_in, k_exude_sp_in, k_exude_phaeo_in, &      
           K_Nit_diatoms_in, K_Nit_sp_in, K_Nit_phaeo_in, &        
           K_Am_diatoms_in, K_Am_sp_in, K_Am_phaeo_in, &          
           K_Sil_diatoms_in, K_Sil_sp_in, K_Sil_phaeo_in, &         
           K_Fe_diatoms_in, K_Fe_sp_in, K_Fe_phaeo_in, &           
           f_don_protein_in, kn_bac_protein_in, f_don_Am_protein_in, & 
           f_doc_s_in, f_doc_l_in, f_exude_s_in, f_exude_l_in, k_bac_s_in, k_bac_l_in, &             
           T_max_in, fsal_in, op_dep_min_in, fr_graze_s_in, fr_graze_e_in, fr_mort2min_in, &        
           fr_dFe_in, k_nitrif_in, t_iron_conv_in, max_loss_in, max_dfe_doc1_in, fr_resp_s_in, &          
           y_sk_DMS_in, t_sk_conv_in, t_sk_ox_in, algaltype_diatoms_in, algaltype_sp_in, &       
           algaltype_phaeo_in, nitratetype_in, ammoniumtype_in, silicatetype_in, &       
           dmspptype_in, dmspdtype_in, humtype_in, doctype_s_in, doctype_l_in, dontype_protein_in, &     
           fedtype_1_in, feptype_1_in, zaerotype_bc1_in, zaerotype_bc2_in, zaerotype_dust1_in, &     
           zaerotype_dust2_in, zaerotype_dust3_in, zaerotype_dust4_in, ratio_C2N_diatoms_in, &   
           ratio_C2N_sp_in, ratio_C2N_phaeo_in, ratio_chl2N_diatoms_in, ratio_chl2N_sp_in, &      
           ratio_chl2N_phaeo_in, F_abs_chl_diatoms_in, F_abs_chl_sp_in, F_abs_chl_phaeo_in, &
           ratio_C2N_proteins_in)
           !restore_bgc_in)

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
             calc_Tsfc_in       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(in), optional :: &
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
             mu_rdg_in, &      ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 
             Cf_in             ! ratio of ridging work to PE change in ridging (kstrength = 1)

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

        logical (kind=log_kind), intent(in), optional :: &
             oceanmixed_ice_in           ! if true, use ocean mixed layer
        
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

   !  character(char_len_long), intent(in), optional :: & 
   !     bgc_data_dir_in   ! directory for biogeochemistry data

     character(char_len), intent(in), optional :: &     
        bgc_flux_type_in    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      
    !    sil_data_type_in  , & ! 'default', 'clim'
    !    nit_data_type_in  , & ! 'default', 'clim'   
    !    fe_data_type_in   , & ! 'default', 'clim'      

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
         grid_o_t_in    , & ! top grid point length scale 
         initbio_frac_in, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav_in , & ! multiple of ocean tracer concentration due to frazil scavenging
         phi_snow_in        ! snow porosity at the ice/snow interface 

      real (kind=dbl_kind), intent(in), optional :: & 
         grid_oS_in     , & ! for bottom flux (zsalinity)
         l_skS_in           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(in), optional :: &
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

     !logical (kind=log_kind), intent(in), optional :: & 
     !   restore_bgc_in      ! if true, restore nitrate

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
        
        if (present(ktherm_in)               ) ktherm        = ktherm_in
        if (present(conduct_in)              ) conduct       = conduct_in
        if (present(fbot_xfer_type_in)       ) fbot_xfer_type    = fbot_xfer_type_in
        if (present(calc_Tsfc_in)            ) calc_Tsfc         = calc_Tsfc_in
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
        if (present(Cf_in)                   ) Cf            = Cf_in
        if (present(atmbndy_in)              ) atmbndy       = atmbndy_in
        if (present(calc_strair_in)          ) calc_strair   = calc_strair_in
        if (present(formdrag_in)             ) formdrag      = formdrag_in
        if (present(highfreq_in)             ) highfreq      = highfreq_in
        if (present(natmiter_in)             ) natmiter      = natmiter_in
        if (present(oceanmixed_ice_in)       ) oceanmixed_ice = oceanmixed_ice_in
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
     !  if (present(bgc_data_dir_in)         ) bgc_data_dir  = bgc_data_dir_in
     !  if (present(sil_data_type_in)        ) sil_data_type = sil_data_type_in
     !  if (present(nit_data_type_in)        ) nit_data_type = nit_data_type_in
     !  if (present(fe_data_type_in)         ) fe_data_type  = fe_data_type_in
        if (present(bgc_flux_type_in)        ) bgc_flux_type = bgc_flux_type_in
        if (present(z_tracers_in)            ) z_tracers     = z_tracers_in
        if (present(scale_bgc_in)            ) scale_bgc     = scale_bgc_in
        if (present(solve_zbgc_in)           ) solve_zbgc    = solve_zbgc_in
        if (present(dEdd_algae_in)           ) dEdd_algae    = dEdd_algae_in
        if (present(skl_bgc_in)              ) skl_bgc       = skl_bgc_in
        if (present(grid_o_in)               ) grid_o        = grid_o_in
        if (present(l_sk_in)                 ) l_sk          = l_sk_in
        if (present(grid_o_t_in)             ) grid_o_t      = grid_o_t_in
        if (present(initbio_frac_in)         ) initbio_frac  = initbio_frac_in
        if (present(frazil_scav_in)          ) frazil_scav   = frazil_scav_in
        if (present(grid_oS_in)              ) grid_oS       = grid_oS_in
        if (present(l_skS_in)                ) l_skS         = l_skS_in
        if (present(phi_snow_in)             ) phi_snow      = phi_snow_in
     !  if (present(restore_bgc_in)          ) restore_bgc   = restore_bgc_in
        if (present(ratio_Si2N_diatoms_in)   ) ratio_Si2N_diatoms = ratio_Si2N_diatoms_in
        if (present(ratio_Si2N_sp_in)        ) ratio_Si2N_sp      = ratio_Si2N_sp_in
        if (present(ratio_Si2N_phaeo_in)     ) ratio_Si2N_phaeo   = ratio_Si2N_phaeo_in
        if (present(ratio_S2N_diatoms_in)    ) ratio_S2N_diatoms  = ratio_S2N_diatoms_in
        if (present(ratio_S2N_sp_in)         ) ratio_S2N_sp       = ratio_S2N_sp_in
        if (present(ratio_S2N_phaeo_in)      ) ratio_S2N_phaeo    = ratio_S2N_phaeo_in
        if (present(ratio_Fe2C_diatoms_in)   ) ratio_Fe2C_diatoms = ratio_Fe2C_diatoms_in
        if (present(ratio_Fe2C_sp_in)        ) ratio_Fe2C_sp      = ratio_Fe2C_sp_in
        if (present(ratio_Fe2C_phaeo_in)     ) ratio_Fe2C_phaeo   = ratio_Fe2C_phaeo_in
        if (present(ratio_Fe2N_diatoms_in)   ) ratio_Fe2N_diatoms = ratio_Fe2N_diatoms_in
        if (present(ratio_Fe2N_sp_in)        ) ratio_Fe2N_sp      = ratio_Fe2N_sp_in
        if (present(ratio_Fe2N_phaeo_in)     ) ratio_Fe2N_phaeo   = ratio_Fe2N_phaeo_in
        if (present(ratio_Fe2DON_in)         ) ratio_Fe2DON       = ratio_Fe2DON_in
        if (present(ratio_Fe2DOC_s_in)       ) ratio_Fe2DOC_s     = ratio_Fe2DOC_s_in
        if (present(ratio_Fe2DOC_l_in)       ) ratio_Fe2DOC_l     = ratio_Fe2DOC_l_in
        if (present(fr_resp_in)              ) fr_resp          = fr_resp_in
        if (present(tau_min_in)              ) tau_min          = tau_min_in
        if (present(tau_max_in)              ) tau_max          = tau_max_in
        if (present(algal_vel_in)            ) algal_vel        = algal_vel_in
        if (present(R_dFe2dust_in)           ) R_dFe2dust       = R_dFe2dust_in
        if (present(dustFe_sol_in)           ) dustFe_sol       = dustFe_sol_in
        if (present(chlabs_diatoms_in)       ) chlabs_diatoms   = chlabs_diatoms_in
        if (present(chlabs_sp_in)            ) chlabs_sp        = chlabs_sp_in
        if (present(chlabs_phaeo_in)         ) chlabs_phaeo     = chlabs_phaeo_in
        if (present(alpha2max_low_diatoms_in)) alpha2max_low_diatoms = alpha2max_low_diatoms_in
        if (present(alpha2max_low_sp_in)     ) alpha2max_low_sp = alpha2max_low_sp_in
        if (present(alpha2max_low_phaeo_in)  ) alpha2max_low_phaeo = alpha2max_low_phaeo_in
        if (present(beta2max_diatoms_in)     ) beta2max_diatoms = beta2max_diatoms_in
        if (present(beta2max_sp_in)          ) beta2max_sp      = beta2max_sp_in
        if (present(beta2max_phaeo_in)       ) beta2max_phaeo   = beta2max_phaeo_in
        if (present(mu_max_diatoms_in)       ) mu_max_diatoms   = mu_max_diatoms_in
        if (present(mu_max_sp_in)            ) mu_max_sp        = mu_max_sp_in
        if (present(mu_max_phaeo_in)         ) mu_max_phaeo     = mu_max_phaeo_in
        if (present(grow_Tdep_diatoms_in)    ) grow_Tdep_diatoms= grow_Tdep_diatoms_in
        if (present(grow_Tdep_sp_in)         ) grow_Tdep_sp     = grow_Tdep_sp_in
        if (present(grow_Tdep_phaeo_in)      ) grow_Tdep_phaeo  = grow_Tdep_phaeo_in
        if (present(fr_graze_diatoms_in)     ) fr_graze_diatoms = fr_graze_diatoms_in
        if (present(fr_graze_sp_in)          ) fr_graze_sp      = fr_graze_sp_in
        if (present(fr_graze_phaeo_in)       ) fr_graze_phaeo   = fr_graze_phaeo_in
        if (present(mort_pre_diatoms_in)     ) mort_pre_diatoms = mort_pre_diatoms_in
        if (present(mort_pre_sp_in)          ) mort_pre_sp      = mort_pre_sp_in
        if (present(mort_pre_phaeo_in)       ) mort_pre_phaeo   = mort_pre_phaeo_in
        if (present(mort_Tdep_diatoms_in)    ) mort_Tdep_diatoms= mort_Tdep_diatoms_in
        if (present(mort_Tdep_sp_in)         ) mort_Tdep_sp     = mort_Tdep_sp_in
        if (present(mort_Tdep_phaeo_in)      ) mort_Tdep_phaeo  = mort_Tdep_phaeo_in
        if (present(k_exude_diatoms_in)      ) k_exude_diatoms  = k_exude_diatoms_in
        if (present(k_exude_sp_in)           ) k_exude_sp       = k_exude_sp_in
        if (present(k_exude_phaeo_in)        ) k_exude_phaeo    = k_exude_phaeo_in
        if (present(K_Nit_diatoms_in)        ) K_Nit_diatoms    = K_Nit_diatoms_in
        if (present(K_Nit_sp_in)             ) K_Nit_sp         = K_Nit_sp_in
        if (present(K_Nit_phaeo_in)          ) K_Nit_phaeo      = K_Nit_phaeo_in
        if (present(K_Am_diatoms_in)         ) K_Am_diatoms     = K_Am_diatoms_in
        if (present(K_Am_sp_in)              ) K_Am_sp          = K_Am_sp_in
        if (present(K_Am_phaeo_in)           ) K_Am_phaeo       = K_Am_phaeo_in
        if (present(K_Sil_diatoms_in)        ) K_Sil_diatoms    = K_Sil_diatoms_in
        if (present(K_Sil_sp_in)             ) K_Sil_sp         = K_Sil_sp_in
        if (present(K_Sil_phaeo_in)          ) K_Sil_phaeo      = K_Sil_phaeo_in
        if (present(K_Fe_diatoms_in)         ) K_Fe_diatoms     = K_Fe_diatoms_in
        if (present(K_Fe_sp_in)              ) K_Fe_sp          = K_Fe_sp_in
        if (present(K_Fe_phaeo_in)           ) K_Fe_phaeo       = K_Fe_phaeo_in
        if (present(f_don_protein_in)        ) f_don_protein    = f_don_protein_in
        if (present(kn_bac_protein_in)       ) kn_bac_protein   = kn_bac_protein_in
        if (present(f_don_Am_protein_in)     ) f_don_Am_protein = f_don_Am_protein_in
        if (present(f_doc_s_in)              ) f_doc_s          = f_doc_s_in
        if (present(f_doc_l_in)              ) f_doc_l          = f_doc_l_in
        if (present(f_exude_s_in)            ) f_exude_s        = f_exude_s_in
        if (present(f_exude_l_in)            ) f_exude_l        = f_exude_l_in
        if (present(k_bac_s_in)              ) k_bac_s          = k_bac_s_in
        if (present(k_bac_l_in)              ) k_bac_l          = k_bac_l_in
        if (present(T_max_in)                ) T_max            = T_max_in
        if (present(fsal_in)                 ) fsal             = fsal_in
        if (present(op_dep_min_in)           ) op_dep_min       = op_dep_min_in
        if (present(fr_graze_s_in)           ) fr_graze_s       = fr_graze_s_in
        if (present(fr_graze_e_in)           ) fr_graze_e       = fr_graze_e_in
        if (present(fr_mort2min_in)          ) fr_mort2min      = fr_mort2min_in
        if (present(fr_dFe_in)               ) fr_dFe           = fr_dFe_in
        if (present(k_nitrif_in)             ) k_nitrif         = k_nitrif_in
        if (present(t_iron_conv_in)          ) t_iron_conv      = t_iron_conv_in
        if (present(max_loss_in)             ) max_loss         = max_loss_in
        if (present(max_dfe_doc1_in)         ) max_dfe_doc1     = max_dfe_doc1_in
        if (present(fr_resp_s_in)            ) fr_resp_s        = fr_resp_s_in
        if (present(y_sk_DMS_in)             ) y_sk_DMS         = y_sk_DMS_in
        if (present(t_sk_conv_in)            ) t_sk_conv        = t_sk_conv_in
        if (present(t_sk_ox_in)              ) t_sk_ox          = t_sk_ox_in
        if (present(algaltype_diatoms_in)    ) algaltype_diatoms  = algaltype_diatoms_in
        if (present(algaltype_sp_in)         ) algaltype_sp       = algaltype_sp_in
        if (present(algaltype_phaeo_in)      ) algaltype_phaeo    = algaltype_phaeo_in
        if (present(nitratetype_in)          ) nitratetype        = nitratetype_in
        if (present(ammoniumtype_in)         ) ammoniumtype       = ammoniumtype_in
        if (present(silicatetype_in)         ) silicatetype       = silicatetype_in
        if (present(dmspptype_in)            ) dmspptype          = dmspptype_in
        if (present(dmspdtype_in)            ) dmspdtype          = dmspdtype_in
        if (present(humtype_in)              ) humtype            = humtype_in
        if (present(doctype_s_in)            ) doctype_s          = doctype_s_in
        if (present(doctype_l_in)            ) doctype_l          = doctype_l_in
        if (present(dontype_protein_in)      ) dontype_protein    = dontype_protein_in
        if (present(fedtype_1_in)            ) fedtype_1          = fedtype_1_in
        if (present(feptype_1_in)            ) feptype_1          = feptype_1_in
        if (present(zaerotype_bc1_in)        ) zaerotype_bc1      = zaerotype_bc1_in
        if (present(zaerotype_bc2_in)        ) zaerotype_bc2      = zaerotype_bc2_in
        if (present(zaerotype_dust1_in)      ) zaerotype_dust1    = zaerotype_dust1_in
        if (present(zaerotype_dust2_in)      ) zaerotype_dust2    = zaerotype_dust2_in
        if (present(zaerotype_dust3_in)      ) zaerotype_dust3    = zaerotype_dust3_in
        if (present(zaerotype_dust4_in)      ) zaerotype_dust4    = zaerotype_dust4_in
        if (present(ratio_C2N_diatoms_in)    ) ratio_C2N_diatoms  = ratio_C2N_diatoms_in
        if (present(ratio_C2N_sp_in)         ) ratio_C2N_sp       = ratio_C2N_sp_in
        if (present(ratio_C2N_phaeo_in)      ) ratio_C2N_phaeo    = ratio_C2N_phaeo_in
        if (present(ratio_chl2N_diatoms_in)  ) ratio_chl2N_diatoms= ratio_chl2N_diatoms_in
        if (present(ratio_chl2N_sp_in)       ) ratio_chl2N_sp     = ratio_chl2N_sp_in
        if (present(ratio_chl2N_phaeo_in)    ) ratio_chl2N_phaeo  = ratio_chl2N_phaeo_in
        if (present(F_abs_chl_diatoms_in)    ) F_abs_chl_diatoms  = F_abs_chl_diatoms_in
        if (present(F_abs_chl_sp_in)         ) F_abs_chl_sp       = F_abs_chl_sp_in
        if (present(F_abs_chl_phaeo_in)      ) F_abs_chl_phaeo    = F_abs_chl_phaeo_in
        if (present(ratio_C2N_proteins_in)   ) ratio_C2N_proteins = ratio_C2N_proteins_in

      end subroutine icepack_init_parameters

!=======================================================================

! subroutine to query the column package internal parameters

      subroutine icepack_query_parameters(   &
           max_algae_out, max_dic_out, max_doc_out, max_don_out, max_fe_out, &
           nmodal1_out, nmodal2_out, max_aero_out, max_nbtrcr_out, &
           ktherm_out, conduct_out, fbot_xfer_type_out, calc_Tsfc_out, ustar_min_out, a_rapid_mode_out, &
           Rac_rapid_mode_out, aspect_rapid_mode_out, dSdt_slow_mode_out, phi_c_slow_mode_out, &
           phi_i_mushy_out, shortwave_out, albedo_type_out, albicev_out, albicei_out, albsnowv_out, &
           albsnowi_out, ahmax_out, R_ice_out, R_pnd_out, R_snw_out, dT_mlt_out, rsnw_mlt_out, &
           kalg_out, kstrength_out, krdg_partic_out, krdg_redist_out, mu_rdg_out, Cf_out, &
           atmbndy_out, calc_strair_out, formdrag_out, highfreq_out, natmiter_out, &
           oceanmixed_ice_out, tfrz_option_out, kitd_out, kcatbound_out, hs0_out, frzpnd_out, &
           dpscale_out, rfracmin_out, rfracmax_out, pndaspect_out, hs1_out, hp1_out, &
         ! bgc_data_dir_out, sil_data_type_out, nit_data_type_out, fe_data_type_out, &
           bgc_flux_type_out, z_tracers_out, scale_bgc_out, solve_zbgc_out, dEdd_algae_out, &
           modal_aero_out, skl_bgc_out, solve_zsal_out, grid_o_out, l_sk_out, &
           grid_o_t_out, initbio_frac_out, frazil_scav_out, grid_oS_out, l_skS_out, &
           phi_snow_out, ratio_Si2N_diatoms_out, ratio_Si2N_sp_out, ratio_Si2N_phaeo_out, &
           ratio_S2N_diatoms_out, ratio_S2N_sp_out, ratio_S2N_phaeo_out, ratio_Fe2C_diatoms_out, & 
           ratio_Fe2C_sp_out, ratio_Fe2C_phaeo_out, ratio_Fe2N_diatoms_out, ratio_Fe2N_sp_out, &     
           ratio_Fe2N_phaeo_out, ratio_Fe2DON_out, ratio_Fe2DOC_s_out, ratio_Fe2DOC_l_out, &     
           fr_resp_out, tau_min_out, tau_max_out, algal_vel_out, R_dFe2dust_out, dustFe_sol_out, &         
           chlabs_diatoms_out, chlabs_sp_out, chlabs_phaeo_out, alpha2max_low_diatoms_out, &
           alpha2max_low_sp_out, alpha2max_low_phaeo_out, beta2max_diatoms_out, beta2max_sp_out, &       
           beta2max_phaeo_out, mu_max_diatoms_out, mu_max_sp_out, mu_max_phaeo_out, &      
           grow_Tdep_diatoms_out, grow_Tdep_sp_out, grow_Tdep_phaeo_out, &   
           fr_graze_diatoms_out, fr_graze_sp_out, fr_graze_phaeo_out, &    
           mort_pre_diatoms_out, mort_pre_sp_out, mort_pre_phaeo_out, &    
           mort_Tdep_diatoms_out, mort_Tdep_sp_out, mort_Tdep_phaeo_out, &    
           k_exude_diatoms_out, k_exude_sp_out, k_exude_phaeo_out, &      
           K_Nit_diatoms_out, K_Nit_sp_out, K_Nit_phaeo_out, &        
           K_Am_diatoms_out, K_Am_sp_out, K_Am_phaeo_out, &          
           K_Sil_diatoms_out, K_Sil_sp_out, K_Sil_phaeo_out, &         
           K_Fe_diatoms_out, K_Fe_sp_out, K_Fe_phaeo_out, &           
           f_don_protein_out, kn_bac_protein_out, f_don_Am_protein_out, & 
           f_doc_s_out, f_doc_l_out, f_exude_s_out, f_exude_l_out, k_bac_s_out, k_bac_l_out, &             
           T_max_out, fsal_out, op_dep_min_out, fr_graze_s_out, fr_graze_e_out, fr_mort2min_out, &        
           fr_dFe_out, k_nitrif_out, t_iron_conv_out, max_loss_out, max_dfe_doc1_out, fr_resp_s_out, &          
           y_sk_DMS_out, t_sk_conv_out, t_sk_ox_out, algaltype_diatoms_out, algaltype_sp_out, &       
           algaltype_phaeo_out, nitratetype_out, ammoniumtype_out, silicatetype_out, &       
           dmspptype_out, dmspdtype_out, humtype_out, doctype_s_out, doctype_l_out, dontype_protein_out, &     
           fedtype_1_out, feptype_1_out, zaerotype_bc1_out, zaerotype_bc2_out, zaerotype_dust1_out, &     
           zaerotype_dust2_out, zaerotype_dust3_out, zaerotype_dust4_out, ratio_C2N_diatoms_out, &   
           ratio_C2N_sp_out, ratio_C2N_phaeo_out, ratio_chl2N_diatoms_out, ratio_chl2N_sp_out, &      
           ratio_chl2N_phaeo_out, F_abs_chl_diatoms_out, F_abs_chl_sp_out, F_abs_chl_phaeo_out, &
           ratio_C2N_proteins_out)
           !restore_bgc_out)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(out), optional :: &
             max_algae_out, & ! maximum number of algal types
             max_dic_out,   & ! maximum number of dissolved inorganic carbon types
             max_doc_out,   & ! maximum number of dissolved organic carbon types
             max_don_out,   & ! maximum number of dissolved organic nitrogen types
             max_fe_out,    & ! maximum number of iron types
             nmodal1_out,   & ! dimension for modal aerosol radiation parameters
             nmodal2_out,   & ! dimension for modal aerosol radiation parameters
             max_aero_out,  & ! maximum number of aerosols
             max_nbtrcr_out   !

        integer (kind=int_kind), intent(out), optional :: &
             ktherm_out         ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(out), optional :: &
             conduct_out, &      ! 'MU71' or 'bubbly'
             fbot_xfer_type_out  ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(out), optional :: &
             calc_Tsfc_out       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(out), optional :: &
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
             mu_rdg_out, &      ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 
             Cf_out             ! ratio of ridging work to PE change in ridging (kstrength = 1)

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

        logical (kind=log_kind), intent(out), optional :: &
             oceanmixed_ice_out           ! if true, use ocean mixed layer
        
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

   !  character(char_len_long), intent(out), optional :: & 
   !     bgc_data_dir_out   ! directory for biogeochemistry data

     character(char_len), intent(out), optional :: &     
        bgc_flux_type_out    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      
    !    sil_data_type_out  , & ! 'default', 'clim'
    !    nit_data_type_out  , & ! 'default', 'clim'   
    !    fe_data_type_out   , & ! 'default', 'clim'      

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
         grid_o_t_out    , & ! top grid point length scale 
         initbio_frac_out, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav_out , & ! multiple of ocean tracer concentration due to frazil scavenging
         phi_snow_out        ! snow porosity at the ice/snow interface 

      real (kind=dbl_kind), intent(out), optional :: & 
         grid_oS_out     , & ! for bottom flux (zsalinity)
         l_skS_out           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(out), optional :: &
         ratio_Si2N_diatoms_out, &   ! algal Si to N (mol/mol)
         ratio_Si2N_sp_out     , &
         ratio_Si2N_phaeo_out  , &
         ratio_S2N_diatoms_out , &   ! algal S  to N (mol/mol)
         ratio_S2N_sp_out      , &
         ratio_S2N_phaeo_out   , &
         ratio_Fe2C_diatoms_out, &   ! algal Fe to C  (umol/mol)
         ratio_Fe2C_sp_out     , &
         ratio_Fe2C_phaeo_out  , &
         ratio_Fe2N_diatoms_out, &   ! algal Fe to N  (umol/mol)
         ratio_Fe2N_sp_out     , &
         ratio_Fe2N_phaeo_out  , &
         ratio_Fe2DON_out      , &   ! Fe to N of DON (nmol/umol)
         ratio_Fe2DOC_s_out    , &   ! Fe to C of DOC (nmol/umol) saccharids
         ratio_Fe2DOC_l_out    , &   ! Fe to C of DOC (nmol/umol) lipids
         fr_resp_out           , &   ! fraction of algal growth lost due to respiration
         tau_min_out           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
         tau_max_out           , &   ! long time mobile to stationary exchanges (s) = 2 days
         algal_vel_out         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_out        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_out        , &   ! solubility fraction
         chlabs_diatoms_out   , & ! chl absorption (1/m/(mg/m^3))
         chlabs_sp_out        , & !
         chlabs_phaeo_out     , & !
         alpha2max_low_diatoms_out , & ! light limitation (1/(W/m^2))  
         alpha2max_low_sp_out      , & 
         alpha2max_low_phaeo_out   , & 
         beta2max_diatoms_out , & ! light inhibition (1/(W/m^2))  
         beta2max_sp_out      , & 
         beta2max_phaeo_out   , & 
         mu_max_diatoms_out   , & ! maximum growth rate (1/day)       
         mu_max_sp_out        , & 
         mu_max_phaeo_out     , & 
         grow_Tdep_diatoms_out, & ! Temperature dependence of growth (1/C)
         grow_Tdep_sp_out     , & 
         grow_Tdep_phaeo_out  , & 
         fr_graze_diatoms_out , & ! Fraction grazed
         fr_graze_sp_out      , & 
         fr_graze_phaeo_out   , & 
         mort_pre_diatoms_out , & ! Mortality (1/day)
         mort_pre_sp_out      , & 
         mort_pre_phaeo_out   , & 
         mort_Tdep_diatoms_out, & ! T dependence of mortality (1/C) 
         mort_Tdep_sp_out     , &  
         mort_Tdep_phaeo_out  , &  
         k_exude_diatoms_out  , & ! algal exudation (1/d)
         k_exude_sp_out       , &  
         k_exude_phaeo_out    , &  
         K_Nit_diatoms_out    , & ! nitrate half saturation (mmol/m^3)
         K_Nit_sp_out         , &  
         K_Nit_phaeo_out      , &  
         K_Am_diatoms_out     , & ! ammonium half saturation (mmol/m^3)
         K_Am_sp_out          , &   
         K_Am_phaeo_out       , &   
         K_Sil_diatoms_out    , & ! silicate half saturation (mmol/m^3)
         K_Sil_sp_out         , &   
         K_Sil_phaeo_out      , &   
         K_Fe_diatoms_out     , & ! iron half saturation (nM)
         K_Fe_sp_out          , &   
         K_Fe_phaeo_out       , &    
         f_don_protein_out    , & ! fraction of spilled grazing to proteins            
         kn_bac_protein_out   , & ! Bacterial degredation of DON (1/d)                  
         f_don_Am_protein_out , & ! fraction of remineralized DON to ammonium          
         f_doc_s_out          , & ! fraction of mortality to DOC 
         f_doc_l_out          , &   
         f_exude_s_out        , & ! fraction of exudation to DOC
         f_exude_l_out        , & 
         k_bac_s_out          , & ! Bacterial degredation of DOC (1/d)
         k_bac_l_out          , & 
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
         t_sk_ox_out          , &   ! DMS oxidation time (d)
         algaltype_diatoms_out  , & ! mobility type
         algaltype_sp_out       , & !
         algaltype_phaeo_out    , & !
         nitratetype_out        , & !
         ammoniumtype_out       , & !
         silicatetype_out       , & !
         dmspptype_out          , & !
         dmspdtype_out          , & !
         humtype_out            , & !
         doctype_s_out          , & !
         doctype_l_out          , & !
         dontype_protein_out    , & !
         fedtype_1_out          , & !
         feptype_1_out          , & !
         zaerotype_bc1_out      , & !
         zaerotype_bc2_out      , & !
         zaerotype_dust1_out    , & !
         zaerotype_dust2_out    , & !
         zaerotype_dust3_out    , & !
         zaerotype_dust4_out    , & !
         ratio_C2N_diatoms_out  , & ! algal C to N ratio (mol/mol)
         ratio_C2N_sp_out       , & !
         ratio_C2N_phaeo_out    , & !
         ratio_chl2N_diatoms_out, & ! algal chlorophyll to N ratio (mg/mmol)
         ratio_chl2N_sp_out     , & !
         ratio_chl2N_phaeo_out  , & !
         F_abs_chl_diatoms_out  , & ! scales absorbed radiation for dEdd
         F_abs_chl_sp_out       , & !
         F_abs_chl_phaeo_out    , & !
         ratio_C2N_proteins_out     ! ratio of C to N in proteins (mol/mol)       

     !logical (kind=log_kind), intent(out), optional :: & 
     !   restore_bgc_out      ! if true, restore nitrate

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
        
        if (present(max_algae_out) ) max_algae_out  = max_algae
        if (present(max_dic_out)   ) max_dic_out    = max_dic
        if (present(max_doc_out)   ) max_doc_out    = max_doc
        if (present(max_don_out)   ) max_don_out    = max_don
        if (present(max_fe_out)    ) max_fe_out     = max_fe
        if (present(nmodal1_out)   ) nmodal1_out    = nmodal1
        if (present(nmodal2_out)   ) nmodal2_out    = nmodal2
        if (present(max_aero_out)  ) max_aero_out   = max_aero
        if (present(max_nbtrcr_out)) max_nbtrcr_out = max_nbtrcr

        if (present(ktherm_out)               ) ktherm_out        = ktherm
        if (present(conduct_out)              ) conduct_out       = conduct
        if (present(fbot_xfer_type_out)       ) fbot_xfer_type_out    = fbot_xfer_type
        if (present(calc_Tsfc_out)            ) calc_Tsfc_out         = calc_Tsfc
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
        if (present(Cf_out)                   ) Cf_out            = Cf
        if (present(atmbndy_out)              ) atmbndy_out       = atmbndy
        if (present(calc_strair_out)          ) calc_strair_out   = calc_strair
        if (present(formdrag_out)             ) formdrag_out      = formdrag
        if (present(highfreq_out)             ) highfreq_out      = highfreq
        if (present(natmiter_out)             ) natmiter_out      = natmiter
        if (present(oceanmixed_ice_out)       ) oceanmixed_ice_out = oceanmixed_ice
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
     !  if (present(bgc_data_dir_out)         ) bgc_data_dir_out  = bgc_data_dir
     !  if (present(sil_data_type_out)        ) sil_data_type_out = sil_data_type
     !  if (present(nit_data_type_out)        ) nit_data_type_out = nit_data_type
     !  if (present(fe_data_type_out)         ) fe_data_type_out  = fe_data_type
        if (present(bgc_flux_type_out)        ) bgc_flux_type_out = bgc_flux_type
        if (present(z_tracers_out)            ) z_tracers_out     = z_tracers
        if (present(scale_bgc_out)            ) scale_bgc_out     = scale_bgc
        if (present(solve_zbgc_out)           ) solve_zbgc_out    = solve_zbgc
        if (present(dEdd_algae_out)           ) dEdd_algae_out    = dEdd_algae
        if (present(skl_bgc_out)              ) skl_bgc_out       = skl_bgc
        if (present(grid_o_out)               ) grid_o_out        = grid_o
        if (present(l_sk_out)                 ) l_sk_out          = l_sk
        if (present(grid_o_t_out)             ) grid_o_t_out      = grid_o_t
        if (present(initbio_frac_out)         ) initbio_frac_out  = initbio_frac
        if (present(frazil_scav_out)          ) frazil_scav_out   = frazil_scav
        if (present(grid_oS_out)              ) grid_oS_out       = grid_oS
        if (present(l_skS_out)                ) l_skS_out         = l_skS
        if (present(phi_snow_out)             ) phi_snow_out      = phi_snow
     !  if (present(restore_bgc_out)          ) restore_bgc_out   = restore_bgc
        if (present(ratio_Si2N_diatoms_out)   ) ratio_Si2N_diatoms_out = ratio_Si2N_diatoms
        if (present(ratio_Si2N_sp_out)        ) ratio_Si2N_sp_out      = ratio_Si2N_sp
        if (present(ratio_Si2N_phaeo_out)     ) ratio_Si2N_phaeo_out   = ratio_Si2N_phaeo
        if (present(ratio_S2N_diatoms_out)    ) ratio_S2N_diatoms_out  = ratio_S2N_diatoms
        if (present(ratio_S2N_sp_out)         ) ratio_S2N_sp_out       = ratio_S2N_sp
        if (present(ratio_S2N_phaeo_out)      ) ratio_S2N_phaeo_out    = ratio_S2N_phaeo
        if (present(ratio_Fe2C_diatoms_out)   ) ratio_Fe2C_diatoms_out = ratio_Fe2C_diatoms
        if (present(ratio_Fe2C_sp_out)        ) ratio_Fe2C_sp_out      = ratio_Fe2C_sp
        if (present(ratio_Fe2C_phaeo_out)     ) ratio_Fe2C_phaeo_out   = ratio_Fe2C_phaeo
        if (present(ratio_Fe2N_diatoms_out)   ) ratio_Fe2N_diatoms_out = ratio_Fe2N_diatoms
        if (present(ratio_Fe2N_sp_out)        ) ratio_Fe2N_sp_out      = ratio_Fe2N_sp
        if (present(ratio_Fe2N_phaeo_out)     ) ratio_Fe2N_phaeo_out   = ratio_Fe2N_phaeo
        if (present(ratio_Fe2DON_out)         ) ratio_Fe2DON_out       = ratio_Fe2DON
        if (present(ratio_Fe2DOC_s_out)       ) ratio_Fe2DOC_s_out     = ratio_Fe2DOC_s
        if (present(ratio_Fe2DOC_l_out)       ) ratio_Fe2DOC_l_out     = ratio_Fe2DOC_l
        if (present(fr_resp_out)              ) fr_resp_out          = fr_resp
        if (present(tau_min_out)              ) tau_min_out          = tau_min
        if (present(tau_max_out)              ) tau_max_out          = tau_max
        if (present(algal_vel_out)            ) algal_vel_out        = algal_vel
        if (present(R_dFe2dust_out)           ) R_dFe2dust_out       = R_dFe2dust
        if (present(dustFe_sol_out)           ) dustFe_sol_out       = dustFe_sol
        if (present(chlabs_diatoms_out)       ) chlabs_diatoms_out   = chlabs_diatoms
        if (present(chlabs_sp_out)            ) chlabs_sp_out        = chlabs_sp
        if (present(chlabs_phaeo_out)         ) chlabs_phaeo_out     = chlabs_phaeo
        if (present(alpha2max_low_diatoms_out)) alpha2max_low_diatoms_out = alpha2max_low_diatoms
        if (present(alpha2max_low_sp_out)     ) alpha2max_low_sp_out = alpha2max_low_sp
        if (present(alpha2max_low_phaeo_out)  ) alpha2max_low_phaeo_out = alpha2max_low_phaeo
        if (present(beta2max_diatoms_out)     ) beta2max_diatoms_out = beta2max_diatoms
        if (present(beta2max_sp_out)          ) beta2max_sp_out      = beta2max_sp
        if (present(beta2max_phaeo_out)       ) beta2max_phaeo_out   = beta2max_phaeo
        if (present(mu_max_diatoms_out)       ) mu_max_diatoms_out   = mu_max_diatoms
        if (present(mu_max_sp_out)            ) mu_max_sp_out        = mu_max_sp
        if (present(mu_max_phaeo_out)         ) mu_max_phaeo_out     = mu_max_phaeo
        if (present(grow_Tdep_diatoms_out)    ) grow_Tdep_diatoms_out = grow_Tdep_diatoms
        if (present(grow_Tdep_sp_out)         ) grow_Tdep_sp_out     = grow_Tdep_sp
        if (present(grow_Tdep_phaeo_out)      ) grow_Tdep_phaeo_out  = grow_Tdep_phaeo
        if (present(fr_graze_diatoms_out)     ) fr_graze_diatoms_out = fr_graze_diatoms
        if (present(fr_graze_sp_out)          ) fr_graze_sp_out      = fr_graze_sp
        if (present(fr_graze_phaeo_out)       ) fr_graze_phaeo_out   = fr_graze_phaeo
        if (present(mort_pre_diatoms_out)     ) mort_pre_diatoms_out = mort_pre_diatoms
        if (present(mort_pre_sp_out)          ) mort_pre_sp_out      = mort_pre_sp
        if (present(mort_pre_phaeo_out)       ) mort_pre_phaeo_out   = mort_pre_phaeo
        if (present(mort_Tdep_diatoms_out)    ) mort_Tdep_diatoms_out = mort_Tdep_diatoms
        if (present(mort_Tdep_sp_out)         ) mort_Tdep_sp_out     = mort_Tdep_sp
        if (present(mort_Tdep_phaeo_out)      ) mort_Tdep_phaeo_out  = mort_Tdep_phaeo
        if (present(k_exude_diatoms_out)      ) k_exude_diatoms_out  = k_exude_diatoms
        if (present(k_exude_sp_out)           ) k_exude_sp_out       = k_exude_sp
        if (present(k_exude_phaeo_out)        ) k_exude_phaeo_out    = k_exude_phaeo
        if (present(K_Nit_diatoms_out)        ) K_Nit_diatoms_out    = K_Nit_diatoms
        if (present(K_Nit_sp_out)             ) K_Nit_sp_out         = K_Nit_sp
        if (present(K_Nit_phaeo_out)          ) K_Nit_phaeo_out      = K_Nit_phaeo
        if (present(K_Am_diatoms_out)         ) K_Am_diatoms_out     = K_Am_diatoms
        if (present(K_Am_sp_out)              ) K_Am_sp_out          = K_Am_sp
        if (present(K_Am_phaeo_out)           ) K_Am_phaeo_out       = K_Am_phaeo
        if (present(K_Sil_diatoms_out)        ) K_Sil_diatoms_out    = K_Sil_diatoms
        if (present(K_Sil_sp_out)             ) K_Sil_sp_out         = K_Sil_sp
        if (present(K_Sil_phaeo_out)          ) K_Sil_phaeo_out      = K_Sil_phaeo
        if (present(K_Fe_diatoms_out)         ) K_Fe_diatoms_out     = K_Fe_diatoms
        if (present(K_Fe_sp_out)              ) K_Fe_sp_out          = K_Fe_sp
        if (present(K_Fe_phaeo_out)           ) K_Fe_phaeo_out       = K_Fe_phaeo
        if (present(f_don_protein_out)        ) f_don_protein_out    = f_don_protein
        if (present(kn_bac_protein_out)       ) kn_bac_protein_out   = kn_bac_protein
        if (present(f_don_Am_protein_out)     ) f_don_Am_protein_out = f_don_Am_protein
        if (present(f_doc_s_out)              ) f_doc_s_out          = f_doc_s
        if (present(f_doc_l_out)              ) f_doc_l_out          = f_doc_l
        if (present(f_exude_s_out)            ) f_exude_s_out        = f_exude_s
        if (present(f_exude_l_out)            ) f_exude_l_out        = f_exude_l
        if (present(k_bac_s_out)              ) k_bac_s_out          = k_bac_s
        if (present(k_bac_l_out)              ) k_bac_l_out          = k_bac_l
        if (present(T_max_out)                ) T_max_out            = T_max
        if (present(fsal_out)                 ) fsal_out             = fsal
        if (present(op_dep_min_out)           ) op_dep_min_out       = op_dep_min
        if (present(fr_graze_s_out)           ) fr_graze_s_out       = fr_graze_s
        if (present(fr_graze_e_out)           ) fr_graze_e_out       = fr_graze_e
        if (present(fr_mort2min_out)          ) fr_mort2min_out      = fr_mort2min
        if (present(fr_dFe_out)               ) fr_dFe_out           = fr_dFe
        if (present(k_nitrif_out)             ) k_nitrif_out         = k_nitrif
        if (present(t_iron_conv_out)          ) t_iron_conv_out      = t_iron_conv
        if (present(max_loss_out)             ) max_loss_out         = max_loss
        if (present(max_dfe_doc1_out)         ) max_dfe_doc1_out     = max_dfe_doc1
        if (present(fr_resp_s_out)            ) fr_resp_s_out        = fr_resp_s
        if (present(y_sk_DMS_out)             ) y_sk_DMS_out         = y_sk_DMS
        if (present(t_sk_conv_out)            ) t_sk_conv_out        = t_sk_conv
        if (present(t_sk_ox_out)              ) t_sk_ox_out          = t_sk_ox
        if (present(algaltype_diatoms_out)    ) algaltype_diatoms_out  = algaltype_diatoms
        if (present(algaltype_sp_out)         ) algaltype_sp_out       = algaltype_sp
        if (present(algaltype_phaeo_out)      ) algaltype_phaeo_out    = algaltype_phaeo
        if (present(nitratetype_out)          ) nitratetype_out        = nitratetype
        if (present(ammoniumtype_out)         ) ammoniumtype_out       = ammoniumtype
        if (present(silicatetype_out)         ) silicatetype_out       = silicatetype
        if (present(dmspptype_out)            ) dmspptype_out          = dmspptype
        if (present(dmspdtype_out)            ) dmspdtype_out          = dmspdtype
        if (present(humtype_out)              ) humtype_out            = humtype
        if (present(doctype_s_out)            ) doctype_s_out          = doctype_s
        if (present(doctype_l_out)            ) doctype_l_out          = doctype_l
        if (present(dontype_protein_out)      ) dontype_protein_out    = dontype_protein
        if (present(fedtype_1_out)            ) fedtype_1_out          = fedtype_1
        if (present(feptype_1_out)            ) feptype_1_out          = feptype_1
        if (present(zaerotype_bc1_out)        ) zaerotype_bc1_out      = zaerotype_bc1
        if (present(zaerotype_bc2_out)        ) zaerotype_bc2_out      = zaerotype_bc2
        if (present(zaerotype_dust1_out)      ) zaerotype_dust1_out    = zaerotype_dust1
        if (present(zaerotype_dust2_out)      ) zaerotype_dust2_out    = zaerotype_dust2
        if (present(zaerotype_dust3_out)      ) zaerotype_dust3_out    = zaerotype_dust3
        if (present(zaerotype_dust4_out)      ) zaerotype_dust4_out    = zaerotype_dust4
        if (present(ratio_C2N_diatoms_out)    ) ratio_C2N_diatoms_out  = ratio_C2N_diatoms
        if (present(ratio_C2N_sp_out)         ) ratio_C2N_sp_out       = ratio_C2N_sp
        if (present(ratio_C2N_phaeo_out)      ) ratio_C2N_phaeo_out    = ratio_C2N_phaeo
        if (present(ratio_chl2N_diatoms_out)  ) ratio_chl2N_diatoms_out = ratio_chl2N_diatoms
        if (present(ratio_chl2N_sp_out)       ) ratio_chl2N_sp_out     = ratio_chl2N_sp
        if (present(ratio_chl2N_phaeo_out)    ) ratio_chl2N_phaeo_out  = ratio_chl2N_phaeo
        if (present(F_abs_chl_diatoms_out)    ) F_abs_chl_diatoms_out  = F_abs_chl_diatoms
        if (present(F_abs_chl_sp_out)         ) F_abs_chl_sp_out       = F_abs_chl_sp
        if (present(F_abs_chl_phaeo_out)      ) F_abs_chl_phaeo_out    = F_abs_chl_phaeo
        if (present(ratio_C2N_proteins_out)   ) ratio_C2N_proteins_out = ratio_C2N_proteins

      end subroutine icepack_query_parameters

!=======================================================================

! subroutine to write the column package internal parameters

      subroutine icepack_write_parameters(iounit)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             iounit   ! unit number for output

        write(iounit,*) "icepack_write_parameters:"
        write(iounit,*) "  max_algae  = ",max_algae
        write(iounit,*) "  max_dic    = ",max_dic
        write(iounit,*) "  max_doc    = ",max_doc
        write(iounit,*) "  max_don    = ",max_don
        write(iounit,*) "  max_fe     = ",max_fe
        write(iounit,*) "  nmodal1    = ",nmodal1
        write(iounit,*) "  nmodal2    = ",nmodal2
        write(iounit,*) "  max_aero   = ",max_aero
        write(iounit,*) "  max_nbtrcr = ",max_nbtrcr
        write(iounit,*) "  ktherm        = ", ktherm
        write(iounit,*) "  conduct       = ", conduct
        write(iounit,*) "  fbot_xfer_type    = ", fbot_xfer_type
        write(iounit,*) "  calc_Tsfc         = ", calc_Tsfc
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
        write(iounit,*) "  Cf            = ", Cf
        write(iounit,*) "  atmbndy       = ", atmbndy
        write(iounit,*) "  calc_strair   = ", calc_strair
        write(iounit,*) "  formdrag      = ", formdrag
        write(iounit,*) "  highfreq      = ", highfreq
        write(iounit,*) "  natmiter      = ", natmiter
        write(iounit,*) "  oceanmixed_ice = ", oceanmixed_ice
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
     !  write(iounit,*) "  bgc_data_dir  = ", bgc_data_dir
     !  write(iounit,*) "  sil_data_type = ", sil_data_type
     !  write(iounit,*) "  nit_data_type = ", nit_data_type
     !  write(iounit,*) "  fe_data_type  = ", fe_data_type
        write(iounit,*) "  bgc_flux_type = ", bgc_flux_type
        write(iounit,*) "  z_tracers     = ", z_tracers
        write(iounit,*) "  scale_bgc     = ", scale_bgc
        write(iounit,*) "  solve_zbgc    = ", solve_zbgc
        write(iounit,*) "  dEdd_algae    = ", dEdd_algae
        write(iounit,*) "  skl_bgc       = ", skl_bgc
        write(iounit,*) "  grid_o        = ", grid_o
        write(iounit,*) "  l_sk          = ", l_sk
        write(iounit,*) "  grid_o_t      = ", grid_o_t
        write(iounit,*) "  initbio_frac  = ", initbio_frac
        write(iounit,*) "  frazil_scav   = ", frazil_scav
        write(iounit,*) "  grid_oS       = ", grid_oS
        write(iounit,*) "  l_skS         = ", l_skS
        write(iounit,*) "  phi_snow      = ", phi_snow
     !  write(iounit,*) "  restore_bgc   = ", restore_bgc
        write(iounit,*) "  ratio_Si2N_diatoms = ", ratio_Si2N_diatoms
        write(iounit,*) "  ratio_Si2N_sp      = ", ratio_Si2N_sp
        write(iounit,*) "  ratio_Si2N_phaeo   = ", ratio_Si2N_phaeo
        write(iounit,*) "  ratio_S2N_diatoms  = ", ratio_S2N_diatoms
        write(iounit,*) "  ratio_S2N_sp       = ", ratio_S2N_sp
        write(iounit,*) "  ratio_S2N_phaeo    = ", ratio_S2N_phaeo
        write(iounit,*) "  ratio_Fe2C_diatoms = ", ratio_Fe2C_diatoms
        write(iounit,*) "  ratio_Fe2C_sp      = ", ratio_Fe2C_sp
        write(iounit,*) "  ratio_Fe2C_phaeo   = ", ratio_Fe2C_phaeo
        write(iounit,*) "  ratio_Fe2N_diatoms = ", ratio_Fe2N_diatoms
        write(iounit,*) "  ratio_Fe2N_sp      = ", ratio_Fe2N_sp
        write(iounit,*) "  ratio_Fe2N_phaeo   = ", ratio_Fe2N_phaeo
        write(iounit,*) "  ratio_Fe2DON       = ", ratio_Fe2DON
        write(iounit,*) "  ratio_Fe2DOC_s     = ", ratio_Fe2DOC_s
        write(iounit,*) "  ratio_Fe2DOC_l     = ", ratio_Fe2DOC_l
        write(iounit,*) "  fr_resp          = ", fr_resp
        write(iounit,*) "  tau_min          = ", tau_min
        write(iounit,*) "  tau_max          = ", tau_max
        write(iounit,*) "  algal_vel        = ", algal_vel
        write(iounit,*) "  R_dFe2dust       = ", R_dFe2dust
        write(iounit,*) "  dustFe_sol       = ", dustFe_sol
        write(iounit,*) "  chlabs_diatoms   = ", chlabs_diatoms
        write(iounit,*) "  chlabs_sp        = ", chlabs_sp
        write(iounit,*) "  chlabs_phaeo     = ", chlabs_phaeo
        write(iounit,*) "  alpha2max_low_diatoms = ", alpha2max_low_diatoms
        write(iounit,*) "  alpha2max_low_sp = ", alpha2max_low_sp
        write(iounit,*) "  alpha2max_low_phaeo = ", alpha2max_low_phaeo
        write(iounit,*) "  beta2max_diatoms = ", beta2max_diatoms
        write(iounit,*) "  beta2max_sp      = ", beta2max_sp
        write(iounit,*) "  beta2max_phaeo   = ", beta2max_phaeo
        write(iounit,*) "  mu_max_diatoms   = ", mu_max_diatoms
        write(iounit,*) "  mu_max_sp        = ", mu_max_sp
        write(iounit,*) "  mu_max_phaeo     = ", mu_max_phaeo
        write(iounit,*) "  grow_Tdep_diatoms= ", grow_Tdep_diatoms
        write(iounit,*) "  grow_Tdep_sp     = ", grow_Tdep_sp
        write(iounit,*) "  grow_Tdep_phaeo  = ", grow_Tdep_phaeo
        write(iounit,*) "  fr_graze_diatoms = ", fr_graze_diatoms
        write(iounit,*) "  fr_graze_sp      = ", fr_graze_sp
        write(iounit,*) "  fr_graze_phaeo   = ", fr_graze_phaeo
        write(iounit,*) "  mort_pre_diatoms = ", mort_pre_diatoms
        write(iounit,*) "  mort_pre_sp      = ", mort_pre_sp
        write(iounit,*) "  mort_pre_phaeo   = ", mort_pre_phaeo
        write(iounit,*) "  mort_Tdep_diatoms= ", mort_Tdep_diatoms
        write(iounit,*) "  mort_Tdep_sp     = ", mort_Tdep_sp
        write(iounit,*) "  mort_Tdep_phaeo  = ", mort_Tdep_phaeo
        write(iounit,*) "  k_exude_diatoms  = ", k_exude_diatoms
        write(iounit,*) "  k_exude_sp       = ", k_exude_sp
        write(iounit,*) "  k_exude_phaeo    = ", k_exude_phaeo
        write(iounit,*) "  K_Nit_diatoms    = ", K_Nit_diatoms
        write(iounit,*) "  K_Nit_sp         = ", K_Nit_sp
        write(iounit,*) "  K_Nit_phaeo      = ", K_Nit_phaeo
        write(iounit,*) "  K_Am_diatoms     = ", K_Am_diatoms
        write(iounit,*) "  K_Am_sp          = ", K_Am_sp
        write(iounit,*) "  K_Am_phaeo       = ", K_Am_phaeo
        write(iounit,*) "  K_Sil_diatoms    = ", K_Sil_diatoms
        write(iounit,*) "  K_Sil_sp         = ", K_Sil_sp
        write(iounit,*) "  K_Sil_phaeo      = ", K_Sil_phaeo
        write(iounit,*) "  K_Fe_diatoms     = ", K_Fe_diatoms
        write(iounit,*) "  K_Fe_sp          = ", K_Fe_sp
        write(iounit,*) "  K_Fe_phaeo       = ", K_Fe_phaeo
        write(iounit,*) "  f_don_protein    = ", f_don_protein
        write(iounit,*) "  kn_bac_protein   = ", kn_bac_protein
        write(iounit,*) "  f_don_Am_protein = ", f_don_Am_protein
        write(iounit,*) "  f_doc_s          = ", f_doc_s
        write(iounit,*) "  f_doc_l          = ", f_doc_l
        write(iounit,*) "  f_exude_s        = ", f_exude_s
        write(iounit,*) "  f_exude_l        = ", f_exude_l
        write(iounit,*) "  k_bac_s          = ", k_bac_s
        write(iounit,*) "  k_bac_l          = ", k_bac_l
        write(iounit,*) "  T_max            = ", T_max
        write(iounit,*) "  fsal             = ", fsal
        write(iounit,*) "  op_dep_min       = ", op_dep_min
        write(iounit,*) "  fr_graze_s       = ", fr_graze_s
        write(iounit,*) "  fr_graze_e       = ", fr_graze_e
        write(iounit,*) "  fr_mort2min      = ", fr_mort2min
        write(iounit,*) "  fr_dFe           = ", fr_dFe
        write(iounit,*) "  k_nitrif         = ", k_nitrif
        write(iounit,*) "  t_iron_conv      = ", t_iron_conv
        write(iounit,*) "  max_loss         = ", max_loss
        write(iounit,*) "  max_dfe_doc1     = ", max_dfe_doc1
        write(iounit,*) "  fr_resp_s        = ", fr_resp_s
        write(iounit,*) "  y_sk_DMS         = ", y_sk_DMS
        write(iounit,*) "  t_sk_conv        = ", t_sk_conv
        write(iounit,*) "  t_sk_ox          = ", t_sk_ox
        write(iounit,*) "  algaltype_diatoms  = ", algaltype_diatoms
        write(iounit,*) "  algaltype_sp       = ", algaltype_sp
        write(iounit,*) "  algaltype_phaeo    = ", algaltype_phaeo
        write(iounit,*) "  nitratetype        = ", nitratetype
        write(iounit,*) "  ammoniumtype       = ", ammoniumtype
        write(iounit,*) "  silicatetype       = ", silicatetype
        write(iounit,*) "  dmspptype          = ", dmspptype
        write(iounit,*) "  dmspdtype          = ", dmspdtype
        write(iounit,*) "  humtype            = ", humtype
        write(iounit,*) "  doctype_s          = ", doctype_s
        write(iounit,*) "  doctype_l          = ", doctype_l
        write(iounit,*) "  dontype_protein    = ", dontype_protein
        write(iounit,*) "  fedtype_1          = ", fedtype_1
        write(iounit,*) "  feptype_1          = ", feptype_1
        write(iounit,*) "  zaerotype_bc1      = ", zaerotype_bc1
        write(iounit,*) "  zaerotype_bc2      = ", zaerotype_bc2
        write(iounit,*) "  zaerotype_dust1    = ", zaerotype_dust1
        write(iounit,*) "  zaerotype_dust2    = ", zaerotype_dust2
        write(iounit,*) "  zaerotype_dust3    = ", zaerotype_dust3
        write(iounit,*) "  zaerotype_dust4    = ", zaerotype_dust4
        write(iounit,*) "  ratio_C2N_diatoms  = ", ratio_C2N_diatoms
        write(iounit,*) "  ratio_C2N_sp       = ", ratio_C2N_sp
        write(iounit,*) "  ratio_C2N_phaeo    = ", ratio_C2N_phaeo
        write(iounit,*) "  ratio_chl2N_diatoms= ", ratio_chl2N_diatoms
        write(iounit,*) "  ratio_chl2N_sp     = ", ratio_chl2N_sp
        write(iounit,*) "  ratio_chl2N_phaeo  = ", ratio_chl2N_phaeo
        write(iounit,*) "  F_abs_chl_diatoms  = ", F_abs_chl_diatoms
        write(iounit,*) "  F_abs_chl_sp       = ", F_abs_chl_sp
        write(iounit,*) "  F_abs_chl_phaeo    = ", F_abs_chl_phaeo
        write(iounit,*) "  ratio_C2N_proteins = ", ratio_C2N_proteins

      end subroutine icepack_write_parameters

!=======================================================================

    end module icepack_parameters

!=======================================================================
