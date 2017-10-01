!  SVN:$Id: icepack_intfc.F90 1227 2017-05-22 22:49:10Z tcraig $
!=========================================================================
!
! flags and interface routines for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_intfc

      use icepack_kinds_mod    ! kinds

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

      use icepack_intfc_shared, only: icepack_init_parameters
      use icepack_intfc_shared, only: icepack_query_parameters
      use icepack_intfc_shared, only: icepack_write_parameters

      use icepack_intfc_tracers, only: icepack_compute_tracers

      implicit none

      public

      ! initialization
      public :: &
           icepack_init_tracer_flags, &
           icepack_init_tracer_indices, &
           icepack_init_tracer_numbers

!=======================================================================

      contains

!=======================================================================
!     Initialization routines
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
