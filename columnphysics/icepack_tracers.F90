!=======================================================================
! Indices and flags associated with the tracer infrastructure. 
! Grid-dependent and max_trcr-dependent arrays are declared in ice_state.F90.
!
! author Elizabeth C. Hunke, LANL

      module icepack_tracers

      use icepack_kinds
      use icepack_parameters, only: c0, c1, puny, Tocnfrz
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: icepack_compute_tracers
      public :: icepack_query_tracer_sizes
      public :: icepack_write_tracer_sizes
      public :: icepack_init_tracer_flags
      public :: icepack_query_tracer_flags
      public :: icepack_write_tracer_flags
      public :: icepack_init_tracer_indices
      public :: icepack_query_tracer_indices
      public :: icepack_write_tracer_indices
      public :: icepack_init_tracer_numbers
      public :: icepack_query_tracer_numbers
      public :: icepack_write_tracer_numbers

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

      integer (kind=int_kind), public :: &
         ntrcr   , & ! number of tracers in use
         ntrcr_o , & ! number of non-bio tracers in use
         n_aero  , & ! number of aerosols in use
         n_zaero , & ! number of z aerosols in use 
         n_algae , & ! number of algae in use 
         n_doc   , & ! number of DOC pools in use
         n_dic   , & ! number of DIC pools in use
         n_don   , & ! number of DON pools in use
         n_fed   , & ! number of Fe  pools in use dissolved Fe
         n_fep       ! number of Fe  pools in use particulate Fe

      integer (kind=int_kind), public :: &
         nt_Tsfc  , & ! ice/snow temperature
         nt_qice  , & ! volume-weighted ice enthalpy (in layers)
         nt_qsno  , & ! volume-weighted snow enthalpy (in layers)
         nt_sice  , & ! volume-weighted ice bulk salinity (CICE grid layers)
         nt_fbri  , & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_iage  , & ! volume-weighted ice age
         nt_FY    , & ! area-weighted first-year ice area
         nt_alvl  , & ! level ice area fraction
         nt_vlvl  , & ! level ice volume fraction
         nt_apnd  , & ! melt pond area fraction
         nt_hpnd  , & ! melt pond depth
         nt_ipnd  , & ! melt pond refrozen lid thickness
         nt_aero  , & ! starting index for aerosols in ice
         nt_bgc_Nit,   & ! nutrients  
         nt_bgc_Am,    & ! 
         nt_bgc_Sil,   & !
         nt_bgc_DMSPp, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd, & ! 
         nt_bgc_DMS,   & ! 
         nt_bgc_PON,   & ! zooplankton and detritus 
         nt_bgc_hum,   & ! humic material 
         nt_zbgc_frac, & ! fraction of tracer in the mobile phase
         nt_bgc_S        ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

      logical (kind=log_kind), public :: &
         tr_iage     , & ! if .true., use age tracer
         tr_FY       , & ! if .true., use first-year area tracer
         tr_lvl      , & ! if .true., use level ice tracer
         tr_pond     , & ! if .true., use melt pond tracer
         tr_pond_cesm, & ! if .true., use cesm pond tracer
         tr_pond_lvl , & ! if .true., use level-ice pond tracer
         tr_pond_topo, & ! if .true., use explicit topography-based ponds
         tr_aero     , & ! if .true., use aerosol tracers
         tr_brine        ! if .true., brine height differs from ice thickness

      !-----------------------------------------------------------------
      !  biogeochemistry
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: & 
         tr_zaero,       & ! if .true., black carbon as tracers  (n_zaero)
         tr_bgc_Nit,     & ! if .true. Nitrate tracer in ice 
         tr_bgc_N,       & ! if .true., algal nitrogen tracers  (n_algae)
         tr_bgc_DON,     & ! if .true., DON pools are tracers  (n_don)
         tr_bgc_C,       & ! if .true., algal carbon tracers + DOC and DIC 
         tr_bgc_chl,     & ! if .true., algal chlorophyll tracers 
         tr_bgc_Am,      & ! if .true., ammonia/um as nutrient tracer 
         tr_bgc_Sil,     & ! if .true., silicon as nutrient tracer 
         tr_bgc_DMS,     & ! if .true., DMS as  tracer 
         tr_bgc_Fe,      & ! if .true., Fe as  tracer 
         tr_bgc_PON,     & ! if .true., PON as tracer 
         tr_bgc_hum        ! if .true., humic material as tracer 

      integer (kind=int_kind), public :: &
         nbtrcr,         & ! number of bgc tracers in use
         nbtrcr_sw,      & ! number of bgc tracers which impact shortwave
         nlt_chl_sw        ! points to total chla in trcrn_sw

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero_sw       ! points to aerosol in trcrn_sw
  
      integer (kind=int_kind), dimension(max_algae), public :: &
         nlt_bgc_N      , & ! algae 
         nlt_bgc_C      , & ! 
         nlt_bgc_chl   

      integer (kind=int_kind), dimension(max_doc), public :: &
         nlt_bgc_DOC        ! disolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: &
         nlt_bgc_DON        !

      integer (kind=int_kind), dimension(max_dic), public :: &
         nlt_bgc_DIC        ! disolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: &
         nlt_bgc_Fed    , & !
         nlt_bgc_Fep        !

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero          ! non-reacting layer aerosols

      integer (kind=int_kind), public :: &
         nlt_bgc_Nit   ,   & ! nutrients  
         nlt_bgc_Am    ,   & ! 
         nlt_bgc_Sil   ,   & !
         nlt_bgc_DMSPp ,   & ! trace gases (skeletal layer)
         nlt_bgc_DMSPd ,   & ! 
         nlt_bgc_DMS   ,   & ! 
         nlt_bgc_PON   ,   & ! zooplankton and detritus
         nlt_bgc_hum         ! humic material

      integer (kind=int_kind), dimension(max_algae), public :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(max_doc), public :: &  
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: & 
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(max_dic), public :: &  
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(max_aero), public :: &  
         nt_zaero       !  black carbon and other aerosols
      
      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index_o         ! relates nlt_bgc_NO to ocean concentration index
                             ! see ocean_bio_all

      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index           ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N 

!=======================================================================

      contains

!=======================================================================
! query tracer sizes

      subroutine icepack_query_tracer_sizes( &
           max_algae_out  , max_dic_out    , max_doc_out    , &
           max_don_out    , max_fe_out     , nmodal1_out    , &
           nmodal2_out    , max_aero_out   , max_nbtrcr_out )

        integer (kind=int_kind), intent(out), optional :: &
             max_algae_out  , & ! maximum number of algal types
             max_dic_out    , & ! maximum number of dissolved inorganic carbon types
             max_doc_out    , & ! maximum number of dissolved organic carbon types
             max_don_out    , & ! maximum number of dissolved organic nitrogen types
             max_fe_out     , & ! maximum number of iron types
             nmodal1_out    , & ! dimension for modal aerosol radiation parameters
             nmodal2_out    , & ! dimension for modal aerosol radiation parameters
             max_aero_out   , & ! maximum number of aerosols
             max_nbtrcr_out     ! algal nitrogen and chlorophyll

        if (present(max_algae_out))  max_algae_out = max_algae
        if (present(max_dic_out))    max_dic_out   = max_dic
        if (present(max_doc_out))    max_doc_out   = max_doc
        if (present(max_don_out))    max_don_out   = max_don
        if (present(max_fe_out))     max_fe_out    = max_fe
        if (present(nmodal1_out))    nmodal1_out   = nmodal1
        if (present(nmodal2_out))    nmodal2_out   = nmodal2
        if (present(max_aero_out))   max_aero_out  = max_aero
        if (present(max_nbtrcr_out)) max_nbtrcr_out= max_nbtrcr

      end subroutine icepack_query_tracer_sizes

!=======================================================================
! write tracer sizes

      subroutine icepack_write_tracer_sizes(iounit)

        integer, intent(in) :: iounit

        write(iounit,*) "icepack_write_tracer_sizes:"
        write(iounit,*) '  max_algae_out =', max_algae
        write(iounit,*) '  max_dic_out   =', max_dic
        write(iounit,*) '  max_doc_out   =', max_doc
        write(iounit,*) '  max_don_out   =', max_don
        write(iounit,*) '  max_fe_out    =', max_fe
        write(iounit,*) '  nmodal1_out   =', nmodal1
        write(iounit,*) '  nmodal2_out   =', nmodal2
        write(iounit,*) '  max_aero_out  =', max_aero
        write(iounit,*) '  max_nbtrcr_out=', max_nbtrcr

      end subroutine icepack_write_tracer_sizes

!=======================================================================
! set tracer active flags

      subroutine icepack_init_tracer_flags(&
           tr_iage_in, tr_FY_in, tr_lvl_in, &
           tr_pond_in, tr_pond_cesm_in, tr_pond_lvl_in, tr_pond_topo_in, &
           tr_aero_in, tr_brine_in, tr_zaero_in, &
           tr_bgc_Nit_in, tr_bgc_N_in, tr_bgc_DON_in, tr_bgc_C_in, tr_bgc_chl_in, &
           tr_bgc_Am_in, tr_bgc_Sil_in, tr_bgc_DMS_in, tr_bgc_Fe_in, tr_bgc_hum_in, &
           tr_bgc_PON_in)

        logical, intent(in), optional :: &
             tr_iage_in      , & ! if .true., use age tracer
             tr_FY_in        , & ! if .true., use first-year area tracer
             tr_lvl_in       , & ! if .true., use level ice tracer
             tr_pond_in      , & ! if .true., use melt pond tracer
             tr_pond_cesm_in , & ! if .true., use cesm pond tracer
             tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
             tr_aero_in      , & ! if .true., use aerosol tracers
             tr_brine_in     , & ! if .true., brine height differs from ice thickness
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

        if (present(tr_iage_in)) tr_iage = tr_iage_in
        if (present(tr_FY_in)  ) tr_FY   = tr_FY_in
        if (present(tr_lvl_in) ) tr_lvl  = tr_lvl_in
        if (present(tr_pond_in)) tr_pond = tr_pond_in
        if (present(tr_pond_cesm_in)) tr_pond_cesm = tr_pond_cesm_in
        if (present(tr_pond_lvl_in) ) tr_pond_lvl  = tr_pond_lvl_in
        if (present(tr_pond_topo_in)) tr_pond_topo = tr_pond_topo_in
        if (present(tr_aero_in)   ) tr_aero    = tr_aero_in
        if (present(tr_brine_in)  ) tr_brine   = tr_brine_in
        if (present(tr_zaero_in)  ) tr_zaero   = tr_zaero_in 
        if (present(tr_bgc_Nit_in)) tr_bgc_Nit = tr_bgc_Nit_in
        if (present(tr_bgc_N_in)  ) tr_bgc_N   = tr_bgc_N_in 
        if (present(tr_bgc_DON_in)) tr_bgc_DON = tr_bgc_DON_in
        if (present(tr_bgc_C_in)  ) tr_bgc_C   = tr_bgc_C_in 
        if (present(tr_bgc_chl_in)) tr_bgc_chl = tr_bgc_chl_in
        if (present(tr_bgc_Am_in) ) tr_bgc_Am  = tr_bgc_Am_in
        if (present(tr_bgc_Sil_in)) tr_bgc_Sil = tr_bgc_Sil_in
        if (present(tr_bgc_DMS_in)) tr_bgc_DMS = tr_bgc_DMS_in
        if (present(tr_bgc_Fe_in )) tr_bgc_Fe  = tr_bgc_Fe_in 
        if (present(tr_bgc_hum_in)) tr_bgc_hum = tr_bgc_hum_in
        if (present(tr_bgc_PON_in)) tr_bgc_PON = tr_bgc_PON_in 

      end subroutine icepack_init_tracer_flags

!=======================================================================
! query tracer active flags

      subroutine icepack_query_tracer_flags(&
           tr_iage_out, tr_FY_out, tr_lvl_out, &
           tr_pond_out, tr_pond_cesm_out, tr_pond_lvl_out, tr_pond_topo_out, &
           tr_aero_out, tr_brine_out, tr_zaero_out, &
           tr_bgc_Nit_out, tr_bgc_N_out, tr_bgc_DON_out, tr_bgc_C_out, tr_bgc_chl_out, &
           tr_bgc_Am_out, tr_bgc_Sil_out, tr_bgc_DMS_out, tr_bgc_Fe_out, tr_bgc_hum_out, &
           tr_bgc_PON_out)

        logical, intent(out), optional :: &
             tr_iage_out      , & ! if .true., use age tracer
             tr_FY_out        , & ! if .true., use first-year area tracer
             tr_lvl_out       , & ! if .true., use level ice tracer
             tr_pond_out      , & ! if .true., use melt pond tracer
             tr_pond_cesm_out , & ! if .true., use cesm pond tracer
             tr_pond_lvl_out  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_out , & ! if .true., use explicit topography-based ponds
             tr_aero_out      , & ! if .true., use aerosol tracers
             tr_brine_out     , & ! if .true., brine height differs from ice thickness
             tr_zaero_out     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit_out   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N_out     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON_out   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C_out     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl_out   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am_out    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil_out   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS_out   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe_out    , & ! if .true., Fe as product tracer 
             tr_bgc_hum_out   , & ! if .true., hum as product tracer 
             tr_bgc_PON_out       ! if .true., PON as product tracer 

        if (present(tr_iage_out)) tr_iage_out = tr_iage
        if (present(tr_FY_out)  ) tr_FY_out   = tr_FY
        if (present(tr_lvl_out) ) tr_lvl_out  = tr_lvl
        if (present(tr_pond_out)) tr_pond_out = tr_pond
        if (present(tr_pond_cesm_out)) tr_pond_cesm_out = tr_pond_cesm
        if (present(tr_pond_lvl_out) ) tr_pond_lvl_out  = tr_pond_lvl
        if (present(tr_pond_topo_out)) tr_pond_topo_out = tr_pond_topo
        if (present(tr_aero_out)   ) tr_aero_out    = tr_aero
        if (present(tr_brine_out)  ) tr_brine_out   = tr_brine
        if (present(tr_zaero_out)  ) tr_zaero_out   = tr_zaero
        if (present(tr_bgc_Nit_out)) tr_bgc_Nit_out = tr_bgc_Nit
        if (present(tr_bgc_N_out)  ) tr_bgc_N_out   = tr_bgc_N
        if (present(tr_bgc_DON_out)) tr_bgc_DON_out = tr_bgc_DON
        if (present(tr_bgc_C_out)  ) tr_bgc_C_out   = tr_bgc_C
        if (present(tr_bgc_chl_out)) tr_bgc_chl_out = tr_bgc_chl
        if (present(tr_bgc_Am_out) ) tr_bgc_Am_out  = tr_bgc_Am
        if (present(tr_bgc_Sil_out)) tr_bgc_Sil_out = tr_bgc_Sil
        if (present(tr_bgc_DMS_out)) tr_bgc_DMS_out = tr_bgc_DMS
        if (present(tr_bgc_Fe_out )) tr_bgc_Fe_out  = tr_bgc_Fe
        if (present(tr_bgc_hum_out)) tr_bgc_hum_out = tr_bgc_hum
        if (present(tr_bgc_PON_out)) tr_bgc_PON_out = tr_bgc_PON

      end subroutine icepack_query_tracer_flags

!=======================================================================
! write tracer active flags

      subroutine icepack_write_tracer_flags(iounit)

        integer, intent(in) :: iounit

        write(iounit,*) "icepack_write_tracer_flags:"
        write(iounit,*) "  tr_iage = ",tr_iage
        write(iounit,*) "  tr_FY   = ",tr_FY  
        write(iounit,*) "  tr_lvl  = ",tr_lvl 
        write(iounit,*) "  tr_pond = ",tr_pond
        write(iounit,*) "  tr_pond_cesm = ",tr_pond_cesm
        write(iounit,*) "  tr_pond_lvl  = ",tr_pond_lvl 
        write(iounit,*) "  tr_pond_topo = ",tr_pond_topo
        write(iounit,*) "  tr_aero    = ",tr_aero   
        write(iounit,*) "  tr_brine   = ",tr_brine  
        write(iounit,*) "  tr_zaero   = ",tr_zaero  
        write(iounit,*) "  tr_bgc_Nit = ",tr_bgc_Nit
        write(iounit,*) "  tr_bgc_N   = ",tr_bgc_N  
        write(iounit,*) "  tr_bgc_DON = ",tr_bgc_DON
        write(iounit,*) "  tr_bgc_C   = ",tr_bgc_C  
        write(iounit,*) "  tr_bgc_chl = ",tr_bgc_chl
        write(iounit,*) "  tr_bgc_Am  = ",tr_bgc_Am 
        write(iounit,*) "  tr_bgc_Sil = ",tr_bgc_Sil
        write(iounit,*) "  tr_bgc_DMS = ",tr_bgc_DMS
        write(iounit,*) "  tr_bgc_Fe  = ",tr_bgc_Fe 
        write(iounit,*) "  tr_bgc_hum = ",tr_bgc_hum
        write(iounit,*) "  tr_bgc_PON = ",tr_bgc_PON

      end subroutine icepack_write_tracer_flags

!=======================================================================
! set the number of column tracer indices

      subroutine icepack_init_tracer_indices(&
           nt_Tsfc_in, nt_qice_in, nt_qsno_in, nt_sice_in, &
           nt_fbri_in, nt_iage_in, nt_FY_in, & 
           nt_alvl_in, nt_vlvl_in, nt_apnd_in, nt_hpnd_in, nt_ipnd_in, &
           nt_aero_in, nt_zaero_in, &
           nt_bgc_N_in, nt_bgc_chl_in, nt_bgc_DOC_in, nt_bgc_DON_in, &
           nt_bgc_DIC_in, nt_bgc_Fed_in, nt_bgc_Fep_in, nt_bgc_Nit_in, nt_bgc_Am_in, &
           nt_bgc_Sil_in, nt_bgc_DMSPp_in, nt_bgc_DMSPd_in, nt_bgc_DMS_in, nt_bgc_hum_in, &
           nt_bgc_PON_in, nlt_zaero_in, nlt_bgc_N_in, nlt_bgc_chl_in, &
           nlt_bgc_DOC_in, nlt_bgc_DON_in, nlt_bgc_DIC_in, nlt_bgc_Fed_in, &
           nlt_bgc_Fep_in, nlt_bgc_Nit_in, nlt_bgc_Am_in, nlt_bgc_Sil_in, &
           nlt_bgc_DMSPp_in, nlt_bgc_DMSPd_in, nlt_bgc_DMS_in, nlt_bgc_hum_in, &
           nlt_bgc_PON_in, nt_zbgc_frac_in, nt_bgc_S_in, nlt_chl_sw_in, &
           nlt_zaero_sw_in, n_algae_in, n_DOC_in, &
           n_DON_in, n_DIC_in, n_fed_in, n_fep_in, n_zaero_in, &
           bio_index_o_in, bio_index_in, nbtrcr_in)

        integer, intent(in), optional :: &
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

       integer, intent(in), optional :: &
             n_algae_in,    & !  Dimensions
             n_DOC_in,      & !
             n_DON_in,      & !
             n_DIC_in,      & !
             n_fed_in,      & !
             n_fep_in,      & ! 
             n_zaero_in,    & !
             nbtrcr_in

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             bio_index_o_in, & 
             bio_index_in  

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
!            nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
!            nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_in   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DOC_in, & !  dissolved organic carbon
             nlt_bgc_DOC_in   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DON_in, & !  dissolved organic nitrogen
             nlt_bgc_DON_in   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DIC_in, & ! dissolved inorganic carbon
             nlt_bgc_DIC_in   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_Fed_in, & !  dissolved iron
             nt_bgc_Fep_in, & !  particulate iron
             nlt_bgc_Fed_in,& !  dissolved iron
             nlt_bgc_Fep_in   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_zaero_in,   & !  black carbon and other aerosols
             nlt_zaero_in,  & !  black carbon and other aerosols
             nlt_zaero_sw_in  ! black carbon and dust in trcrn_sw

        ! local
        integer (kind=int_kind) :: k

        if (present(nt_Tsfc_in)) nt_Tsfc = nt_Tsfc_in
        if (present(nt_qice_in)) nt_qice = nt_qice_in
        if (present(nt_qsno_in)) nt_qsno = nt_qsno_in
        if (present(nt_sice_in)) nt_sice = nt_sice_in
        if (present(nt_fbri_in)) nt_fbri = nt_fbri_in
        if (present(nt_iage_in)) nt_iage = nt_iage_in
        if (present(nt_FY_in)  ) nt_FY   = nt_FY_in
        if (present(nt_alvl_in)) nt_alvl = nt_alvl_in
        if (present(nt_vlvl_in)) nt_vlvl = nt_vlvl_in
        if (present(nt_apnd_in)) nt_apnd = nt_apnd_in
        if (present(nt_hpnd_in)) nt_hpnd = nt_hpnd_in
        if (present(nt_ipnd_in)) nt_ipnd = nt_ipnd_in
        if (present(nt_aero_in)) nt_aero = nt_aero_in
        if (present(nt_bgc_Nit_in)   ) nt_bgc_Nit    = nt_bgc_Nit_in
        if (present(nt_bgc_Am_in)    ) nt_bgc_Am     = nt_bgc_Am_in
        if (present(nt_bgc_Sil_in)   ) nt_bgc_Sil    = nt_bgc_Sil_in
        if (present(nt_bgc_DMSPp_in) ) nt_bgc_DMSPp  = nt_bgc_DMSPp_in
        if (present(nt_bgc_DMSPd_in) ) nt_bgc_DMSPd  = nt_bgc_DMSPd_in
        if (present(nt_bgc_DMS_in)   ) nt_bgc_DMS    = nt_bgc_DMS_in
        if (present(nt_bgc_hum_in)   ) nt_bgc_hum    = nt_bgc_hum_in
        if (present(nt_bgc_PON_in)   ) nt_bgc_PON    = nt_bgc_PON_in
        if (present(nlt_bgc_Nit_in)  ) nlt_bgc_Nit   = nlt_bgc_Nit_in
        if (present(nlt_bgc_Am_in)   ) nlt_bgc_Am    = nlt_bgc_Am_in
        if (present(nlt_bgc_Sil_in)  ) nlt_bgc_Sil   = nlt_bgc_Sil_in
        if (present(nlt_bgc_DMSPp_in)) nlt_bgc_DMSPp = nlt_bgc_DMSPp_in
        if (present(nlt_bgc_DMSPd_in)) nlt_bgc_DMSPd = nlt_bgc_DMSPd_in
        if (present(nlt_bgc_DMS_in)  ) nlt_bgc_DMS   = nlt_bgc_DMS_in
        if (present(nlt_bgc_hum_in)  ) nlt_bgc_hum   = nlt_bgc_hum_in
        if (present(nlt_bgc_PON_in)  ) nlt_bgc_PON   = nlt_bgc_PON_in
        if (present(nlt_chl_sw_in)   ) nlt_chl_sw    = nlt_chl_sw_in
        if (present(nt_zbgc_frac_in) ) nt_zbgc_frac  = nt_zbgc_frac_in
        if (present(nt_bgc_S_in)     ) nt_bgc_S      = nt_bgc_S_in

        if (present(nbtrcr_in)) then
           nbtrcr = nbtrcr_in
           do k = 1, nbtrcr_in
              if (present(bio_index_o_in)) bio_index_o(k)= bio_index_o_in(k)
              if (present(bio_index_in)  ) bio_index(k)  = bio_index_in(k)
           enddo
        endif

        if (present(n_algae_in)) then
           n_algae = n_algae_in
           do k = 1, n_algae_in
              if (present(nt_bgc_N_in) ) nt_bgc_N(k) = nt_bgc_N_in(k) 
              if (present(nlt_bgc_N_in)) nlt_bgc_N(k)= nlt_bgc_N_in(k) 
              if (present(nt_bgc_chl_in) ) nt_bgc_chl(k) = nt_bgc_chl_in(k) 
              if (present(nlt_bgc_chl_in)) nlt_bgc_chl(k)= nlt_bgc_chl_in(k) 
           enddo
        endif

! algal C is not yet distinct from algal N
!        if (present(n_algalC_in)) then
!           n_algalC = n_algalC_in
!           do k = 1, n_algalC_in
!              if (present(nt_bgc_C_in) ) nt_bgc_C(k) = nt_bgc_C_in(k) 
!              if (present(nlt_bgc_C_in)) nlt_bgc_C(k)= nlt_bgc_C_in(k) 
!           enddo
!        endif

        if (present(n_DOC_in)) then
           n_DOC = n_DOC_in
           do k = 1, n_DOC_in
              if (present(nt_bgc_DOC_in) ) nt_bgc_DOC(k) = nt_bgc_DOC_in(k) 
              if (present(nlt_bgc_DOC_in)) nlt_bgc_DOC(k)= nlt_bgc_DOC_in(k) 
           enddo
        endif

        if (present(n_DON_in)) then
           n_DON = n_DON_in
           do k = 1, n_DON_in
              if (present(nt_bgc_DON_in) ) nt_bgc_DON(k) = nt_bgc_DON_in(k) 
              if (present(nlt_bgc_DON_in)) nlt_bgc_DON(k)= nlt_bgc_DON_in(k) 
           enddo
        endif

        if (present(n_DIC_in)) then
           n_DIC = n_DIC_in
           do k = 1, n_DIC_in
              if (present(nt_bgc_DIC_in) ) nt_bgc_DIC(k) = nt_bgc_DIC_in(k) 
              if (present(nlt_bgc_DIC_in)) nlt_bgc_DIC(k)= nlt_bgc_DIC_in(k) 
           enddo
        endif

        if (present(n_fed_in)) then
           n_fed = n_fed_in
           do k = 1, n_fed_in
              if (present(nt_bgc_Fed_in) ) nt_bgc_Fed(k) = nt_bgc_Fed_in(k) 
              if (present(nlt_bgc_Fed_in)) nlt_bgc_Fed(k)= nlt_bgc_Fed_in(k) 
           enddo
        endif

        if (present(n_fep_in)) then
           n_fed = n_fep_in
           do k = 1, n_fep_in
              if (present(nt_bgc_Fep_in) ) nt_bgc_Fep(k) = nt_bgc_Fep_in(k) 
              if (present(nlt_bgc_Fep_in)) nlt_bgc_Fep(k)= nlt_bgc_Fep_in(k) 
           enddo
        endif

        if (present(n_zaero_in)) then
           n_zaero = n_zaero_in
           do k = 1, n_zaero_in
              if (present(nt_zaero_in)    ) nt_zaero(k)    = nt_zaero_in(k)   
              if (present(nlt_zaero_in)   ) nlt_zaero(k)   = nlt_zaero_in(k)   
              if (present(nlt_zaero_sw_in)) nlt_zaero_sw(k)= nlt_zaero_sw_in(k)   
           enddo
        endif

      end subroutine icepack_init_tracer_indices

!=======================================================================
! query the number of column tracer indices

      subroutine icepack_query_tracer_indices(&
           nt_Tsfc_out, nt_qice_out, nt_qsno_out, nt_sice_out, &
           nt_fbri_out, nt_iage_out, nt_FY_out, & 
           nt_alvl_out, nt_vlvl_out, nt_apnd_out, nt_hpnd_out, nt_ipnd_out, &
           nt_aero_out, nt_zaero_out, &
           nt_bgc_N_out, nt_bgc_C_out, nt_bgc_chl_out, nt_bgc_DOC_out, nt_bgc_DON_out, &
           nt_bgc_DIC_out, nt_bgc_Fed_out, nt_bgc_Fep_out, nt_bgc_Nit_out, nt_bgc_Am_out, &
           nt_bgc_Sil_out, nt_bgc_DMSPp_out, nt_bgc_DMSPd_out, nt_bgc_DMS_out, nt_bgc_hum_out, &
           nt_bgc_PON_out, nlt_zaero_out, nlt_bgc_N_out, nlt_bgc_C_out, nlt_bgc_chl_out, &
           nlt_bgc_DOC_out, nlt_bgc_DON_out, nlt_bgc_DIC_out, nlt_bgc_Fed_out, &
           nlt_bgc_Fep_out, nlt_bgc_Nit_out, nlt_bgc_Am_out, nlt_bgc_Sil_out, &
           nlt_bgc_DMSPp_out, nlt_bgc_DMSPd_out, nlt_bgc_DMS_out, nlt_bgc_hum_out, &
           nlt_bgc_PON_out, nt_zbgc_frac_out, nt_bgc_S_out, nlt_chl_sw_out, &
           nlt_zaero_sw_out, &
           bio_index_o_out, bio_index_out)

        integer, intent(out), optional :: &
             nt_Tsfc_out, & ! ice/snow temperature
             nt_qice_out, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_out, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_out, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_out, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_out, & ! volume-weighted ice age
             nt_FY_out, & ! area-weighted first-year ice area
             nt_alvl_out, & ! level ice area fraction
             nt_vlvl_out, & ! level ice volume fraction
             nt_apnd_out, & ! melt pond area fraction
             nt_hpnd_out, & ! melt pond depth
             nt_ipnd_out, & ! melt pond refrozen lid thickness
             nt_aero_out, & ! starting index for aerosols in ice
             nt_bgc_Nit_out, & ! nutrients  
             nt_bgc_Am_out,  & ! 
             nt_bgc_Sil_out, & !
             nt_bgc_DMSPp_out,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_out,&! 
             nt_bgc_DMS_out, & ! 
             nt_bgc_hum_out, & ! 
             nt_bgc_PON_out, & ! zooplankton and detritus   
             nlt_bgc_Nit_out,& ! nutrients  
             nlt_bgc_Am_out, & ! 
             nlt_bgc_Sil_out,& !
             nlt_bgc_DMSPp_out,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_out,&! 
             nlt_bgc_DMS_out,& ! 
             nlt_bgc_hum_out,& ! 
             nlt_bgc_PON_out,& ! zooplankton and detritus  
             nt_zbgc_frac_out,&! fraction of tracer in the mobile phase
             nt_bgc_S_out,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_out    ! points to total chla in trcrn_sw

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             bio_index_o_out, & 
             bio_index_out  

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_N_out ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C_out ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_out, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_out , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C_out , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_out   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_DOC_out, & !  dissolved organic carbon
             nlt_bgc_DOC_out   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_DON_out, & !  dissolved organic nitrogen
             nlt_bgc_DON_out   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_DIC_out, & ! dissolved inorganic carbon
             nlt_bgc_DIC_out   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_Fed_out, & !  dissolved iron
             nt_bgc_Fep_out, & !  particulate iron
             nlt_bgc_Fed_out,& !  dissolved iron
             nlt_bgc_Fep_out   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_zaero_out,   & !  black carbon and other aerosols
             nlt_zaero_out,  & !  black carbon and other aerosols
             nlt_zaero_sw_out  ! black carbon and dust in trcrn_sw

        if (present(nt_Tsfc_out)) nt_Tsfc_out = nt_Tsfc
        if (present(nt_qice_out)) nt_qice_out = nt_qice
        if (present(nt_qsno_out)) nt_qsno_out = nt_qsno
        if (present(nt_sice_out)) nt_sice_out = nt_sice
        if (present(nt_fbri_out)) nt_fbri_out = nt_fbri
        if (present(nt_iage_out)) nt_iage_out = nt_iage
        if (present(nt_FY_out)  ) nt_FY_out   = nt_FY
        if (present(nt_alvl_out)) nt_alvl_out = nt_alvl
        if (present(nt_vlvl_out)) nt_vlvl_out = nt_vlvl
        if (present(nt_apnd_out)) nt_apnd_out = nt_apnd
        if (present(nt_hpnd_out)) nt_hpnd_out = nt_hpnd
        if (present(nt_ipnd_out)) nt_ipnd_out = nt_ipnd
        if (present(nt_aero_out)) nt_aero_out = nt_aero
        if (present(nt_bgc_Nit_out)   ) nt_bgc_Nit_out    = nt_bgc_Nit
        if (present(nt_bgc_Am_out)    ) nt_bgc_Am_out     = nt_bgc_Am
        if (present(nt_bgc_Sil_out)   ) nt_bgc_Sil_out    = nt_bgc_Sil
        if (present(nt_bgc_DMSPp_out) ) nt_bgc_DMSPp_out  = nt_bgc_DMSPp
        if (present(nt_bgc_DMSPd_out) ) nt_bgc_DMSPd_out  = nt_bgc_DMSPd
        if (present(nt_bgc_DMS_out)   ) nt_bgc_DMS_out    = nt_bgc_DMS
        if (present(nt_bgc_hum_out)   ) nt_bgc_hum_out    = nt_bgc_hum
        if (present(nt_bgc_PON_out)   ) nt_bgc_PON_out    = nt_bgc_PON
        if (present(nlt_bgc_Nit_out)  ) nlt_bgc_Nit_out   = nlt_bgc_Nit
        if (present(nlt_bgc_Am_out)   ) nlt_bgc_Am_out    = nlt_bgc_Am
        if (present(nlt_bgc_Sil_out)  ) nlt_bgc_Sil_out   = nlt_bgc_Sil
        if (present(nlt_bgc_DMSPp_out)) nlt_bgc_DMSPp_out = nlt_bgc_DMSPp
        if (present(nlt_bgc_DMSPd_out)) nlt_bgc_DMSPd_out = nlt_bgc_DMSPd
        if (present(nlt_bgc_DMS_out)  ) nlt_bgc_DMS_out   = nlt_bgc_DMS
        if (present(nlt_bgc_hum_out)  ) nlt_bgc_hum_out   = nlt_bgc_hum
        if (present(nlt_bgc_PON_out)  ) nlt_bgc_PON_out   = nlt_bgc_PON
        if (present(nlt_chl_sw_out)   ) nlt_chl_sw_out    = nlt_chl_sw
        if (present(nt_zbgc_frac_out) ) nt_zbgc_frac_out  = nt_zbgc_frac
        if (present(nt_bgc_S_out)     ) nt_bgc_S_out      = nt_bgc_S

        if (present(bio_index_o_out) ) bio_index_o_out  = bio_index_o
        if (present(bio_index_out)   ) bio_index_out    = bio_index
        if (present(nt_bgc_N_out)    ) nt_bgc_N_out     = nt_bgc_N 
        if (present(nlt_bgc_N_out)   ) nlt_bgc_N_out    = nlt_bgc_N 
        if (present(nt_bgc_C_out)    ) nt_bgc_C_out     = nt_bgc_C 
        if (present(nlt_bgc_C_out)   ) nlt_bgc_C_out    = nlt_bgc_C 
        if (present(nt_bgc_chl_out)  ) nt_bgc_chl_out   = nt_bgc_chl 
        if (present(nlt_bgc_chl_out) ) nlt_bgc_chl_out  = nlt_bgc_chl 
        if (present(nt_bgc_DOC_out)  ) nt_bgc_DOC_out   = nt_bgc_DOC 
        if (present(nlt_bgc_DOC_out) ) nlt_bgc_DOC_out  = nlt_bgc_DOC 
        if (present(nt_bgc_DON_out)  ) nt_bgc_DON_out   = nt_bgc_DON 
        if (present(nlt_bgc_DON_out) ) nlt_bgc_DON_out  = nlt_bgc_DON 
        if (present(nt_bgc_DIC_out)  ) nt_bgc_DIC_out   = nt_bgc_DIC 
        if (present(nlt_bgc_DIC_out) ) nlt_bgc_DIC_out  = nlt_bgc_DIC 
        if (present(nt_bgc_Fed_out)  ) nt_bgc_Fed_out   = nt_bgc_Fed 
        if (present(nlt_bgc_Fed_out) ) nlt_bgc_Fed_out  = nlt_bgc_Fed 
        if (present(nt_bgc_Fep_out)  ) nt_bgc_Fep_out   = nt_bgc_Fep 
        if (present(nlt_bgc_Fep_out) ) nlt_bgc_Fep_out  = nlt_bgc_Fep 
        if (present(nt_zaero_out)    ) nt_zaero_out     = nt_zaero   
        if (present(nlt_zaero_out)   ) nlt_zaero_out    = nlt_zaero   
        if (present(nlt_zaero_sw_out)) nlt_zaero_sw_out = nlt_zaero_sw   

      end subroutine icepack_query_tracer_indices

!=======================================================================
! write the number of column tracer indices

      subroutine icepack_write_tracer_indices(iounit)

        integer, intent(in), optional :: iounit 

        ! local
        integer (kind=int_kind) :: k

        write(iounit,*) "icepack_write_tracer_indices:"
        write(iounit,*) "  nt_Tsfc = ",nt_Tsfc
        write(iounit,*) "  nt_qice = ",nt_qice
        write(iounit,*) "  nt_qsno = ",nt_qsno
        write(iounit,*) "  nt_sice = ",nt_sice
        write(iounit,*) "  nt_fbri = ",nt_fbri
        write(iounit,*) "  nt_iage = ",nt_iage
        write(iounit,*) "  nt_FY   = ",nt_FY  
        write(iounit,*) "  nt_alvl = ",nt_alvl
        write(iounit,*) "  nt_vlvl = ",nt_vlvl
        write(iounit,*) "  nt_apnd = ",nt_apnd
        write(iounit,*) "  nt_hpnd = ",nt_hpnd
        write(iounit,*) "  nt_ipnd = ",nt_ipnd
        write(iounit,*) "  nt_aero = ",nt_aero
        write(iounit,*) "  nt_bgc_Nit    = ",nt_bgc_Nit   
        write(iounit,*) "  nt_bgc_Am     = ",nt_bgc_Am    
        write(iounit,*) "  nt_bgc_Sil    = ",nt_bgc_Sil   
        write(iounit,*) "  nt_bgc_DMSPp  = ",nt_bgc_DMSPp 
        write(iounit,*) "  nt_bgc_DMSPd  = ",nt_bgc_DMSPd 
        write(iounit,*) "  nt_bgc_DMS    = ",nt_bgc_DMS   
        write(iounit,*) "  nt_bgc_hum    = ",nt_bgc_hum   
        write(iounit,*) "  nt_bgc_PON    = ",nt_bgc_PON   
        write(iounit,*) "  nlt_bgc_Nit   = ",nlt_bgc_Nit  
        write(iounit,*) "  nlt_bgc_Am    = ",nlt_bgc_Am   
        write(iounit,*) "  nlt_bgc_Sil   = ",nlt_bgc_Sil  
        write(iounit,*) "  nlt_bgc_DMSPp = ",nlt_bgc_DMSPp
        write(iounit,*) "  nlt_bgc_DMSPd = ",nlt_bgc_DMSPd
        write(iounit,*) "  nlt_bgc_DMS   = ",nlt_bgc_DMS  
        write(iounit,*) "  nlt_bgc_hum   = ",nlt_bgc_hum  
        write(iounit,*) "  nlt_bgc_PON   = ",nlt_bgc_PON  
        write(iounit,*) "  nlt_chl_sw    = ",nlt_chl_sw   
        write(iounit,*) "  nt_zbgc_frac  = ",nt_zbgc_frac 
        write(iounit,*) "  nt_bgc_S      = ",nt_bgc_S     

        write(iounit,*) "  max_nbtrcr = ",max_nbtrcr
        do k = 1, max_nbtrcr
           write(iounit,*) "  bio_index_o(k) = ",k,bio_index_o(k)
           write(iounit,*) "  bio_index(k)   = ",k,bio_index(k)  
        enddo

        write(iounit,*) "  max_algae = ",max_algae
        do k = 1, max_algae
           write(iounit,*) "  nt_bgc_N(k)  = ",k,nt_bgc_N(k)
           write(iounit,*) "  nlt_bgc_N(k) = ",k,nlt_bgc_N(k)
!           write(iounit,*) "  nt_bgc_C(k)  = ",k,nt_bgc_C(k)
!           write(iounit,*) "  nlt_bgc_C(k) = ",k,nlt_bgc_C(k)
           write(iounit,*) "  nt_bgc_chl(k)  = ",k,nt_bgc_chl(k) 
           write(iounit,*) "  nlt_bgc_chl(k) = ",k,nlt_bgc_chl(k)
        enddo

        write(iounit,*) "  max_DOC = ",max_DOC
        do k = 1, max_DOC
           write(iounit,*) "  nt_bgc_DOC(k)  = ",k,nt_bgc_DOC(k) 
           write(iounit,*) "  nlt_bgc_DOC(k) = ",k,nlt_bgc_DOC(k)
        enddo

        write(iounit,*) "  max_DON = ",max_DON
        do k = 1, max_DON
           write(iounit,*) "  nt_bgc_DON(k)  = ",k,nt_bgc_DON(k) 
           write(iounit,*) "  nlt_bgc_DON(k) = ",k,nlt_bgc_DON(k)
        enddo

        write(iounit,*) "  max_DIC = ",max_DIC
        do k = 1, max_DIC
           write(iounit,*) "  nt_bgc_DIC(k)  = ",k,nt_bgc_DIC(k) 
           write(iounit,*) "  nlt_bgc_DIC(k) = ",k,nlt_bgc_DIC(k)
        enddo

        write(iounit,*) "  max_fe = ",max_fe
        do k = 1, max_fe
           write(iounit,*) "  nt_bgc_Fed(k)  = ",k,nt_bgc_Fed(k) 
           write(iounit,*) "  nlt_bgc_Fed(k) = ",k,nlt_bgc_Fed(k)
           write(iounit,*) "  nt_bgc_Fep(k)  = ",k,nt_bgc_Fep(k) 
           write(iounit,*) "  nlt_bgc_Fep(k) = ",k,nlt_bgc_Fep(k)
        enddo

        write(iounit,*) "  max_aero = ",max_aero
        do k = 1, max_aero
           write(iounit,*) "  nt_zaero(k)     = ",k,nt_zaero(k)    
           write(iounit,*) "  nlt_zaero(k)    = ",k,nlt_zaero(k)   
           write(iounit,*) "  nlt_zaero_sw(k) = ",k,nlt_zaero_sw(k)
        enddo

      end subroutine icepack_write_tracer_indices

!=======================================================================
! set the number of column tracers

      subroutine icepack_init_tracer_numbers(&
         ntrcr_in, ntrcr_o_in, nbtrcr_in, nbtrcr_sw_in)

      integer (kind=int_kind), intent(in), optional :: &
         ntrcr_in  , &! number of tracers in use
         ntrcr_o_in, &! number of non-bio tracers in use
         nbtrcr_in , &! number of bio tracers in use
         nbtrcr_sw_in ! number of shortwave bio tracers in use

        if (present(ntrcr_in)    ) ntrcr     = ntrcr_in
        if (present(ntrcr_o_in)  ) ntrcr_o   = ntrcr_o_in
        if (present(nbtrcr_in)   ) nbtrcr    = nbtrcr_in
        if (present(nbtrcr_sw_in)) nbtrcr_sw = nbtrcr_sw_in

      end subroutine icepack_init_tracer_numbers

!=======================================================================
! query the number of column tracers

      subroutine icepack_query_tracer_numbers(&
         ntrcr_out, ntrcr_o_out, nbtrcr_out, nbtrcr_sw_out)

      integer (kind=int_kind), intent(out), optional :: &
         ntrcr_out  , &! number of tracers in use
         ntrcr_o_out, &! number of non-bio tracers in use
         nbtrcr_out , &! number of bio tracers in use
         nbtrcr_sw_out ! number of shortwave bio tracers in use

        if (present(ntrcr_out)    ) ntrcr_out     = ntrcr
        if (present(ntrcr_o_out)  ) ntrcr_o_out   = ntrcr_o
        if (present(nbtrcr_out)   ) nbtrcr_out    = nbtrcr
        if (present(nbtrcr_sw_out)) nbtrcr_sw_out = nbtrcr_sw

      end subroutine icepack_query_tracer_numbers

!=======================================================================
! write the number of column tracers

      subroutine icepack_write_tracer_numbers(iounit)

      integer (kind=int_kind), intent(in) :: iounit

        write(iounit,*) "icepack_write_tracer_numbers:"
        write(iounit,*) "  ntrcr     = ",ntrcr    
        write(iounit,*) "  nbtrcr    = ",nbtrcr   
        write(iounit,*) "  nbtrcr_sw = ",nbtrcr_sw

      end subroutine icepack_write_tracer_numbers

!=======================================================================

! Compute tracer fields.
! Given atrcrn = aicen*trcrn (or vicen*trcrn, vsnon*trcrn), compute trcrn.
!
! author: William H. Lipscomb, LANL

      subroutine icepack_compute_tracers (ntrcr,     trcr_depend,    &
                                         atrcrn,    aicen,          &
                                         vicen,     vsnon,          &
                                         trcr_base, n_trcr_strata,  &
                                         nt_strata, trcrn)

      integer (kind=int_kind), intent(in) :: &
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (ntrcr), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         atrcrn    ! aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (ntrcr), intent(out) :: &
         trcrn     ! ice tracers

      ! local variables

      integer (kind=int_kind) :: &
         it,     & ! tracer index
         itl,    & ! tracer index
         ntr,    & ! tracer index
         k         ! loop index

      real (kind=dbl_kind), dimension(3) :: &
         divisor   ! base quantity on which tracers are carried

      real (kind=dbl_kind) :: &
         work      ! temporary scalar

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do it = 1, ntrcr
         divisor(1) = trcr_base(it,1)*aicen
         divisor(2) = trcr_base(it,2)*vicen
         divisor(3) = trcr_base(it,3)*vsnon

         if (trcr_depend(it) == 0) then ! ice area tracers
            if (aicen > puny) then  
               trcrn(it) = atrcrn(it) / aicen
            else
               trcrn(it) = c0
               if (it == nt_Tsfc) trcrn(it) = Tocnfrz  ! surface temperature
            endif

         else

            work = c0
            do k = 1, 3
               if (divisor(k) > c0) then
                  work = atrcrn(it) / divisor(k)
               endif
            enddo
            trcrn(it) = work                ! save it
            if (n_trcr_strata(it) > 0) then          ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  if (trcrn(ntr) > c0) then
                      trcrn(it) = trcrn(it) / trcrn(ntr)
                  else
                      trcrn(it) = c0
                  endif
               enddo
            endif
            if (vicen <= c0 .and. it == nt_fbri) trcrn(it) = c1

         endif ! trcr_depend=0

      enddo

    end subroutine icepack_compute_tracers

!=======================================================================

      end module icepack_tracers

!=======================================================================
