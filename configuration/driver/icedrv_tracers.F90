!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model 
!
! author Elizabeth C. Hunke, LANL

      module icedrv_tracers

      use icedrv_kinds
      use icepack_tracers, only: ntrcr, ntrcr_o, nbtrcr, nbtrcr_sw

      use icepack_tracers, only: tr_brine, tr_iage, tr_fy, tr_lvl
      use icepack_tracers, only: tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_pond, tr_aero
      use icepack_tracers, only: tr_bgc_DMS, tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C
      use icepack_tracers, only: tr_bgc_DON, tr_bgc_hum, tr_bgc_Fe
      use icepack_tracers, only: tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil
      use icepack_tracers, only: tr_zaero, tr_bgc_n, tr_bgc_chl

      use icepack_tracers, only: nt_iage, nt_fbri, nt_Tsfc, nt_qice, nt_qsno, nt_sice
      use icepack_tracers, only: nt_fy, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd
      use icepack_tracers, only: nt_aero, nt_bgc_s, nt_zaero
      use icepack_tracers, only: nt_bgc_n, nt_bgc_doc, nt_bgc_don, nt_bgc_fed, nt_bgc_fep
      use icepack_tracers, only: nt_bgc_nit, nt_bgc_sil, nt_bgc_am, nt_bgc_hum, nt_bgc_pon
      use icepack_tracers, only: nt_bgc_dmspp, nt_bgc_dmspd, nt_bgc_dms

      use icepack_tracers, only: nlt_chl_sw, nlt_zaero, nlt_zaero_sw
      use icepack_tracers, only: nlt_bgc_n, nlt_bgc_doc, nlt_bgc_don, nlt_bgc_doc
      use icepack_tracers, only: nlt_bgc_dms, nlt_bgc_dmspp
      use icepack_tracers, only: nlt_bgc_Fed, nlt_bgc_Fep, nlt_bgc_hum
      use icepack_tracers, only: nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil
      
      use icepack_tracers, only: nt_bgc_chl, nt_bgc_C, nt_bgc_DOC, nt_bgc_DIC, nt_bgc_AM
      use icepack_tracers, only: nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_DMS
      use icepack_tracers, only: nt_bgc_PON, nt_bgc_DON, nt_bgc_Fed, nt_bgc_Fep
      use icepack_tracers, only: nt_zbgc_frac, nt_zaero

      use icepack_tracers, only: bio_index_o

!=======================================================================

      end module icedrv_tracers

!=======================================================================
