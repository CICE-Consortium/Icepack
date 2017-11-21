!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model 
!
! author Elizabeth C. Hunke, LANL

      module icepack_drv_parameters

      use icepack_drv_kinds
      use icepack_parameters, only: max_nbtrcr, max_algae, max_aero
      use icepack_parameters, only: nmodal1, nmodal2, modal_aero
      use icepack_parameters, only: max_doc, max_don, max_dic, max_fe
      use icepack_parameters, only: formdrag, calc_Tsfc, ktherm, calc_strair
      use icepack_parameters, only: ustar_min, oceanmixed_ice
      use icepack_parameters, only: albicev, albicei, albsnowv, albsnowi
      use icepack_parameters, only: ahmax, shortwave, albedo_type
      use icepack_parameters, only: R_ice, R_pnd, R_snw
      use icepack_parameters, only: dT_mlt, rsnw_mlt, rhosi
      use icepack_parameters, only: kstrength, krdg_partic, krdg_redist, mu_rdg
      use icepack_parameters, only: atmbndy, highfreq
      use icepack_parameters, only: natmiter, kitd, kcatbound
      use icepack_parameters, only: hs0, dpscale, frzpnd, pndaspect
      use icepack_parameters, only: rfracmin, rfracmax
      use icepack_parameters, only: hs1, hp1, cf, heat_capacity
      use icepack_parameters, only: conduct, tfrz_option, kalg, fbot_xfer_type
      use icepack_parameters, only: a_rapid_mode, rac_rapid_mode, aspect_rapid_mode
      use icepack_parameters, only: dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy
      use icepack_parameters, only: dEdd_algae, solve_zsal, solve_zbgc, phi_snow
      use icepack_parameters, only: nit_data_type, sil_data_type, bgc_data_dir, max_dic
      use icepack_parameters, only: fe_data_type

!=======================================================================

      end module icepack_drv_parameters

!=======================================================================
