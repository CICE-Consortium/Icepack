!=======================================================================
!
! Biogeochemistry driver
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module icepack_zbgc

      use icepack_kinds
      use icepack_parameters

      !use icepack_parameters, only: c0, c1, c2, p001, p1, p5, puny
      !use icepack_parameters, only: depressT, rhosi, min_salin, salt_loss
      !use icepack_parameters, only: fr_resp, algal_vel, R_dFe2dust, dustFe_sol, T_max
      !use icepack_parameters, only: op_dep_min, fr_graze_s, fr_graze_e, fr_mort2min, fr_dFe
      !use icepack_parameters, only: k_nitrif, t_iron_conv, max_loss, max_dfe_doc1
      !use icepack_parameters, only: fr_resp_s, y_sk_DMS, t_sk_conv, t_sk_ox
      !use icepack_parameters, only: scale_bgc, ktherm, skl_bgc
      !use icepack_parameters, only: z_tracers, fsal, conserv_check

      use icepack_tracers, only: max_algae, max_dic, max_doc, max_don, max_fe
      use icepack_tracers, only: max_aero, max_nbtrcr
      use icepack_tracers, only: ncat, nilyr, nslyr, nblyr, nbtrcr, ntrcr, ntrcr_o
      use icepack_tracers, only: tr_brine, nt_fbri, nt_qice, nt_Tsfc
      use icepack_tracers, only: tr_zaero, tr_bgc_Nit, tr_bgc_N
      use icepack_tracers, only: tr_bgc_DON, tr_bgc_C, tr_bgc_chl
      use icepack_tracers, only: tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS
      use icepack_tracers, only: tr_bgc_Fe, tr_bgc_hum, tr_bgc_PON
      use icepack_tracers, only: nt_bgc_Nit, nlt_bgc_Nit
      use icepack_tracers, only: nt_bgc_N, nlt_bgc_N, nt_bgc_Am, nlt_bgc_Am
      use icepack_tracers, only: nt_bgc_DMSPp, nlt_bgc_DMSPp, nt_bgc_Sil, nlt_bgc_Sil
      use icepack_tracers, only: nt_bgc_DMSPd, nlt_bgc_DMSPd, nt_bgc_DMS, nlt_bgc_DMS
      use icepack_tracers, only: nt_bgc_hum, nlt_bgc_hum, nt_bgc_PON, nlt_bgc_PON
      use icepack_tracers, only: nt_bgc_C, nlt_bgc_C, nt_bgc_chl, nlt_bgc_chl
      use icepack_tracers, only: nt_bgc_DOC, nlt_bgc_DOC, nt_bgc_DON, nlt_bgc_DON
      use icepack_tracers, only: nt_bgc_DIC, nlt_bgc_DIC, nt_bgc_Fed, nlt_bgc_Fed
      use icepack_tracers, only: nt_sice, bio_index, bio_index_o
      use icepack_tracers, only: nt_zaero, nlt_zaero, nt_bgc_Fep, nlt_bgc_Fep
      use icepack_tracers, only: nlt_zaero_sw, nlt_chl_sw, nt_zbgc_frac
      use icepack_tracers, only: n_algae, n_doc, n_dic, n_don, n_fed, n_fep, n_zaero

      use icepack_zbgc_shared, only: zbgc_init_frac
      use icepack_zbgc_shared, only: zbgc_frac_init
      use icepack_zbgc_shared, only: bgrid, cgrid, igrid, icgrid
      use icepack_zbgc_shared, only: bgc_tracer_type, remap_zbgc
      use icepack_zbgc_shared, only: R_S2N, R_Si2N, R_Fe2C, R_Fe2N, R_Fe2DON, R_Fe2DOC
      use icepack_zbgc_shared, only: chlabs, alpha2max_low, beta2max
      use icepack_zbgc_shared, only: mu_max, grow_Tdep, fr_graze
      use icepack_zbgc_shared, only: mort_pre, mort_Tdep, k_exude
      use icepack_zbgc_shared, only: K_Nit, K_Am, K_Sil, K_Fe
      use icepack_zbgc_shared, only: f_don, kn_bac, f_don_Am, f_doc
      use icepack_zbgc_shared, only: f_exude, k_bac
      use icepack_zbgc_shared, only: tau_ret, tau_rel
      use icepack_zbgc_shared, only: R_C2N, R_chl2N, f_abs_chl, R_C2N_DON
      use icepack_zbgc_shared, only: doc_pool_fractions
      use icepack_zbgc_shared, only: algaltype, doctype, dictype
      use icepack_zbgc_shared, only: dontype, fedtype, feptype, zaerotype

      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_brine, only: preflushing_changes, compute_microS_mushy
      use icepack_brine, only: update_hbrine
      use icepack_algae, only: zbio, sklbio
      use icepack_therm_shared, only: calculate_Tin_from_qin
      use icepack_itd, only: column_sum, column_conservation_check

      implicit none

      private
      public :: add_new_ice_bgc, &
                lateral_melt_bgc, &
                icepack_init_bgc, &
                icepack_init_zbgc, &
                icepack_biogeochemistry, &
                icepack_load_ocean_bio_array, &
                icepack_init_ocean_bio

!=======================================================================

      contains

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (dt,         ncats,                &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vin0new,  &
                                  trcrn,                            &
                                  ocean_bio,  flux_bio,   hsurp,    &
                                  d_an_tot)

      integer (kind=int_kind), intent(in) :: &
         ncats       ! 1 without floe size distribution or ncat

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! initial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         d_an_tot

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         vin0new          ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         vi0_init        ! volume of new ice added to cat 1 (intial)

      real (kind=dbl_kind), intent(in) :: &
         hsurp           ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio       ! ocean concentration of biological tracer

! local

      integer (kind=int_kind) :: &
         location    , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         m           , & ! bio index
         n           , & ! ice category index
         k           , & ! ice layer index
         nbiolayer

      real (kind=dbl_kind) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final      ! brine volume summed over categories

      real (kind=dbl_kind) :: &
         vsurp       , & ! volume of new ice added to each cat
         vtmp            ! total volume of new and old ice

      real (kind=dbl_kind), dimension (ncat) :: &
         vbrin           ! trcrn(nt_fbri,n)*vicen(n)

      real (kind=dbl_kind) :: &
         vice_new    , & ! vicen_init + vsurp
         bio0new         ! ocean_bio * zbgc_init_fac

      real (kind=dbl_kind) :: &
         Tmlts       ! melting temperature (oC)

      character (len=char_len) :: &
         fieldid         ! field identifier

      character(len=*),parameter :: subname='(add_new_ice_bgc)'

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing
      !-----------------------------------------------------------------
      ! brine
      !-----------------------------------------------------------------

      zspace(:)  = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = p5*zspace(2)
      zspace(nblyr+1) = zspace(1)
      vbrin(:) = c0
      do n = 1, ncat
         vbrin(n) = vicen_init(n)
         if (tr_brine) vbrin(n) =  trcrn(nt_fbri,n)*vicen_init(n)
      enddo

      call column_sum (ncat,  vbrin,  vbri_init)
      if (icepack_warnings_aborted(subname)) return

      vbri_init = vbri_init + vi0_init
      do k = 1, nbtrcr
         flux_bio(k) = flux_bio(k) &
                            - vi0_init/dt*ocean_bio(k)*zbgc_init_frac(k)
      enddo
      !-----------------------------------------------------------------
      ! Distribute bgc in new ice volume among all ice categories by
      ! increasing ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area does not change
         ! add salt to the bottom , location = 1

       vsurp = c0
       vtmp = c0

      do n = 1,ncat

      if (hsurp > c0) then

         vtmp = vbrin(n)
         vsurp = hsurp * aicen_init(n)
         vbrin(n) = vbrin(n) + vsurp
         vice_new = vicen_init(n) + vsurp

         if (tr_brine .and. vice_new > puny) then !c0) then
            trcrn(nt_fbri,n) = vbrin(n)/vice_new
         elseif (tr_brine .and. vicen(n) <= c0) then
            trcrn(nt_fbri,n) = c1
         endif

         if (nbtrcr > 0) then
            do m = 1, nbtrcr
               bio0new = ocean_bio(m)*zbgc_init_frac(m)
               nbiolayer = nblyr+1
               call update_vertical_bio_tracers(nbiolayer, trcrn(bio_index(m):bio_index(m) + nblyr,n), &
                    vtmp, vbrin(n), bio0new,zspace(:))
            enddo !nbtrcr
            if (icepack_warnings_aborted(subname)) return
         endif       ! nbtrcr
      endif          ! hsurp > 0
      enddo          ! n

      !-----------------------------------------------------------------
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------
      do n = 1, ncats
      if (vin0new(n) > c0 .and. d_an_tot(n) > c0) then

         vbri1    = vbrin(n)
         vbrin(n) = vbrin(n) + vin0new(n)
         if (tr_brine .and. vicen(n) > puny) then
            trcrn(nt_fbri,n) = vbrin(n)/vicen(n)
         elseif (tr_brine .and. vicen(n) <= puny) then
            trcrn(nt_fbri,n) = c1
         endif

      ! Diffuse_bio handles concentration changes from ice growth/melt
      ! ice area changes
      ! add salt throughout, location = 0

         if (nbtrcr > 0 .and. vbrin(n) > puny) then
            do m = 1, nbtrcr
               bio0new = ocean_bio(m)*zbgc_init_frac(m)
               do k = 1, nblyr+1
                  trcrn(bio_index(m) + k-1,n) = &
                       (trcrn(bio_index(m) + k-1,n)*vbri1 + bio0new * vin0new(n))/vbrin(n)
               enddo
            enddo

            if (icepack_warnings_aborted(subname)) return

         endif           ! nbtrcr > 0
      endif              ! vin0new(n) > 0
      enddo              ! n = 1,ncats

      if (tr_brine .and. conserv_check) then
         call column_sum (ncat,   vbrin,  vbri_final)
         if (icepack_warnings_aborted(subname)) return

         fieldid = subname//':vbrin'
         call column_conservation_check (fieldid,                  &
                                         vbri_init, vbri_final,    &
                                         puny)
         if (icepack_warnings_aborted(subname)) return

      endif   ! conserv_check

      end subroutine add_new_ice_bgc

!=======================================================================

! When sea ice melts laterally, flux bgc to ocean

      subroutine lateral_melt_bgc (dt,                  &
                                   rsiden,   vicen_init,&
                                   trcrn,    flux_bio)

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vicen_init! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcrn     ! tracer array

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         rsiden     ! fraction of ice that melts laterally

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

      ! local variables

      real (kind=dbl_kind) :: &
         total_bio_initial, & ! initial column tracer concentration (mmol/m2)
         total_bio_final      ! final column tracer concentration (mmol/m20

      integer (kind=int_kind) :: &
         k     , & ! layer index
         m     , & !
         n         ! category index

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing

      character(len=*),parameter :: subname='(lateral_melt_bgc)'

      zspace(:)       = c1/real(nblyr,kind=dbl_kind)
      zspace(1)       = p5*zspace(2)
      zspace(nblyr+1) = zspace(1)

      do m = 1, nbtrcr
         do n = 1, ncat
         do k = 1, nblyr+1
            flux_bio(m) = flux_bio(m) + trcrn(nt_fbri,n) &
                        * vicen_init(n)*zspace(k)*trcrn(bio_index(m)+k-1,n) &
                        * rsiden(n)/dt
         enddo
         enddo
      enddo

      end subroutine lateral_melt_bgc

!=======================================================================
!autodocument_start icepack_init_bgc
!
      subroutine icepack_init_bgc( &
         sicen, trcrn, sss, ocean_bio_all, DOCPoolFractions)

      real (kind=dbl_kind), dimension(nilyr, ncat), intent(in) :: &
         sicen     ! salinity on the cice grid

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! subset of tracer array (only bgc)

      real (kind=dbl_kind), intent(in) :: &
         sss       ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      real (kind=dbl_kind), dimension (:), optional, intent(out) :: &
         DOCPoolFractions   ! Fraction of DOC in polysacharids, lipids, and proteins

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index
         n     , & ! category index
         mm        ! bio tracer index

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp     ! temporary, remapped tracers

      character(len=*),parameter :: subname='(icepack_init_bgc)'

      !-----------------------------------------------------------------------------
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------
      if (.not. restartbgc) then

         if (skl_bgc) then

            do  n = 1,ncat
            do mm = 1,nbtrcr
               ! bulk concentration (mmol or mg per m^3, or 10^-3 mmol/m^3)
               trcrn(bio_index(mm)-ntrcr_o, n) = ocean_bio_all(bio_index_o(mm))
            enddo       ! nbtrcr
            enddo       ! n

      !-----------------------------------------------------------------------------
      !    zbgc Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The vertical layer model uses prognosed porosity and layer depth
      !-----------------------------------------------------------------------------

         else   ! not skl_bgc

            if (scale_bgc .and. ktherm == 2) then
               trtmp(:) = c0
               do n = 1,ncat
                  call remap_zbgc(nilyr,    &
                                  1,                          &
                                  sicen(:,n),       trtmp,    &
                                  0,                nblyr+1,  &
                                  c1,               c1,       &
                                  cgrid(2:nilyr+1),           &
                                  igrid(1:nblyr+1),           &
                                  sicen(1,n)                  )
                  if (icepack_warnings_aborted(subname)) return

                  do mm = 1,nbtrcr
                  do k = 1, nblyr + 1
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) =   &
                          (trtmp(k)/sss*ocean_bio_all(bio_index_o(mm)))
                     trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
                  enddo  ! k
                  enddo  ! mm
               enddo     ! n

            elseif (nbtrcr > 0 .and. nt_fbri > 0) then ! not scale_bgc

               do n = 1,ncat
               do mm = 1,nbtrcr
               do k = 1, nblyr+1
                  trcrn(bio_index(mm)+k-1-ntrcr_o,n) = ocean_bio_all(bio_index_o(mm)) &
                                             * zbgc_init_frac(mm)
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo    ! k
               trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
               enddo    ! mm
               enddo    ! n

            endif  ! scale_bgc
         endif     ! skl_bgc
      endif       ! restart

      if (present(DOCPoolFractions)) then
         DOCPoolFractions(:) = c1
         if (.not. use_macromolecules) then
           do mm = 1,max_doc
              DOCPoolFractions(mm) = doc_pool_fractions(mm)
           end do
         end if
      endif

      end subroutine icepack_init_bgc

!=======================================================================
!autodocument_start icepack_init_zbgc
!

      subroutine icepack_init_zbgc (&
           zbgc_frac_init_in, zbgc_init_frac_in, tau_ret_in, tau_rel_in, bgc_tracer_type_in)

      real (kind=dbl_kind), dimension (:), intent(in), optional :: zbgc_frac_init_in(:)  ! initializes mobile fraction
      real (kind=dbl_kind), dimension (:), intent(in), optional :: bgc_tracer_type_in(:) ! described tracer in mobile or stationary phases
      real (kind=dbl_kind), dimension (:), intent(in), optional :: zbgc_init_frac_in(:)  ! fraction of ocean tracer  concentration in new ice
      real (kind=dbl_kind), dimension (:), intent(in), optional :: tau_ret_in(:)         ! retention timescale  (s), mobile to stationary phase
      real (kind=dbl_kind), dimension (:), intent(in), optional :: tau_rel_in(:)         ! release timescale    (s), stationary to mobile phase

!autodocument_end

      character(len=*),parameter :: subname='(icepack_init_zbgc)'

      !--------

      if (present(zbgc_frac_init_in))  zbgc_frac_init(:)  = zbgc_frac_init_in(:)
      if (present(bgc_tracer_type_in)) bgc_tracer_type(:) = bgc_tracer_type_in(:)
      if (present(zbgc_init_frac_in))  zbgc_init_frac(:)  =  zbgc_init_frac_in(:)
      if (present(tau_ret_in)) tau_ret(:) = tau_ret_in(:)
      if (present(tau_rel_in)) tau_rel(:) = tau_rel_in(:)

      R_Si2N(1) = ratio_Si2N_diatoms
      R_Si2N(2) = ratio_Si2N_sp
      R_Si2N(3) = ratio_Si2N_phaeo

      R_S2N(1) = ratio_S2N_diatoms
      R_S2N(2) = ratio_S2N_sp
      R_S2N(3) = ratio_S2N_phaeo

      R_Fe2C(1) = ratio_Fe2C_diatoms
      R_Fe2C(2) = ratio_Fe2C_sp
      R_Fe2C(3) = ratio_Fe2C_phaeo

      R_Fe2N(1) = ratio_Fe2N_diatoms
      R_Fe2N(2) = ratio_Fe2N_sp
      R_Fe2N(3) = ratio_Fe2N_phaeo

      R_C2N(1) = ratio_C2N_diatoms
      R_C2N(2) = ratio_C2N_sp
      R_C2N(3) = ratio_C2N_phaeo

      R_chl2N(1) = ratio_chl2N_diatoms
      R_chl2N(2) = ratio_chl2N_sp
      R_chl2N(3) = ratio_chl2N_phaeo

      F_abs_chl(1) = F_abs_chl_diatoms
      F_abs_chl(2) = F_abs_chl_sp
      F_abs_chl(3) = F_abs_chl_phaeo

      R_Fe2DON(1) = ratio_Fe2DON
      R_C2N_DON(1) = ratio_C2N_proteins

      R_Fe2DOC(1) = ratio_Fe2DOC_s
      R_Fe2DOC(2) = ratio_Fe2DOC_l
      R_Fe2DOC(3) = c0

      chlabs(1) = chlabs_diatoms
      chlabs(2) = chlabs_sp
      chlabs(3) = chlabs_phaeo

      alpha2max_low(1) = alpha2max_low_diatoms
      alpha2max_low(2) = alpha2max_low_sp
      alpha2max_low(3) = alpha2max_low_phaeo

      beta2max(1) = beta2max_diatoms
      beta2max(2) = beta2max_sp
      beta2max(3) = beta2max_phaeo

      mu_max(1) = mu_max_diatoms
      mu_max(2) = mu_max_sp
      mu_max(3) = mu_max_phaeo

      grow_Tdep(1) = grow_Tdep_diatoms
      grow_Tdep(2) = grow_Tdep_sp
      grow_Tdep(3) = grow_Tdep_phaeo

      fr_graze(1) = fr_graze_diatoms
      fr_graze(2) = fr_graze_sp
      fr_graze(3) = fr_graze_phaeo

      mort_pre(1) = mort_pre_diatoms
      mort_pre(2) = mort_pre_sp
      mort_pre(3) = mort_pre_phaeo

      mort_Tdep(1) = mort_Tdep_diatoms
      mort_Tdep(2) = mort_Tdep_sp
      mort_Tdep(3) = mort_Tdep_phaeo

      k_exude(1) = k_exude_diatoms
      k_exude(2) = k_exude_sp
      k_exude(3) = k_exude_phaeo

      K_Nit(1) = K_Nit_diatoms
      K_Nit(2) = K_Nit_sp
      K_Nit(3) = K_Nit_phaeo

      K_Am(1) = K_Am_diatoms
      K_Am(2) = K_Am_sp
      K_Am(3) = K_Am_phaeo

      K_Sil(1) = K_Sil_diatoms
      K_Sil(2) = K_Sil_sp
      K_Sil(3) = K_Sil_phaeo

      K_Fe(1) = K_Fe_diatoms
      K_Fe(2) = K_Fe_sp
      K_Fe(3) = K_Fe_phaeo

      f_doc(:) = c0
      f_doc(1) = f_doc_s
      f_doc(2) = f_doc_l

      f_don(1) = f_don_protein
      kn_bac(1) = kn_bac_protein
      f_don_Am(1) = f_don_Am_protein

      f_exude(:) = c0
      f_exude(1) = f_exude_s
      f_exude(2) = f_exude_l

      k_bac(:) = c0
      k_bac(1) = k_bac_s
      k_bac(2) = k_bac_l

      algaltype(1) = algaltype_diatoms
      algaltype(2) = algaltype_sp
      algaltype(3) = algaltype_phaeo

      doctype(:) = c0
      doctype(1) = doctype_s
      doctype(2) = doctype_l

      dictype(:) = dictype_1

      dontype(:) = dontype_protein

      fedtype(:) = fedtype_1
      feptype(:) = feptype_1

      zaerotype(1) = zaerotype_bc1
      zaerotype(2) = zaerotype_bc2
      zaerotype(3) = zaerotype_dust1
      zaerotype(4) = zaerotype_dust2
      zaerotype(5) = zaerotype_dust3
      zaerotype(6) = zaerotype_dust4

      end subroutine icepack_init_zbgc

!=======================================================================
!autodocument_start icepack_biogeochemistry
!
      subroutine icepack_biogeochemistry(dt, &
                           upNO, upNH, iDi, iki, zfswin, &
                           darcy_V, grow_net,  &
                           PP_net, hbri,dhbr_bot, dhbr_top, Zoo,&
                           fbio_snoice, fbio_atmice, ocean_bio_dh, ocean_bio, &
                           first_ice, fswpenln, bphi, bTiz, ice_bio_net,  &
                           snow_bio_net, totalChla, fswthrun, &
                           meltbn, melttn, congeln, snoicen, &
                           sst, sss, Tf, fsnow, meltsn, & !hmix, &
                           hin_old, flux_bio, flux_bio_atm, &
                           aicen_init, vicen_init, aicen, vicen, vsnon, &
                           aice0, trcrn, vsnon_init, &
                           flux_bion, bioPorosityIceCell, &
                           bioSalinityIceCell, bioTemperatureIceCell)

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         fbio_snoice   , &  ! fluxes from snow to ice
         fbio_atmice   , &  ! fluxes from atm to ice
         dhbr_top      , &  ! brine top change
         dhbr_bot      , &  ! brine bottom change
         darcy_V       , &  ! darcy velocity positive up (m/s)
         hin_old       , &  ! old ice thickness
         ice_bio_net   , &  ! depth integrated tracer (mmol/m^2)
         snow_bio_net  , &  ! depth integrated snow tracer (mmol/m^2)
         flux_bio           ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio          ! contains the ocean bgc tracer concentrations in use (mmol/m^3)

      real (kind=dbl_kind), optional, dimension (:), intent(out) :: &
         ocean_bio_dh       ! The ocean bgc tracer concentrations in use * brine thickness * phi (mmol/m^2)

      logical (kind=log_kind), dimension (:), intent(inout) :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

      real (kind=dbl_kind), optional, dimension (:,:), intent(out) :: &
         flux_bion      ! per categeory ice to ocean biogeochemistry flux (mmol/m2/s)

      real (kind=dbl_kind), optional, dimension (:), intent(inout) :: &
         bioPorosityIceCell, & ! category average porosity on the interface bio grid
         bioSalinityIceCell, & ! (ppt) category average porosity on the interface bio grid
         bioTemperatureIceCell ! (oC) category average porosity on the interface bio grid

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         Zoo            , & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                            ! mmol/m^3
         bphi           , & ! porosity of layers
         bTiz           , & ! layer temperatures interpolated on bio grid (C)
         zfswin         , & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi            , & ! igrid Diffusivity (m^2/s)
         iki            , & ! Ice permeability (m^2)
         trcrn     ! tracers

      real (kind=dbl_kind), intent(inout) :: &
         grow_net       , & ! Specific growth rate (/s) per grid cell
         PP_net         , & ! Total production (mg C/m^2/s) per grid cell
         hbri           , & ! brine height, area-averaged for comparison with hi (m)
         upNO           , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH               ! ammonium uptake rate (mmol/m^2/d) times aice

      real (kind=dbl_kind), optional, intent(inout) :: &
         totalChla          ! ice integrated chla and summed over all algal groups (mg/m^2)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fswpenln        ! visible SW entering ice layers (W m-2)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         meltsn      , & ! snow melt in category n (m)
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen     , & ! snow-ice formation in category n (m)
         flux_bio_atm, & ! all bio fluxes to ice from atmosphere
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init  , & ! initial snow volume (m), for aerosol
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(in) :: &
         aice0   , & ! open water area fraction
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         !hmix    , & ! mixed layer depth (m)
         Tf      , & ! basal freezing temperature (C)
         fsnow       ! snowfall rate (kg/m^2 s)

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k           , & ! vertical index
         n, mm           ! thickness category index

      real (kind=dbl_kind) :: &
         hin         , & ! new ice thickness
         hsn         , & ! snow thickness  (m)
         hbr_old     , & ! old brine thickness before growh/melt
         dhice       , & ! change due to sublimation/condensation (m)
         kavg        , & ! average ice permeability (m^2)
         bphi_o      , & ! surface ice porosity
         hbrin       , & ! brine height
         dh_direct       ! surface flooding or runoff

      real (kind=dbl_kind), dimension (nblyr+2) :: &
      ! Defined on Bio Grid points
         bSin        , & ! salinity on the bio grid  (ppt)
         brine_sal   , & ! brine salinity (ppt)
         brine_rho       ! brine_density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr+1) :: &
      ! Defined on Bio Grid interfaces
         iphin       , & ! porosity
         ibrine_sal  , & ! brine salinity  (ppt)
         ibrine_rho  , & ! brine_density (kg/m^3)
         iSin        , & ! Salinity on the interface grid (ppt)
         iTin            ! Temperature on the interface grid (oC)

      real (kind=dbl_kind) :: &
         sloss           ! brine flux contribution from surface runoff (g/m^2)

      real (kind=dbl_kind), dimension (nbtrcr) :: &
         flux_bion_n     ! temporary

      ! for bgc sk
      real (kind=dbl_kind) :: &
         dh_bot_chl  , & ! Chlorophyll may or may not flush
         dh_top_chl  , & ! Chlorophyll may or may not flush
         darcy_V_chl

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing

      character(len=*),parameter :: subname='(icepack_biogeochemistry)'

      zspace(:)       = c1/real(nblyr,kind=dbl_kind)
      zspace(1)       = p5*zspace(2)
      zspace(nblyr+1) = zspace(1)

      if (present(bioPorosityIceCell)) bioPorosityIceCell(:) = c0
      if (present(bioSalinityIceCell)) bioSalinityIceCell(:) = c0
      if (present(bioTemperatureIceCell)) bioTemperatureIceCell(:) = c0

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------
         flux_bion_n(:) = c0
         hin_old(n) = c0
         if (aicen_init(n) > puny) then
            hin_old(n) = vicen_init(n) &
                                / aicen_init(n)
         else
            first_ice(n) = .true.
            if (tr_brine) trcrn(nt_fbri,n) = c1
            if (z_tracers) then
               do mm = 1,nbtrcr
                  trcrn(nt_zbgc_frac-1+mm,n) = zbgc_frac_init(mm)
               enddo
            endif
         endif

         if (aicen(n) > puny) then

            dh_top_chl = c0
            dh_bot_chl = c0
            darcy_V_chl= c0
            bSin(:)    = c0
            hsn        = c0
            hin        = c0
            hbrin      = c0
            kavg       = c0
            bphi_o     = c0
            sloss      = c0

      !-----------------------------------------------------------------
      ! brine dynamics
      !-----------------------------------------------------------------

            dhbr_top(n) = c0
            dhbr_bot(n) = c0

            if (tr_brine) then
               if (trcrn(nt_fbri,n) .le. c0) trcrn(nt_fbri,n) = c1

               dhice = c0
               call preflushing_changes  (aicen  (n),   &
                                 vicen   (n), vsnon  (n),   &
                                 meltbn  (n), melttn (n),   &
                                 congeln (n), snoicen(n),   &
                                 hin_old (n), dhice,        &
                                 trcrn(nt_fbri,n),          &
                                 dhbr_top(n), dhbr_bot(n),  &
                                 hbr_old,     hin,          &
                                 hsn)
               if (icepack_warnings_aborted(subname)) return

               ! Requires the average ice permeability = kavg(:)
               ! and the surface ice porosity = zphi_o(:)
               ! computed in "compute_microS" or from "thermosaline_vertical"

               iDi(:,n) = c0

               call compute_microS_mushy ( &
                           trcrn(:,n),    hin_old(n),    hbr_old,     &
                           sss,           sst,           bTiz(:,n),   &
                           iTin(:),       bphi(:,n),     kavg,        &
                           bphi_o,        bSin(:),     &
                           brine_sal(:),  brine_rho(:),  iphin(:),    &
                           ibrine_rho(:), ibrine_sal(:),              &
                           iDi(:,n)     , iSin(:))
               if (icepack_warnings_aborted(subname)) return

               call update_hbrine (melttn(n),   &
                                   meltsn  (n), dt,          &
                                   hin,         hsn,         &
                                   hin_old (n), hbrin,       &
                                   hbr_old,     &
                                   trcrn(nt_fbri,n),         &
                                   dhbr_top(n), dhbr_bot(n), &
                                   dh_top_chl,  dh_bot_chl,  &
                                   kavg,        bphi_o,      &
                                   darcy_V (n), darcy_V_chl, &
                                   bphi(2,n),   aice0,       &
                                   dh_direct)
               if (icepack_warnings_aborted(subname)) return

               hbri = hbri + hbrin * aicen(n)

            endif ! tr_brine

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

            if (z_tracers) then

               call zbio (dt,                                            &
                          melttn(n),                                     &
                          meltsn(n),             meltbn  (n),            &
                          congeln(n),            snoicen(n),             &
                          fsnow,                 trcrn(1:ntrcr,n),       &
                          bio_index(1:nbtrcr),   bio_index_o(:),         &
                          aicen_init(n),                                 &
                          vicen_init(n),         vsnon_init(n),          &
                          vicen(n),              vsnon(n),               &
                          aicen(n),              flux_bio_atm(1:nbtrcr), &
                          n,                     first_ice(n),           &
                          hin_old(n),            ocean_bio(1:nbtrcr),    &
                          ocean_bio_dh,                                  &
                          bphi(:,n),             iphin,                  &
                          iDi(:,n),                                      &
                          fswpenln(:,n),                                 &
                          dhbr_top(n),           dhbr_bot(n),            &
                          zfswin(:,n),                                   &
                          hbrin,                 hbr_old,                &
                          darcy_V(n),                                    &
                          bphi_o,                                        &
                          iTin,                                          &
                          Zoo(:,n),                                      &
                          flux_bio(1:nbtrcr),    dh_direct,              &
                          upNO,                  upNH,                   &
                          fbio_snoice,           fbio_atmice,            &
                          PP_net,                ice_bio_net (1:nbtrcr), &
                          snow_bio_net(1:nbtrcr),grow_net,               &
                          totalChla,                                     &
                          flux_bion_n(1:nbtrcr), iSin,                   &
                          bioPorosityIceCell, bioSalinityIceCell,        &
                          bioTemperatureIceCell                          )
               if (icepack_warnings_aborted(subname)) return

               if (present(flux_bion)) then
                  do mm = 1, nbtrcr
                     flux_bion(mm,n) = flux_bion_n(mm)
                  enddo
               endif

            elseif (skl_bgc) then

               call sklbio (dt,                      Tf,                  &
                            flux_bio (1:nbtrcr),     ocean_bio(1:nbtrcr), &
                            aicen    (n),        &
                            meltbn   (n),            congeln  (n),        &
                            fswthrun (n),            first_ice(n),        &
                            trcrn    (1:ntrcr,n), &
                            PP_net,                  upNO,                &
                            upNH,                    grow_net             )
               if (icepack_warnings_aborted(subname)) return

            endif  ! skl_bgc

            first_ice(n) = .false.

         else
            if (z_tracers) then
            do mm = 1, nbtrcr
               do k  = 1, nblyr+1
                  if (present(flux_bion)) then
                     flux_bion(mm,n) = flux_bion(mm,n) + trcrn(bio_index(mm) + k-1,n) *  &
                        hin_old(n) * zspace(k)/dt * trcrn(nt_fbri,n)
                  endif
                  flux_bio(mm) = flux_bio(mm) + trcrn(bio_index(mm) + k-1,n) * &
                     vicen_init(n) * zspace(k)/dt * trcrn(nt_fbri,n)
                  trcrn(bio_index(mm) + k-1,n) = c0
                enddo
            enddo
            endif
         endif             ! aicen > puny
      enddo                ! ncat

      end subroutine icepack_biogeochemistry

!=======================================================================
!autodocument_start icepack_load_ocean_bio_array
! basic initialization for ocean_bio_all

      subroutine icepack_load_ocean_bio_array(&
          nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, ocean_bio_all, hum)

      real (kind=dbl_kind), intent(in) :: &
         nit         , & ! ocean nitrate (mmol/m^3)
         amm         , & ! ammonia/um (mmol/m^3)
         sil         , & ! silicate (mmol/m^3)
         dmsp        , & ! dmsp (mmol/m^3)
         dms         , & ! dms (mmol/m^3)
         hum             ! humic material (mmol/m^3)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         algalN          ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         don             ! ocean don (mmol/m^3)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         dic             ! ocean dic (mmol/m^3)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fed, fep        ! ocean disolved and particulate fe (nM)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zaeros          ! ocean aerosols (mmol/m^3)

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k, ks           ! tracer indices

      character(len=*),parameter :: subname='(icepack_load_ocean_bio_array)'

      ocean_bio_all(:) = c0

      do k = 1, max_algae
         ocean_bio_all(k)      = algalN(k)           ! N
         ks = max_algae + max_doc + max_dic + 1
         ocean_bio_all(ks + k) = R_chl2N(k)*algalN(k)!chl
      enddo

      ks = max_algae + 1
      do k = 1, max_doc
         ocean_bio_all(ks + k) = doc(k)              ! doc
      enddo
      ks = ks + max_doc
      do k = 1, max_dic
         ocean_bio_all(ks + k) = dic(k)              ! dic
      enddo

      ks = 2*max_algae + max_doc + max_dic + 7
      do k = 1, max_don
         ocean_bio_all(ks + k) = don(k)              ! don
      enddo

      ks = max_algae + 1
      ocean_bio_all(ks) = nit                        ! nit

      ks = 2*max_algae + max_doc + 2 + max_dic
      ocean_bio_all(ks) = amm                        ! Am
      ks = ks + 1
      ocean_bio_all(ks) = sil                        ! Sil
      ks = ks + 1
      ocean_bio_all(ks) =  R_S2N(1)*algalN(1) &      ! DMSPp
                        +  R_S2N(2)*algalN(2) &
                        +  R_S2N(3)*algalN(3)
      ks = ks + 1
      ocean_bio_all(ks) = dmsp                       ! DMSPd
      ks = ks + 1
      ocean_bio_all(ks) = dms                        ! DMS
      ks = ks + 1
      ocean_bio_all(ks) = nit                        ! PON
      ks = 2*max_algae + max_doc + 7 + max_dic + max_don
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fed(k)              ! fed
      enddo
      ks = ks + max_fe
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fep(k)              ! fep
      enddo
      ks = ks + max_fe
      do k = 1, max_aero
         ocean_bio_all(ks+k) = zaeros(k)             ! zaero
      enddo
      ks = ks + max_aero + 1
      ocean_bio_all(ks)  = hum                       ! humics

      end subroutine icepack_load_ocean_bio_array

!=======================================================================
!autodocument_start icepack_init_ocean_bio
!  Initialize ocean concentration

      subroutine icepack_init_ocean_bio (amm, dmsp, dms, algalN, doc, dic, don, &
             fed, fep, hum, nit, sil, zaeros,CToN, CToN_DON)

      real (kind=dbl_kind), intent(out), optional:: &
       amm      , & ! ammonium
       dmsp     , & ! DMSPp
       dms      , & ! DMS
       hum      , & ! humic material
       nit      , & ! nitrate
       sil          ! silicate

      real (kind=dbl_kind), dimension(:), intent(out), optional:: &
       algalN   , & ! algae
       doc      , & ! DOC
       dic      , & ! DIC
       don      , & ! DON
       fed      , & ! Dissolved Iron
       fep      , & ! Particulate Iron
       zaeros       ! BC and dust

      real (kind=dbl_kind), dimension(:), intent(out), optional :: &
       CToN     , & ! carbon to nitrogen ratio for algae
       CToN_DON     ! nitrogen to carbon ratio for proteins

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
        k

      character(len=*),parameter :: subname='(icepack_init_ocean_bio)'

       if (present(CToN)) then
          CToN(:) = c0
          CToN(1) = R_C2N(1)
          CToN(2) = R_C2N(2)
          CToN(3) = R_C2N(3)
       endif

       if (present(CToN_DON)) then
          CToN_DON(:) = c0
          CToN_DON(1) = R_C2N_DON(1)
       endif

       if (present(amm)) &
          amm  = c1 ! ISPOL < 1 mmol/m^3
       if (present(dmsp)) &
          dmsp = p1
       if (present(dms)) &
          dms  = p1
       if (present(algalN)) then
          algalN(:) = c0
          algalN(1) = c1  !0.0026_dbl_kind ! ISPOL, Lannuzel 2013(pennate)
          algalN(2) = 0.0057_dbl_kind ! ISPOL, Lannuzel 2013(small plankton)
          algalN(3) = 0.0027_dbl_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
       endif
       if (present(doc))then
          doc(:) = c0
          doc(1) = 16.2_dbl_kind ! 18% saccharides
          doc(2) = 9.0_dbl_kind  ! lipids
          doc(3) = c1 !
       endif
       if (present(dic)) then
          dic(:) = c0
          dic(1) = 1950.0_dbl_kind ! 1950-2260 mmol C/m3 (Tynan et al. 2015)
       endif
       if (present(don)) then
          don(:) = c0
          don(1) = 12.9_dbl_kind ! 64.3_dbl_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
       endif
       if (present(fed)) then
          fed(:) = c0
          fed(1) = 0.4_dbl_kind ! c1 (nM) Lannuzel2007 DFe,! range 0.14-2.6 (nM) van der Merwe 2011
           ! Tagliabue 2012 (0.4 nM)
       endif
       if (present(fep)) then
          fep(:) = c0
          fep(1) = c2 ! (nM) van der Merwe 2011
          ! (0.6 to 2.9 nM ocean)
       endif
       if (present(hum)) &
            hum  = c1        ! mmol C/m^3
       if (present(nit)) &
            nit  = 12.0_dbl_kind
       if (present(sil)) &
            sil  = 25.0_dbl_kind
       if (present(zaeros)) &
          zaeros(:) = c0

      end subroutine icepack_init_ocean_bio
!
!=======================================================================
!
! Given some added new ice to the base of the existing ice, recalculate
! vertical bio tracer so that new grid cells are all the same size.
!
! author: N. Jeffery, LANL
! date  : Mar 3, 2024
!
! Based on update_vertical_tracers modified for vertical biogeochemistry
!
      subroutine update_vertical_bio_tracers(nbiolyr, trc, h1, h2, trc0, zspace)

      integer (kind=int_kind), intent(in) :: &
         nbiolyr ! number of bio layers nblyr+1

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
           trc ! vertical tracer

      real (kind=dbl_kind), intent(in) :: &
         h1, & ! old thickness
         h2, & ! new thickness
         trc0  ! tracer value of added ice on ice bottom

      real (kind=dbl_kind), dimension(nbiolyr), intent(in) :: &
         zspace

      ! local variables

      real(kind=dbl_kind), dimension(nbiolyr) :: trc2 ! updated tracer temporary

      ! vertical indices for old and new grid
      integer :: k1, k2

      real (kind=dbl_kind) :: &
         z1a, z1b, & ! upper, lower boundary of old cell/added new ice at bottom
         z2a, z2b, & ! upper, lower boundary of new cell
         overlap , & ! overlap between old and new cell
         rnilyr

        !rnilyr = real(nilyr,dbl_kind)
        z2a = c0
        z2b = c0
        if (h2 > puny) then
        ! loop over new grid cells
        do k2 = 1, nbiolyr

           ! initialize new tracer
           trc2(k2) = c0

           ! calculate upper and lower boundary of new cell
           z2a = z2b  !((k2 - 1) * h2) * zspace(k2)+z2b ! / rnilyr
           z2b = z2b + h2 * zspace(k2) !(k2       * h2) * zspace(k2)+z2a !/ rnilyr

           z1a = c0
           z1b = c0
           ! loop over old grid cells
           do k1 = 1, nbiolyr

              ! calculate upper and lower boundary of old cell
              z1a = z1b !((k1 - 1) * h1) * zspace(k1)+z1b !/ rnilyr
              z1b = z1b + h1 * zspace(k1) !(k1       * h1) * zspace(k1)+z1a !/ rnilyr

              ! calculate overlap between old and new cell
              overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)

              ! aggregate old grid cell contribution to new cell
              trc2(k2) = trc2(k2) + overlap * trc(k1)

           enddo ! k1

           ! calculate upper and lower boundary of added new ice at bottom
           z1a = h1
           z1b = h2

           ! calculate overlap between added ice and new cell
           overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
           ! aggregate added ice contribution to new cell
           trc2(k2) = trc2(k2) + overlap * trc0
           ! renormalize new grid cell
           trc2(k2) = trc2(k2)/zspace(k2)/h2 !(rnilyr * trc2(k2)) / h2

        enddo ! k2
     else
        trc2 = trc
     endif
        ! update vertical tracer array with the adjusted tracer
        trc = trc2

      end subroutine update_vertical_bio_tracers

!=======================================================================

      end module icepack_zbgc

!=======================================================================
