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

      use icepack_tracers, only: nilyr, nslyr, nblyr
      use icepack_tracers, only: max_algae, max_dic, max_doc, max_don, max_fe
      use icepack_tracers, only: max_aero, max_nbtrcr
      use icepack_tracers, only: nt_sice, bio_index, bio_index_o
      use icepack_tracers, only: tr_brine, nt_fbri, nt_qice, nt_Tsfc
      use icepack_tracers, only: tr_zaero, tr_bgc_Nit, tr_bgc_N
      use icepack_tracers, only: tr_bgc_DON, tr_bgc_C, tr_bgc_chl
      use icepack_tracers, only: tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS
      use icepack_tracers, only: tr_bgc_Fe, tr_bgc_hum, tr_bgc_PON
      use icepack_tracers, only: ntrcr, ntrcr_o, nt_bgc_Nit, nlt_bgc_Nit
      use icepack_tracers, only: nt_bgc_N, nlt_bgc_N, nt_bgc_Am, nlt_bgc_Am
      use icepack_tracers, only: nt_bgc_DMSPp, nlt_bgc_DMSPp, nt_bgc_Sil, nlt_bgc_Sil
      use icepack_tracers, only: nt_bgc_DMSPd, nlt_bgc_DMSPd, nt_bgc_DMS, nlt_bgc_DMS
      use icepack_tracers, only: nt_bgc_hum, nlt_bgc_hum, nt_bgc_PON, nlt_bgc_PON
      use icepack_tracers, only: nt_bgc_C, nlt_bgc_C, nt_bgc_chl, nlt_bgc_chl
      use icepack_tracers, only: nt_bgc_DOC, nlt_bgc_DOC, nt_bgc_DON, nlt_bgc_DON
      use icepack_tracers, only: nt_bgc_DIC, nlt_bgc_DIC, nt_bgc_Fed, nlt_bgc_Fed
      use icepack_tracers, only: nt_zaero, nlt_zaero, nt_bgc_Fep, nlt_bgc_Fep
      use icepack_tracers, only: nlt_zaero_sw, nlt_chl_sw, nt_zbgc_frac
      use icepack_tracers, only: ntrcr, ntrcr_o, nbtrcr, nbtrcr_sw
      use icepack_tracers, only: n_algae, n_doc, n_dic, n_don, n_fed, n_fep, n_zaero

      use icepack_zbgc_shared, only: zbgc_init_frac
      use icepack_zbgc_shared, only: zbgc_frac_init
      use icepack_zbgc_shared, only: bgc_tracer_type, remap_zbgc
      use icepack_zbgc_shared, only: regrid_stationary
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
                icepack_init_zbgc_tracer_indices, &
                icepack_biogeochemistry, &
                icepack_load_ocean_bio_array, &
                icepack_init_ocean_bio

!=======================================================================

      contains

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (dt,         nblyr,      ncats,    &
                                  ncat,       nilyr,                &
                                  bgrid,      cgrid,      igrid,    &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vin0new,  &
                                  ntrcr,      trcrn,      nbtrcr,   &
                                  ocean_bio,  flux_bio,   hsurp,    &
                                  d_an_tot)

      integer (kind=int_kind), intent(in) :: &
         nblyr   , & ! number of bio layers
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nbtrcr  , & ! number of biology tracers
         ntrcr   , & ! number of tracers in use
         ncats       ! 1 without floe size distribution or ncat

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate

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
                                   ncat,     nblyr,     &
                                   rsiden,   vicen_init,&
                                   trcrn,    flux_bio,  &
                                   nbltrcr)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nblyr , & ! number of bio layers
         nbltrcr   ! number of biology tracers

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

      do m = 1, nbltrcr
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
!
! Add new ice tracers to the ice bottom and adjust the vertical profile
!
! author: Nicole Jeffery, LANL

!      subroutine adjust_tracer_profile (nbtrcr, dt, ntrcr, &
!                                        aicen,      vbrin, &
!                                        vicen,      trcrn, &
!                                        vtmp,              &
!                                        vsurp,             &
!                                        nilyr,      nblyr, &
!                                        bgrid,             &
!                                        cgrid,      ocean_bio, &
!                                        igrid,      location)!!!

!      integer (kind=int_kind), intent(in) :: &
!         location          , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
!         ntrcr             , & ! number of tracers in use
!         nilyr             , & ! number of ice layers
!         nbtrcr            , & ! number of biology tracers
!         nblyr                 ! number of biology layers

!      real (kind=dbl_kind), intent(in) :: &
!         dt              ! timestep (s)!

!      real (kind=dbl_kind), intent(in) :: &
!         aicen   , & ! concentration of ice
!         vicen   , & ! volume of ice
!         vsurp   , & ! volume of new ice added to each cat
!         vtmp        ! total volume of new and old ice

!      real (kind=dbl_kind), dimension (nbtrcr), intent(in) :: &
!         ocean_bio

!      real (kind=dbl_kind), intent(in) :: &
!         vbrin       ! fbri*volume per unit area of ice  (m)

!      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
!         igrid       ! zbio grid

!      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
!         bgrid       ! zsal grid

!      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
!         cgrid       ! CICE grid

!      real (kind=dbl_kind), dimension (ntrcr), intent(inout) :: &
!         trcrn       ! ice tracers

      ! local variables

!      real (kind=dbl_kind), dimension (ntrcr+2) :: &
!         trtmp0, &      ! temporary, remapped tracers
!         trtmp          ! temporary, remapped tracers

!      integer (kind=int_kind) :: &
!         k, m

!      real (kind=dbl_kind), dimension (nblyr+1) ::  &
!         C_stationary      ! stationary bulk concentration*h (mmol/m^2)

!      real(kind=dbl_kind) :: &
!         top_conc     , & ! salinity or bgc ocean concentration of frazil
!         fluxb        , & ! needed for regrid (set to zero here)
!         hbri_old     , & ! previous timestep brine height
!         hbri             ! brine height

!      character(len=*),parameter :: subname='(adjust_tracer_profile)'

!      trtmp0(:) = c0
!      trtmp(:) = c0
!      fluxb = c0

 !     if (location == 1 .and. vbrin > c0) then  ! add frazil to bottom

!         hbri     = vbrin
!         hbri_old = vtmp

!         do m = 1, nbtrcr
!            top_conc = ocean_bio(m)*zbgc_init_frac(m)
!            do k = 1, nblyr+1
!               C_stationary(k) = trcrn(bio_index(m) + k-1)* hbri_old
!            enddo !k
!            call regrid_stationary (C_stationary, hbri_old, &
!                                    hbri,         dt,       &
!                                    ntrcr,                  &
!                                    nblyr,        top_conc, &
!                                    igrid,        fluxb )
!            if (icepack_warnings_aborted(subname)) return
!            do k = 1, nblyr+1
!               trcrn(bio_index(m) + k-1) =  C_stationary(k)/hbri
!            enddo !k
!         enddo !m

!      elseif (vbrin > c0) then   ! add frazil throughout  location == 0 .and.

!         do k = 1, nblyr+1
 !           do m = 1, nbtrcr
!               trcrn(bio_index(m) + k-1) = (trcrn(bio_index(m) + k-1) * vtmp &
!                         + ocean_bio(m)*zbgc_init_frac(m) * vsurp) / vbrin
!            enddo
!         enddo

!      endif     ! location

!      end subroutine adjust_tracer_profile

!=======================================================================
!autodocument_start icepack_init_bgc
!
      subroutine icepack_init_bgc(ncat, nblyr, nilyr, &
         cgrid, igrid, sicen, trcrn, sss, ocean_bio_all, &
         DOCPoolFractions)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr     ! number of bio layers

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         igrid     ! biology vertical interface points

      real (kind=dbl_kind), dimension (nilyr+1), intent(inout) :: &
         cgrid     ! CICE vertical coordinate

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

      subroutine icepack_init_zbgc ( &
                 R_Si2N_in, R_S2N_in, R_Fe2C_in, R_Fe2N_in, R_C2N_in, R_C2N_DON_in, &
                 R_chl2N_in, F_abs_chl_in, R_Fe2DON_in, R_Fe2DOC_in, chlabs_in, &
                 alpha2max_low_in, beta2max_in, mu_max_in, fr_graze_in, mort_pre_in, &
                 mort_Tdep_in, k_exude_in, K_Nit_in, K_Am_in, K_sil_in, K_Fe_in, &
                 f_don_in, kn_bac_in, f_don_Am_in, f_doc_in, f_exude_in, k_bac_in, &
                 grow_Tdep_in, zbgc_frac_init_in, &
                 zbgc_init_frac_in, tau_ret_in, tau_rel_in, bgc_tracer_type_in, &
                 fr_resp_in, algal_vel_in, R_dFe2dust_in, dustFe_sol_in, T_max_in, &
                 op_dep_min_in, fr_graze_s_in, fr_graze_e_in, fr_mort2min_in, fr_dFe_in, &
                 k_nitrif_in, t_iron_conv_in, max_loss_in, max_dfe_doc1_in, &
                 fr_resp_s_in, y_sk_DMS_in, t_sk_conv_in, t_sk_ox_in, fsal_in)

      real (kind=dbl_kind), optional :: R_C2N_in(:)        ! algal C to N (mole/mole)
      real (kind=dbl_kind), optional :: R_chl2N_in(:)      ! 3 algal chlorophyll to N (mg/mmol)
      real (kind=dbl_kind), optional :: F_abs_chl_in(:)    ! to scale absorption in Dedd
      real (kind=dbl_kind), optional :: R_C2N_DON_in(:)    ! increase compare to algal R_Fe2C
      real (kind=dbl_kind), optional :: R_Si2N_in(:)       ! algal Sil to N (mole/mole)
      real (kind=dbl_kind), optional :: R_S2N_in(:)        ! algal S to N (mole/mole)
      real (kind=dbl_kind), optional :: R_Fe2C_in(:)       ! algal Fe to carbon (umol/mmol)
      real (kind=dbl_kind), optional :: R_Fe2N_in(:)       ! algal Fe to N (umol/mmol)
      real (kind=dbl_kind), optional :: R_Fe2DON_in(:)     ! Fe to N of DON (nmol/umol)
      real (kind=dbl_kind), optional :: R_Fe2DOC_in(:)     ! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), optional :: fr_resp_in         ! frac of algal growth lost due to respiration
      real (kind=dbl_kind), optional :: algal_vel_in       ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
      real (kind=dbl_kind), optional :: R_dFe2dust_in      ! g/g (3.5% content) Tagliabue 2009
      real (kind=dbl_kind), optional :: dustFe_sol_in      ! solubility fraction
      real (kind=dbl_kind), optional :: T_max_in           ! maximum temperature (C)
      real (kind=dbl_kind), optional :: op_dep_min_in      ! Light attenuates for optical depths exceeding min
      real (kind=dbl_kind), optional :: fr_graze_s_in      ! fraction of grazing spilled or slopped
      real (kind=dbl_kind), optional :: fr_graze_e_in      ! fraction of assimilation excreted
      real (kind=dbl_kind), optional :: fr_mort2min_in     ! fractionation of mortality to Am
      real (kind=dbl_kind), optional :: fr_dFe_in          ! fraction of remineralized nitrogen
                                                           ! (in units of algal iron)
      real (kind=dbl_kind), optional :: k_nitrif_in        ! nitrification rate (1/day)
      real (kind=dbl_kind), optional :: t_iron_conv_in     ! desorption loss pFe to dFe (day)
      real (kind=dbl_kind), optional :: max_loss_in        ! restrict uptake to % of remaining value
      real (kind=dbl_kind), optional :: max_dfe_doc1_in    ! max ratio of dFe to saccharides in the ice (nM Fe/muM C)
      real (kind=dbl_kind), optional :: fr_resp_s_in       ! DMSPd fraction of respiration loss as DMSPd
      real (kind=dbl_kind), optional :: y_sk_DMS_in        ! fraction conversion given high yield
      real (kind=dbl_kind), optional :: t_sk_conv_in       ! Stefels conversion time (d)
      real (kind=dbl_kind), optional :: t_sk_ox_in         ! DMS oxidation time (d)
      real (kind=dbl_kind), optional :: fsal_in            ! salinity limitation factor (1)

      real (kind=dbl_kind), optional :: chlabs_in(:)       ! chla absorption 1/m/(mg/m^3)
      real (kind=dbl_kind), optional :: alpha2max_low_in(:)  ! light limitation (1/(W/m^2))
      real (kind=dbl_kind), optional :: beta2max_in(:)     ! light inhibition (1/(W/m^2))
      real (kind=dbl_kind), optional :: mu_max_in(:)       ! maximum growth rate (1/d)
      real (kind=dbl_kind), optional :: grow_Tdep_in(:)    ! T dependence of growth (1/C)
      real (kind=dbl_kind), optional :: fr_graze_in(:)     ! fraction of algae grazed
      real (kind=dbl_kind), optional :: mort_pre_in(:)     ! mortality (1/day)
      real (kind=dbl_kind), optional :: mort_Tdep_in(:)    ! T dependence of mortality (1/C)
      real (kind=dbl_kind), optional :: k_exude_in(:)      ! algal carbon  exudation rate (1/d)
      real (kind=dbl_kind), optional :: K_Nit_in(:)        ! nitrate half saturation (mmol/m^3)
      real (kind=dbl_kind), optional :: K_Am_in(:)         ! ammonium half saturation (mmol/m^3)
      real (kind=dbl_kind), optional :: K_Sil_in(:)        ! silicon half saturation (mmol/m^3)
      real (kind=dbl_kind), optional :: K_Fe_in(:)         ! iron half saturation  or micromol/m^3
      real (kind=dbl_kind), optional :: f_don_in(:)        ! fraction of spilled grazing to DON
      real (kind=dbl_kind), optional :: kn_bac_in(:)       ! Bacterial degredation of DON (1/d)
      real (kind=dbl_kind), optional :: f_don_Am_in(:)     ! fraction of remineralized DON to Am
      real (kind=dbl_kind), optional :: f_doc_in(:)        ! fraction of mort_N that goes to each doc pool
      real (kind=dbl_kind), optional :: f_exude_in(:)      ! fraction of exuded carbon to each DOC pool
      real (kind=dbl_kind), optional :: k_bac_in(:)        ! Bacterial degredation of DOC (1/d)

      real (kind=dbl_kind), optional :: zbgc_frac_init_in(:)  ! initializes mobile fraction
      real (kind=dbl_kind), optional :: bgc_tracer_type_in(:) ! described tracer in mobile or stationary phases
      real (kind=dbl_kind), optional :: zbgc_init_frac_in(:)  ! fraction of ocean tracer  concentration in new ice
      real (kind=dbl_kind), optional :: tau_ret_in(:)         ! retention timescale  (s), mobile to stationary phase
      real (kind=dbl_kind), optional :: tau_rel_in(:)         ! release timescale    (s), stationary to mobile phase

!autodocument_end

      character(len=*),parameter :: subname='(icepack_init_zbgc)'

      !--------

      if (present(R_C2N_in)) then
         R_C2N(:)     = R_C2N_in(:)
      else
         R_C2N(1) = ratio_C2N_diatoms
         R_C2N(2) = ratio_C2N_sp
         R_C2N(3) = ratio_C2N_phaeo
      endif
      if (present(R_chl2N_in)) then
         R_chl2N(:)   = R_chl2N_in(:)
      else
         R_chl2N(1) = ratio_chl2N_diatoms
         R_chl2N(2) = ratio_chl2N_sp
         R_chl2N(3) = ratio_chl2N_phaeo
      endif
      if (present(F_abs_chl_in)) then
         F_abs_chl(:) = F_abs_chl_in(:)
      else
         F_abs_chl(1) = F_abs_chl_diatoms
         F_abs_chl(2) = F_abs_chl_sp
         F_abs_chl(3) = F_abs_chl_phaeo
      endif
      if (present(R_C2N_DON_in)) then
         R_C2N_DON(:) = R_C2N_DON_in(:)
      else
         R_C2N_DON(1) = ratio_C2N_proteins
      endif
      if (present(R_Si2N_in)) then
         R_Si2N(:)    = R_Si2N_in(:)
      else
         R_Si2N(1) = ratio_Si2N_diatoms
         R_Si2N(2) = ratio_Si2N_sp
         R_Si2N(3) = ratio_Si2N_phaeo
      endif
      if (present(R_S2N_in)) then
         R_S2N(:)     = R_S2N_in(:)
      else
         R_S2N(1) = ratio_S2N_diatoms
         R_S2N(2) = ratio_S2N_sp
         R_S2N(3) = ratio_S2N_phaeo
      endif
      if (present(R_Fe2C_in)) then
         R_Fe2C(:)    = R_Fe2C_in(:)
      else
         R_Fe2C(1) = ratio_Fe2C_diatoms
         R_Fe2C(2) = ratio_Fe2C_sp
         R_Fe2C(3) = ratio_Fe2C_phaeo
      endif
      if (present(R_Fe2N_in)) then
         R_Fe2N(:)    = R_Fe2N_in(:)
      else
         R_Fe2N(1) = ratio_Fe2N_diatoms
         R_Fe2N(2) = ratio_Fe2N_sp
         R_Fe2N(3) = ratio_Fe2N_phaeo
      endif
      if (present(R_Fe2DON_in)) then
         R_Fe2DON(:)  = R_Fe2DON_in(:)
      else
         R_Fe2DON(1) = ratio_Fe2DON
      endif

      if (present(R_Fe2DOC_in)) then
         R_Fe2DOC(:)  = R_Fe2DOC_in(:)
      else
         R_Fe2DOC(1) = ratio_Fe2DOC_s
         R_Fe2DOC(2) = ratio_Fe2DOC_l
         R_Fe2DOC(3) = c0
      endif

      if (present(chlabs_in)) then
         chlabs(:)    = chlabs_in(:)
      else
         chlabs(1) = chlabs_diatoms
         chlabs(2) = chlabs_sp
         chlabs(3) = chlabs_phaeo
      endif

      if (present(alpha2max_low_in)) then
         alpha2max_low(:) = alpha2max_low_in(:)
      else
         alpha2max_low(1) = alpha2max_low_diatoms
         alpha2max_low(2) = alpha2max_low_sp
         alpha2max_low(3) = alpha2max_low_phaeo
      endif
      if (present(beta2max_in)) then
         beta2max(:)  = beta2max_in(:)
      else
         beta2max(1) = beta2max_diatoms
         beta2max(2) = beta2max_sp
         beta2max(3) = beta2max_phaeo
      endif
      if (present(mu_max_in)) then
         mu_max(:)    = mu_max_in(:)
      else
         mu_max(1) = mu_max_diatoms
         mu_max(2) = mu_max_sp
         mu_max(3) = mu_max_phaeo
      endif
      if (present(grow_Tdep_in)) then
         grow_Tdep(:) = grow_Tdep_in(:)
      else
         grow_Tdep(1) = grow_Tdep_diatoms
         grow_Tdep(2) = grow_Tdep_sp
         grow_Tdep(3) = grow_Tdep_phaeo
      endif
      if (present(fr_graze_in)) then
         fr_graze(:)  = fr_graze_in(:)
      else
         fr_graze(1) = fr_graze_diatoms
         fr_graze(2) = fr_graze_sp
         fr_graze(3) = fr_graze_phaeo
      endif
      if (present(mort_pre_in)) then
         mort_pre(:)  = mort_pre_in(:)
      else
         mort_pre(1) = mort_pre_diatoms
         mort_pre(2) = mort_pre_sp
         mort_pre(3) = mort_pre_phaeo
      endif
      if (present(mort_Tdep_in)) then
         mort_Tdep(:) = mort_Tdep_in(:)
      else
         mort_Tdep(1) = mort_Tdep_diatoms
         mort_Tdep(2) = mort_Tdep_sp
         mort_Tdep(3) = mort_Tdep_phaeo
      endif
      if (present(k_exude_in)) then
         k_exude(:)   = k_exude_in(:)
      else
         k_exude(1) = k_exude_diatoms
         k_exude(2) = k_exude_sp
         k_exude(3) = k_exude_phaeo
      endif
      if (present(K_Nit_in)) then
         K_Nit(:)     = K_Nit_in(:)
      else
         K_Nit(1) = K_Nit_diatoms
         K_Nit(2) = K_Nit_sp
         K_Nit(3) = K_Nit_phaeo
      endif
      if (present(K_Am_in)) then
         K_Am(:)      = K_Am_in(:)
      else
         K_Am(1) = K_Am_diatoms
         K_Am(2) = K_Am_sp
         K_Am(3) = K_Am_phaeo
      endif
      if (present(K_Sil_in)) then
         K_Sil(:)     = K_Sil_in(:)
      else
         K_Sil(1) = K_Sil_diatoms
         K_Sil(2) = K_Sil_sp
         K_Sil(3) = K_Sil_phaeo
      endif
      if (present(K_Fe_in)) then
         K_Fe(:)      = K_Fe_in(:)
      else
         K_Fe(1) = K_Fe_diatoms
         K_Fe(2) = K_Fe_sp
         K_Fe(3) = K_Fe_phaeo
      endif
      if (present(f_don_in)) then
         f_don(:)     = f_don_in(:)
      else
         f_don(1) = f_don_protein
      endif
      if (present(kn_bac_in)) then
         kn_bac(:)    = kn_bac_in(:)
      else
         kn_bac(1) = kn_bac_protein
      endif
      if (present(f_don_Am_in)) then
         f_don_Am(:)  = f_don_Am_in(:)
      else
         f_don_Am(1) = f_don_Am_protein
      endif
      if (present(f_doc_in)) then
         f_doc(:)     = f_doc_in(:)
      else
         f_doc(1) = f_doc_s
         f_doc(2) = f_doc_l
      endif
      if (present(f_exude_in)) then
         f_exude(:)   = f_exude_in(:)
      else
         f_exude(1) = f_exude_s
         f_exude(2) = f_exude_l
      endif
      if (present(k_bac_in)) then
         k_bac(:)     = k_bac_in(:)
      else
         k_bac(1) = k_bac_s
         k_bac(2) = k_bac_l
      endif

      ! should already be defined with tracer indices
      if (present(zbgc_frac_init_in))  zbgc_frac_init(:)  = zbgc_frac_init_in(:)
      if (present(bgc_tracer_type_in)) bgc_tracer_type(:) = bgc_tracer_type_in(:)
      if (present(zbgc_init_frac_in))  zbgc_init_frac(:)  =  zbgc_init_frac_in(:)
      if (present(tau_ret_in)) tau_ret(:) = tau_ret_in(:)
      if (present(tau_rel_in)) tau_rel(:) = tau_rel_in(:)

      ! defined in init_parameters (not needed here)
      if (present(fr_resp_in))      fr_resp      = fr_resp_in
      if (present(algal_vel_in))    algal_vel    = algal_vel_in
      if (present(R_dFe2dust_in))   R_dFe2dust   = R_dFe2dust_in
      if (present(dustFe_sol_in))   dustFe_sol   = dustFe_sol_in
      if (present(T_max_in))        T_max        = T_max_in
      if (present(op_dep_min_in))   op_dep_min   = op_dep_min_in
      if (present(fr_graze_s_in))   fr_graze_s   = fr_graze_s_in
      if (present(fr_graze_e_in))   fr_graze_e   = fr_graze_e_in
      if (present(fr_mort2min_in))  fr_mort2min  = fr_mort2min_in
      if (present(fr_dFe_in))       fr_dFe       = fr_dFe_in
      if (present(k_nitrif_in))     k_nitrif     = k_nitrif_in
      if (present(t_iron_conv_in))  t_iron_conv  = t_iron_conv_in
      if (present(max_loss_in))     max_loss     = max_loss_in
      if (present(max_dfe_doc1_in)) max_dfe_doc1 = max_dfe_doc1_in
      if (present(fr_resp_s_in))    fr_resp_s    = fr_resp_s_in
      if (present(y_sk_DMS_in))     y_sk_DMS     = y_sk_DMS_in
      if (present(t_sk_conv_in))    t_sk_conv    = t_sk_conv_in
      if (present(t_sk_ox_in))      t_sk_ox      = t_sk_ox_in
      if (present(fsal_in))         fsal         = fsal_in


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
                           bgrid, igrid, icgrid, cgrid,  &
                           nblyr, nilyr, nslyr, ncat, &
                           meltbn, melttn, congeln, snoicen, &
                           sst, sss, Tf, fsnow, meltsn, & !hmix, &
                           hin_old, flux_bio, flux_bio_atm, &
                           aicen_init, vicen_init, aicen, vicen, vsnon, &
                           aice0, trcrn, vsnon_init, &
                           flux_bion, bioPorosityIceCell, &
                           bioSalinityIceCell, bioTemperatureIceCell)

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat, &
         nilyr, &
         nslyr, &
         nblyr

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         bgrid         , &  ! biology nondimensional vertical grid points
         igrid         , &  ! biology vertical interface points
         cgrid         , &  ! CICE vertical coordinate
         icgrid        , &  ! interface grid for CICE (shortwave variable)
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
            !if (z_tracers) then
               do mm = 1,nbtrcr
                  trcrn(nt_zbgc_frac-1+mm,n) = zbgc_frac_init(mm)
               enddo
            !endif
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

               call compute_microS_mushy (nilyr,         nblyr,       &
                           bgrid,         cgrid,         igrid,       &
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

               call zbio (dt,                    nblyr,                  &
                          nslyr,                 nilyr,                  &
                          melttn(n),                                     &
                          meltsn(n),             meltbn  (n),            &
                          congeln(n),            snoicen(n),             &
                          nbtrcr,                fsnow,                  &
                          ntrcr,                 trcrn(1:ntrcr,n),       &
                          bio_index(1:nbtrcr),   bio_index_o(:),         &
                          aicen_init(n),                                 &
                          vicen_init(n),         vsnon_init(n),          &
                          vicen(n),              vsnon(n),               &
                          aicen(n),              flux_bio_atm(1:nbtrcr), &
                          n,                     n_algae,                &
                          n_doc,                 n_dic,                  &
                          n_don,                                         &
                          n_fed,                 n_fep,                  &
                          n_zaero,               first_ice(n),           &
                          hin_old(n),            ocean_bio(1:nbtrcr),    &
                          ocean_bio_dh,                                  &
                          bphi(:,n),             iphin,                  &
                          iDi(:,n),                                      &
                          fswpenln(:,n),                                 &
                          dhbr_top(n),           dhbr_bot(n),            &
                          zfswin(:,n),                                   &
                          hbrin,                 hbr_old,                &
                          darcy_V(n),                                    &
                          bgrid,                                         &
                          igrid,                 icgrid,                 &
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
                          bioPorosityIceCell(:), bioSalinityIceCell(:),  &
                          bioTemperatureIceCell(:)                       )

               if (icepack_warnings_aborted(subname)) return

               if (present(flux_bion)) then
                  do mm = 1, nbtrcr
                     flux_bion(mm,n) = flux_bion_n(mm)
                  enddo
               endif

            elseif (skl_bgc) then

               call sklbio (dt,                      Tf,                  &
                            ntrcr,                                        &
                            nbtrcr,                  n_algae,             &
                            n_doc,               &
                            n_dic,                   n_don,               &
                            n_fed,                   n_fep,               &
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

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
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

      real (kind=dbl_kind), intent(out):: &
       amm      , & ! ammonium
       dmsp     , & ! DMSPp
       dms      , & ! DMS
       hum      , & ! humic material
       nit      , & ! nitrate
       sil          ! silicate

      real (kind=dbl_kind), dimension(:), intent(out):: &
       algalN   , & ! algae
       doc      , & ! DOC
       dic      , & ! DIC
       don      , & ! DON
       fed      , & ! Dissolved Iron
       fep      , & ! Particulate Iron
       zaeros       ! BC and dust

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
       CToN     , & ! carbon to nitrogen ratio for algae
       CToN_DON     ! nitrogen to carbon ratio for proteins

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
        k

      character(len=*),parameter :: subname='(icepack_init_ocean_bio)'

       if (present(CToN)) then
         CToN(1) = R_C2N(1)
         CToN(2) = R_C2N(2)
         CToN(3) = R_C2N(3)
       endif

       if (present(CToN_DON)) then
         CToN_DON(1) = R_C2N_DON(1)
       endif

       amm  = c1 ! ISPOL < 1 mmol/m^3
       dmsp = p1
       dms  = p1
       algalN(1) = c1  !0.0026_dbl_kind ! ISPOL, Lannuzel 2013(pennate)
       algalN(2) = 0.0057_dbl_kind ! ISPOL, Lannuzel 2013(small plankton)
       algalN(3) = 0.0027_dbl_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
                                     ! 0.024_dbl_kind ! 5% of 1 mgchl/m^3
       doc(1) = 16.2_dbl_kind ! 18% saccharides
       doc(2) = 9.0_dbl_kind  ! lipids
       doc(3) = c1 !
       do k = 1, max_dic
            dic(k) = 1950.0_dbl_kind ! 1950-2260 mmol C/m3 (Tynan et al. 2015)
       enddo
       do k = 1, max_don
            don(k) = 12.9_dbl_kind
            ! 64.3_dbl_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
       enddo
       !ki = 1
       !if (trim(fe_data_type) == 'clim') ki = 2
       do k = 1, max_fe ! ki, max_fe
            fed(k) = 0.4_dbl_kind ! c1 (nM) Lannuzel2007 DFe,
                                  ! range 0.14-2.6 (nM) van der Merwe 2011
                                  ! Tagliabue 2012 (0.4 nM)
            fep(k) = c2 ! (nM) van der Merwe 2011
                        ! (0.6 to 2.9 nM ocean)
       enddo
       hum  = c1        ! mmol C/m^3
       nit  = 12.0_dbl_kind
       sil  = 25.0_dbl_kind
       do k = 1, max_aero
         zaeros(k) = c0
       enddo


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
!
!  Initialize zbgc tracer indices if not done in the driver

      subroutine icepack_init_zbgc_tracer_indices(&
                 nilyr_in, nslyr_in, nblyr_in, &
                 n_algae_in, n_zaero_in, n_doc_in, n_dic_in, n_don_in, n_fed_in, n_fep_in, &
                 trcr_base, trcr_depend, n_trcr_strata, nt_strata, nbtrcr_sw_out, tr_brine_in,&
                 nt_fbri_out, ntrcr_out,ntrcr_o_out, nbtrcr_out, nt_bgc_Nit_out, nt_bgc_Am_out, &
                 nt_bgc_Sil_out, nt_bgc_DMS_out, nt_bgc_PON_out, nt_bgc_N_out, &
                 nt_bgc_C_out, nt_bgc_chl_out, nt_bgc_DOC_out, nt_bgc_DON_out, nt_bgc_DIC_out, &
                 nt_zaero_out, nt_bgc_DMSPp_out, nt_bgc_DMSPd_out, nt_bgc_Fed_out, nt_bgc_Fep_out, &
                 nt_zbgc_frac_out, tr_bgc_Nit_in, tr_bgc_Am_in, tr_bgc_Sil_in, tr_bgc_DMS_in, &
                 tr_bgc_PON_in, tr_bgc_N_in, tr_bgc_C_in, tr_bgc_chl_in, &
                 tr_bgc_DON_in, tr_bgc_Fe_in, tr_zaero_in, nlt_zaero_sw_out, nlt_chl_sw_out, &
                 nlt_bgc_N_out, nlt_bgc_Nit_out, nlt_bgc_Am_out, nlt_bgc_Sil_out, &
                 nlt_bgc_DMS_out, nlt_bgc_DMSPp_out, nlt_bgc_DMSPd_out, &
                 nlt_bgc_C_out, nlt_bgc_chl_out, nlt_bgc_DIC_out, nlt_bgc_DOC_out, &
                 nlt_bgc_PON_out, nlt_bgc_DON_out, nlt_bgc_Fed_out, nlt_bgc_Fep_out, &
                 nlt_zaero_out, &
                 nt_bgc_hum_out, nlt_bgc_hum_out, tr_bgc_hum_in, &
                 skl_bgc_in, z_tracers_in, dEdd_algae_in, solve_zbgc_in, &
                 bio_index_o_out, bio_index_out, printdiags)

      integer (kind=int_kind), optional, intent(in) ::&
         nilyr_in     , & ! number of ice layers
         nslyr_in     , & ! number of snow layers
         nblyr_in     , & ! number of biology layers
         n_zaero_in   , & ! number of z aerosols in use
         n_algae_in   , & ! number of algae in use
         n_doc_in     , & ! number of DOC pools in use
         n_dic_in     , & ! number of DIC pools in use
         n_don_in     , & ! number of DON pools in use
         n_fed_in     , & ! number of Fe  pools in use dissolved Fe
         n_fep_in         ! number of Fe  pools in use particulate Fe

      integer (kind=int_kind), optional, intent(inout) :: &
         ntrcr_out,       & ! number of tracers
         ntrcr_o_out,     & ! number of non-bio tracers in use
         nbtrcr_out,      & ! number of bgc tracers in use
         nbtrcr_sw_out      ! size of shorwave tracer vector

      integer (kind=int_kind), dimension (:), intent(inout) :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), dimension (:), intent(inout) :: &
         n_trcr_strata ! number of underlying tracer layers

      integer (kind=int_kind), dimension (:,:), intent(inout) :: &
         nt_strata     ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcr_base     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno

      logical (kind=log_kind), optional, intent(in) :: &
         tr_brine_in,       & ! if .true., brine height differs from ice thickness
         tr_zaero_in,       & ! if .true., black carbon is tracers  (n_zaero)
         tr_bgc_Nit_in,     & ! if .true. Nitrate tracer in ice
         tr_bgc_N_in,       & ! if .true., algal nitrogen tracers  (n_algae)
         tr_bgc_DON_in,     & ! if .true., DON pools are tracers  (n_don)
         tr_bgc_C_in,       & ! if .true., algal carbon tracers + DOC and DIC
         tr_bgc_chl_in,     & ! if .true., algal chlorophyll tracers
         tr_bgc_Am_in,      & ! if .true., ammonia/um as nutrient tracer
         tr_bgc_Sil_in,     & ! if .true., silicon as nutrient tracer
         tr_bgc_DMS_in,     & ! if .true., DMS as  tracer
         tr_bgc_Fe_in,      & ! if .true., Fe as  tracer
         tr_bgc_PON_in,     & ! if .true., PON as tracer
         tr_bgc_hum_in,     & ! if .true., humic material as tracer
         z_tracers_in,      & ! if .true., bgc or aerosol tracers are vertically resolved
         solve_zbgc_in,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_in,     & ! if .true., algal absorption of Shortwave is computed in the
         skl_bgc_in           ! if true, solve skeletal biochemistry

       integer (kind=int_kind), optional, intent(out) :: &
         nt_fbri_out,      & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_bgc_Nit_out,   & ! nutrients
         nt_bgc_Am_out,    & !
         nt_bgc_Sil_out,   & !
         nt_bgc_DMSPp_out, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd_out, & !
         nt_bgc_DMS_out,   & !
         nt_bgc_PON_out,   & ! zooplankton and detritus
         nt_bgc_hum_out,   & ! humic material
                         ! bio layer indicess
         nlt_bgc_Nit_out,  & ! nutrients
         nlt_bgc_Am_out,   & !
         nlt_bgc_Sil_out,  & !
         nlt_bgc_DMSPp_out,& ! trace gases (skeletal layer)
         nlt_bgc_DMSPd_out,& !
         nlt_bgc_DMS_out,  & !
         nlt_bgc_PON_out,  & ! zooplankton and detritus
         nlt_bgc_hum_out,  & ! humic material
         nlt_chl_sw_out,   & ! points to total chla in trcrn_sw
         nt_zbgc_frac_out    ! fraction of tracer in the mobile phase

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         nt_bgc_N_out , & ! diatoms, phaeocystis, pico/small
         nt_bgc_C_out , & ! diatoms, phaeocystis, pico/small
         nt_bgc_chl_out,& ! diatoms, phaeocystis, pico/small
         nlt_bgc_N_out ,& ! diatoms, phaeocystis, pico/small
         nlt_bgc_C_out ,& ! diatoms, phaeocystis, pico/small
         nlt_bgc_chl_out   ! diatoms, phaeocystis, pico/small

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         nt_bgc_DOC_out,   & !  dissolved organic carbon
         nlt_bgc_DOC_out     !  dissolved organic carbon

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         nt_bgc_DON_out,   & !  dissolved organic nitrogen
         nlt_bgc_DON_out     !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         nt_bgc_DIC_out,    & !  dissolved inorganic carbon
         nlt_bgc_DIC_out      !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         nt_bgc_Fed_out,     & !  dissolved iron
         nt_bgc_Fep_out,     & !  particulate iron
         nlt_bgc_Fed_out,    & !  dissolved iron
         nlt_bgc_Fep_out       !  particulate iron

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         nt_zaero_out,    & !  black carbon and other aerosols
         nlt_zaero_out,   & !  black carbon and other aerosols
         nlt_zaero_sw_out

      integer (kind=int_kind), dimension(:), optional, intent(out) :: &
         bio_index_o_out , & ! nlt  to appropriate value in ocean data array
         bio_index_out       ! nlt to nt

      logical (kind=log_kind), optional, intent(in) :: &
         printdiags    ! print diagnostics to warning package

!autodocument_end

      ! local variables

      real (kind=dbl_kind), dimension (:), allocatable :: &
        algaltype, &
        doctype, &
        dictype, &
        dontype, &
        fedtype, &
        feptype, &
        zaerotype

      integer (kind=int_kind) :: &
        k, mm    , & ! loop index
        ntd      , & ! for tracer dependency calculation
        nk       , & !
        nt_depend

      character(len=*),parameter :: subname='(icepack_init_zbgc_tracer_indices)'

      !------------

      ! BGC Tracer Dimensions
      if (present(nilyr_in)    ) nilyr     = nilyr_in
      if (present(nslyr_in)    ) nslyr     = nslyr_in
      if (present(nblyr_in)    ) nblyr     = nblyr_in
      if (present(n_algae_in)  ) n_algae   = n_algae_in
      if (present(n_zaero_in)  ) n_zaero   =  n_zaero_in
      if (present(n_doc_in)    ) n_doc     = n_doc_in
      if (present(n_dic_in)    ) n_dic     = n_dic_in
      if (present(n_don_in)    ) n_don     = n_don_in
      if (present(n_fed_in)    ) n_fed     = n_fed_in
      if (present(n_fep_in)    ) n_fep     = n_fep_in

      ! BGC Tracer Flags
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

      !BGC Config Flags
      if (present(z_tracers_in)         ) z_tracers        = z_tracers_in
      if (present(solve_zbgc_in)        ) solve_zbgc       = solve_zbgc_in
      if (present(dEdd_algae_in)        ) dEdd_algae       = dEdd_algae_in
      if (present(skl_bgc_in)           ) skl_bgc          = skl_bgc_in

      !  Total tracer counts.  Initialize first
      if (present(ntrcr_out)                ) ntrcr            = ntrcr_out

      ntrcr_o = ntrcr
      nt_fbri = 0
      if (tr_brine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
          trcr_depend(nt_fbri)   = 1   ! volume-weighted
          trcr_base  (nt_fbri,1) = c0  ! volume-weighted
          trcr_base  (nt_fbri,2) = c1  ! volume-weighted
          trcr_base  (nt_fbri,3) = c0  ! volume-weighted
          n_trcr_strata(nt_fbri) = 0
          nt_strata  (nt_fbri,1) = 0
          nt_strata  (nt_fbri,2) = 0
      endif

      ntd = 0                    ! if nt_fbri /= 0 then use fbri dependency
      if (nt_fbri == 0) ntd = -1 ! otherwise make tracers depend on ice volume

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      nbtrcr = 0
      nbtrcr_sw = 0

      ! vectors of size max_algae
      nlt_bgc_N(:) = 0
      nlt_bgc_C(:) = 0
      nlt_bgc_chl(:) = 0
      nt_bgc_N(:) = 0
      nt_bgc_C(:) = 0
      nt_bgc_chl(:) = 0

      ! vectors of size max_dic
      nlt_bgc_DIC(:) = 0
      nt_bgc_DIC(:) = 0

      ! vectors of size max_doc
      nlt_bgc_DOC(:) = 0
      nt_bgc_DOC(:) = 0

      ! vectors of size max_don
      nlt_bgc_DON(:) = 0
      nt_bgc_DON(:) = 0

      ! vectors of size max_fe
      nlt_bgc_Fed(:) = 0
      nlt_bgc_Fep(:) = 0
      nt_bgc_Fed(:) = 0
      nt_bgc_Fep(:) = 0

      ! vectors of size max_aero
      nlt_zaero(:) = 0
      nlt_zaero_sw(:) = 0
      nt_zaero(:) = 0

      nlt_bgc_Nit    = 0
      nlt_bgc_Am     = 0
      nlt_bgc_Sil    = 0
      nlt_bgc_DMSPp  = 0
      nlt_bgc_DMSPd  = 0
      nlt_bgc_DMS    = 0
      nlt_bgc_PON    = 0
      nlt_bgc_hum    = 0
      nlt_chl_sw     = 0
      bio_index(:)   = 0
      bio_index_o(:) = 0

      nt_bgc_Nit    = 0
      nt_bgc_Am     = 0
      nt_bgc_Sil    = 0
      nt_bgc_DMSPp  = 0
      nt_bgc_DMSPd  = 0
      nt_bgc_DMS    = 0
      nt_bgc_PON    = 0
      nt_bgc_hum    = 0

      allocate(algaltype(max_algae))
      allocate(doctype(max_doc))
      allocate(dictype(max_dic))
      allocate(dontype(max_don))
      allocate(fedtype(max_fe))
      allocate(feptype(max_fe))
      allocate(zaerotype(max_aero))

      algaltype(1) = algaltype_diatoms
      algaltype(2) = algaltype_sp
      algaltype(3) = algaltype_phaeo

      doctype(1) = doctype_s
      doctype(2) = doctype_l

      dictype(1) = dictype_1

      dontype(1) = dontype_protein

      fedtype(1) = fedtype_1
      feptype(1) = feptype_1

      zaerotype(1) = zaerotype_bc1
      zaerotype(2) = zaerotype_bc2
      zaerotype(3) = zaerotype_dust1
      zaerotype(4) = zaerotype_dust2
      zaerotype(5) = zaerotype_dust3
      zaerotype(6) = zaerotype_dust4

      if (skl_bgc) then

         nk = 1
         nt_depend = 0

         if (dEdd_algae) then
           nlt_chl_sw = 1
           nbtrcr_sw = nilyr+nslyr+2  ! only the bottom layer
                                      ! will be nonzero
          endif
      elseif (z_tracers) then ! defined on nblyr+1 in ice
                              ! and 2 snow layers (snow surface + interior)
         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd

         if (tr_bgc_N) then
            if (dEdd_algae) then
               nlt_chl_sw = 1
               nbtrcr_sw =  nilyr+nslyr+2
            endif
         endif ! tr_bgc_N

      endif ! skl_bgc or z_tracers


      if (skl_bgc .or. z_tracers) then

      !-----------------------------------------------------------------
      ! assign tracer indices and dependencies
      ! bgc_tracer_type: < 0  purely mobile , >= 0 stationary
      !------------------------------------------------------------------

      if (tr_bgc_N) then
         do mm = 1, n_algae
            call icepack_init_bgc_trcr(&
                                      nk,             nt_fbri,       &
                                      nt_bgc_N(mm),    nlt_bgc_N(mm), &
                                      algaltype(mm),   nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_N(mm)) = mm
         enddo   ! mm
         if (icepack_warnings_aborted(subname)) return
      endif ! tr_bgc_N

      if (tr_bgc_Nit) then
         call icepack_init_bgc_trcr(&
                                      nk,             nt_fbri,       &
                                      nt_bgc_Nit,      nlt_bgc_Nit,   &
                                      nitratetype,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_Nit) = max_algae + 1
      endif ! tr_bgc_Nit

      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do mm = 1,n_algae
       !     call icepack_init_bgc_trcr(nk,             nt_fbri,       &
       !                               nt_bgc_C(mm),    nlt_bgc_C(mm), &
       !                               algaltype(mm),   nt_depend,     &
       !                               ntrcr,           nbtrcr,        &
       !                               bgc_tracer_type, trcr_depend,   &
       !                               trcr_base,       n_trcr_strata, &
       !                               nt_strata,       bio_index)
       !     bio_index_o(nlt_bgc_C(mm)) = max_algae + 1 + mm
       !  enddo   ! mm

         do mm = 1, n_doc
            call icepack_init_bgc_trcr(&
                                      nk,               nt_fbri,       &
                                      nt_bgc_DOC(mm),    nlt_bgc_DOC(mm), &
                                      doctype(mm),       nt_depend,     &
                                      ntrcr,             nbtrcr,        &
                                      bgc_tracer_type,   trcr_depend,   &
                                      trcr_base,         n_trcr_strata, &
                                      nt_strata,         bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_DOC(mm)) = max_algae + 1 + mm
         enddo   ! mm
         do mm = 1, n_dic
            call icepack_init_bgc_trcr(&
                                      nk,               nt_fbri,       &
                                      nt_bgc_DIC(mm),    nlt_bgc_DIC(mm), &
                                      dictype(mm),       nt_depend,     &
                                      ntrcr,             nbtrcr,        &
                                      bgc_tracer_type,   trcr_depend,   &
                                      trcr_base,         n_trcr_strata, &
                                      nt_strata,         bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_DIC(mm)) = max_algae + max_doc + 1 + mm
         enddo   ! mm
      endif      ! tr_bgc_C

      if (tr_bgc_chl) then
         do mm = 1, n_algae
            call icepack_init_bgc_trcr(&
                                      nk,               nt_fbri,       &
                                      nt_bgc_chl(mm),    nlt_bgc_chl(mm), &
                                      algaltype(mm),     nt_depend,     &
                                      ntrcr,             nbtrcr,        &
                                      bgc_tracer_type,   trcr_depend,   &
                                      trcr_base,         n_trcr_strata, &
                                      nt_strata,         bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_chl(mm)) = max_algae + 1 + max_doc + max_dic + mm
         enddo   ! mm
      endif      ! tr_bgc_chl

      if (tr_bgc_Am) then
         call icepack_init_bgc_trcr(&
                                      nk,             nt_fbri,       &
                                      nt_bgc_Am,       nlt_bgc_Am,    &
                                      ammoniumtype,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_Am) = 2*max_algae + max_doc + max_dic + 2
      endif
      if (tr_bgc_Sil) then
         call icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc_Sil,      nlt_bgc_Sil,   &
                                      silicatetype,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_Sil) = 2*max_algae + max_doc + max_dic + 3
      endif
      if (tr_bgc_DMS) then   ! all together
         call icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc_DMSPp,    nlt_bgc_DMSPp, &
                                      dmspptype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_DMSPp) = 2*max_algae + max_doc + max_dic + 4

            call icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc_DMSPd,    nlt_bgc_DMSPd, &
                                      dmspdtype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_DMSPd) = 2*max_algae + max_doc + max_dic + 5

            call icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc_DMS,      nlt_bgc_DMS,   &
                                      dmspdtype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_DMS) = 2*max_algae + max_doc + max_dic + 6
      endif
      if (tr_bgc_PON) then
         call icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc_PON,      nlt_bgc_PON, &
                                      nitratetype,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_PON) =  2*max_algae + max_doc + max_dic + 7
      endif
      if (tr_bgc_DON) then
         do mm = 1, n_don
            call icepack_init_bgc_trcr(&
                                      nk,                nt_fbri,       &
                                      nt_bgc_DON(mm),    nlt_bgc_DON(mm), &
                                      dontype(mm),       nt_depend,     &
                                      ntrcr,             nbtrcr,        &
                                      bgc_tracer_type,   trcr_depend,   &
                                      trcr_base,         n_trcr_strata, &
                                      nt_strata,         bio_index)
            bio_index_o(nlt_bgc_DON(mm)) = 2*max_algae + max_doc + max_dic + 7 + mm
         enddo   ! mm
         if (icepack_warnings_aborted(subname)) return
      endif      ! tr_bgc_DON
      if (tr_bgc_Fe) then
         do mm = 1, n_fed
            call icepack_init_bgc_trcr(&
                                      nk,                nt_fbri,       &
                                      nt_bgc_Fed(mm),    nlt_bgc_Fed(mm), &
                                      fedtype(mm),       nt_depend,     &
                                      ntrcr,             nbtrcr,        &
                                      bgc_tracer_type,   trcr_depend,   &
                                      trcr_base,         n_trcr_strata, &
                                      nt_strata,         bio_index)
            bio_index_o(nlt_bgc_Fed(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + 7 + mm
         enddo   ! mm
         if (icepack_warnings_aborted(subname)) return
         do mm = 1, n_fep
            call icepack_init_bgc_trcr(&
                                      nk,                nt_fbri,       &
                                      nt_bgc_Fep(mm),    nlt_bgc_Fep(mm), &
                                      feptype(mm),       nt_depend,     &
                                      ntrcr,             nbtrcr,        &
                                      bgc_tracer_type,   trcr_depend,   &
                                      trcr_base,         n_trcr_strata, &
                                      nt_strata,         bio_index)
            bio_index_o(nlt_bgc_Fep(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + max_fe + 7 + mm
         enddo   ! mm
         if (icepack_warnings_aborted(subname)) return
      endif      ! tr_bgc_Fe

      if (tr_bgc_hum) then
         call icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc_hum,      nlt_bgc_hum,   &
                                      humtype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            if (icepack_warnings_aborted(subname)) return
            bio_index_o(nlt_bgc_hum) =   2*max_algae + max_doc + 8 + max_dic &
                                         + max_don + 2*max_fe + max_aero
      endif
      endif  ! skl_bgc or z_tracers

      if (z_tracers) then ! defined on nblyr+1 in ice
                              ! and 2 snow layers (snow surface + interior)

         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd

         ! z layer aerosols
         if (tr_zaero) then
            do mm = 1, n_zaero
               if (dEdd_algae) then
                  nlt_zaero_sw(mm) = nbtrcr_sw + 1
                  nbtrcr_sw = nbtrcr_sw + nilyr + nslyr+2
               endif
               call icepack_init_bgc_trcr(&
                                         nk,             nt_fbri,       &
                                         nt_zaero(mm),    nlt_zaero(mm), &
                                         zaerotype(mm),   nt_depend,     &
                                         ntrcr,           nbtrcr,        &
                                         bgc_tracer_type, trcr_depend,   &
                                         trcr_base,       n_trcr_strata, &
                                         nt_strata,       bio_index)
               bio_index_o(nlt_zaero(mm)) = 2*max_algae + max_doc + max_dic &
                                          + max_don + 2*max_fe + 7 + mm
            enddo   ! mm
            if (icepack_warnings_aborted(subname)) return
         endif      ! tr_zaero

         nt_zbgc_frac = 0
         if (nbtrcr > 0) then
            nt_zbgc_frac = ntrcr + 1
            ntrcr = ntrcr + nbtrcr
            do k = 1,nbtrcr
               zbgc_frac_init(k) = c1
               trcr_depend(nt_zbgc_frac+k-1) =  2+nt_fbri
               trcr_base(nt_zbgc_frac+ k - 1,1)  = c0
               trcr_base(nt_zbgc_frac+ k - 1,2)  = c1
               trcr_base(nt_zbgc_frac+ k - 1,3)  = c0
               n_trcr_strata(nt_zbgc_frac+ k - 1)= 1
               nt_strata(nt_zbgc_frac+ k - 1,1)  = nt_fbri
               nt_strata(nt_zbgc_frac+ k - 1,2)  = 0
               tau_ret(k) = c1
               tau_rel(k) = c1
               if (bgc_tracer_type(k) >=  c0 .and. bgc_tracer_type(k) < p5) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= p5 .and. bgc_tracer_type(k) < c1) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c1 .and. bgc_tracer_type(k) < c2) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c2 ) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               endif
            enddo
         endif

      endif ! z_tracers

      do k = 1, nbtrcr
         zbgc_init_frac(k) = frazil_scav
         if (bgc_tracer_type(k) < c0)  zbgc_init_frac(k) = initbio_frac
      enddo

      !-----------------------------------------------------------------
      ! final consistency checks
      !-----------------------------------------------------------------

      if (nbtrcr > max_nbtrcr) then
         write (warnstr,'(a,2i6)') subname//'ERROR:  nbtrcr > max_nbtrcr:',nbtrcr, max_nbtrcr
         call icepack_warnings_add(warnstr)
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
      endif
      if (.NOT. dEdd_algae) nbtrcr_sw = 1

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------

      if (present(printdiags)) then
      if (printdiags) then
         write(warnstr,1030) subname//' output:'
         if (skl_bgc) then
            write(warnstr,1010) ' skl_bgc                = ', skl_bgc
            write(warnstr,1030) ' bgc_flux_type          = ', bgc_flux_type
            write(warnstr,1020) ' number of bio tracers  = ', nbtrcr
            write(warnstr,1020) ' number of Isw tracers  = ', nbtrcr_sw
            write(warnstr,1020) ' number of autotrophs   = ', n_algae
            write(warnstr,1020) ' number of doc          = ', n_doc
            write(warnstr,1020) ' number of dic          = ', n_dic
            write(warnstr,1020) ' number of don          = ', n_don
            write(warnstr,1020) ' number of fed          = ', n_fed
            write(warnstr,1020) ' number of fep          = ', n_fep
            write(warnstr,1010) ' tr_bgc_N               = ', tr_bgc_N
            write(warnstr,1010) ' tr_bgc_C               = ', tr_bgc_C
            write(warnstr,1010) ' tr_bgc_chl             = ', tr_bgc_chl
            write(warnstr,1010) ' tr_bgc_Nit             = ', tr_bgc_Nit
            write(warnstr,1010) ' tr_bgc_Am              = ', tr_bgc_Am
            write(warnstr,1010) ' tr_bgc_Sil             = ', tr_bgc_Sil
            write(warnstr,1010) ' tr_bgc_hum             = ', tr_bgc_hum
            write(warnstr,1010) ' tr_bgc_DMS             = ', tr_bgc_DMS
            write(warnstr,1010) ' tr_bgc_PON             = ', tr_bgc_PON
            write(warnstr,1010) ' tr_bgc_DON             = ', tr_bgc_DON
            write(warnstr,1010) ' tr_bgc_Fe              = ', tr_bgc_Fe
         elseif (z_tracers) then
            write(warnstr,1010) ' dEdd_algae             = ', dEdd_algae
            write(warnstr,1010) ' modal_aero             = ', modal_aero
            write(warnstr,1010) ' scale_bgc              = ', scale_bgc
            write(warnstr,1010) ' solve_zbgc             = ', solve_zbgc
            write(warnstr,1020) ' number of ztracers     = ', nbtrcr
            write(warnstr,1020) ' number of Isw tracers  = ', nbtrcr_sw
            write(warnstr,1020) ' number of autotrophs   = ', n_algae
            write(warnstr,1020) ' number of doc          = ', n_doc
            write(warnstr,1020) ' number of dic          = ', n_dic
            write(warnstr,1020) ' number of fed          = ', n_fed
            write(warnstr,1020) ' number of fep          = ', n_fep
            write(warnstr,1020) ' number of aerosols     = ', n_zaero
            write(warnstr,1010) ' tr_zaero               = ', tr_zaero
            write(warnstr,1010) ' tr_bgc_Nit             = ', tr_bgc_Nit
            write(warnstr,1010) ' tr_bgc_N               = ', tr_bgc_N
            write(warnstr,1010) ' tr_bgc_Am              = ', tr_bgc_Am
            write(warnstr,1010) ' tr_bgc_C               = ', tr_bgc_C
            write(warnstr,1010) ' tr_bgc_Sil             = ', tr_bgc_Sil
            write(warnstr,1010) ' tr_bgc_hum             = ', tr_bgc_hum
            write(warnstr,1010) ' tr_bgc_chl             = ', tr_bgc_chl
            write(warnstr,1010) ' tr_bgc_DMS             = ', tr_bgc_DMS
            write(warnstr,1010) ' tr_bgc_PON             = ', tr_bgc_PON
            write(warnstr,1010) ' tr_bgc_DON             = ', tr_bgc_DON
            write(warnstr,1010) ' tr_bgc_Fe              = ', tr_bgc_Fe
            write(warnstr,1000) ' grid_o                 = ', grid_o
            write(warnstr,1005) ' l_sk                   = ', l_sk
            write(warnstr,1000) ' initbio_frac           = ', initbio_frac
            write(warnstr,1000) ' frazil_scav            = ', frazil_scav
      endif
      endif
      endif

      ! BGC Indices
      if (present(bio_index_out)    ) bio_index_out    = bio_index
      if (present(bio_index_o_out)  ) bio_index_o_out  = bio_index_o
      if (present(nt_fbri_out)      ) nt_fbri_out      = nt_fbri
      if (present(nt_bgc_Nit_out)   ) nt_bgc_Nit_out   = nt_bgc_Nit
      if (present(nlt_bgc_Nit_out)  ) nlt_bgc_Nit_out  = nlt_bgc_Nit
      if (present(nt_bgc_Am_out)    ) nt_bgc_Am_out    = nt_bgc_Am
      if (present(nlt_bgc_Am_out)   ) nlt_bgc_Am_out   = nlt_bgc_Am
      if (present(nt_bgc_Sil_out)   ) nt_bgc_Sil_out   = nt_bgc_Sil
      if (present(nlt_bgc_Sil_out)  ) nlt_bgc_Sil_out  = nlt_bgc_Sil
      if (present(nt_bgc_DMSPp_out) ) nt_bgc_DMSPp_out = nt_bgc_DMSPp
      if (present(nlt_bgc_DMSPp_out)) nlt_bgc_DMSPp_out= nlt_bgc_DMSPp
      if (present(nt_bgc_DMSPd_out) ) nt_bgc_DMSPd_out = nt_bgc_DMSPd
      if (present(nlt_bgc_DMSPd_out)) nlt_bgc_DMSPd_out= nlt_bgc_DMSPd
      if (present(nt_bgc_DMS_out)   ) nt_bgc_DMS_out   = nt_bgc_DMS
      if (present(nlt_bgc_DMS_out)  ) nlt_bgc_DMS_out  = nlt_bgc_DMS
      if (present(nt_bgc_hum_out)   ) nt_bgc_hum_out   = nt_bgc_hum
      if (present(nlt_bgc_hum_out)  ) nlt_bgc_hum_out  = nlt_bgc_hum
      if (present(nt_bgc_PON_out)   ) nt_bgc_PON_out   = nt_bgc_PON
      if (present(nlt_bgc_PON_out)  ) nlt_bgc_PON_out  = nlt_bgc_PON
      if (present(nt_bgc_N_out)     ) nt_bgc_N_out    = nt_bgc_N
      if (present(nlt_bgc_N_out)    ) nlt_bgc_N_out   = nlt_bgc_N
      if (present(nt_bgc_C_out)     ) nt_bgc_C_out    = nt_bgc_C
      if (present(nlt_bgc_C_out)    ) nlt_bgc_C_out   = nlt_bgc_C
      if (present(nt_bgc_chl_out)   ) nt_bgc_chl_out  = nt_bgc_chl
      if (present(nlt_bgc_chl_out)  ) nlt_bgc_chl_out = nlt_bgc_chl
      if (present(nt_bgc_DOC_out)   ) nt_bgc_DOC_out  = nt_bgc_DOC
      if (present(nlt_bgc_DOC_out)  ) nlt_bgc_DOC_out = nlt_bgc_DOC
      if (present(nt_bgc_DON_out)   ) nt_bgc_DON_out  = nt_bgc_DON
      if (present(nlt_bgc_DON_out)  ) nlt_bgc_DON_out = nlt_bgc_DON
      if (present(nt_bgc_DIC_out)   ) nt_bgc_DIC_out  = nt_bgc_DIC
      if (present(nlt_bgc_DIC_out)  ) nlt_bgc_DIC_out = nlt_bgc_DIC
      if (present(nt_bgc_Fed_out)   ) nt_bgc_Fed_out  = nt_bgc_Fed
      if (present(nlt_bgc_Fed_out)  ) nlt_bgc_Fed_out = nlt_bgc_Fed
      if (present(nt_bgc_Fep_out)   ) nt_bgc_Fep_out  = nt_bgc_Fep
      if (present(nlt_bgc_Fep_out)  ) nlt_bgc_Fep_out = nlt_bgc_Fep
      if (present(nt_zaero_out)     ) nt_zaero_out    = nt_zaero
      if (present(nlt_zaero_out)    ) nlt_zaero_out   = nlt_zaero
      if (present(nlt_zaero_sw_out) ) nlt_zaero_sw_out = nlt_zaero_sw
      if (present(nlt_chl_sw_out)   ) nlt_chl_sw_out   = nlt_chl_sw
      if (present(nt_zbgc_frac_out) ) nt_zbgc_frac_out = nt_zbgc_frac

      if (present(ntrcr_out)    ) ntrcr_out     = ntrcr
      if (present(ntrcr_o_out)  ) ntrcr_o_out   = ntrcr_o
      if (present(nbtrcr_out)   ) nbtrcr_out    = nbtrcr
      if (present(nbtrcr_sw_out)) nbtrcr_sw_out = nbtrcr_sw

      deallocate(algaltype)
      deallocate(doctype)
      deallocate(dictype)
      deallocate(dontype)
      deallocate(fedtype)
      deallocate(feptype)
      deallocate(zaerotype)

 1000    format (a,2x,f9.2)  ! float
 1005    format (a,2x,f9.6)  ! float
 1010    format (a,2x,l6)    ! logical
 1020    format (a,2x,i6)    ! integer
 1030    format (a,   a8)    ! character

      end subroutine icepack_init_zbgc_tracer_indices

!=======================================================================

      subroutine icepack_init_bgc_trcr(&
                                      nk,              nt_fbri,       &
                                      nt_bgc,          nlt_bgc,       &
                                      bgctype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)


      integer (kind=int_kind), intent(in) :: &
         nk           , & ! counter
         nt_depend    , & ! tracer dependency index
         nt_fbri

      integer (kind=int_kind), intent(inout) :: &
         ntrcr        , & ! number of tracers
         nbtrcr       , & ! number of bio tracers
         nt_bgc       , & ! tracer index
         nlt_bgc          ! bio tracer index

      integer (kind=int_kind), dimension(:), intent(inout) :: &
         trcr_depend  , & ! tracer dependencies
         n_trcr_strata, & ! number of underlying tracer layers
         bio_index        !

      integer (kind=int_kind), dimension(:,:), intent(inout) :: &
         nt_strata        ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcr_base        ! = 0 or 1 depending on tracer dependency
                          ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=dbl_kind), intent(in) :: &
         bgctype          ! bio tracer transport type (mobile vs stationary)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         bgc_tracer_type  ! bio tracer transport type array

      ! local variables

      integer (kind=int_kind) :: &
         k         , & ! loop index
         n_strata  , & ! temporary values
         nt_strata1, & !
         nt_strata2

      real (kind=dbl_kind) :: &
         trcr_base1, & ! temporary values
         trcr_base2, &
         trcr_base3

      character(len=*),parameter :: subname='(icepack_init_bgc_trcr)'

         nt_bgc = ntrcr + 1
         nbtrcr = nbtrcr + 1
         nlt_bgc = nbtrcr
         bgc_tracer_type(nbtrcr) = bgctype

         if (nk > 1) then
            ! include vertical bgc in snow
            do k = nk, nk+1
               ntrcr = ntrcr + 1
               trcr_depend  (nt_bgc + k  ) = 2 ! snow volume
               trcr_base    (nt_bgc + k,1) = c0
               trcr_base    (nt_bgc + k,2) = c0
               trcr_base    (nt_bgc + k,3) = c1
               n_trcr_strata(nt_bgc + k  ) = 0
               nt_strata    (nt_bgc + k,1) = 0
               nt_strata    (nt_bgc + k,2) = 0
            enddo

            trcr_base1 = c0
            trcr_base2 = c1
            trcr_base3 = c0
            n_strata = 1
            nt_strata1 = nt_fbri
            nt_strata2 = 0
         else  ! nk = 1
            trcr_base1 = c1
            trcr_base2 = c0
            trcr_base3 = c0
            n_strata = 0
            nt_strata1 = 0
            nt_strata2 = 0
         endif ! nk

         do k = 1, nk     !in ice
            ntrcr = ntrcr + 1
            trcr_depend  (nt_bgc + k - 1  ) = nt_depend
            trcr_base    (nt_bgc + k - 1,1) = trcr_base1
            trcr_base    (nt_bgc + k - 1,2) = trcr_base2
            trcr_base    (nt_bgc + k - 1,3) = trcr_base3
            n_trcr_strata(nt_bgc + k - 1  ) = n_strata
            nt_strata    (nt_bgc + k - 1,1) = nt_strata1
            nt_strata    (nt_bgc + k - 1,2) = nt_strata2
         enddo

         bio_index (nlt_bgc) = nt_bgc

       end subroutine icepack_init_bgc_trcr

!=======================================================================

      end module icepack_zbgc

!=======================================================================
