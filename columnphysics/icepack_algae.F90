!=======================================================================
!
! Compute sea ice biogeochemistry (vertical or skeletal layer)
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module icepack_algae

      use icepack_kinds

      use icepack_parameters, only: p05, p5, c0, c1, c2, c6, c10, p1
      use icepack_parameters, only: pi, secday, puny
      use icepack_parameters, only: hs_ssl, sk_l

      use icepack_parameters, only: dEdd_algae, solve_zbgc, use_atm_dust_iron
      use icepack_parameters, only: R_dFe2dust, dustFe_sol, algal_vel
      use icepack_parameters, only: bgc_flux_type
      use icepack_parameters, only: grid_o
      use icepack_parameters, only: T_max, fsal      , fr_resp
      use icepack_parameters, only: op_dep_min       , fr_graze_s
      use icepack_parameters, only: fr_graze_e       , fr_mort2min
      use icepack_parameters, only: fr_dFe           , k_nitrif
      use icepack_parameters, only: t_iron_conv      , max_loss
      use icepack_parameters, only: max_dfe_doc1     , fr_resp_s
      use icepack_parameters, only: y_sk_DMS         , t_sk_conv
      use icepack_parameters, only: t_sk_ox

      use icepack_tracers, only: nblyr, nilyr, nslyr, ntrcr, nbtrcr
      use icepack_tracers, only: n_algae, n_doc, n_dic, n_don, n_fed, n_fep, n_zaero
      use icepack_tracers, only: bio_index, bio_index_o
      use icepack_tracers, only: nt_bgc_N, nt_fbri, nt_zbgc_frac
      use icepack_tracers, only: tr_brine
      use icepack_tracers, only: nt_bgc_DON, nt_bgc_hum, nt_bgc_DOC
      use icepack_tracers, only: tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil
      use icepack_tracers, only: tr_bgc_DMS,    tr_bgc_PON,   tr_bgc_hum
      use icepack_tracers, only: tr_bgc_N,      tr_bgc_C,     tr_bgc_chl
      use icepack_tracers, only: tr_bgc_DON,    tr_bgc_Fe,    tr_zaero
      use icepack_tracers, only: nlt_bgc_Nit,   nlt_bgc_Am,   nlt_bgc_Sil
      use icepack_tracers, only: nlt_bgc_DMS,   nlt_bgc_PON
      use icepack_tracers, only: nlt_bgc_N,     nlt_bgc_C,    nlt_bgc_chl
      use icepack_tracers, only: nlt_bgc_DOC,   nlt_bgc_DON,  nlt_bgc_DIC
      use icepack_tracers, only: nlt_zaero  ,   nlt_bgc_DMSPp,nlt_bgc_DMSPd
      use icepack_tracers, only: nlt_bgc_Fed,   nlt_bgc_Fep,  nlt_bgc_hum

      use icepack_zbgc_shared, only: remap_zbgc, regrid_stationary
      use icepack_zbgc_shared, only: merge_bgc_fluxes
      use icepack_zbgc_shared, only: merge_bgc_fluxes_skl
      use icepack_zbgc_shared, only: bgrid, cgrid, igrid, icgrid
      use icepack_zbgc_shared, only: phi_sk, bgc_tracer_type
      use icepack_zbgc_shared, only: zbgc_init_frac
      use icepack_zbgc_shared, only: zbgc_frac_init
      use icepack_zbgc_shared, only: tau_rel, tau_ret, thinS
      use icepack_zbgc_shared, only: r_Si2N, R_Fe2N, R_S2N, R_C2N_DON
      use icepack_zbgc_shared, only: chlabs, alpha2max_low, beta2max, mu_max
      use icepack_zbgc_shared, only: k_exude, K_Nit, K_Am, K_Sil, K_Fe
      use icepack_zbgc_shared, only: grow_Tdep, fr_graze, mort_pre, mort_Tdep
      use icepack_zbgc_shared, only: f_don, kn_bac, f_don_Am
      use icepack_zbgc_shared, only: f_doc, f_exude, k_bac, R_chl2N, R_C2N
      use icepack_zbgc_shared, only: graze_exponent, graze_conc, large_bgc

      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_aerosol, only: update_snow_bgc

      implicit none

      private
      public :: zbio, sklbio

      real (kind=dbl_kind), parameter :: &
           exp_argmax = c10,           & ! maximum argument of exponential
           accuracy = 1.0e-13_dbl_kind   ! accuracy parameter for bgc tracers
!=======================================================================

      contains

!=======================================================================

      subroutine zbio   (dt,                        &
                         meltt,        melts,       &
                         meltb,        congel,      &
                         snoice,       fsnow,       &
                         trcrn,        bio_index,   &
                         bio_index_o,  aice_old,    &
                         vice_old,     vsno_old,    &
                         vicen,        vsnon,       &
                         aicen,        flux_bio_atm,&
                         n_cat,        first_ice,   &
                         hice_old,     ocean_bio,   &
                         ocean_bio_dh,              &
                         bphin,        iphin,       &
                         iDin,                      &
                         fswthrul,                  &
                         dh_top,       dh_bot,      &
                         zfswin,                    &
                         hbri,         hbri_old,    &
                         darcy_V,                   &
                         bphi_min,                  &
                         iTin,                      &
                         Zoo,                       &
                         flux_bio,     dh_direct,   &
                         upNO,         upNH,        &
                         fbio_snoice,  fbio_atmice, &
                         PP_net,       ice_bio_net, &
                         snow_bio_net, grow_net,    &
                         totalChla,                 &
                         flux_bion,    iSin,        &
                         bioPorosityIceCell,        &
                         bioSalinityIceCell,        &
                         bioTemperatureIceCell      )

      integer (kind=int_kind), intent(in) :: &
         n_cat        ! category number

      integer (kind=int_kind), dimension (nbtrcr), intent(in) :: &
         bio_index, & ! references index of bio tracer (nbtrcr) to tracer array (ntrcr)
         bio_index_o  ! references index of data arrays (eg. kscavz)

      real (kind=dbl_kind), intent(in) :: &
         dt,       &  ! time step
         hbri,     &  ! brine height  (m)
         bphi_min, &  ! surface porosity
         meltt,    &  ! thermodynamic melt/growth rates in dt (m)
         melts,    &
         meltb,    &
         congel,   &
         snoice,   &
         fsnow,    & ! snowfall rate (kg/m^2 s)
         hice_old, & ! ice height (m)
         vicen,    & ! ice volume (m)
         vsnon,    & ! snow volume (m)
         aicen,    & ! ice area fraction
         aice_old, & ! values prior to thermodynamic changes
         vice_old, &
         vsno_old, &
         darcy_V,  & ! darcy velocity
!        darcy_V_chl,& ! darcy velocity for algae
         dh_bot,     & ! change in brine bottom (m)
         dh_top,     & ! change in brine top (m)
         dh_direct     ! surface flooding or surface runoff (m)

      real (kind=dbl_kind), dimension (nbtrcr), intent(inout) :: &
         snow_bio_net,& ! net bio tracer in snow (mmol/m^2)
         ice_bio_net, & ! net bio tracer in ice (mmol/m^2)
         fbio_atmice, & ! bio flux from atm to ice (mmol/m^2/s)
         fbio_snoice, & ! bio flux from snow to ice  (mmol/m^2/s)
         flux_bio,    & ! total ocean tracer flux (mmol/m^2/s)
         flux_bion      ! category ocean tracer flux (mmol/m^2/s)

      real (kind=dbl_kind), intent(in) :: &
         hbri_old       ! brine height  (m)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         iTin       , & ! salinity vertical interface points
         iphin      , & ! Porosity on the igrid
         iDin       , & ! Diffusivity/h on the igrid (1/s)
         iSin           ! Salinity on vertical interface points (ppt)

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         fswthrul       ! visible short wave radiation on icgrid (W/m^2)

      real (kind=dbl_kind), dimension(nbtrcr), intent(in) :: &
         flux_bio_atm   ! aerosol/bgc deposition rate (mmol/m^2 s)

      real (kind=dbl_kind), dimension(ntrcr), intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         zfswin         ! visible Short wave flux on igrid (W/m^2)

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         Zoo            ! N losses to the system from reaction terms
                        ! (ie. zooplankton/bacteria) (mmol/m^3)

      real (kind=dbl_kind), dimension (nbtrcr), intent(in) :: &
         !change to  inout when updating ocean fields
         ocean_bio      ! ocean concentrations (mmol/m^3)

      real (kind=dbl_kind), optional, dimension (:), intent(out) :: &
         !change to  inout when updating ocean fields
         ocean_bio_dh   ! ocean concentrations * hbrine * phi (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bphin          ! Porosity on the bgrid

      real (kind=dbl_kind), intent(inout):: &
         PP_net     , & ! net PP (mg C/m^2/d)  times aice
         grow_net   , & ! net specific growth (m/d) times vice
         upNO       , & ! tot nitrate uptake rate (mmol/m^2/d) times aice
         upNH           ! tot ammonium uptake rate (mmol/m^2/d) times aice

      real (kind=dbl_kind), optional, intent(inout):: &
         totalChla      ! total chla (mg chla/m^2)

      real (kind=dbl_kind), optional, dimension (:), intent(inout):: &  ! diagnostics (nblyr+1)
         bioPorosityIceCell , & ! porosity on vertical interface points
         bioSalinityIceCell , & ! salinity on vertical interface points (ppt)
         bioTemperatureIceCell  ! temperature on vertical interface points (oC)

      logical (kind=log_kind), intent(in) :: &
         first_ice      ! initialized values should be used

      ! local variables

      integer (kind=int_kind) :: &
         mm              ! thickness category index

      real (kind=dbl_kind), dimension (nblyr+1,n_algae) :: &
         upNOn      , & ! algal nitrate uptake rate  (mmol/m^3/s)
         upNHn      , & ! algal ammonium uptake rate (mmol/m^3/s)
         grow_alg       ! algal growth rate          (mmol/m^3/s)

      real (kind=dbl_kind),dimension(nbtrcr) :: &
         zbgc_snown, & ! aerosol contribution from snow to ice
         zbgc_atmn     ! and atm to ice concentration * volume (mmol/m^3*m)

      real (kind=dbl_kind), dimension(nbtrcr) :: &
         Tot_BGC_i, & ! initial column sum, ice + snow,  of BGC tracer (mmol/m^2)
         Tot_BGC_f, & ! final column sum
         flux_bio_sno !

      real (kind=dbl_kind) :: &
         hsnow_i,  & ! initial snow thickness (m)
         hsnow_f, & ! final snow thickness (m)
         carbonError ! carbon conservation error (mmol/m2)

      real (kind=dbl_kind) :: &
         carbonInitial, & ! initial carbon content (mmol/m2)
         carbonFinal,   & ! final carbon content (mmol/m2)
         carbonFlux       ! carbon flux (mmol/m2/s)

      logical (kind=log_kind) :: &
         write_flux_diag, &
         write_carbon_errors

      character(len=*),parameter :: subname='(zbio)'

      zbgc_snown(:) = c0
      zbgc_atmn (:) = c0
      flux_bion (:) = c0
      flux_bio_sno(:) = c0
      Tot_BGC_i (:) = c0
      Tot_BGC_f (:) = c0
      Zoo (:) = c0
      hsnow_i = c0
      hsnow_f = c0
      write_flux_diag = .false.
      write_carbon_errors = .true.
      if (.not. tr_bgc_C) write_carbon_errors = .false.

      call bgc_carbon_sum(hbri_old, trcrn(:), carbonInitial)
      if (icepack_warnings_aborted(subname)) return

      if (aice_old > puny) then
          hsnow_i = vsno_old/aice_old
          do  mm = 1,nbtrcr
             call bgc_column_sum (hsnow_i, hbri_old, &
                              trcrn(bio_index(mm):bio_index(mm)+nblyr+2), &
                              Tot_BGC_i(mm))
             if (icepack_warnings_aborted(subname)) return
          enddo
      endif

      call update_snow_bgc     (dt,                      &
                                meltt,     melts,        &
                                meltb,     congel,       &
                                snoice,    fsnow,        &
                                trcrn,     bio_index,    &
                                aice_old,  zbgc_snown,   &
                                vice_old,  vsno_old,     &
                                vicen,     vsnon,        &
                                aicen,     flux_bio_atm, &
                                zbgc_atmn, flux_bio_sno, &
                                bio_index_o)

      if (icepack_warnings_aborted(subname)) return

      call z_biogeochemistry   (n_cat,        dt,        &
                                first_ice, &
                                aicen,        vicen,     &
                                hice_old,     ocean_bio, &
                                ocean_bio_dh,            &
                                flux_bion,    bphin,     &
                                iphin,        trcrn,     &
                                iDin,                    &
                                fswthrul,     grow_alg,  &
                                upNOn,        upNHn,     &
                                dh_top,       dh_bot,    &
                                zfswin,       hbri,      &
                                hbri_old,     darcy_V,   &
                                bphi_min,     zbgc_snown,&
                                zbgc_atmn,               &
                                iTin,         dh_direct, &
                                Zoo,          meltb,     &
                                congel                   )
      if (icepack_warnings_aborted(subname)) return

      do mm = 1,nbtrcr
         flux_bion(mm) = flux_bion(mm) + flux_bio_sno(mm)
      enddo

      call bgc_carbon_sum(hbri, trcrn(:), carbonFinal)
      if (icepack_warnings_aborted(subname)) return
      call bgc_carbon_flux(flux_bio_atm,flux_bion,carbonFlux)
      if (icepack_warnings_aborted(subname)) return

      carbonError = carbonInitial-carbonFlux*dt-carbonFinal

      if (abs(carbonError) > max(puny,accuracy * maxval ((/carbonInitial, carbonFinal/))) .and. write_carbon_errors) then
            write(warnstr,*) subname, 'carbonError:', carbonError
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'carbonInitial:', carbonInitial
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'carbonFinal:', carbonFinal
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'carbonFlux (positive into ocean):', carbonFlux
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'accuracy * maxval ((/carbonInitial, carbonFinal/:)', accuracy * maxval ((/carbonInitial, carbonFinal/))
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'puny', puny
            call icepack_warnings_add(warnstr)
            if (aicen > c0) then
            hsnow_f = vsnon/aicen
            write(warnstr,*) subname, 'after z_biogeochemistry'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'Remaining carbon after algal_dyn: Zoo'
            call icepack_warnings_add(warnstr)
            do mm = 1,nblyr+1
               write(warnstr,*) subname, 'layer mm, Zoo(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, mm,Zoo(mm)
               call icepack_warnings_add(warnstr)
            end do
            do mm = 1,nbtrcr
               call bgc_column_sum (hsnow_f, hbri, &
                              trcrn(bio_index(mm):bio_index(mm)+nblyr+2), &
                              Tot_BGC_f(mm))
               write(warnstr,*) subname, 'mm, Tot_BGC_i(mm), Tot_BGC_f(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  mm, Tot_BGC_i(mm), Tot_BGC_f(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'flux_bion(mm), flux_bio_atm(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  flux_bion(mm), flux_bio_atm(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zbgc_snown(mm),zbgc_atmn(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  zbgc_snown(mm),zbgc_atmn(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Tot_BGC_i(mm) + flux_bio_atm(mm)*dt - flux_bion(mm)*dt'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  Tot_BGC_i(mm) + flux_bio_atm(mm)*dt - flux_bion(mm)*dt
               call icepack_warnings_add(warnstr)
            enddo
         endif
         !call icepack_warnings_add(warnstr)
         !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//" zbio: Carbon conservation failure after z_biogeochemistry")
      endif
      if (icepack_warnings_aborted(subname)) return

      if (write_flux_diag) then
         if (aicen > c0) then
            hsnow_f = vsnon/aicen
            do mm = 1,nbtrcr
               call bgc_column_sum (hsnow_f, hbri, &
                              trcrn(bio_index(mm):bio_index(mm)+nblyr+2), &
                              Tot_BGC_f(mm))
               write(warnstr,*) subname, 'mm, Tot_BGC_i(mm), Tot_BGC_f(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  mm, Tot_BGC_i(mm), Tot_BGC_f(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'flux_bion(mm), flux_bio_atm(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  flux_bion(mm), flux_bio_atm(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zbgc_snown(mm),zbgc_atmn(mm)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  zbgc_snown(mm),zbgc_atmn(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Tot_BGC_i(mm) + flux_bio_atm(mm)*dt - flux_bion(mm)*dt'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  Tot_BGC_i(mm) + flux_bio_atm(mm)*dt - flux_bion(mm)*dt
               call icepack_warnings_add(warnstr)
            enddo
         endif
      endif
      if (icepack_warnings_aborted(subname)) return

      call merge_bgc_fluxes   (dt,                       &
                               bio_index,                &
                               aicen,                    &
                               vicen,        vsnon,      &
                               iphin,                    &! ntrcr
                               trcrn,        aice_old,   &!aice_old
                               flux_bion,    flux_bio,   &
                               upNOn,        upNHn,      &
                               upNO,         upNH,       &
                               zbgc_snown,   zbgc_atmn,  &
                               fbio_snoice,  fbio_atmice,&
                               PP_net,       ice_bio_net,&
                               snow_bio_net, grow_alg,   &
                               grow_net,     totalChla,  &
                               iTin,         iSin,       &
                               bioPorosityIceCell,       &
                               bioSalinityIceCell,       &
                               bioTemperatureIceCell)
      if (icepack_warnings_aborted(subname)) return

      if (write_flux_diag) then
         if (aicen > c0) then
            write(warnstr,*) subname, 'after merge_bgc_fluxes, n_cat:', n_cat
            call icepack_warnings_add(warnstr)
            do mm = 1,nbtrcr
               write(warnstr,*) subname,  'mm, flux_bio(mm):',mm,flux_bio(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'fbio_snoice(mm)',fbio_snoice(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'fbio_atmice(mm)',fbio_atmice(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  'flux_bio_atm(mm)', flux_bio_atm(mm)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,  'flux_bio_atm(mm)*aicen', flux_bio_atm(mm)*aicen
               call icepack_warnings_add(warnstr)
            enddo
         endif
      endif
      if (icepack_warnings_aborted(subname)) return

      end subroutine zbio

!=======================================================================

      subroutine sklbio       (dt,       Tf,         &
                               flux_bio, ocean_bio,  &
                               aicen,      &
                               meltb,    congel,     &
                               fswthru,  first_ice,  &
                               trcrn,  &
                               PP_net,   upNO,       &
                               upNH,     grow_net    )

      logical (kind=log_kind), intent(in) :: &
         first_ice      ! initialized values should be used

      real (kind=dbl_kind), intent(in) :: &
         dt,       &  ! time step
         Tf,       &  ! basal freezing temperature (C)
!        hmix,     &  ! mixed layer depth (m)
         aicen,    &  ! ice area fraction
         meltb,    &  ! bottom melt (m)
         congel,   &  ! bottom growth (m)
         fswthru      ! visible shortwave passing to ocean(W/m^2)

      real (kind=dbl_kind), dimension(ntrcr), intent(inout) :: &
         trcrn      ! bulk concentration per m^3

      real (kind=dbl_kind), dimension (nbtrcr), intent(inout) :: &
         flux_bio   ! ocean tracer flux (mmol/m^2/s) positive into ocean

      real (kind=dbl_kind), dimension (nbtrcr), intent(in) :: &
         ocean_bio  ! ocean tracer concentration (mmol/m^3)

      ! history output
      real (kind=dbl_kind), intent(inout):: &
         PP_net  , & ! Bulk net PP (mg C/m^2/s)
         grow_net, & ! net specific growth (/s)
         upNO    , & ! tot nitrate uptake rate (mmol/m^2/s)
         upNH        ! tot ammonium uptake rate (mmol/m^2/s)

      ! local variables

      real (kind=dbl_kind), dimension (n_algae) :: &
         upNOn      , & ! algal nitrate uptake rate  (mmol/m^3/s)
         upNHn      , & ! algal ammonium uptake rate (mmol/m^3/s)
         grow_alg       ! algal growth rate          (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nbtrcr) :: &
         flux_bion       !tracer flux to ocean

      character(len=*),parameter :: subname='(sklbio)'

      flux_bion (:) = c0
      upNOn     (:) = c0
      upNHn     (:) = c0
      grow_alg  (:) = c0

      call skl_biogeochemistry       (dt, &
                                      flux_bion, ocean_bio, &
!                                     hmix,      aicen,     &
                                      meltb,     congel,    &
                                      fswthru,   first_ice, &
                                      trcrn,     upNOn,     &
                                      upNHn,     grow_alg,  &
                                      Tf)
      if (icepack_warnings_aborted(subname)) return

      call merge_bgc_fluxes_skl    ( &
                                    aicen,     trcrn,       &
                                    flux_bion, flux_bio,    &
                                    PP_net,    upNOn,       &
                                    upNHn,     upNO,        &
                                    upNH,      grow_net,    &
                                    grow_alg)
      if (icepack_warnings_aborted(subname)) return

      end subroutine sklbio

!=======================================================================
!
! skeletal layer biochemistry
!
      subroutine skl_biogeochemistry (dt, &
                                      flux_bio,   ocean_bio,    &
!                                     hmix,       aicen,        &
                                      meltb,      congel,       &
                                      fswthru,    first_ice,    &
                                      trcrn,      upNOn,        &
                                      upNHn,      grow_alg_skl, &
                                      Tf)

      real (kind=dbl_kind), intent(in) :: &
         dt     , & ! time step
!        hmix   , & ! mixed layer depth
!        aicen  , & ! ice area
         meltb  , & ! bottom ice melt
         congel , & ! bottom ice growth
         Tf     , & ! bottom freezing temperature
         fswthru    ! shortwave passing through ice to ocean

      logical (kind=log_kind), intent(in) :: &
         first_ice  ! initialized values should be used

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         trcrn      ! bulk concentration per m^3

      ! history variables

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         flux_bio   ! ocean tracer flux (mmol/m^2/s) positive into ocean

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio  ! ocean tracer concentration (mmol/m^3)

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         grow_alg_skl, & ! tot algal growth rate (mmol/m^3/s)
         upNOn       , & !  algal NO uptake rate (mmol/m^3/s)
         upNHn           !  algal NH uptake rate (mmol/m^3/s)

      ! local variables

      integer (kind=int_kind) :: nn

      real (kind=dbl_kind), dimension(nbtrcr):: &
         react        , & ! biological sources and sinks (mmol/m^3)
         cinit        , & ! initial brine concentration*sk_l (mmol/m^2)
         cinit_v      , & ! initial brine concentration (mmol/m^3)
         congel_alg   , & ! congelation flux contribution to ice algae (mmol/m^2 s)
                          ! (used as initialization)
         f_meltn      , & ! vertical melt fraction of skeletal layer in dt
         flux_bio_temp, & ! tracer flux to ocean (mmol/m^2 s)
         PVflag       , & ! 1 for tracers that flow with the brine, 0 otherwise
         cling            ! 1 for tracers that cling, 0 otherwise

      real (kind=dbl_kind) :: &
         Zoo_skl      , & ! N losses from zooplankton/bacteria ... (mmol/m^3)
         iTin         , &
         PVt          , & ! type 'Jin2006' piston velocity (m/s)
         ice_growth   , & ! Jin2006 definition: either congel rate or bottom melt rate  (m/s)
         grow_val     , & ! (m/x)
         rphi_sk      , & ! 1 / skeletal layer porosity
         cinit_tmp    , & ! temporary variable for concentration (mmol/m^2)
         Cerror       , & ! change in total carbon from reactions (mmol/m^3)
         nitrification    ! nitrate from nitrification (mmol/m^3)

      real (kind=dbl_kind), parameter :: &
         PVc = 1.e-6_dbl_kind           , & ! type 'constant' piston velocity for interface (m/s)
         PV_scale_growth = p5           , & ! scale factor in Jin code PV during ice growth
         PV_scale_melt = p05            , & ! scale factor in Jin code PV during ice melt
         growth_max = 1.85e-5_dbl_kind , & ! PVt function reaches maximum here.  (m/s)
         MJ1 = 9.667e-9_dbl_kind        , & ! (m/s) coefficients in Jin2008
         MJ2 = 38.8_dbl_kind            , & ! (1) from:4.49e-4_dbl_kind*secday
         MJ3 = 1.04e7_dbl_kind          , & ! 1/(m/s) from: 1.39e-3_dbl_kind*secday^2
         PV_frac_max = 0.9_dbl_kind         ! Maximum Piston velocity is 90% of skeletal layer/dt

      logical (kind=log_kind) :: conserve_C

      character(len=*),parameter :: subname='(skl_biogeochemistry)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      conserve_C = .true.
      Zoo_skl    = c0
      rphi_sk    = c1/phi_sk
      PVt        = c0
      iTin       = Tf
      ice_growth = (congel-meltb)/dt

      do nn = 1, nbtrcr
         cinit     (nn) = c0
         cinit_v   (nn) = c0
         congel_alg(nn) = c0
         f_meltn   (nn) = c0
         react     (nn) = c0
         PVflag    (nn) = c1
         cling     (nn) = c0

      !-----------------------------------------------------------------
      ! only the dominant tracer_type affects behavior
      !  < 0 is purely mobile:  > 0 stationary behavior
      ! NOTE: retention times are not used in skl model
      !-----------------------------------------------------------------

         if (bgc_tracer_type(nn) >= c0) then
            PVflag(nn) = c0
            cling (nn) = c1
         endif

         cinit  (nn) = trcrn(bio_index(nn)) * sk_l * rphi_sk
         cinit_v(nn) = cinit(nn)/sk_l
         if (cinit(nn) < c0) then
            write(warnstr,*) subname,'initial sk_bgc < 0, nn,nbtrcr,cinit(nn)', &
                 nn,nbtrcr,cinit(nn)
            call icepack_warnings_add(warnstr)
            call icepack_warnings_add(subname//' cinit < c0')
            !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            return
         endif
      enddo     ! nbtrcr

      if (icepack_warnings_aborted(subname)) return

      if (trim(bgc_flux_type) == 'Jin2006') then

      !-----------------------------------------------------------------
      ! 'Jin2006':
      ! 1. congel/melt dependent piston velocity (PV) for growth and melt
      ! 2. If congel > melt use 'congel'; if melt > congel use 'melt'
      ! 3. For algal N, PV for ice growth only provides a seeding concentration
      ! 4. Melt affects nutrients and algae in the same manner through PV(melt)
      !-----------------------------------------------------------------

         if (ice_growth > c0) then  ! ice_growth = congel/dt
            grow_val = min(ice_growth,growth_max)
            PVt = -min(abs(PV_scale_growth*(MJ1 + MJ2*grow_val &
                                                - MJ3*grow_val**2)), &
                           PV_frac_max*sk_l/dt)
         else                       ! ice_growth = -meltb/dt
            PVt =  min(abs(PV_scale_melt  *(      MJ2*ice_growth &
                                                - MJ3*ice_growth**2)), &
                           PV_frac_max*sk_l/dt)
         endif
         do nn = 1, nbtrcr
            if (bgc_tracer_type(nn) >= c0) then
               if (ice_growth < c0) then ! flux from ice to ocean
                  ! Algae and clinging tracers melt like nutrients
                  f_meltn(nn) = PVt*cinit_v(nn) ! for algae only
               elseif (ice_growth > c0 .AND. &
                  cinit(nn) < ocean_bio(nn)*sk_l/phi_sk) then
                  ! Growth only contributes to seeding from ocean
                  congel_alg(nn) = (ocean_bio(nn)*sk_l/phi_sk - cinit(nn))/dt
               endif ! PVt > c0
            endif
         enddo

      else   ! bgc_flux_type = 'constant'

      !-----------------------------------------------------------------
      ! 'constant':
      ! 1. Constant PV for congel > melt
      ! 2. For algae, PV for ice growth only provides a seeding concentration
      ! 3. Melt loss (f_meltn) affects algae only and is proportional to melt
      !-----------------------------------------------------------------

         if (ice_growth > c0) PVt = -PVc
         do nn = 1, nbtrcr
            if (bgc_tracer_type(nn) >= c0 ) then
               if (ice_growth >= c0 .AND. cinit_v(nn) < ocean_bio(nn)/phi_sk) then
                  congel_alg(nn) = (ocean_bio(nn)*sk_l/phi_sk - cinit(nn))/dt
               elseif (ice_growth < c0) then
                  f_meltn(nn) = min(c1, meltb/sk_l)*cinit(nn)/dt
               endif
            endif
         enddo ! nn

      endif  ! bgc_flux_type

      !-----------------------------------------------------------------
      ! begin building biogeochemistry terms
      !-----------------------------------------------------------------

      react(:) = c0
      grow_alg_skl(:) = c0

      call algal_dyn (dt,                         &
                      dEdd_algae,                 &
                      fswthru,         react,     &
                      cinit_v,                    &
                      grow_alg_skl(:),            &
                      iTin,                       &
                      upNOn(:),        upNHn(:),  &
                      Zoo_skl,                    &
                      Cerror,          conserve_C,&
                      nitrification)
      if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! compute new tracer concencentrations
      !-----------------------------------------------------------------

      do nn = 1, nbtrcr

      !-----------------------------------------------------------------
      ! if PVt > 0, ie melt, then ocean_bio term drops out (MJ2006)
      ! Combine boundary fluxes
      !-----------------------------------------------------------------

         PVflag(nn) = SIGN(PVflag(nn),PVt)
         cinit_tmp = max(c0, cinit_v(nn) + react(nn))
         flux_bio_temp(nn) = (PVflag(nn)*PVt*cinit_tmp &
                           -  PVflag(nn)*min(c0,PVt)*ocean_bio(nn)) &
                           + f_meltn(nn)*cling(nn) - congel_alg(nn)

         if (cinit_tmp*sk_l < flux_bio_temp(nn)*dt) then
            flux_bio_temp(nn) = cinit_tmp*sk_l/dt*(c1-puny)
         endif

         cinit(nn) = cinit_tmp*sk_l - flux_bio_temp(nn)*dt
         flux_bio(nn) = flux_bio(nn) + flux_bio_temp(nn)*phi_sk

         ! Uncomment to update ocean concentration
         ! Currently not coupled with ocean biogeochemistry
!         ocean_bio(nn) = ocean_bio(nn) + flux_bio(nn)/hmix*aicen

         if (.not. conserve_C) then
              write(warnstr,*) subname, 'C not conserved in skl_bgc, Cerror:',Cerror
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'sk_bgc < 0 after algal fluxes, nn,cinit,flux_bio',&
                               nn,cinit(nn),flux_bio(nn)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'cinit_tmp,flux_bio_temp,f_meltn,congel_alg,PVt,PVflag: '
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, cinit_tmp,flux_bio_temp(nn),f_meltn(nn), &
                               congel_alg(nn),PVt,PVflag(nn)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'congel, meltb: ',congel,meltb
              call icepack_warnings_add(warnstr)
              call icepack_warnings_add(subname//' N not conserved in skl_bgc')
              !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (cinit(nn) < c0) then
              write(warnstr,*) subname, 'sk_bgc < 0 after algal fluxes, nn,cinit,flux_bio',&
                               nn,cinit(nn),flux_bio(nn)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'cinit_tmp,flux_bio_temp,f_meltn,congel_alg,PVt,PVflag: '
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, cinit_tmp,flux_bio_temp(nn),f_meltn(nn), &
                               congel_alg(nn),PVt,PVflag(nn)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'congel, meltb: ',congel,meltb
              call icepack_warnings_add(warnstr)
              call icepack_warnings_add(subname//'sk_bgc < 0 after algal fluxes')
              !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         endif

         if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! reload tracer array:  Bulk tracer concentration (mmol or mg per m^3)
      !-----------------------------------------------------------------

         trcrn(bio_index(nn)) = cinit(nn) * phi_sk/sk_l

       enddo  !nbtrcr

      end subroutine skl_biogeochemistry

!=======================================================================
!
! Solve the scalar vertical diffusion equation implicitly using
! tridiag_solver. Calculate the diffusivity from temperature and salinity.
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with
! dynamic salinity or the height ratio == hinS/vicen*aicen, where hinS is the
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds).
!

      subroutine z_biogeochemistry (n_cat,        dt,        &
                                    first_ice, &
                                    aicen,        vicen,     &
                                    hice_old,     ocean_bio, &
                                    ocean_bio_dh,            &
                                    flux_bio,     bphin,     &
                                    iphin,        trcrn,     &
                                    iDin,                    &
                                    fswthrul,     grow_alg,  &
                                    upNOn,        upNHn,     &
                                    dh_top,       dh_bot,    &
                                    zfswin,       hbri,      &
                                    hbri_old,     darcy_V,   &
                                    bphi_min,     zbgc_snow, &
                                    zbgc_atm,                &
                                    iTin,         dh_direct, &
                                    Zoo,          meltb,     &
                                    congel                   )

      integer (kind=int_kind), intent(in) :: &
         n_cat          ! category number

      logical (kind=log_kind), intent(in) :: &
         first_ice      ! initialized values should be used

      real (kind=dbl_kind), intent(in) :: &
         dt         , & ! time step
         hbri       , & ! brine height  (m)
         bphi_min   , & ! surface porosity
         aicen      , & ! concentration of ice
         vicen      , & ! volume per unit area of ice  (m)
         hice_old   , & ! ice height (m)
         meltb      , & ! bottom melt in dt (m)
         congel     , & ! bottom growth in dt (m)
         darcy_V    , & ! darcy velocity
         dh_bot     , & ! change in brine bottom (m)
         dh_top     , & ! change in brine top (m)
         dh_direct      ! surface flooding or runoff (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         flux_bio   , & ! total ocean tracer flux (mmol/m^2/s)
         zfswin     , & ! visible Short wave flux on igrid (W/m^2)
         Zoo        , & ! N losses to the system from reaction terms
                        ! (ie. zooplankton/bacteria) (mmol/m^3)
         trcrn          ! bulk tracer concentration (mmol/m^3)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         iTin       , & ! salinity vertical interface points
         iphin      , & ! Porosity on the igrid
         iDin       , & ! Diffusivity/h on the igrid (1/s)
         fswthrul   , & ! visible short wave radiation on icgrid (W/m^2)
         zbgc_snow  , & ! tracer input from snow (mmol/m^3*m)
         zbgc_atm   , & ! tracer input from  atm (mmol/m^3 *m)
         ocean_bio  , & ! ocean concentrations (mmol/m^3)
         bphin          ! Porosity on the bgrid

      real (kind=dbl_kind), optional, dimension (:), intent(out) :: &
         ocean_bio_dh   ! ocean concentrations * hbrine * phi (mmol/m^2)

      real (kind=dbl_kind), intent(in) :: &
         hbri_old       ! brine height  (m)

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         upNOn      , & ! algal nitrate uptake rate  (mmol/m^3/s)
         upNHn      , & ! algal ammonium uptake rate (mmol/m^3/s)
         grow_alg       ! algal growth rate          (mmol/m^3/s)

      !-----------------------------------------------------------------------------
      ! algae absorption coefficient for 0.5 m thick layer
      ! Grenfell (1991): SA = specific absorption coefficient= 0.004 m^2/mg Chla
      ! generalizing kalg_bio(k) = SA*\sum R_chl2N(m)*trcrn(i,j,nt_bgc_N(m)+k-1)
      ! output kalg on the icgrid
      !-----------------------------------------------------------------------------

      ! local variables

      integer (kind=int_kind) :: &
         k, m, mm        ! vertical biology layer index

      real (kind=dbl_kind) :: &
         hin         , & ! ice thickness (m)
         hin_old     , & ! ice thickness before current melt/growth (m)
         ice_conc    , & ! algal concentration in ice above hin > hinS
         sum_initial , & !
         sum_old     , & !
         sum_new     , & !
         sum_tot     , & !
         zspace      , & ! 1/nblyr
         darcyV      , & !
         dhtop       , & !
         dhbot       , & !
         dhflood         ! >=0 (m) surface flooding from the ocean

      real (kind=dbl_kind), dimension (nblyr+2) :: &
         bphin_N         ! porosity for tracer model has minimum
                         ! bphin_N >= bphimin

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         iphin_N      , & ! tracer porosity on the igrid
         sbdiagz      , & ! sub-diagonal matrix elements
         diagz        , & ! diagonal matrix elements
         spdiagz      , & ! super-diagonal matrix elements
         rhsz         , & ! rhs of tri-diagonal matrix equation
         ML_diag      , & ! lumped mass matrix
         D_spdiag     , & ! artificial diffusion matrix
         D_sbdiag     , & ! artificial diffusion matrix
         biomat_low   , & ! Low order solution
         Cerror           ! Change in N after reactions

      real (kind=dbl_kind), dimension(nblyr+1,nbtrcr):: &
         react            ! biological sources and sinks for equation matrix

      real (kind=dbl_kind), dimension(nblyr+1,nbtrcr):: &
         in_init_cons , & ! Initial bulk concentration*h (mmol/m^2)
         biomat_cons  , & ! Matrix output of (mmol/m^2)
         biomat_brine     ! brine concentration (mmol/m^3)

      real (kind=dbl_kind), dimension(nbtrcr):: &
         C_top,         & ! bulk top tracer source: h phi C(meltwater) (mmol/m^2)
         C_bot,         & ! bulk bottom tracer source: h phi C(ocean) (mmol/m^2)
         Source_top,    & ! For cons: (+) top tracer source into ice (mmol/m^2/s)
         Source_bot,    & ! For cons: (+) bottom tracer source into ice (mmol/m^2/s)
         Sink_bot,      & ! For cons: (+ or -) remaining bottom flux into ice(mmol/m^2/s)
         Sink_top,      & ! For cons: (+ or -) remaining bottom flux into ice(mmol/m^2/s)
         exp_ret,       & ! exp dt/retention frequency
         exp_rel,       & ! exp dt/release frequency
         atm_add_cons , & ! zbgc_snow+zbgc_atm (mmol/m^3*m)
         dust_Fe      , & ! contribution of dust surface flux to dFe (umol/m*3*m)
         source           ! mmol/m^2 surface input from snow/atmosphere

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0       , & ! temporary, remapped tracers
         trtmp            ! temporary, remapped tracers

      logical (kind=log_kind), dimension(nblyr+1) :: &
         conserve_C

      real (kind=dbl_kind), dimension(nblyr+1):: &  ! temporary variables for
         Diff         , & ! diffusivity
         initcons     , & ! initial concentration
         biocons      , & !  new concentration
         dmobile      , & !
         initcons_mobile,&!
         initcons_stationary, &
         dz           , & ! normalized vertical grid spacing
         nitrification    ! nitrate produced from nitrification (mmol/m3)

      real (kind=dbl_kind):: &
         top_conc         ! 1% (min_bgc) of surface concentration
                          ! when hin > hbri:  just used in sw calculation

      real (kind=dbl_kind):: &
         bio_tmp, &       ! temporary variable
         exp_min          ! temporary exp var

      real (kind=dbl_kind):: &
         Sat_conc   , & ! adsorbing saturation concentration  (mmols/m^3)
         phi_max    , & ! maximum porosity
         S_col      , & ! surface area of collector (um^2)
         P_b        , & ! projected area of diatoms (um^2)
         V_c        , & ! volume of collector  (um^3)
         V_alg          ! volume of algae (um^3)

      real (kind=dbl_kind), dimension(nbtrcr) :: &
         mobile     , & ! c1 if mobile, c0 otherwise
         flux_bio_tmp

      ! local parameters

      real (kind=dbl_kind), parameter :: &
         r_c  = 3.0e3_dbl_kind     , & ! ice crystal radius (um)
         r_bac= 4.7_dbl_kind     , & ! diatom large radius (um)
         r_alg= 10.0_dbl_kind    , & ! diatom small radius (um)
         Nquota_A = 0.88_dbl_kind, & ! slope in Nitrogen quota to cell volume fit
                                     ! (Lomas et al. 2019, Edwards et al. 2012)
         Nquota_I = 0.0408_dbl_kind, & ! Intercept in N quota to cell volume fit
         f_s = p1, & ! fracton of sites available for saturation
         f_a = 0.3_dbl_kind, & !c1 , &  ! fraction of collector available for attachment
         f_v = 0.7854_dbl_kind ! fraction of algal coverage on area availabel for attachment
                       ! 4(pi r^2)/(4r)^2  [Johnson et al, 1995, water res. research]

      integer, parameter :: &
         nt_zfswin = 1    ! for interpolation of short wave to bgrid

      character(len=*),parameter :: subname='(z_biogeochemistry)'

      logical (kind=log_kind) :: &
         write_carbon_errors

  !-------------------------------------
  ! Initialize
  !-----------------------------------

      write_carbon_errors = .true.
      if (.not. tr_bgc_C) write_carbon_errors = .false.

      zspace = c1/real(nblyr,kind=dbl_kind)
      dz(:) = zspace
      dz(1) = zspace/c2
      dz(nblyr+1) = zspace/c2
      in_init_cons(:,:) = c0
      atm_add_cons(:) = c0
      dhtop = c0
      dhbot = c0
      darcyV = c0
      C_top(:) = c0
      C_bot(:) = c0
      if (present(ocean_bio_dh)) ocean_bio_dh(:) = c0
      mobile(:) = c0
      conserve_C(:) = .true.
      nitrification(:) = c0
      flux_bio_tmp(:) = c0

      do m = 1, nbtrcr
         do k  = 1, nblyr+1

            bphin_N(nblyr+2) =c1
            bphin_N(k) = bphin(k)
            iphin_N(k) = iphin(k)
            bphin_N(1) = bphi_min

            if (abs(trcrn(bio_index(m) + k-1)) < accuracy) then
               flux_bio_tmp(m) = trcrn(bio_index(m) + k-1)* hbri_old * dz(k)/dt
               trcrn(bio_index(m) + k-1) = c0
               in_init_cons(k,m) = c0
            else
               in_init_cons(k,m) = trcrn(bio_index(m) + k-1)* hbri_old
            endif

            if (trcrn(bio_index(m) + k-1) < c0  ) then
               write(warnstr,*) subname,'zbgc initialization error, first ice = ', first_ice
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,'Category,m:',n_cat,m
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,'hbri,hbri_old'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, hbri,hbri_old
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname,'trcrn(bio_index(m) + k-1)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, trcrn(bio_index(m) + k-1)
               call icepack_warnings_add(warnstr)
               call icepack_warnings_add(subname//' zbgc initialization error')
               !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            endif
            if (icepack_warnings_aborted(subname)) return
        enddo         !k
      enddo           !m

      !-----------------------------------------------------------------
      !     boundary conditions
      !-----------------------------------------------------------------

      ice_conc = c0
      hin = vicen/aicen
      hin_old = hice_old

      !-----------------------------------------------------------------
      !    calculate the saturation concentration for attachment: Sat_conc
      !-----------------------------------------------------------------

      phi_max = maxval(bphin_N(2:nblyr+1))
      S_col   = 4.0_dbl_kind*pi*r_c**2
      P_b     = pi*r_bac**2    !*10-6 for colloids
      V_c     = 4.0_dbl_kind*pi*r_c**3/3.0_dbl_kind  ! (m^3) sphere
      V_alg   = pi/6.0_dbl_kind*r_bac*r_alg**2       ! prolate spheroid (*10-9 for colloids)
      Sat_conc= f_s*f_a*f_v*(c1-phi_max)/V_c*S_col/P_b*(V_alg)**Nquota_A*Nquota_I * 1.0e9_dbl_kind
      !mmol/m^3 (algae, don, hum...) and umols/m^3 for colloids
      !-----------------------------------------------------------------
      !    convert surface dust flux (n_zaero > 2) to dFe(1) flux
      !-----------------------------------------------------------------

      dust_Fe(:) = c0

      if (tr_zaero .and. n_zaero > 2 .and. tr_bgc_Fe .and. use_atm_dust_iron) then
       do m = 3,n_zaero
         dust_Fe(nlt_bgc_Fed(1)) = dust_Fe(nlt_bgc_Fed(1)) + &
              (zbgc_snow(nlt_zaero(m)) + zbgc_atm(nlt_zaero(m))) * &
               R_dFe2dust * dustFe_sol
        ! dust_Fe(nlt_zaero(m)) = -(zbgc_snow(nlt_zaero(m)) + zbgc_atm(nlt_zaero(m))) * &
        !       dustFe_sol
       enddo
      endif

      do m = 1,nbtrcr
      !-----------------------------------------------------------------
      !   time constants for mobile/stationary phase changes
      !-----------------------------------------------------------------

         exp_rel(m) = c0
         exp_ret(m) = c0
         if (tau_ret(m) > c0) then
            exp_min = min(dt/tau_ret(m),exp_argmax)
            exp_ret(m) = exp(-exp_min)
         endif
         if (tau_rel(m) > c0) then
            exp_min = min(dt/tau_rel(m),exp_argmax)
            exp_rel(m) = exp(-exp_min)
         endif
         if (m .ne. nlt_bgc_N(1)) then
            if (hin_old  > hin) then  !melting
               exp_ret(m) = c1
            else                              !not melting
               exp_rel(m) = c1
            endif
         elseif (tr_bgc_N .and. hin_old > hin + algal_vel*dt) then
               exp_ret(m) = c1
         elseif (tr_bgc_N) then
               exp_rel(m) = c1
         endif

         dhtop      = dh_top
         dhbot      = dh_bot
         darcyV     = darcy_V
         C_top(m)   = in_init_cons(1,m)*trcrn(nt_zbgc_frac+m-1)!mobile fraction
         source(m)  = abs(zbgc_snow(m) + zbgc_atm(m) + dust_Fe(m))
         dhflood  = max(c0,-dh_direct)                              ! ocean water flooding surface

         if (dhtop+darcyV/bphin_N(1)*dt < -puny) then !snow/top ice melt
             C_top(m) = (zbgc_snow(m)+zbgc_atm(m) + dust_Fe(m))/abs(dhtop &
                        + darcyV/bphin_N(1)*dt + puny)*hbri_old
         elseif (dhtop+darcyV/bphin_N(1)*dt >= -puny .and. &
                        abs((zbgc_snow(m)+zbgc_atm(m) + dust_Fe(m)) + &
                        ocean_bio(m)*bphin_N(1)*dhflood) >  puny) then
              atm_add_cons(m) =  abs(zbgc_snow(m) + zbgc_atm(m)+ dust_Fe(m)) + &
                                      ocean_bio(m)*bphin_N(1)*dhflood
         else   ! only positive fluxes
              atm_add_cons(m) =  abs(zbgc_snow(m) + zbgc_atm(m)+ dust_Fe(m))
         endif

         C_bot(m) = ocean_bio(m)*hbri_old*iphin_N(nblyr+1)
         if (present(ocean_bio_dh)) ocean_bio_dh(m) = C_bot(m)

      enddo             ! m

      !-----------------------------------------------------------------
      ! Interpolate shortwave flux, fswthrul (defined at top to bottom with nilyr+1
      !  evenly spaced  with spacing = (1/nilyr) to grid variable zfswin:
      !-----------------------------------------------------------------

      trtmp(:) = c0
      trtmp0(:)= c0
      zfswin(:) = c0

      do k = 1, nilyr+1
         ! contains cice values (fswthrul(1) is surface value)
         ! and fwsthrul(nilyr+1) is output
         trtmp0(nt_zfswin+k-1) = fswthrul(k)
      enddo   !k

      call remap_zbgc(nilyr+1,  &
                      nt_zfswin,                  &
                      trtmp0(1:ntrcr),  trtmp(1:ntrcr+2), &
                      0,                nblyr+1,  &
                      hin,              hbri,     &
                      icgrid(1:nilyr+1),          &
                      igrid(1:nblyr+1),ice_conc   )

      if (icepack_warnings_aborted(subname)) return

      do k = 1,nblyr+1
         zfswin(k) = trtmp(nt_zfswin+k-1)
      enddo

      !-----------------------------------------------------------------
      ! Initialze Biology
      !-----------------------------------------------------------------

      do mm = 1, nbtrcr
         mobile(mm) = c0
         if (bgc_tracer_type(mm) .GE. c0) mobile(mm) = c1

         do k = 1, nblyr+1
            biomat_cons(k,mm) = in_init_cons(k,mm)
         enddo  !k
      enddo  !mm

      !-----------------------------------------------------------------
      ! Compute FCT
      !-----------------------------------------------------------------

      do mm = 1, nbtrcr

         if (hbri_old > thinS .and. hbri > thinS) then
            do k = 1,nblyr+1
               initcons_mobile(k) = in_init_cons(k,mm)*trcrn(nt_zbgc_frac+mm-1)
               initcons_stationary(k) = max(c0,in_init_cons(k,mm)-initcons_mobile(k))
               ! Allow release of Nitrate/silicate to mobile phase, but not adsorption
               dmobile(k) = mobile(mm)*(initcons_mobile(k)*(exp_ret(mm)-c1) + &
                    initcons_stationary(k)*(c1-exp_rel(mm))) + &
                    (1-mobile(mm))*initcons_stationary(k)*(c1-exp_rel(mm))
               initcons_mobile(k) = max(c0,initcons_mobile(k) + dmobile(k))
               initcons_stationary(k) = max(c0,initcons_stationary(k) - dmobile(k))
               if (initcons_stationary(k)/hbri_old > Sat_conc) then
                  initcons_mobile(k) = initcons_mobile(k) + initcons_stationary(k) - Sat_conc*hbri_old
                  initcons_stationary(k) = Sat_conc*hbri_old
               endif

               Diff(k) = iDin(k)
               initcons(k) = initcons_mobile(k)
               biocons(k) =  initcons_mobile(k)
            enddo

            call compute_FCT_matrix &
                                (initcons,sbdiagz, dt,         &
                                diagz, spdiagz, rhsz,          &
                                darcyV,    dhtop,              &
                                dhbot,   iphin_N,              &
                                Diff, hbri_old,                &
                                atm_add_cons(mm), bphin_N,     &
                                C_top(mm), C_bot(mm),          &
                                Source_bot(mm), Source_top(mm),&
                                Sink_bot(mm),Sink_top(mm),     &
                                D_sbdiag, D_spdiag, ML_diag)
            if (icepack_warnings_aborted(subname)) return

            call tridiag_solverz &
                               (nblyr+1, sbdiagz,               &
                                diagz,   spdiagz,               &
                                rhsz,    biocons)
            if (icepack_warnings_aborted(subname)) return

            call check_conservation_FCT &
                               (initcons,    &
                                biocons,     &
                                biomat_low,               &
                                Source_top(mm),        &
                                Source_bot(mm),        &
                                Sink_bot(mm),          &
                                Sink_top(mm),          &
                                dt, flux_bio(mm),     &
                                source(mm))
            if (icepack_warnings_aborted(subname)) return

            call compute_FCT_corr &
                                (initcons, biocons, dt, &
                                 D_sbdiag, D_spdiag, ML_diag)
            if (icepack_warnings_aborted(subname)) return

            top_conc = c0        ! or frazil ice concentration

            ! assume diatoms actively maintain there relative position in the ice

            if (mm .ne. nlt_bgc_N(1)) then

               call regrid_stationary &
                                (initcons_stationary,    hbri_old,    &
                                 hbri,                   dt,          &
                                 top_conc,                            &
                                 igrid,                  flux_bio(mm),&
                                 meltb,                  congel)
               if (icepack_warnings_aborted(subname)) return

            elseif (tr_bgc_N .and. mm .eq. nlt_bgc_N(1)) then
               if (meltb > algal_vel*dt .or. aicen < 0.001_dbl_kind) then

                  call regrid_stationary &
                                (initcons_stationary,    hbri_old,    &
                                 hbri,                   dt,          &
                                 top_conc,                            &
                                 igrid,                  flux_bio(mm),&
                                 meltb,                  congel)
                  if (icepack_warnings_aborted(subname)) return

               endif
            endif

            biomat_cons(:,mm) =  biocons(:) +  initcons_stationary(:)

            sum_initial = (in_init_cons(1,mm) + in_init_cons(nblyr+1,mm))*zspace/c2
            sum_old = (biomat_low(1) + biomat_low(nblyr+1))*zspace/c2
            sum_new = (biocons(1)+ biocons(nblyr+1))*zspace/c2
            sum_tot = (biomat_cons(1,mm) + biomat_cons(nblyr+1,mm))*zspace/c2
            do k = 2,nblyr
               sum_initial = sum_initial + in_init_cons(k,mm)*zspace
               sum_old = sum_old + biomat_low(k)*zspace
               sum_new = sum_new + biocons(k)*zspace
               sum_tot = sum_tot + biomat_cons(k,mm)*zspace
            enddo
            trcrn(nt_zbgc_frac+mm-1) = zbgc_frac_init(mm)
            if (sum_tot > c0) trcrn(nt_zbgc_frac+mm-1) = sum_new/sum_tot

            if ((abs((sum_initial-sum_tot+source(mm))/dt-flux_bio(mm)) > max(puny, accuracy*abs(flux_bio(mm)))) &
                .or. (minval(biocons(:)) < c0)  .or. (minval(initcons_stationary(:)) < c0)  &
                .or. icepack_warnings_aborted()) then
                write(warnstr,*) subname,'zbgc FCT tracer solution failed, mm:', mm
                call icepack_warnings_add(warnstr)
                write(warnstr,*)'sum_new,sum_tot,sum_initial,flux_bio(mm),source(mm):'
                call icepack_warnings_add(warnstr)
                write(warnstr,*)sum_new,sum_tot,sum_initial,flux_bio(mm),source(mm)
                call icepack_warnings_add(warnstr)
                write(warnstr,*)'error = (sum_initial-sum_tot+source(mm))/dt-flux_bio(mm)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*)(sum_initial-sum_tot+source(mm))/dt-flux_bio(mm)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,'sum_new,sum_old:',sum_new,sum_old
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,'mobile(mm):',mobile(mm)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'trcrn(nt_zbgc_frac+mm-1):',trcrn(nt_zbgc_frac+mm-1)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'exp_ret( mm),exp_rel( mm)',exp_ret( mm),exp_rel( mm)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,'darcyV,dhtop,dhbot'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,darcyV,dhtop,dhbot
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,'Category,mm:',n_cat,mm
                call icepack_warnings_add(warnstr)
                call icepack_warnings_add(subname//'zbgc FCT tracer solution failed')
                !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            endif
            if (icepack_warnings_aborted(subname)) return

         else

            call thin_ice_flux(hbri,hbri_old,biomat_cons(:,mm), &
                               flux_bio(mm),source(mm), &
                               dt, ocean_bio(mm))
            if (icepack_warnings_aborted(subname)) return

         endif ! thin or not

         do k = 1,nblyr+1
            biomat_brine(k,mm) =  biomat_cons(k,mm)/hbri/iphin_N(k)
         enddo ! k
      enddo ! mm

      react(:,:) = c0
      grow_alg(:,:) = c0

      if (solve_zbgc) then
         do k = 1, nblyr+1
            call algal_dyn (dt,                            &
                         dEdd_algae,                       &
                         zfswin(k),        react(k,:),     &
                         biomat_brine(k,:),                &
                         grow_alg(k,:),                    &
                         iTin(k),                          &
                         upNOn(k,:),       upNHn(k,:),     &
                         Zoo(k),                           &
                         Cerror(k),        conserve_C(k),  &
                         nitrification(k))
            if (icepack_warnings_aborted(subname)) return

         enddo ! k
      endif    ! solve_zbgc

      !-----------------------------------------------------------------
      ! Update the tracer variable
      !-----------------------------------------------------------------

      sum_new = c0
      sum_tot = c0

      do m = 1,nbtrcr
         do k = 1,nblyr+1                  ! back to bulk quantity
            bio_tmp = (biomat_brine(k,m) + react(k,m))*iphin_N(k)

            if (tr_bgc_C .and. m .eq. nlt_bgc_DIC(1) .and. bio_tmp .le. -accuracy) then  ! satisfy DIC demands from ocean
                !Uncomment for additional diagnostics
                !write(warnstr,*) subname, 'DIC demand from ocean'
                !call icepack_warnings_add(warnstr)
                !write(warnstr,*) subname, 'm, k, nlt_bgc_DIC(1), bio_tmp, react(k,m):'
                !call icepack_warnings_add(warnstr)
                !write(warnstr,*) subname, m, k, nlt_bgc_DIC(1), bio_tmp, react(k,m)
                !call icepack_warnings_add(warnstr)
                !write(warnstr,*) subname, 'flux_bio(m), hbri, hbri_old:'
                !call icepack_warnings_add(warnstr)
                !write(warnstr,*) subname, flux_bio(m), hbri, hbri_old
                !call icepack_warnings_add(warnstr)
                flux_bio(m) = flux_bio(m) + bio_tmp*dz(k)*hbri/dt
                bio_tmp = c0
                !write(warnstr,*) subname, 'flux_bio(m) Final:'
                !call icepack_warnings_add(warnstr)
                !write(warnstr,*) subname, flux_bio(m)
                !call icepack_warnings_add(warnstr)
            end if
            if (m .eq. nlt_bgc_Nit) then
               initcons_mobile(k) = max(c0,(biomat_brine(k,m)-nitrification(k) + &
                  react(k,m))*iphin_N(k)*trcrn(nt_zbgc_frac+m-1))
               initcons_stationary(k) = max(c0,((c1-trcrn(nt_zbgc_frac+m-1))*(biomat_brine(k,m)- &
                  nitrification(k) + react(k,m)) + nitrification(k))*iphin_N(k))

               sum_new = sum_new + initcons_mobile(k)*dz(k)
               sum_tot = sum_tot + (initcons_mobile(k) + initcons_stationary(k))*dz(k)

            end if  ! m .eq. nlt_bgc_Nit
            if (.not. conserve_C(k) .and. write_carbon_errors) then
                write(warnstr,*) subname, 'C in algal_dyn not conserved'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'Cerror(k):', Cerror(k)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'k,m,hbri,hbri_old,bio_tmp,biomat_cons(k,m),ocean_bio(m)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  k,m,hbri,hbri_old,bio_tmp,biomat_cons(k,m),ocean_bio(m)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'react(k,m),iphin_N(k),biomat_brine(k,m)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  react(k,m),iphin_N(k),biomat_brine(k,m)
                call icepack_warnings_add(warnstr)
                call icepack_warnings_add(subname//' C in algal_dyn not conserved')
                !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            elseif (abs(bio_tmp) < accuracy) then
               flux_bio(m) = flux_bio(m) + bio_tmp*dz(k)*hbri/dt
               bio_tmp = c0
            elseif (bio_tmp > large_bgc) then
                write(warnstr,*) subname, 'very large bgc value'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'k,m,hbri,hbri_old,bio_tmp,biomat_cons(k,m),ocean_bio(m)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  k,m,hbri,hbri_old,bio_tmp,biomat_cons(k,m),ocean_bio(m)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'react(k,m),iphin_N(k),biomat_brine(k,m)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  react(k,m),iphin_N(k),biomat_brine(k,m)
                call icepack_warnings_add(warnstr)
                call icepack_warnings_add(subname//' very large bgc value')
                !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            elseif (bio_tmp < c0) then
                write(warnstr,*) subname, 'negative bgc'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'k,m,nlt_bgc_Nit,hbri,hbri_old'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  k,m,nlt_bgc_Nit,hbri,hbri_old
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'bio_tmp,biomat_cons(k,m),ocean_bio(m)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  bio_tmp,biomat_cons(k,m),ocean_bio(m)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'react(k,m),iphin_N(k),biomat_brine(k,m)'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  react(k,m),iphin_N(k),biomat_brine(k,m)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'exp_ret( m),exp_ret( m)',exp_ret( m),exp_ret( m)
                call icepack_warnings_add(warnstr)
                call icepack_warnings_add(subname//'negative bgc')
                !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            endif
            trcrn(bio_index(m)+k-1) = max(c0, bio_tmp)
            if (icepack_warnings_aborted()) then
                write(warnstr,*) subname, 'trcrn(nt_zbgc_frac+m-1):',trcrn(nt_zbgc_frac+m-1)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'in_init_cons(k,m):',in_init_cons(k,m)
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'trcrn(bio_index(m) + k-1), bio_tmp'
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname,  trcrn(bio_index(m) + k-1), bio_tmp
                call icepack_warnings_add(warnstr)
                write(warnstr,*) subname, 'Category,m:',n_cat,m
                call icepack_warnings_add(warnstr)
                return
            endif
         enddo        ! k
         if (m .eq. nlt_bgc_Nit .and. MAXVAL(nitrification) > c0) then
            trcrn(nt_zbgc_frac+m-1) = zbgc_frac_init(m)
            if (sum_tot > c0) trcrn(nt_zbgc_frac+m-1) = sum_new/sum_tot
         end if
         flux_bio(m) = flux_bio(m) + flux_bio_tmp(m)
      enddo        ! m

      end subroutine z_biogeochemistry

!=======================================================================
!
! Do biogeochemistry from subroutine algal_dynamics
! authors: Scott Elliott, LANL
!          Nicole Jeffery, LANL

      subroutine algal_dyn (dt,                         &
                            dEdd_algae,                 &
                            fswthru,      reactb,       &
                            ltrcrn,                     &
                            grow_alg,                   &
                            T_bot,                      &
                            upNOn,        upNHn,        &
                            Zoo,                        &
                            Cerror,       conserve_C,   &
                            nitrification)

      real (kind=dbl_kind), intent(in) :: &
         dt      , & ! time step
         T_bot   , & ! ice temperature (oC)
         fswthru     ! average shortwave passing through current ice layer (W/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         Zoo,     & ! N losses from zooplankton/bacteria... (mmol/m^3)
         Cerror,  & ! Change in C after reactions (mmol/m^3)
         nitrification ! nitrate produced through nitrification (mmol/m3)

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         grow_alg,& !  algal growth rate   (mmol/m^3/s)
         upNOn,   & !  algal NO uptake rate   (mmol/m^3/s)
         upNHn      !  algal NH uptake rate   (mmol/m^3/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         reactb     ! biological reaction terms (mmol/m3)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         ltrcrn     ! brine concentrations  in layer (mmol/m^3)

      logical (kind=log_kind), intent(inout) :: &
         conserve_C

      logical (kind=log_kind), intent(in) :: &
         dEdd_algae  ! .true.  chla impact on shortwave computed in dEdd

      !  local variables
      !------------------------------------------------------------------------------------
      !            3 possible autotrophs nt_bgc_N(1:3):  diatoms, flagellates, phaeocystis
      !                2 types of dissolved organic carbon nt_bgc_DOC(1:2):
      !                        polysaccharids, lipids
      !                1 DON (proteins)
      !                1 particulate iron (nt_bgc_Fe) n_fep
      !                1 dossp;ved orpm m+fed
      ! Limiting macro/micro nutrients: nt_bgc_Nit -> nitrate, nt_bgc_NH -> ammonium,
      !                        nt_bgc_Sil -> silicate, nt_bgc_Fe -> dissolved iron
      ! --------------------------------------------------------------------------------------

      integer (kind=int_kind) :: k, n

      real (kind=dbl_kind), dimension(n_algae) :: &
         Nin        , &     ! algal nitrogen concentration on volume (mmol/m^3)
         Cin        , &     ! algal carbon concentration on volume (mmol/m^3)
         chlin              ! algal chlorophyll concentration on volume (mg/m^3)

      real (kind=dbl_kind), dimension(n_doc) :: &
         DOCin              ! dissolved organic carbon concentration on volume (mmolC/m^3)

      real (kind=dbl_kind), dimension(n_dic) :: &
         DICin              ! dissolved inorganic carbon concentration on volume (mmolC/m^3)

      real (kind=dbl_kind), dimension(n_don) :: &  !proteins
         DONin              ! dissolved organic nitrogen concentration on volume (mmolN/m^3)

      real (kind=dbl_kind), dimension(n_fed) :: &  !iron
         Fedin              ! dissolved iron concentration on volume (umol/m^3)

      real (kind=dbl_kind), dimension(n_fep) :: &  !iron
         Fepin              ! algal nitrogen concentration on volume (umol/m^3)

      real (kind=dbl_kind) :: &
         Nitin      , &     ! nitrate concentration on volume (mmol/m^3)
         Amin       , &     ! ammonia/um concentration on volume (mmol/m^3)
         Silin      , &     ! silicon concentration on volume (mmol/m^3)
!        DMSPpin    , &     ! DMSPp concentration on volume (mmol/m^3)
         DMSPdin    , &     ! DMSPd concentration on volume (mmol/m^3)
         DMSin      , &     ! DMS concentration on volume (mmol/m^3)
 !        PONin      , &     ! PON concentration on volume (mmol/m^3)
         op_dep     , &     ! bottom layer attenuation exponent (optical depth)
         Iavg_loc           ! bottom layer attenuated Fswthru (W/m^2)

      real (kind=dbl_kind), dimension(n_algae) :: &
         L_lim    , &  ! overall light limitation
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonium limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         Fe_lim   , &  ! overall iron limitation
         fr_Nit   , &  ! fraction of local ecological growth as nitrate
         fr_Am    , &  ! fraction of local ecological growth as ammonia
         growmax_N, &  ! maximum growth rate in N currency (mmol/m^3/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^3/s)
         potU_Am  , &  ! potential ammonium uptake (mmol/m^3/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^3/s)
         U_Am     , &  ! actual ammonium uptake (mmol/m^3/s)
         U_Sil    , &  ! actual silicon uptake (mmol/m^3/s)
         U_Fe     , &  ! actual iron uptake   (umol/m^3/s)
         U_Nit_f  , &  ! fraction of Nit uptake due to each algal species
         U_Am_f   , &  ! fraction of Am uptake due to each algal species
         U_Sil_f  , &  ! fraction of Sil uptake due to each algal species
         U_Fe_f        ! fraction of Fe uptake due to each algal species

      real (kind=dbl_kind) :: &
         dTemp        , &  ! sea ice temperature minus sst (oC) < 0
         U_Nit_tot    , &  ! actual nitrate uptake (mmol/m^3/s)
         U_Am_tot     , &  ! actual ammonium uptake (mmol/m^3/s)
         U_Sil_tot    , &  ! actual silicon uptake (mmol/m^3/s)
         U_Fe_tot     , &  ! actual iron uptake   (umol/m^3/s)
         nitrif       , &  ! nitrification (mmol/m^3/s)
         mort_N       , &  ! total algal mortality (mmol N/m^3)
         mort_C       , &  ! total algal mortality (mmol C/m^3)
         graze_N      , &  ! total algae grazed (mmol N/m^3)
         graze_C      , &  ! total algae grazed (mmol C/m^3)
         exude_C      , &  ! total carbon exuded by algae (mmol C/m^3)
         resp_N       , &  ! total N in respiration (mmol N/m^3)
         growth_N          ! total algal growth (mmol N/m^3)

      real (kind=dbl_kind), dimension(n_algae) :: &
         resp     , &  ! respiration (mmol/m^3/s)
         graze    , &  ! grazing (mmol/m^3/s)
         mort          ! sum of mortality and excretion (mmol/m^3/s)

!  source terms underscore s, removal underscore r

      real (kind=dbl_kind), dimension(n_algae) :: &
         N_s       , &  ! net algal nitrogen sources (mmol/m^3)
         N_r            ! net algal nitrogen removal (mmol/m^3)

      real (kind=dbl_kind), dimension(n_doc) :: &
         DOC_r      , &  ! net DOC removal (mmol/m^3)
         DOC_s           ! net DOC sources (mmol/m^3)

      real (kind=dbl_kind), dimension(n_dic) :: &
         DIC_r      , &  ! net DIC removal (mmol/m^3)
         DIC_s           ! net DIC sources (mmol/m^3)

      real (kind=dbl_kind), dimension(n_don) :: &
         DON_r      , &  ! net DON removal (mmol/m^3)
         DON_s           ! net DON sources (mmol/m^3)

      real (kind=dbl_kind), dimension(n_fed) :: &
         Fed_r_l     , &  ! removal due to loss of binding saccharids (nM)
         Fed_r       , &  ! net Fed removal (nM)
         Fed_s       , &  ! net Fed sources (nM)
         rFed             ! ratio of dissolved Fe to tot Fed

      real (kind=dbl_kind), dimension(n_fep) :: &
         Fep_r       , &  ! net Fep removal (nM)
         Fep_s       , &  ! net Fep sources (nM)
         rFep             ! ratio of particulate Fe to tot Fep

      real (kind=dbl_kind) :: &
         dN        , &  ! change in N (mmol/m^3)
         dC        , &  ! change in Carbon (mmol C/m^3)
         N_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^3)
         N_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^3)
         N_r_r     , &  ! algal nitrogen losses to respiration (mmol/m^3)
         N_r_mo    , &  ! algal nitrogen losses to mortality (mmol/m^3)
         Nit_s_n   , &  ! nitrate from nitrification (mmol/m^3)
         Nit_r_p   , &  ! nitrate uptake by algae (mmol/m^3)
         Nit_s     , &  ! net nitrate sources (mmol/m^3)
         Nit_r     , &  ! net nitrate removal (mmol/m^3)
         Am_s_e    , &  ! ammonium source from excretion (mmol/m^3)
         Am_s_r    , &  ! ammonium source from respiration (mmol/m^3)
         Am_s_mo   , &  ! ammonium source from mort/remin (mmol/m^3)
         Am_r_p    , &  ! ammonium uptake by algae (mmol/m^3)
         Am_s      , &  ! net ammonium sources (mmol/m^3)
         Am_r      , &  ! net ammonium removal (mmol/m^3)
         Sil_r_p   , &  ! silicon uptake by algae (mmol/m^3)
         Sil_r     , &  ! net silicon removal (mmol/m^3)
         Fe_r_p         ! iron uptake by algae  (nM)
!        DOC_r_c   , &  ! net doc removal from bacterial consumption (mmol/m^3)
!        doc_s_m   , &  ! protein source due to algal mortality (mmol/m^3)
!        doc_s_g        ! protein source due to grazing (mmol/m^3)

      real (kind=dbl_kind) :: &
         DMSPd_s_r , &  ! skl dissolved DMSP from respiration (mmol/m^3)
         DMSPd_s_mo, &  ! skl dissolved DMSP from MBJ algal mortexc (mmol/m^3)
         DMSPd_r   , &  ! skl dissolved DMSP conversion (mmol/m^3) DMSPD_sk_r
         DMSPd_s   , &  ! net skl dissolved DMSP sources (mmol/m^3)
         DMS_s_c   , &  ! skl DMS source via conversion (mmol/m^3)
         DMS_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^3)
         DMS_s     , &  ! net skl DMS sources (mmol/m^3)
         DMS_r     , &  ! net skl DMS removal (mmol/m^3)
         Fed_tot   , &  ! total dissolved iron from all sources (nM)
         Fed_tot_r , &  ! total dissolved iron losses (nM)
         Fed_tot_s , &  ! total dissolved iron sources (nM)
         Fep_tot   , &  ! total particulate iron from all sources (nM)
!        Fep_tot_r , &  ! total particulate iron losses (nM)
         Fep_tot_s , &  ! total particulate iron sources (nM)
         Zoo_s_a   , &  ! N Losses due to zooplankton assimilation (mmol/m^3)
         Zoo_s_s   , &  ! N Losses due to grazing spillage (mmol/m^3)
         Zoo_s_m   , &  ! N Losses due to algal mortality (mmol/m^3)
         Zoo_s_b        ! N losses due to bacterial recycling of DON (mmol/m^3)

      character(len=*),parameter :: subname='(algal_dyn)'

      !-----------------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------------

       conserve_C = .true.
       Nin(:)     = c0
       Cin(:)     = c0
       chlin(:)   = c0
       DOCin(:)   = c0
       DICin(:)   = c0
       DONin(:)   = c0
       Fedin(:)   = c0
       Fepin(:)   = c0
       Nitin      = c0
       Amin       = c0
       Silin      = c0
       DMSPdin    = c0
       DMSin      = c0
       U_Am_tot   = c0
       U_Nit_tot  = c0
       U_Sil_tot  = c0
       U_Fe_tot   = c0
       U_Am_f(:)  = c0
       U_Nit_f(:) = c0
       U_Sil_f(:) = c0
       U_Fe_f(:)  = c0
       DOC_s(:)   = c0
       DOC_r(:)   = c0
       nitrif     = c0
       mort_N     = c0
       mort_C     = c0
       graze_N    = c0
       graze_C    = c0
       exude_C    = c0
       resp_N     = c0
       growth_N   = c0
       Nit_r      = c0
       Am_s       = c0
       Am_r       = c0
       Sil_r      = c0
       Fed_r(:)   = c0
       Fed_s(:)   = c0
       Fep_r(:)   = c0
       Fep_s(:)   = c0
       DMSPd_s    = c0
       dTemp      = min(T_bot-T_max,c0)
       Fed_tot    = c0
       Fed_tot_r  = c0
       Fed_tot_s  = c0
       rFed(:)    = c0
       Fep_tot    = c0
       Fep_tot_s  = c0
       rFep(:)    = c0
       DIC_r(:)   = c0
       DIC_s(:)   = c0
       Zoo        = c0

       Nitin     = ltrcrn(nlt_bgc_Nit)
       op_dep = c0
       do k = 1, n_algae
          Nin(k)   = ltrcrn(nlt_bgc_N(k))
          chlin(k) = R_chl2N(k)* Nin(k)
          op_dep = op_dep + chlabs(k)*chlin(k)
       enddo
       if (tr_bgc_C)   then
        ! do k = 1, n_algae
        !     Cin(k)=  ltrcrn(nlt_bgc_C(k))
        ! enddo
         do k = 1, n_doc
             DOCin(k)= ltrcrn(nlt_bgc_DOC(k))
         enddo
         do k = 1, n_dic
             DICin(k)= ltrcrn(nlt_bgc_DIC(k))
         enddo
       endif
       if (tr_bgc_Am)  Amin     = ltrcrn(nlt_bgc_Am)
       if (tr_bgc_Sil) Silin    = ltrcrn(nlt_bgc_Sil)
       if (tr_bgc_DMS) then
       !      DMSPpin  = ltrcrn(nlt_bgc_DMSPp)
             DMSPdin  = ltrcrn(nlt_bgc_DMSPd)
             DMSin    = ltrcrn(nlt_bgc_DMS)
       endif
!       if (tr_bgc_PON) then
!         PONin    = c0
!         PONin    = ltrcrn(nlt_bgc_PON)
!       endif
       if (tr_bgc_DON) then
         do k = 1, n_don
             DONin(k) = ltrcrn(nlt_bgc_DON(k))
         enddo
       endif
       if (tr_bgc_Fe ) then
         do k = 1, n_fed
             Fedin(k) = ltrcrn(nlt_bgc_Fed(k))
         enddo
         do k = 1, n_fep
             Fepin(k) = ltrcrn(nlt_bgc_Fep(k))
         enddo
       endif

      !-----------------------------------------------------------------------
      ! Total iron from all pools
      !-----------------------------------------------------------------------

       do k = 1,n_fed
         Fed_tot = Fed_tot + Fedin(k)
       enddo
       do k = 1,n_fep
         Fep_tot = Fep_tot + Fepin(k)
       enddo
       if (Fed_tot > puny) then
       do k = 1,n_fed
         rFed(k) = Fedin(k)/Fed_tot
       enddo
       endif
       if (Fep_tot > puny) then
       do k = 1,n_fep
         rFep(k) = Fepin(k)/Fep_tot
       enddo
       endif

      !-----------------------------------------------------------------------
      ! Light limitation  (op_dep) defined above
      !-----------------------------------------------------------------------

       if (op_dep > op_dep_min .and. .not. dEdd_algae) then
         Iavg_loc = fswthru * (c1 - exp(-op_dep)) / op_dep
       else
         Iavg_loc = fswthru
       endif

       do k = 1, n_algae
          ! With light inhibition ! Maybe include light inhibition for diatoms but phaeocystis
          L_lim = (c1 - exp(-alpha2max_low(k)*Iavg_loc)) * exp(-beta2max(k)*Iavg_loc)

          ! Without light inhibition
          ! L_lim(k) = (c1 - exp(-alpha2max_low(k)*Iavg_loc))

      !-----------------------------------------------------------------------
      ! Nutrient limitation
      !-----------------------------------------------------------------------

          Nit_lim(k) = Nitin/(Nitin + K_Nit(k))
          Am_lim(k)  = c0
          N_lim(k) = Nit_lim(k)
          if (tr_bgc_Am) then
             Am_lim(k) = Amin/(Amin + K_Am(k))
             N_lim(k)  = min(c1, Nit_lim(k) + Am_lim(k))
          endif
          Sil_lim(k) = c1
          if (tr_bgc_Sil .and. K_Sil(k) > c0) Sil_lim(k) = Silin/(Silin + K_Sil(k))

      !-----------------------------------------------------------------------
      ! Iron limitation
      !-----------------------------------------------------------------------

          Fe_lim(k) = c1
          if (tr_bgc_Fe  .and. K_Fe (k) > c0) Fe_lim (k) = Fed_tot/(Fed_tot + K_Fe(k))

      !----------------------------------------------------------------------------
      ! Growth and uptake computed within the bottom layer
      ! Note here per A93 discussions and MBJ model, salinity is a universal
      ! restriction.  Comparison with available column nutrients inserted
      ! but in tests had no effect.
      ! Primary production reverts to SE form, see MBJ below and be careful
      !----------------------------------------------------------------------------

          growmax_N(k) = mu_max(k) / secday * exp(grow_Tdep(k) * dTemp)* Nin(k) *fsal
          if (n_fed == 0) then
             grow_N(k)    = min(L_lim(k), N_lim(k), Sil_lim(k)) * growmax_N(k)
          else
             grow_N(k)    = min(L_lim(k), N_lim(k), Sil_lim(k), Fe_lim(k)) * growmax_N(k)
          endif
!         potU_Nit(k)  = Nit_lim(k)* growmax_N(k)
          potU_Am(k)   = Am_lim(k)* growmax_N(k)
          U_Am(k)      = min(grow_N(k), potU_Am(k))
          U_Nit(k)     = grow_N(k) - U_Am(k)
          U_Sil(k)     = R_Si2N(k) * grow_N(k)
          U_Fe (k)     = R_Fe2N(k) * grow_N(k)

          U_Am_tot     = U_Am_tot  + U_Am(k)
          U_Nit_tot    = U_Nit_tot + U_Nit(k)
          U_Sil_tot    = U_Sil_tot + U_Sil(k)
          U_Fe_tot     = U_Fe_tot  + U_Fe(k)
       enddo
       do k = 1, n_algae
          if (U_Am_tot > c0) U_Am_f(k) = U_Am(k)/U_Am_tot
          if (U_Nit_tot > c0) U_Nit_f(k) = U_Nit(k)/U_Nit_tot
          if (U_Sil_tot > c0) U_Sil_f(k) = U_Sil(k)/U_Sil_tot
          if (U_Fe_tot > c0) U_Fe_f(k) = U_Fe(k)/U_Fe_tot
       enddo

       if (tr_bgc_Sil) U_Sil_tot = min(U_Sil_tot, max_loss * Silin/dt)
       if (tr_bgc_Fe)  U_Fe_tot  = min(U_Fe_tot, max_loss * Fed_tot/dt)
       U_Nit_tot = min(U_Nit_tot, max_loss * Nitin/dt)
       U_Am_tot  = min(U_Am_tot,  max_loss * Amin/dt)

       do k = 1, n_algae
          U_Am(k)  = U_Am_f(k)*U_Am_tot
          U_Nit(k) = U_Nit_f(k)*U_Nit_tot
          U_Sil(k) = U_Sil_f(k)*U_Sil_tot
          U_Fe(k)  = U_Fe_f(k)*U_Fe_tot
          if (R_Si2N(k) > c0) then
             grow_N(k) = min(U_Sil(k)/R_Si2N(k),U_Nit(k) + U_Am(k), U_Fe(k)/R_Fe2N(k))
          else
             grow_N(k) = min(U_Nit(k) + U_Am(k),U_Fe(k)/R_Fe2N(k))
          endif

          fr_Am(k) = c0
          if (tr_bgc_Am) then
             fr_Am(k) = p5
             if (grow_N(k) > c0) fr_Am(k) = min(U_Am(k)/grow_N(k), c1)
          endif
          fr_Nit(k) = c1 - fr_Am(k)
          U_Nit(k)  = fr_Nit(k) * grow_N(k)
          U_Am(k)   = fr_Am(k)  * grow_N(k)
          U_Sil(k)  = R_Si2N(k) * grow_N(k)
          U_Fe (k)  = R_Fe2N(k) * grow_N(k)

      !-----------------------------------------------------------------------
      ! Define reaction terms
      !-----------------------------------------------------------------------

      ! Since the framework remains incomplete at this point,
      ! it is assumed as a starting expedient that
      ! DMSP loss to melting results in 10% conversion to DMS
      ! which is then given a ten day removal constant.
      ! Grazing losses are channeled into rough spillage and assimilation
      ! then following ammonia there is some recycling.

      !--------------------------------------------------------------------
      ! Algal reaction term
      ! v1: N_react = (grow_N*(c1 - fr_graze-fr_resp) - mort)*dt
      ! v2: N_react = (grow_N*(c1 - fr_graze * (N/graze_conc)**graze_exp-fr_resp) - mort)*dt
      !  with maximum grazing less than max_loss * Nin(k)/dt
      !--------------------------------------------------------------------

          resp(k)   = fr_resp  * grow_N(k)
          graze(k)  = min(max_loss * Nin(k)/dt, grow_N(k) * fr_graze(k) * (Nin(k)/graze_conc)**graze_exponent(k))
          mort(k)   = min(max_loss * Nin(k)/dt, &
                          mort_pre(k)*exp(mort_Tdep(k)*dTemp) * Nin(k)/secday)

        ! history variables
          grow_alg(k) = grow_N(k)
          upNOn(k) = U_Nit(k)
          upNHn(k) = U_Am(k)

          N_s_p  = grow_N(k) * dt
          N_r_g  = graze(k)  * dt
          N_r_r  = resp(k)   * dt
          N_r_mo = mort(k)   * dt
          N_s(k)    = N_s_p
          N_r(k)    = N_r_g + N_r_mo + N_r_r
          graze_N   = graze_N + graze(k)
          graze_C   = graze_C + R_C2N(k)*graze(k)
          mort_N    = mort_N + mort(k)
          mort_C    = mort_C + R_C2N(k)*mort(k)
          resp_N    = resp_N + resp(k)
          growth_N  = growth_N + grow_N(k)

      enddo ! n_algae
      !--------------------------------------------------------------------
      ! Ammonium source: algal grazing, respiration, and mortality
      !--------------------------------------------------------------------

      Am_s_e  = graze_N*(c1-fr_graze_s)*fr_graze_e*dt
      Am_s_mo = mort_N*fr_mort2min*dt
      Am_s_r  = resp_N*dt
      Am_s    = Am_s_r + Am_s_e + Am_s_mo

      !--------------------------------------------------------------------
      ! Nutrient net loss terms: algal uptake
      !--------------------------------------------------------------------

      do k = 1, n_algae
         Am_r_p  = U_Am(k)   * dt
         Am_r    = Am_r + Am_r_p
         Nit_r_p = U_Nit(k)  * dt
         Nit_r   = Nit_r + Nit_r_p
         Sil_r_p = U_Sil(k) * dt
         Sil_r   = Sil_r + Sil_r_p
         Fe_r_p  = U_Fe (k) * dt
         Fed_tot_r = Fed_tot_r + Fe_r_p
         exude_C = exude_C + k_exude(k)* R_C2N(k)*Nin(k) / secday
         DIC_r(1) = DIC_r(1) + (c1-fr_resp)*grow_N(k) * R_C2N(k) * dt
      enddo

      !--------------------------------------------------------------------
      ! nitrification
      !--------------------------------------------------------------------

       nitrification = c0
       nitrif  = k_nitrif /secday * Amin
       Am_r    = Am_r +  nitrif*dt
       Nit_s_n = nitrif * dt  ! source from NH4
       Nit_s   = Nit_s_n

      !--------------------------------------------------------------------
      ! PON:  currently using PON to shadow nitrate
      !
      ! N Losses are counted in Zoo.  These arise from mortality not
      ! remineralized (Zoo_s_m), assimilated grazing not excreted (Zoo_s_a),
      !spilled N not going to DON (Zoo_s_s) and  bacterial recycling
      ! of DON (Zoo_s_b).
      !--------------------------------------------------------------------

      if (tr_bgc_Am) then
         Zoo_s_a = graze_N*(c1-fr_graze_e)*(c1-fr_graze_s) *dt
         Zoo_s_s = graze_N*fr_graze_s*dt
         Zoo_s_m = mort_N*dt  -  Am_s_mo
      else
         Zoo_s_a = graze_N*dt*(c1-fr_graze_s)
         Zoo_s_s = graze_N*fr_graze_s*dt
         Zoo_s_m = mort_N*dt
      endif

      Zoo_s_b = c0

      !--------------------------------------------------------------------
      ! DON (n_don = 1)
      ! Proteins
      !--------------------------------------------------------------------

      DON_r(:) = c0
      DON_s(:) = c0

      if (tr_bgc_DON) then
      do n = 1, n_don
         DON_r(n) = kn_bac(n)/secday * DONin(n) * dt
         DON_s(n) =  graze_N*dt - Am_s_e + mort_N*dt - Am_s_mo
         Zoo_s_s  = Zoo_s_s - DON_s(n)
         Zoo_s_b  = Zoo_s_b + DON_r(n)*(c1-f_don_Am(n))
         Am_s     = Am_s + DON_r(n)*f_don_Am(n)
         DIC_s(1) = DIC_s(1) + DON_r(n) * R_C2N_DON(n)
      enddo
      endif

      Zoo = Zoo_s_a + Zoo_s_s + Zoo_s_m + Zoo_s_b

      !--------------------------------------------------------------------
      ! DOC
      ! polysaccharids, lipids
      !--------------------------------------------------------------------

      do n = 1, n_doc

         DOC_r(n) = k_bac(n)/secday * DOCin(n) * dt
!         DOC_s(n) = f_doc(n)*(fr_graze_s *graze_C + mort_C)*dt &
!                  + f_exude(n)*exude_C
         DOC_s(n) =  f_doc(n) * (graze_C*dt + mort_C*dt - DON_s(1) * R_C2N_DON(1))
         DIC_s(1) = DIC_s(1) + DOC_r(n)
      enddo

      !--------------------------------------------------------------------
      ! Iron sources from remineralization  (follows ammonium but reduced)
      ! only Fed_s(1)  has remineralized sources
      !--------------------------------------------------------------------

      Fed_s(1) = Fed_s(1) + Am_s * R_Fe2N(1) * fr_dFe   ! remineralization source

      !--------------------------------------------------------------------
      !  Conversion to dissolved Fe from Particulate requires DOC(1)
      !  Otherwise the only source of dFe is from remineralization
      !--------------------------------------------------------------------

      if (tr_bgc_C .and. tr_bgc_Fe) then
        if (DOCin(1) > c0) then
        !if (Fed_tot/DOCin(1) > max_dfe_doc1) then
        !  do n = 1,n_fed                                  ! low saccharid:dFe ratio leads to
        !     Fed_r_l(n)  = Fedin(n)/t_iron_conv*dt/secday ! loss of bioavailable Fe to particulate fraction
        !     Fep_tot_s   = Fep_tot_s + Fed_r_l(n)
        !     Fed_r(n)    = Fed_r_l(n)                     ! removal due to particulate scavenging
        !  enddo
        !  do n = 1,n_fep
        !     Fep_s(n) = rFep(n)* Fep_tot_s                ! source from dissolved Fe
        !  enddo
        if (Fed_tot/DOCin(1) < max_dfe_doc1) then
          do n = 1,n_fep                                  ! high saccharid:dFe ratio leads to
             Fep_r(n)  = Fepin(n)/t_iron_conv*dt/secday   ! gain of bioavailable Fe from particulate fraction
             Fed_tot_s = Fed_tot_s + Fep_r(n)
          enddo
          do n = 1,n_fed
             Fed_s(n) = Fed_s(n) + rFed(n)* Fed_tot_s     ! source from particulate Fe
          enddo
        endif
        endif !Docin(1) > c0
      endif
      if (tr_bgc_Fe) then
        do n = 1,n_fed
           Fed_r(n) = Fed_r(n) + rFed(n)*Fed_tot_r        ! scavenging + uptake
        enddo

      ! source from algal mortality/grazing and fraction of remineralized nitrogen that does
      ! not become immediately bioavailable

         do n = 1,n_fep
            Fep_s(n) = Fep_s(n) + rFep(n)* (Am_s * R_Fe2N(1) * (c1-fr_dFe))
         enddo ! losses not direct to Fed
      endif

      !--------------------------------------------------------------------
      ! Sulfur cycle begins here
      !--------------------------------------------------------------------
      ! Grazing losses are channeled into rough spillage and assimilation
      ! then onward and the MBJ mortality channel is included
      ! It is assumed as a starting expedient that
      ! DMSP loss to melting gives partial conversion to DMS in product layer
      ! which then undergoes Stefels removal.

      !--------------------------------------------------------------------
      ! DMSPd  reaction term  (DMSPd conversion is outside the algal loop)
      ! DMSPd_react = R_S2N*((fr_graze_s+fr_excrt_2S*fr_graze_e*fr_graze_a)
      !                      *fr_graze*grow_N + fr_mort2min*mort)*dt
      !             - [\DMSPd]/t_sk_conv*dt
      !--------------------------------------------------------------------
      do k = 1,n_algae
         DMSPd_s_r = fr_resp_s  * R_S2N(k) * resp(k)   * dt  !respiration fraction to DMSPd
         DMSPd_s_mo= fr_mort2min * R_S2N(k)* mort(k)   * dt  !mortality and extracellular excretion
         DMSPd_s = DMSPd_s + DMSPd_s_r + DMSPd_s_mo
      enddo
      DMSPd_r = (c1/t_sk_conv) * (c1/secday)  * (DMSPdin) * dt

      !--------------------------------------------------------------------
      ! DMS reaction term + DMSPd loss term
      ! DMS_react = ([\DMSPd]*y_sk_DMS/t_sk_conv - c1/t_sk_ox *[\DMS])*dt
      !--------------------------------------------------------------------

      DMS_s_c = y_sk_DMS * DMSPd_r
      DMS_r_o = DMSin * dt / (t_sk_ox * secday)
      DMS_s   = DMS_s_c
      DMS_r   = DMS_r_o

      !-----------------------------------------------------------------------
      ! Load reaction array
      !-----------------------------------------------------------------------

      dN = c0
      dC = c0
      do k = 1,n_algae
         reactb(nlt_bgc_N(k))  = N_s(k) - N_r(k)
         dN = dN + reactb(nlt_bgc_N(k))
         dC = dC + reactb(nlt_bgc_N(k)) * R_C2N(k)
      enddo
      if (tr_bgc_C) then
       ! do k = 1,n_algae
       !      reactb(nlt_bgc_C(k))  = R_C2N(k)*reactb(nlt_bgc_N(k))
       ! enddo
         do k = 1,n_doc
            reactb(nlt_bgc_DOC(k))= DOC_s(k) - DOC_r(k)
            dC = dC + reactb(nlt_bgc_DOC(k))
         enddo
         do k = 1,n_dic
            reactb(nlt_bgc_DIC(k))= DIC_s(k) - DIC_r(k)
            dC = dC + reactb(nlt_bgc_DIC(k))
         enddo
      endif
      reactb(nlt_bgc_Nit)   = Nit_s   - Nit_r
      nitrification = Nit_s_n
      dN = dN + reactb(nlt_bgc_Nit)
      if (tr_bgc_Am)  then
            reactb(nlt_bgc_Am)    = Am_s    - Am_r
            dN = dN + reactb(nlt_bgc_Am)
      endif
      if (tr_bgc_Sil) then
            reactb(nlt_bgc_Sil)   =  - Sil_r
      endif
      if (tr_bgc_DON) then
         do k = 1,n_don
            reactb(nlt_bgc_DON(k))= DON_s(k) - DON_r(k)
            dN = dN + reactb(nlt_bgc_DON(k))
            dC = dC + reactb(nlt_bgc_DON(k)) * R_C2N_DON(k)
         enddo
      endif
      Cerror = dC
      if (tr_bgc_Fe ) then
         do k = 1,n_fed
            reactb(nlt_bgc_Fed(k))= Fed_s (k) - Fed_r (k)
         enddo
         do k = 1,n_fep
            reactb(nlt_bgc_Fep(k))= Fep_s (k) - Fep_r (k)
         enddo
      endif
      if (tr_bgc_DMS) then
         reactb(nlt_bgc_DMSPd) = DMSPd_s - DMSPd_r
         reactb(nlt_bgc_DMS)   = DMS_s   - DMS_r
      endif

      if (tr_bgc_C) then
         if (abs(dC) > max(puny,maxval(abs(reactb(:)))*accuracy) .or. &
            abs(dN) > max(puny,maxval(abs(reactb(:)))*accuracy)) then
            conserve_C = .false.
            write(warnstr,*) subname, 'Conservation error!'
            call icepack_warnings_add(warnstr)
            if (tr_bgc_DON) then
               write(warnstr,*) subname, 'dN,DONin(1), kn_bac(1),secday,dt,n_doc'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, dN, DONin(1),kn_bac(1),secday,dt,n_doc
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'reactb(nlt_bgc_DON(1)), DON_r(1),DON_s(1)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, reactb(nlt_bgc_DON(1)),DON_r(1),DON_s(1)
               call icepack_warnings_add(warnstr)
            end if
            write(warnstr,*) subname, 'dN,secday,dt,n_doc'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, dN,secday,dt,n_doc
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'reactb(nlt_bgc_Nit),reactb(nlt_bgc_N(n_algae))'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, reactb(nlt_bgc_Nit),reactb(nlt_bgc_N(n_algae))
            call icepack_warnings_add(warnstr)
            if (tr_bgc_Am) then
               write(warnstr,*) subname, 'reactb(nlt_bgc_Am),Am_r, Am_s'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, reactb(nlt_bgc_Am),Am_r, Am_s
               call icepack_warnings_add(warnstr)
            end if
            write(warnstr,*) subname, 'dC'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, dC
            call icepack_warnings_add(warnstr)
            do k = 1,n_doc
               write(warnstr,*) subname, 'DOCin'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, DOCin(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'reactb(nlt_bgc_DOC)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, reactb(nlt_bgc_DOC(k))
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'DOC_r(k),DOC_s(k),k'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, DOC_r(k),DOC_s(k),k
               call icepack_warnings_add(warnstr)
             end do
             do k = 1,n_dic
               write(warnstr,*) subname, 'DICin'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, DICin(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'reactb(nlt_bgc_DIC)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, reactb(nlt_bgc_DIC(k))
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'DIC_r(k),DIC_s(k),k'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, DIC_r(k),DIC_s(k),k
               call icepack_warnings_add(warnstr)
            end do
            write(warnstr,*) subname, 'Zoo'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, Zoo
            call icepack_warnings_add(warnstr)
         endif
      endif

      end subroutine algal_dyn

!=======================================================================
!
! Find ice-ocean flux when ice is thin and internal dynamics/reactions are
!             assumed to be zero
!
! authors     Nicole Jeffery, LANL

      subroutine thin_ice_flux (hin, hin_old, Cin, flux_o_tot, &
                                source, dt, ocean_bio)

      real (kind=dbl_kind), dimension(nblyr+1), intent(inout) :: &
         Cin               ! initial concentration*hin_old*phin

      real (kind=dbl_kind), intent(in) :: &
         hin_old   , &     ! brine thickness (m)
         hin       , &     ! new brine thickness (m)
         dt        , &     ! time step
         source    , &     ! atm, ocean, dust flux (mmol/m^2)
         ocean_bio         ! ocean tracer concentration (mmol/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         flux_o_tot        ! tracer flux, gravity+molecular drainage flux ,
                           ! and boundary flux to ocean (mmol/m^2/s)
                           ! positive into the ocean

     ! local variables

     integer (kind=int_kind) :: &
         k                 ! vertical biology layer index

     real (kind=dbl_kind) :: &
         sum_bio,   & ! initial bio mass (mmol/m^2)
         zspace,    & ! 1/nblyr
         dC,        & ! added ocean bio mass (mmol/m^2)
         dh           ! change in thickness  (m)

     character(len=*),parameter :: subname='(thin_ice_flux)'

     zspace = c1/real(nblyr,kind=dbl_kind)

     dC = c0
     sum_bio = c0
     dh = hin-hin_old

     if (dh .le. c0 .and. hin_old > puny) then  ! keep the brine concentration fixed
       sum_bio = (Cin(1)+Cin(nblyr+1))/hin_old*zspace*p5
       Cin(1) = Cin(1)/hin_old*hin
       Cin(nblyr+1) = Cin(nblyr+1)/hin_old*hin
       do k = 2, nblyr
        sum_bio = sum_bio + Cin(k)/hin_old*zspace
        Cin(k) = Cin(k)/hin_old*hin + dC
       enddo
     else ! spread evenly in ice layers
       dC = dh*ocean_bio
       do k = 1, nblyr+1
         Cin(k) = Cin(k) + dC
       enddo
     endif

     flux_o_tot = - dh*sum_bio/dt - dC/dt + source/dt

     end subroutine thin_ice_flux

!=======================================================================
!
! Compute matrix elements for the low order solution of FEM-FCT scheme
! Predictor step
!
! July, 2014 by N. Jeffery
!
      subroutine compute_FCT_matrix  (C_in, sbdiag, dt,             &
                                      diag, spdiag, rhs,            &
                                      darcyV, dhtop, dhbot,         &
                                      iphin_N, iDin, hbri_old,      &
                                      atm_add, bphin_N,             &
                                      C_top, C_bot, Qbot, Qtop,     &
                                      Sink_bot, Sink_top,           &
                                      D_sbdiag, D_spdiag, ML)

      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         C_in            ! Initial (bulk) concentration*hbri_old (mmol/m^2)
                         ! conserved quantity on igrid

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         iDin            ! Diffusivity on the igrid (1/s)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         iphin_N         ! Porosity with min condition on igrid

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bphin_N         ! Porosity with min condition on igrid

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs         , & ! rhs of tri-diagonal matrix eqn.
         ML,           & ! lumped mass matrix
         D_spdiag, D_sbdiag     ! artificial diffusion matrix

      real (kind=dbl_kind), intent(in) :: &
         dhtop         , & ! Change in brine top (m)
         dhbot         , & ! Change in brine bottom (m)
         hbri_old      , & ! brine height (m)
         atm_add       , & ! atm-ice flux
         C_top         , & ! bulk surface source (mmol/m^2)
         C_bot         , & ! bulk bottom source (mmol/m^2)
         darcyV            ! Darcy velocity (m/s

      real (kind=dbl_kind), intent(inout) :: &   ! positive into ice
         Qtop         , & ! top  flux source (mmol/m^2/s)
         Qbot         , & ! bottom flux  source (mmol/m^2/s)
         Sink_bot     , & ! rest of bottom flux (sink or source) (mmol/m^2/s)
         Sink_top         ! rest oftop flux (sink or source) (mmol/m^2/s)

      ! local variables

      real (kind=dbl_kind) :: &
         vel, vel2, dphi_dx, zspace

      integer (kind=int_kind) :: &
         k                ! vertical index

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         Q_top,  Q_bot,              & ! surface and bottom source
         K_diag, K_spdiag, K_sbdiag, & ! advection matrix
         S_diag, S_spdiag, S_sbdiag, & ! diffusion matrix
         D_diag, iDin_phi

      real (kind=dbl_kind), dimension (nblyr) :: &
         kvector1, kvectorn1

      character(len=*),parameter :: subname='(compute_FCT_matrix)'

!---------------------------------------------------------------------
!  Diag (jj) solve for j = 1:nblyr+1
!  spdiag(j) == (j,j+1) solve for j = 1:nblyr otherwise 0
!  sbdiag(j) == (j,j-1) solve for j = 2:nblyr+1 otherwise 0
!---------------------------------------------------------------------
      kvector1(:) = c0
      kvector1(1) = c1
      kvectorn1(:) = c1
      kvectorn1(1) = c0

      zspace = c1/real(nblyr,kind=dbl_kind)
      Qbot = c0
      Qtop = c0
      Sink_bot = c0
      Sink_top = c0

! compute the lumped mass matrix

      ML(:) = zspace
      ML(1) = zspace/c2
      ML(nblyr+1) = zspace/c2

! compute matrix K: K_diag , K_sbdiag, K_spdiag
! compute matrix S: S_diag, S_sbdiag, S_spdiag

      K_diag(:) = c0
      D_diag(:) = c0
      D_spdiag(:) = c0
      D_sbdiag(:) = c0
      K_spdiag(:) = c0
      K_sbdiag(:) = c0
      S_diag(:) = c0
      S_spdiag(:) = c0
      S_sbdiag(:) = c0
      iDin_phi(:) = c0

      iDin_phi(1) = c0  !element 1
      iDin_phi(nblyr+1) = iDin(nblyr+1)/iphin_N(nblyr+1)  !outside ice at bottom
      iDin_phi(nblyr) = p5*(iDin(nblyr+1)/iphin_N(nblyr+1)+iDin(nblyr)/iphin_N(nblyr))

      vel = (bgrid(2)*dhbot - (bgrid(2)-c1)*dhtop)/dt+darcyV/bphin_N(2)
      K_diag(1) = p5*vel/hbri_old
      dphi_dx = (iphin_N(nblyr+1) - iphin_N(nblyr))/(zspace)
      vel = (bgrid(nblyr+1)*dhbot - (bgrid(nblyr+1)-c1)*dhtop)/dt  +darcyV/bphin_N(nblyr+1)
      vel = vel/hbri_old
      vel2 = (dhbot/hbri_old/dt +darcyV/hbri_old)
      K_diag(nblyr+1) =   min(c0, vel2) - iDin_phi(nblyr+1)/(zspace+ grid_o/hbri_old) &
                   + p5*(-vel + iDin_phi(nblyr)/bphin_N(nblyr+1)*dphi_dx)

      do k = 1, nblyr-1
         vel = (bgrid(k+1)*dhbot - (bgrid(k+1)-c1)*dhtop)/dt+darcyV/bphin_N(k+1)
         iDin_phi(k+1) = p5*(iDin(k+2)/iphin_N(k+2) + iDin(k+1)/iphin_N(k+1))
         dphi_dx =  (iphin_N(k+1) - iphin_N(k))/(zspace)
         K_spdiag(k)= p5*(vel/hbri_old - &
                         iDin_phi(k)/(bphin_N(k+1))*dphi_dx)

         vel = (bgrid(k+1)*dhbot - (bgrid(k+1)-c1)*dhtop)/dt  +darcyV/bphin_N(k+1)
         dphi_dx = c0
         dphi_dx = kvectorn1(k)*(iphin_N(k+1) - iphin_N(k))/(zspace)
         K_sbdiag(k+1)= -p5*(vel/hbri_old- &
                         iDin_phi(k)/bphin_N(k+1)*dphi_dx)
         K_diag(k) = K_diag(1)*kvector1(k) + (K_spdiag(k) + K_sbdiag(k))*kvectorn1(k)

         S_diag(k+1) =   -(iDin_phi(k)+ iDin_phi(k+1))/zspace
         S_spdiag(k)   = iDin_phi(k)/zspace
         S_sbdiag(k+1) = iDin_phi(k)/zspace
      enddo

      ! k = nblyr

      vel = (bgrid(nblyr+1)*dhbot - (bgrid(nblyr+1)-c1)*dhtop)/dt+darcyV/bphin_N(nblyr+1)
      dphi_dx =  (iphin_N(nblyr+1) - iphin_N(nblyr))/(zspace)
      K_spdiag(nblyr)= p5*(vel/hbri_old - &
                         iDin_phi(nblyr)/(bphin_N(nblyr+1))*dphi_dx)
      vel = (bgrid(nblyr+1)*dhbot - (bgrid(nblyr+1)-c1)*dhtop)/dt  +darcyV/bphin_N(nblyr+1)
      dphi_dx = kvectorn1(nblyr)*(iphin_N(nblyr+1) - iphin_N(nblyr))/(zspace)
      K_sbdiag(nblyr+1)= -p5*(vel/hbri_old- &
                         iDin_phi(nblyr)/bphin_N(nblyr+1)*dphi_dx)
      K_diag(nblyr) = K_spdiag(nblyr) + K_sbdiag(nblyr)
      S_diag(nblyr+1) = -iDin_phi(nblyr)/zspace
      S_spdiag(nblyr)   = iDin_phi(nblyr)/zspace
      S_sbdiag(nblyr+1) = iDin_phi(nblyr)/zspace

! compute matrix artificial D: D_spdiag, D_diag  (D_spdiag(k) = D_sbdiag(k+1))

      do k = 1,nblyr
         D_spdiag(k)    = max(-K_spdiag(k), c0, -K_sbdiag(k+1))
         D_sbdiag(k+1)  = D_spdiag(k)
      enddo
      do k = 1,nblyr+1
         D_diag(k) = D_diag(k) - D_spdiag(k) - D_sbdiag(k)
      enddo

! compute Q_top and Q_bot: top and bottom sources

      vel2 = -(dhtop/hbri_old/dt +darcyV/bphin_N(1)/hbri_old)

      Q_top(:) = c0
      Q_top(1) = max(c0,vel2*C_top + atm_add/dt)
      Qtop = Q_top(1)

      vel = (dhbot/hbri_old/dt +darcyV/hbri_old)  ! going from iphin_N(nblyr+1) to c1 makes a difference

      Q_bot(:) = c0
      Q_bot(nblyr+1) = max(c0,vel*C_bot) + iDin_phi(nblyr+1)*C_bot&
                      /(zspace + grid_o/hbri_old)

      Qbot = Q_bot(nblyr+1)

      Sink_bot = K_diag(nblyr+1) +  K_spdiag(nblyr)
      Sink_top = K_diag(1) + K_sbdiag(2)

!compute matrix elements  (1 to nblyr+1)

     spdiag = -dt * (D_spdiag + K_spdiag + S_spdiag)
     sbdiag = -dt * (D_sbdiag + K_sbdiag + S_sbdiag)
     diag = ML - dt *  (D_diag + K_diag + S_diag)
     rhs = ML * C_in + dt * Q_top + dt* Q_bot

     end subroutine compute_FCT_matrix

!=======================================================================
!
! Compute matrices for final solution FCT for passive tracers
! Corrector step
!
! July, 2014 by N. Jeffery
!
      subroutine compute_FCT_corr(C_in, C_low, dt,  &
                                  D_sbdiag, D_spdiag, ML)

      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         C_in            ! Initial (bulk) concentration*hbri_old (mmol/m^2)
                         ! conserved quantity on igrid

      real (kind=dbl_kind), dimension(nblyr+1), intent(inout) :: &
         C_low           ! Low order solution (mmol/m^2) corrected

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         D_sbdiag      , & ! sub-diagonal artificial diffusion matrix elements
         ML            , & ! Lumped mass diagonal matrix elements
         D_spdiag          ! super-diagonal artificial diffusion matrix elements

      ! local variables

      real (kind=dbl_kind) :: &
        zspace

      integer (kind=int_kind) :: &
         k                ! vertical index

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         M_spdiag, M_sbdiag,         & ! mass matrix
         F_diag, F_spdiag, F_sbdiag, & ! anti-diffusive matrix
         Pplus, Pminus,              & !
         Qplus, Qminus,              & !
         Rplus, Rminus,              & !
         a_spdiag, a_sbdiag            ! weightings of F

      character(len=*),parameter :: subname='(compute_FCT_corr)'

!---------------------------------------------------------------------
!  Diag (jj) solve for j = 1:nblyr+1
!  spdiag(j) == (j,j+1) solve for j = 1:nblyr otherwise 0
!  sbdiag(j) == (j,j-1) solve for j = 2:nblyr+1 otherwise 0
!---------------------------------------------------------------------

      zspace = c1/real(nblyr,kind=dbl_kind)

! compute the mass matrix

      M_spdiag(:) = zspace/c6
      M_spdiag(nblyr+1) = c0
      M_sbdiag(:) = zspace/c6
      M_sbdiag(1) = c0

! compute off matrix F:

      F_diag(:) = c0
      F_spdiag(:) = c0
      F_sbdiag(:) = c0

      do k = 1, nblyr
         F_spdiag(k) = M_spdiag(k)*(C_low(k)-C_in(k) - C_low(k+1)+ C_in(k+1))/dt &
                     + D_spdiag(k)*(C_low(k)-C_low(k+1))
         F_sbdiag(k+1) =  M_sbdiag(k+1)*(C_low(k+1)-C_in(k+1) - C_low(k)+ C_in(k))/dt &
                       + D_sbdiag(k+1)*(C_low(k+1)-C_low(k))

         if (F_spdiag(k)*(C_low(k) - C_low(k+1)) > c0) F_spdiag(k) = c0
         if (F_sbdiag(k+1)*(C_low(k+1) - C_low(k)) > c0) F_sbdiag(k+1) = c0
      enddo

      if (maxval(abs(F_spdiag)) > c0) then

! compute the weighting factors: a_spdiag, a_sbdiag

      a_spdiag(:) = c0
      a_sbdiag(:) = c0

      Pplus(1)  = max(c0, F_spdiag(1))
      Pminus(1) = min(c0, F_spdiag(1))
      Pplus(nblyr+1)  = max(c0, F_sbdiag(nblyr+1))
      Pminus(nblyr+1) = min(c0, F_sbdiag(nblyr+1))
      Qplus(1) = max(c0,C_low(2)-C_low(1))
      Qminus(1)= min(c0,C_low(2)-C_low(1))
      Qplus(nblyr+1) = max(c0,C_low(nblyr)-C_low(nblyr+1))
      Qminus(nblyr+1)= min(c0,C_low(nblyr)-C_low(nblyr+1))
      Rplus(1)  = min(c1, ML(1)*Qplus(1)/dt/(Pplus(1)+puny))
      Rminus(1) = min(c1, ML(1)*Qminus(1)/dt/(Pminus(1)-puny))
      Rplus(nblyr+1)  = min(c1, ML(nblyr+1)*Qplus(nblyr+1)/dt/(Pplus(nblyr+1)+puny))
      Rminus(nblyr+1) = min(c1, ML(nblyr+1)*Qminus(nblyr+1)/dt/(Pminus(nblyr+1)-puny))
      do k = 2,nblyr
         Pplus(k)  = max(c0,F_spdiag(k)) + max(c0,F_sbdiag(k))
         Pminus(k) = min(c0,F_spdiag(k)) + min(c0,F_sbdiag(k))
         Qplus(k)  = max(c0, max(C_low(k+1)-C_low(k),C_low(k-1)-C_low(k)))
         Qminus(k) = min(c0, min(C_low(k+1)-C_low(k),C_low(k-1)-C_low(k)))
         Rplus(k)  = min(c1, ML(k)*Qplus(k)/dt/(Pplus(k)+puny))
         Rminus(k) = min(c1, ML(k)*Qminus(k)/dt/(Pminus(k)-puny))
      enddo

      do k = 1, nblyr
         a_spdiag(k) = min(Rminus(k),Rplus(k+1))
         if (F_spdiag(k) > c0) a_spdiag(k) = min(Rplus(k),Rminus(k+1))
         a_sbdiag(k+1) = min(Rminus(k+1),Rplus(k))
         if (F_sbdiag(k+1) > c0) a_sbdiag(k+1) = min(Rplus(k+1),Rminus(k))
      enddo

!compute F_diag:

      F_diag(1) = a_spdiag(1)*F_spdiag(1)
      F_diag(nblyr+1) = a_sbdiag(nblyr+1)* F_sbdiag(nblyr+1)
      C_low(1) = C_low(1) + dt*F_diag(1)/ML(1)
      C_low(nblyr+1) = C_low(nblyr+1) + dt*F_diag(nblyr+1)/ML(nblyr+1)

      do k = 2,nblyr
         F_diag(k) = a_spdiag(k)*F_spdiag(k) + a_sbdiag(k)*F_sbdiag(k)
         C_low(k) = C_low(k) + dt*F_diag(k)/ML(k)
      enddo

      endif  !F_spdiag is nonzero

      end subroutine compute_FCT_corr

!=======================================================================
!
! Tridiagonal matrix solver-- for salinity
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
      subroutine tridiag_solverz(nmat,      sbdiag,   &
                                 diag,      spdiag,   &
                                 rhs,       xout)

      integer (kind=int_kind), intent(in) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (nmat), intent(in) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (nmat), intent(inout) :: &
         xout            ! solution vector

      ! local variables

      integer (kind=int_kind) :: &
         k               ! row counter

      real (kind=dbl_kind) :: &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(nmat):: &
         wgamma          ! temporary matrix variable

      character(len=*),parameter :: subname='(tridiag_solverz)'

      wbeta = diag(1)
      xout(1) = rhs(1) / wbeta

      do k = 2, nmat
         wgamma(k) = spdiag(k-1) / wbeta
         wbeta = diag(k) - sbdiag(k)*wgamma(k)
         xout(k) = (rhs(k) - sbdiag(k)*xout(k-1)) / wbeta
      enddo                     ! k

      do k = nmat-1, 1, -1
         xout(k) = xout(k) - wgamma(k+1)*xout(k+1)
      enddo                     ! k

      end subroutine tridiag_solverz

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine check_conservation_FCT (C_init, C_new, C_low, S_top, &
                                         S_bot, L_bot, L_top, dt,     &
                                         fluxbio, source)

      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         C_init        , & ! initial bulk concentration * h_old (mmol/m^2)
         C_new             ! new bulk concentration * h_new (mmol/m^2)

      real (kind=dbl_kind), dimension(nblyr+1), intent(out) :: &
         C_low             ! define low order solution = C_new

      real (kind=dbl_kind),  intent(in) :: &
         S_top         , & ! surface flux into ice (mmol/m^2/s)
         S_bot         , & ! bottom flux into ice (mmol/m^2/s)
         L_bot         , & ! remaining  bottom flux into ice (mmol/m^2/s)
         L_top         , & ! remaining  top  flux into ice (mmol/m^2/s)
         dt            , &
         source            ! nutrient source from snow and atmosphere (mmol/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         fluxbio            ! (mmol/m^2/s)  positive down (into the ocean)

      ! local variables

      integer (kind=int_kind) :: &
         k

      real (kind=dbl_kind) :: &
         diff_dt     , &
         C_init_tot  , &
         C_new_tot   , &
         zspace      , &  !1/nblyr
         accuracyC   , &  ! centered difference is Order(zspace^2)
         var_tmp          ! temporary variable

      character(len=*),parameter :: subname='(check_conservation_FCT)'

      zspace = c1/real(nblyr,kind=dbl_kind)

      !-------------------------------------
      !  Ocean flux: positive into the ocean
      !-------------------------------------
      C_init_tot = (C_init(1) + C_init(nblyr+1))*zspace*p5
      C_new_tot = (C_new(1) + C_new(nblyr+1))*zspace*p5
      C_low(1) = C_new(1)
      C_low(nblyr+1) = C_new(nblyr+1)

      do k = 2, nblyr
         C_init_tot = C_init_tot + C_init(k)*zspace
         C_new_tot = C_new_tot + C_new(k)*zspace
         C_low(k) = C_new(k)
      enddo

      accuracyC = 1.0e-11_dbl_kind*max(c1, C_init_tot, C_new_tot)
      fluxbio = fluxbio + (C_init_tot - C_new_tot + source)/dt
      diff_dt =C_new_tot - C_init_tot - (S_top+S_bot+L_bot*C_new(nblyr+1)+L_top*C_new(1))*dt

      if (minval(C_low) < c0) then
         write(warnstr,*) subname, 'Positivity of zbgc low order solution failed: C_low:',C_low
         call icepack_warnings_add(warnstr)
         !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
      endif

      if (abs(diff_dt) > accuracyC ) then
         write(warnstr,*) subname, ''
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Conservation of zbgc low order solution failed: diff_dt:'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, diff_dt
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Total initial tracer'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, C_init_tot
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Total final1  tracer'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, C_new_tot
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'bottom final tracer'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, C_new(nblyr+1)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'top final tracer'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, C_new(1)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Near bottom final tracer'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, C_new(nblyr)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Near top final tracer'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, C_new(2)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Top flux*dt into ice:'
         call icepack_warnings_add(warnstr)
         var_tmp = S_top*dt
         write(warnstr,*) subname, var_tmp
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Bottom flux*dt into ice:'
         call icepack_warnings_add(warnstr)
         var_tmp = S_bot*dt
         write(warnstr,*) subname, var_tmp
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Remaining bot flux*dt into ice:'
         call icepack_warnings_add(warnstr)
         var_tmp = L_bot*C_new(nblyr+1)*dt
         write(warnstr,*) subname, var_tmp
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'S_bot*dt + L_bot*C_new(nblyr+1)*dt'
         call icepack_warnings_add(warnstr)
         var_tmp = S_bot*dt + L_bot*C_new(nblyr+1)*dt
         write(warnstr,*) subname,  var_tmp
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'fluxbio*dt:'
         call icepack_warnings_add(warnstr)
         var_tmp = fluxbio*dt
         write(warnstr,*) subname, var_tmp
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'fluxbio:'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, fluxbio
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Remaining top flux*dt into ice:'
         call icepack_warnings_add(warnstr)
         var_tmp = L_top*C_new(1)*dt
         write(warnstr,*) subname, var_tmp
         call icepack_warnings_add(warnstr)
         !call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
      endif

      end subroutine check_conservation_FCT

!=======================================================================

! For each grid cell, sum field over all ice and snow layers
!
! author: Nicole Jeffery, LANL

      subroutine bgc_column_sum (hsnow, hbrine, xin, xout)

      real (kind=dbl_kind), dimension(nblyr+3), intent(in) :: &
         xin              ! input field

      real (kind=dbl_kind), intent(in) :: &
         hsnow, &         ! snow thickness
         hbrine           ! brine height

      real (kind=dbl_kind), intent(out) :: &
         xout             ! output field

      ! local variables

      real (kind=dbl_kind) :: &
         dzssl, &        ! snow surface layer (m)
         dzint, &        ! snow interior depth (m)
         hslyr, &        ! snow layer thickness (m)
         zspace          ! brine layer thickness/hbrine

      integer (kind=int_kind) :: &
         n                ! category/layer index

      character(len=*),parameter :: subname='(bgc_column_sum)'

      hslyr      = hsnow/real(nslyr,kind=dbl_kind)
      dzssl      = hslyr*p5
      dzint      = max(c0,hsnow - dzssl)
      zspace     = c1/real(nblyr,kind=dbl_kind)

      xout = c0
      xout = (xin(1) + xin(nblyr+1))*hbrine*p5*zspace
      do n = 2, nblyr
         xout = xout + xin(n)*zspace*hbrine
      enddo                 ! n
      xout = xout + dzssl*xin(nblyr+2) + dzint*xin(nblyr+3)

      end subroutine bgc_column_sum

!=======================================================================

! Find the total carbon concentration by summing the appropriate
! biogeochemical tracers in units of mmol C/m2
!
! author: Nicole Jeffery, LANL

      subroutine bgc_carbon_sum (hbrine, xin, xout)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         xin              ! input field, all tracers and column

      real (kind=dbl_kind), intent(in) :: &
         hbrine           ! brine height

      real (kind=dbl_kind), intent(out) :: &
         xout             ! output field  mmol/m2 carbon

      ! local variables

      real (kind=dbl_kind), dimension(nblyr+1) :: &
         zspace          ! brine layer thickness/hbrine

      integer (kind=int_kind) :: &
         n, m, iBioCount, iLayer, nBGC        ! category/layer index

      character(len=*),parameter :: subname='(bgc_carbon_sum)'

      zspace(:)  = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = p5*zspace(1)
      zspace(nblyr+1) = zspace(1)

      xout = c0

      if (tr_bgc_N) then
         iBioCount = c0
         do m = 1, n_algae
            nBGC = nt_bgc_N(1)
            do n = 1, nblyr+1
               iLayer = iBioCount + n-1
               xout = xout + xin(nBGC+iLayer)*zspace(n)*hbrine*R_C2N(m)
            enddo
            iBioCount = iBioCount + nblyr+3
         enddo
      endif
      if (tr_bgc_C) then
         iBioCount = c0
         nBGC = nt_bgc_DOC(1)
         do m = 1, n_doc
            do n = 1, nblyr+1
               iLayer = iBioCount + n-1
               xout = xout + xin(nBGC+iLayer)*zspace(n)*hbrine
            enddo
            iBioCount = iBioCount + nblyr+3
         enddo
         do m = 1, n_dic
            do n = 1, nblyr+1
               iLayer = iBioCount + n-1
               xout = xout + xin(nBGC+iLayer)*zspace(n)*hbrine
            enddo
            iBioCount = iBioCount + nblyr+3
         enddo
      endif

      if (tr_bgc_DON) then
         iBioCount = c0
         do m = 1, n_don
            nBGC = nt_bgc_DON(1)
            do n = 1, nblyr+1
               iLayer = iBioCount + n-1
               xout = xout + xin(nBGC+iLayer)*zspace(n)*hbrine*R_C2N_DON(m)
            enddo
            iBioCount = iBioCount + nblyr+3
         enddo
      endif
      if (tr_bgc_hum) then
         nBGC = nt_bgc_hum
         do n = 1, nblyr+1
            iLayer = n-1
            xout = xout + xin(nBGC+iLayer)*zspace(n)*hbrine
         enddo
      endif

      end subroutine bgc_carbon_sum

!=======================================================================

! Find the total carbon flux by summing the fluxes for the appropriate
! biogeochemical  each grid cell, sum field over all ice and snow layers
!
! author: Nicole Jeffery, LANL

      subroutine bgc_carbon_flux (flux_bio_atm, flux_bion, Tot_Carbon_flux)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         flux_bio_atm, &              ! input field, all tracers and column
         flux_bion

      real (kind=dbl_kind), intent(out) :: &
         Tot_Carbon_flux             ! output field  mmol/m2/s carbon

      ! local variables
      integer (kind=int_kind) :: &
         m        ! biology index

      character(len=*),parameter :: subname='(bgc_carbon_flux)'

      Tot_Carbon_flux = c0

      if (tr_bgc_N) then
      do m = 1, n_algae
         Tot_Carbon_flux = Tot_Carbon_flux  - (flux_bio_atm(nlt_bgc_N(m)) - flux_bion(nlt_bgc_N(m)))*R_C2N(m)
      enddo
      endif
      if (tr_bgc_C) then
      do m = 1, n_doc
         Tot_Carbon_flux = Tot_Carbon_flux - flux_bio_atm(nlt_bgc_DOC(m)) + flux_bion(nlt_bgc_DOC(m))
      enddo
      do m = 1, n_dic
         Tot_Carbon_flux = Tot_Carbon_flux - flux_bio_atm(nlt_bgc_DIC(m)) + flux_bion(nlt_bgc_DIC(m))
      enddo
      endif
      if (tr_bgc_DON) then
      do m = 1, n_don
         Tot_Carbon_flux = Tot_Carbon_flux - (flux_bio_atm(nlt_bgc_DON(m)) - flux_bion(nlt_bgc_DON(m)))*R_C2N_DON(m)
      enddo
      endif
      if (tr_bgc_hum) &
         Tot_Carbon_flux = Tot_Carbon_flux - flux_bio_atm(nlt_bgc_hum) + flux_bion(nlt_bgc_hum)

      end subroutine bgc_carbon_flux

!=======================================================================

      end module icepack_algae

!=======================================================================
