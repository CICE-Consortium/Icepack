!=========================================================================
!
! Initialization routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
      module icedrv_init_column

      use icedrv_kinds
      use icedrv_domain_size, only: ncat, nx, ncat
      use icedrv_domain_size, only: max_ntrcr, nblyr, nilyr, nslyr
      use icedrv_domain_size, only: n_algae, n_zaero, n_doc, n_dic, n_don
      use icedrv_domain_size, only: n_fed, n_fep, max_nsw, n_bgc, n_aero
      use icedrv_constants, only: c1, c2, p5, c0, p1
      use icedrv_constants, only: nu_diag, nu_nml
      use icepack_intfc, only: icepack_max_don, icepack_max_doc, icepack_max_dic
      use icepack_intfc, only: icepack_max_algae, icepack_max_aero, icepack_max_fe
      use icepack_intfc, only: icepack_max_nbtrcr
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_init_tracer_sizes, icepack_init_tracer_flags
      use icepack_intfc, only: icepack_init_tracer_indices
      use icepack_intfc, only: icepack_init_parameters
      use icepack_intfc, only: icepack_query_tracer_sizes, icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_init_zbgc
      use icepack_intfc, only: icepack_init_thermo
      use icepack_intfc, only: icepack_step_radiation, icepack_init_orbit
      use icepack_intfc, only: icepack_init_bgc, icepack_init_zsalinity
      use icepack_intfc, only: icepack_init_ocean_bio, icepack_load_ocean_bio_array
      use icepack_intfc, only: icepack_init_hbrine
      use icedrv_system, only: icedrv_system_abort

      implicit none

      private
      public :: init_thermo_vertical, init_shortwave, &
                init_bgc, init_hbrine, init_zbgc

!=======================================================================

      contains

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_thermo_vertical

      use icedrv_flux, only: salinz, Tmltz

      integer (kind=int_kind) :: &
         i,          &  ! horizontal indices
         k              ! ice layer index

      real (kind=dbl_kind), dimension(nilyr+1) :: &
         sprofile                         ! vertical salinity profile

      real (kind=dbl_kind) :: &
         depressT

      character(len=*), parameter :: subname='(init_thermo_vertical)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(depressT_out=depressT)
      call icepack_init_thermo(nilyr=nilyr, sprofile=sprofile)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      do i = 1, nx
         do k = 1, nilyr+1
            salinz(i,k) = sprofile(k)
            Tmltz (i,k) = -salinz(i,k)*depressT
         enddo ! k
      enddo    ! i

      end subroutine init_thermo_vertical

!=======================================================================
!
!  Initialize shortwave

      subroutine init_shortwave

      use icedrv_arrays_column, only: fswpenln, Iswabsn, Sswabsn, albicen
      use icedrv_arrays_column, only: albsnon, alvdrn, alidrn, alvdfn, alidfn
      use icedrv_arrays_column, only: fswsfcn, fswthrun, ffracn, snowfracn
      use icedrv_arrays_column, only: fswintn, albpndn, apeffn, trcrn_sw, dhsn
      use icedrv_arrays_column, only: kaer_tab, waer_tab, gaer_tab
      use icedrv_arrays_column, only: kaer_bc_tab, waer_bc_tab, gaer_bc_tab
      use icedrv_arrays_column, only: swgrid, igrid, bcenh
      use icedrv_calendar, only: istep1, dt, calendar_type
      use icedrv_calendar, only:    days_per_year, nextsw_cday, yday, sec
      use icedrv_system, only: icedrv_system_abort
      use icedrv_flux, only: alvdf, alidf, alvdr, alidr
      use icedrv_flux, only: alvdr_ai, alidr_ai, alvdf_ai, alidf_ai
      use icedrv_flux, only: swvdr, swvdf, swidr, swidf, scale_factor, snowfrac
      use icedrv_flux, only: albice, albsno, albpnd, apeff_ai, coszen, fsnow
      use icedrv_init, only: tlat, tlon, tmask
      use icedrv_restart_shared, only: restart
      use icedrv_state, only: aicen, vicen, vsnon, trcrn

      integer (kind=int_kind) :: &
         i, k         , & ! horizontal indices
         n                ! thickness category index

      real (kind=dbl_kind) :: &
         netsw            ! flag for shortwave radiation presence

      logical (kind=log_kind) :: &
         l_print_point, & ! flag to print designated grid point diagnostics
         dEdd_algae,    & ! from icepack
         modal_aero       ! from icepack

      character (len=char_len) :: &
         shortwave        ! from icepack

      real (kind=dbl_kind), dimension(ncat) :: &
         fbri             ! brine height to ice thickness

      real (kind=dbl_kind), allocatable, dimension(:,:) :: &
         ztrcr_sw

      logical (kind=log_kind) :: tr_brine, tr_zaero, tr_bgc_N
      integer (kind=int_kind) :: nt_alvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero, &
         nt_fbri, nt_tsfc, ntrcr, nbtrcr_sw, nlt_chl_sw
      integer (kind=int_kind), dimension(icepack_max_aero) :: nlt_zaero_sw
      integer (kind=int_kind), dimension(icepack_max_aero) :: nt_zaero
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_N
      real (kind=dbl_kind) :: puny

      character(len=*), parameter :: subname='(init_shortwave)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_parameters(puny_out=puny)
         call icepack_query_parameters(shortwave_out=shortwave)
         call icepack_query_parameters(dEdd_algae_out=dEdd_algae)
         call icepack_query_parameters(modal_aero_out=modal_aero)
         call icepack_query_tracer_sizes(ntrcr_out=ntrcr, &
              nbtrcr_sw_out=nbtrcr_sw)
         call icepack_query_tracer_flags(tr_brine_out=tr_brine, &
              tr_zaero_out=tr_zaero, tr_bgc_N_out=tr_bgc_N)
         call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, &
              nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
              nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, &
              nt_fbri_out=nt_fbri, nt_tsfc_out=nt_tsfc, &
              nt_bgc_N_out=nt_bgc_N, nt_zaero_out=nt_zaero, &
              nlt_chl_sw_out=nlt_chl_sw, nlt_zaero_sw_out=nlt_zaero_sw)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         allocate(ztrcr_sw(nbtrcr_sw, ncat))

         fswpenln(:,:,:) = c0
         Iswabsn(:,:,:) = c0
         Sswabsn(:,:,:) = c0

         do i = 1, nx

            l_print_point = .false.

            alvdf(i) = c0
            alidf(i) = c0
            alvdr(i) = c0
            alidr(i) = c0
            alvdr_ai(i) = c0
            alidr_ai(i) = c0
            alvdf_ai(i) = c0
            alidf_ai(i) = c0
            albice(i) = c0
            albsno(i) = c0
            albpnd(i) = c0
            snowfrac(i) = c0
            apeff_ai(i) = c0

            do n = 1, ncat
               alvdrn(i,n) = c0
               alidrn(i,n) = c0
               alvdfn(i,n) = c0
               alidfn(i,n) = c0
               fswsfcn(i,n) = c0
               fswintn(i,n) = c0
               fswthrun(i,n) = c0
            enddo   ! ncat

         enddo

         do i = 1, nx

            if (trim(shortwave) == 'dEdd') then ! delta Eddington

               ! initialize orbital parameters
               ! These come from the driver in the coupled model.
               call icepack_warnings_flush(nu_diag)
               call icepack_init_orbit()
               call icepack_warnings_flush(nu_diag)
               if (icepack_warnings_aborted()) &
                  call icedrv_system_abort(i, istep1, subname, __FILE__, __LINE__)
            endif

            fbri(:) = c0
            ztrcr_sw(:,:) = c0
            do n = 1, ncat
               if (tr_brine)  fbri(n) = trcrn(i,nt_fbri,n)
            enddo

            if (tmask(i)) then
            call icepack_step_radiation (              &
                         dt=dt,           ncat=ncat,           &
                         nblyr=nblyr,                          &
                         nilyr=nilyr,     nslyr=nslyr,         &
                         dEdd_algae=dEdd_algae,                &
                         swgrid=swgrid(:),                     &
                         igrid=igrid(:),                       &
                         fbri=fbri(:),                         &
                         aicen=aicen(i,:),                     &
                         vicen=vicen(i,:),                     &
                         vsnon=vsnon(i,:),                     &
                         Tsfcn=trcrn(i,nt_Tsfc,:),             &
                         alvln=trcrn(i,nt_alvl,:),             &
                         apndn=trcrn(i,nt_apnd,:),             &
                         hpndn=trcrn(i,nt_hpnd,:),             &
                         ipndn=trcrn(i,nt_ipnd,:),             &
                         aeron=trcrn(i,nt_aero:nt_aero+4*n_aero-1,:), &
                         bgcNn=trcrn(i,nt_bgc_N(1):nt_bgc_N(1)+n_algae*(nblyr+3)-1,:), &
                         zaeron=trcrn(i,nt_zaero(1):nt_zaero(1)+n_zaero*(nblyr+3)-1,:), &
                         trcrn_bgcsw=ztrcr_sw,                 &
                         TLAT=TLAT(i), TLON=TLON(i),           &
                         calendar_type=calendar_type,          &
                         days_per_year=days_per_year,          &
                         nextsw_cday=nextsw_cday, yday=yday, sec=sec,      &
                         kaer_tab=kaer_tab, kaer_bc_tab=kaer_bc_tab(:,:),  &
                         waer_tab=waer_tab, waer_bc_tab=waer_bc_tab(:,:),  &
                         gaer_tab=gaer_tab, gaer_bc_tab=gaer_bc_tab(:,:),  &
                         bcenh=bcenh(:,:,:),                               &
                         modal_aero=modal_aero,                            &
                         swvdr=swvdr(i),         swvdf=swvdf(i),           &
                         swidr=swidr(i),         swidf=swidf(i),           &
                         coszen=coszen(i),       fsnow=fsnow(i),           &
                         alvdrn=alvdrn(i,:),     alvdfn=alvdfn(i,:),       &
                         alidrn=alidrn(i,:),     alidfn=alidfn(i,:),       &
                         fswsfcn=fswsfcn(i,:),   fswintn=fswintn(i,:),     &
                         fswthrun=fswthrun(i,:), fswpenln=fswpenln(i,:,:), &
                         Sswabsn=Sswabsn(i,:,:), Iswabsn=Iswabsn(i,:,:),   &
                         albicen=albicen(i,:),   albsnon=albsnon(i,:),     &
                         albpndn=albpndn(i,:),   apeffn=apeffn(i,:),       &
                         snowfracn=snowfracn(i,:),                         &
                         dhsn=dhsn(i,:),         ffracn=ffracn(i,:),       &
                         l_print_point=l_print_point,                      &
                         initonly = .true.)
            endif
            call icepack_warnings_flush(nu_diag)
            if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
               file=__FILE__, line=__LINE__)
         
      !-----------------------------------------------------------------
      ! Define aerosol tracer on shortwave grid
      !-----------------------------------------------------------------

            if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
               do n = 1, ncat
                  do k = 1, nbtrcr_sw
                     trcrn_sw(i,k,n) = ztrcr_sw(k,n)
                  enddo
               enddo
            endif

         enddo

      !-----------------------------------------------------------------
      ! Aggregate albedos 
      !-----------------------------------------------------------------

         do n = 1, ncat
            do i = 1, nx
               
               if (aicen(i,n) > puny) then
                  
                  alvdf(i) = alvdf(i) + alvdfn(i,n)*aicen(i,n)
                  alidf(i) = alidf(i) + alidfn(i,n)*aicen(i,n)
                  alvdr(i) = alvdr(i) + alvdrn(i,n)*aicen(i,n)
                  alidr(i) = alidr(i) + alidrn(i,n)*aicen(i,n)
                  
                  netsw = swvdr(i) + swidr(i) + swvdf(i) + swidf(i)
                  if (netsw > puny) then ! sun above horizon
                     albice(i) = albice(i) + albicen(i,n)*aicen(i,n)
                     albsno(i) = albsno(i) + albsnon(i,n)*aicen(i,n)
                     albpnd(i) = albpnd(i) + albpndn(i,n)*aicen(i,n)
                  endif
                  
                  apeff_ai(i) = apeff_ai(i) + apeffn(i,n)*aicen(i,n)
                  snowfrac(i) = snowfrac(i) + snowfracn(i,n)*aicen(i,n)
               
               endif ! aicen > puny
            enddo  ! i
         enddo   ! ncat

         do i = 1, nx

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i) = alvdf  (i)
            alidf_ai  (i) = alidf  (i)
            alvdr_ai  (i) = alvdr  (i)
            alidr_ai  (i) = alidr  (i)
            
      !----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !----------------------------------------------------------------
            if (.not. restart) then
               scale_factor(i) = &
                       swvdr(i)*(c1 - alvdr_ai(i)) &
                     + swvdf(i)*(c1 - alvdf_ai(i)) &
                     + swidr(i)*(c1 - alidr_ai(i)) &
                     + swidf(i)*(c1 - alidf_ai(i))
            endif

         enddo ! i

         deallocate(ztrcr_sw)

      end subroutine init_shortwave

!=======================================================================

!  Initialize vertical profile for biogeochemistry

      subroutine init_bgc() 

      use icedrv_arrays_column, only: zfswin, trcrn_sw
      use icedrv_arrays_column, only: ocean_bio_all, ice_bio_net, snow_bio_net
      use icedrv_arrays_column, only: cgrid, igrid, bphi, iDi, bTiz, iki
      use icedrv_arrays_column, only: Rayleigh_criteria, Rayleigh_real
      use icedrv_calendar,  only: istep1
      use icedrv_system, only: icedrv_system_abort
      use icedrv_flux, only: sss, nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use icedrv_forcing_bgc, only:  get_forcing_bgc
      use icedrv_state, only: trcrn

      ! local variables

      integer (kind=int_kind) :: &
         i                , & ! horizontal indices
         k                , & ! vertical index
         n                    ! category index

      integer (kind=int_kind) :: &
         max_nbtrcr, max_algae, max_don, max_doc, max_dic, max_aero, max_fe

      logical (kind=log_kind) :: &
         RayleighC, &
         solve_zsal

      real(kind=dbl_kind) :: &
         RayleighR

      real(kind=dbl_kind), dimension(max_ntrcr,ncat) :: &
         trcrn_bgc 
      
      real(kind=dbl_kind), dimension(nilyr,ncat) :: &
         sicen    

      integer (kind=int_kind) :: &
         nbtrcr, ntrcr, ntrcr_o, &
         nt_sice, nt_bgc_S

      character(len=*), parameter :: subname='(init_bgc)'

      ! Initialize

      bphi(:,:,:) = c0   ! initial porosity for no ice 
      iDi (:,:,:) = c0   ! interface diffusivity
      bTiz(:,:,:) = c0   ! initial bio grid ice temperature
      iki (:,:,:) = c0   ! permeability

      ocean_bio_all(:,:)   = c0
      ice_bio_net  (:,:)   = c0 ! integrated ice tracer conc (mmol/m^2 or mg/m^2) 
      snow_bio_net (:,:)   = c0 ! integrated snow tracer conc (mmol/m^2 or mg/m^2)
      zfswin       (:,:,:) = c0 ! shortwave flux on bio grid
      trcrn_sw     (:,:,:) = c0 ! tracers active in the shortwave calculation
      trcrn_bgc    (:,:)   = c0

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(solve_zsal_out=solve_zsal)
      call icepack_query_tracer_sizes(max_nbtrcr_out=max_nbtrcr, &
           max_algae_out=max_algae, max_don_out=max_don, max_doc_out=max_doc, &
           max_dic_out=max_dic, max_aero_out=max_aero, max_fe_out=max_fe)
      call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr, ntrcr_out=ntrcr, &
           ntrcr_o_out=ntrcr_o)
      call icepack_query_tracer_indices(nt_sice_out=nt_sice, nt_bgc_S_out=nt_bgc_S)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! zsalinity initialization
      !-----------------------------------------------------------------
      
      if (solve_zsal) then
         do i = 1, nx
            call icepack_init_zsalinity(ncat=ncat, nblyr=nblyr,        &
                                        ntrcr_o           = ntrcr_o,   &
                                        Rayleigh_criteria = RayleighC, &
                                        Rayleigh_real     = RayleighR, &
                                        trcrn_bgc         = trcrn_bgc, &
                                        nt_bgc_S          = nt_bgc_S,  &
                                        sss               = sss(i))
            Rayleigh_real    (i) = RayleighR
            Rayleigh_criteria(i) = RayleighC
            do n = 1,ncat
               do k  = 1, nblyr
                  trcrn(i,nt_bgc_S+k-1,n) = trcrn_bgc(nt_bgc_S-1+k-ntrcr_o,n)
               enddo
            enddo
         enddo      ! i
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__, line=__LINE__)
      endif ! solve_zsal

      !-----------------------------------------------------------------
      ! biogeochemistry initialization
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Initial Ocean Values if not coupled to the ocean bgc
      !-----------------------------------------------------------------
      do i = 1, nx
         call icepack_init_ocean_bio ( &
                      amm=amm(i),   dmsp=dmsp(i),  dms=dms(i),   &
                      doc=doc(i,:), dic =dic(i,:), don=don(i,:), &
                      fed=fed(i,:), fep =fep(i,:), hum=hum(i),   &
                      nit=nit(i),   sil =sil(i),                 &
                      zaeros=zaeros(i,:),                        &
                      algalN=algalN(i,:),                        &
                      max_dic=max_dic, max_don =max_don,         &
                      max_fe =max_fe,  max_aero=max_aero)
      enddo  ! i
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      call get_forcing_bgc                          ! defines nit and sil

      do i = 1, nx

         do n = 1, ncat
            do k = 1, nilyr
               sicen(k,n) = trcrn(i,nt_sice+k-1,n)
            enddo
            do k = ntrcr_o+1, ntrcr
               trcrn_bgc(k-ntrcr_o,n) = trcrn(i,k,n)
            enddo
         enddo

         call icepack_load_ocean_bio_array(max_nbtrcr=max_nbtrcr,             &
                      max_algae=max_algae, max_don=max_don,  max_doc=max_doc, &
                      max_aero =max_aero,  max_dic=max_dic,  max_fe =max_fe,  &
                      nit =nit(i),   amm=amm(i),   sil   =sil(i),             &
                      dmsp=dmsp(i),  dms=dms(i),   algalN=algalN(i,:),        &
                      doc =doc(i,:), don=don(i,:), dic   =dic(i,:),           &  
                      fed =fed(i,:), fep=fep(i,:), zaeros=zaeros(i,:),        &
                      ocean_bio_all=ocean_bio_all(i,:),  hum=hum(i))
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(i, istep1, subname, &
             __FILE__, __LINE__)

      enddo  ! i

      do i = 1, nx
         call icepack_init_bgc(ncat=ncat, nblyr=nblyr, nilyr=nilyr, ntrcr_o=ntrcr_o, &
                      cgrid=cgrid, igrid=igrid, ntrcr=ntrcr, nbtrcr=nbtrcr,          &
                      sicen=sicen(:,:), trcrn=trcrn_bgc(:,:),                        &
                      sss=sss(i), ocean_bio_all=ocean_bio_all(i,:))
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(i, istep1, subname, &
             __FILE__, __LINE__)
      enddo  ! i

      end subroutine init_bgc

!=======================================================================

!  Initialize brine height tracer

      subroutine init_hbrine()

      use icedrv_arrays_column, only: first_ice, bgrid, igrid, cgrid
      use icedrv_arrays_column, only: icgrid, swgrid
      use icedrv_state, only: trcrn

      real (kind=dbl_kind) :: phi_snow
      integer (kind=int_kind) :: nt_fbri
      logical (kind=log_kind) :: tr_brine
      character(len=*), parameter :: subname='(init_hbrine)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(phi_snow_out=phi_snow)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------

      call icepack_init_hbrine(bgrid=bgrid, igrid=igrid, cgrid=cgrid, icgrid=icgrid, &
           swgrid=swgrid, nblyr=nblyr, nilyr=nilyr, phi_snow=phi_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_init_parameters(phi_snow_in=phi_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      first_ice(:,:) = .true.            
      if (tr_brine) trcrn(:,nt_fbri,:) = c1

      end subroutine init_hbrine

!=======================================================================

! Namelist variables, set to default values; may be altered at run time
! 
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine init_zbgc

      use icedrv_state, only: trcr_base, trcr_depend, n_trcr_strata
      use icedrv_state, only: nt_strata
      use icedrv_forcing, only: bgc_data_type

      character (len=char_len) :: &
         shortwave        ! from icepack

      integer (kind=int_kind) :: &
         ntrcr,         nbtrcr,       nbtrcr_sw,    &
         ntrcr_o,       nt_fbri,      &
         nt_bgc_Nit,    nt_bgc_Am,    nt_bgc_Sil,   &
         nt_bgc_DMS,    nt_bgc_PON,   nt_bgc_S,     &
         nt_bgc_DMSPp,  nt_bgc_DMSPd, &
         nt_zbgc_frac,  nlt_chl_sw,   &
         nlt_bgc_Nit,   nlt_bgc_Am,   nlt_bgc_Sil,  &
         nlt_bgc_DMS,   nlt_bgc_DMSPp,nlt_bgc_DMSPd,&
         nlt_bgc_PON,   nt_bgc_hum,   nlt_bgc_hum

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero_sw       ! points to aerosol in trcrn_sw

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nlt_bgc_N      , & ! algae
         nlt_bgc_C      , & !
         nlt_bgc_chl

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nlt_bgc_DOC        ! disolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nlt_bgc_DON        !

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nlt_bgc_DIC        ! disolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nlt_bgc_Fed    , & !
         nlt_bgc_Fep        !

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero          ! non-reacting layer aerosols

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nt_bgc_DOC         !  dissolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nt_bgc_Fed     , & !  dissolved iron
         nt_bgc_Fep         !  particulate iron

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero       , & !  black carbon and other aerosols
         nt_zaero_sw        !  for shortwave calculation

      integer (kind=int_kind), dimension(icepack_max_nbtrcr) :: &
         bio_index_o        ! relates nlt_bgc_NO to ocean concentration index

      integer (kind=int_kind), dimension(icepack_max_nbtrcr) :: &
         bio_index          ! relates bio indices nlt_bgc_N to nt_bgc_N

      logical (kind=log_kind) :: &
          tr_brine, &
          tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
          tr_bgc_DMS,    tr_bgc_PON,                 &
          tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
          tr_bgc_DON,    tr_bgc_Fe,    tr_zaero,     &
          tr_bgc_hum,    tr_aero
 
      integer (kind=int_kind) :: &
          ktherm

      logical (kind=log_kind) :: &
          solve_zsal, skl_bgc, z_tracers, scale_bgc, solve_zbgc, dEdd_algae, &
          modal_aero, restore_bgc

      character (char_len) :: &
          bgc_flux_type

      real (kind=dbl_kind) :: &
          grid_o, l_sk, initbio_frac, &
          frazil_scav, grid_oS, l_skS, &
          phi_snow, &
          ratio_Si2N_diatoms , ratio_Si2N_sp      , ratio_Si2N_phaeo   ,  &
          ratio_S2N_diatoms  , ratio_S2N_sp       , ratio_S2N_phaeo    ,  &
          ratio_Fe2C_diatoms , ratio_Fe2C_sp      , ratio_Fe2C_phaeo   ,  &
          ratio_Fe2N_diatoms , ratio_Fe2N_sp      , ratio_Fe2N_phaeo   ,  &
          ratio_Fe2DON       , ratio_Fe2DOC_s     , ratio_Fe2DOC_l     ,  &
          fr_resp            , tau_min            , tau_max            ,  &
          algal_vel          , R_dFe2dust         , dustFe_sol         ,  &
          chlabs_diatoms     , chlabs_sp          , chlabs_phaeo       ,  &
          alpha2max_low_diatoms,alpha2max_low_sp  , alpha2max_low_phaeo,  &
          beta2max_diatoms   , beta2max_sp        , beta2max_phaeo     ,  &
          mu_max_diatoms     , mu_max_sp          , mu_max_phaeo       ,  &
          grow_Tdep_diatoms  , grow_Tdep_sp       , grow_Tdep_phaeo    ,  &
          fr_graze_diatoms   , fr_graze_sp        , fr_graze_phaeo     ,  &
          mort_pre_diatoms   , mort_pre_sp        , mort_pre_phaeo     ,  &
          mort_Tdep_diatoms  , mort_Tdep_sp       , mort_Tdep_phaeo    ,  &
          k_exude_diatoms    , k_exude_sp         , k_exude_phaeo      ,  &
          K_Nit_diatoms      , K_Nit_sp           , K_Nit_phaeo        ,  &
          K_Am_diatoms       , K_Am_sp            , K_Am_phaeo         ,  &
          K_Sil_diatoms      , K_Sil_sp           , K_Sil_phaeo        ,  &
          K_Fe_diatoms       , K_Fe_sp            , K_Fe_phaeo         ,  &
          f_don_protein      , kn_bac_protein     , f_don_Am_protein   ,  &
          f_doc_s            , f_doc_l            , f_exude_s          ,  &
          f_exude_l          , k_bac_s            , k_bac_l            ,  &
          T_max              , fsal               , op_dep_min         ,  &
          fr_graze_s         , fr_graze_e         , fr_mort2min        ,  &
          fr_dFe             , k_nitrif           , t_iron_conv        ,  &
          max_loss           , max_dfe_doc1       , fr_resp_s          ,  &
          y_sk_DMS           , t_sk_conv          , t_sk_ox            ,  &
          algaltype_diatoms  , algaltype_sp       , algaltype_phaeo    ,  &
          nitratetype        , ammoniumtype       , silicatetype       ,  &
          dmspptype          , dmspdtype          , humtype            ,  &
          doctype_s          , doctype_l          , dontype_protein    ,  &
          fedtype_1          , feptype_1          , zaerotype_bc1      ,  &
          zaerotype_bc2      , zaerotype_dust1    , zaerotype_dust2    ,  &
          zaerotype_dust3    , zaerotype_dust4    , ratio_C2N_diatoms  ,  &
          ratio_C2N_sp       , ratio_C2N_phaeo    , ratio_chl2N_diatoms,  & 
          ratio_chl2N_sp     , ratio_chl2N_phaeo  , F_abs_chl_diatoms  ,  &
          F_abs_chl_sp       , F_abs_chl_phaeo    , ratio_C2N_proteins 

      real (kind=dbl_kind), dimension(icepack_max_dic) :: &
         dictype

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         algaltype   ! tau_min for both retention and release

      real (kind=dbl_kind), dimension(icepack_max_doc) :: &
         doctype

      real (kind=dbl_kind), dimension(icepack_max_don) :: &
         dontype

      real (kind=dbl_kind), dimension(icepack_max_fe) :: &
         fedtype

      real (kind=dbl_kind), dimension(icepack_max_fe) :: &
         feptype

      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         zaerotype

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         R_C2N     ,      & ! algal C to N (mole/mole)
         R_chl2N   ,      & ! 3 algal chlorophyll to N (mg/mmol)
         F_abs_chl          ! to scale absorption in Dedd

      real (kind=dbl_kind), dimension(icepack_max_don) :: &  ! increase compare to algal R_Fe2C
         R_C2N_DON

       real (kind=dbl_kind),  dimension(icepack_max_algae) :: &
         R_Si2N     , & ! algal Sil to N (mole/mole) 
         R_S2N      , & ! algal S to N (mole/mole)
         ! Marchetti et al 2006, 3 umol Fe/mol C for iron limited Pseudo-nitzschia
         R_Fe2C     , & ! algal Fe to carbon (umol/mmol)
         R_Fe2N         ! algal Fe to N (umol/mmol)

      real (kind=dbl_kind), dimension(icepack_max_don) :: & 
         R_Fe2DON       ! Fe to N of DON (nmol/umol)

      real (kind=dbl_kind), dimension(icepack_max_doc) :: &  
         R_Fe2DOC       ! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         chlabs           , & ! chla absorption 1/m/(mg/m^3)
         alpha2max_low    , & ! light limitation (1/(W/m^2))
         beta2max         , & ! light inhibition (1/(W/m^2))
         mu_max           , & ! maximum growth rate (1/d)
         grow_Tdep        , & ! T dependence of growth (1/C)
         fr_graze         , & ! fraction of algae grazed
         mort_pre         , & ! mortality (1/day)
         mort_Tdep        , & ! T dependence of mortality (1/C)
         k_exude          , & ! algal carbon  exudation rate (1/d)
         K_Nit            , & ! nitrate half saturation (mmol/m^3) 
         K_Am             , & ! ammonium half saturation (mmol/m^3) 
         K_Sil            , & ! silicon half saturation (mmol/m^3)
         K_Fe                 ! iron half saturation  or micromol/m^3
            
      real (kind=dbl_kind), dimension(icepack_max_DON) :: &
         f_don            , & ! fraction of spilled grazing to DON
         kn_bac           , & ! Bacterial degredation of DON (1/d)
         f_don_Am             ! fraction of remineralized DON to Am

      real (kind=dbl_kind), dimension(icepack_max_DOC) :: &
         f_exude          , & ! fraction of exuded carbon to each DOC pool
         k_bac                ! Bacterial degredation of DOC (1/d)    

      real (kind=dbl_kind), dimension(icepack_max_nbtrcr) :: &
         zbgc_frac_init,&! initializes mobile fraction
         bgc_tracer_type ! described tracer in mobile or stationary phases      
                         ! < 0 is purely mobile (eg. nitrate)
                         ! > 0 has timescales for transitions between 
                         ! phases based on whether the ice is melting or growing

     real (kind=dbl_kind), dimension(icepack_max_nbtrcr) :: &
         zbgc_init_frac, &   ! fraction of ocean tracer  concentration in new ice
         tau_ret,        &   ! retention timescale  (s), mobile to stationary phase
         tau_rel             ! release timescale    (s), stationary to mobile phase

      character (32) :: &
         nml_filename = 'icepack_in' ! namelist input file name

      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        k, mm    , & ! loop index
        ntd      , & ! for tracer dependency calculation
        nk       , & !
        nt_depend

      character(len=*), parameter :: subname='(init_zbgc)'

      !------------------------------------------------------------
      !        Tracers have mobile and  stationary phases. 
      ! ice growth allows for retention, ice melt facilitates mobility
      ! bgc_tracer_type defines the exchange timescales between these phases
      ! -1 : entirely in the mobile phase, no exchange  (this is the default)
      !  0 : retention time scale is tau_min, release time scale is tau_max
      !  1 : retention time scale is tau_max, release time scale is tau_min
      ! 0.5: retention time scale is tau_min, release time scale is tau_min
      !  2 : retention time scale is tau_max, release time scale is tau_max
      !------------------------------------------------------------

      !-----------------------------------------------------------------
      ! namelist variables
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
        tr_brine, tr_zaero, modal_aero, skl_bgc, &
        z_tracers, dEdd_algae, solve_zbgc, bgc_flux_type, &
        restore_bgc, scale_bgc, solve_zsal, bgc_data_type, &
        tr_bgc_Nit, tr_bgc_C, tr_bgc_chl, tr_bgc_Am, tr_bgc_Sil, &
        tr_bgc_DMS, tr_bgc_PON, tr_bgc_hum, tr_bgc_DON, tr_bgc_Fe, &
        grid_o, l_sk, grid_oS, &
        l_skS, phi_snow,  initbio_frac, frazil_scav, &
        ratio_Si2N_diatoms , ratio_Si2N_sp      , ratio_Si2N_phaeo   ,  &
        ratio_S2N_diatoms  , ratio_S2N_sp       , ratio_S2N_phaeo    ,  &
        ratio_Fe2C_diatoms , ratio_Fe2C_sp      , ratio_Fe2C_phaeo   ,  &
        ratio_Fe2N_diatoms , ratio_Fe2N_sp      , ratio_Fe2N_phaeo   ,  &
        ratio_Fe2DON       , ratio_Fe2DOC_s     , ratio_Fe2DOC_l     ,  &
        fr_resp            , tau_min            , tau_max            ,  &
        algal_vel          , R_dFe2dust         , dustFe_sol         ,  &
        chlabs_diatoms     , chlabs_sp          , chlabs_phaeo       ,  &
        alpha2max_low_diatoms,alpha2max_low_sp  , alpha2max_low_phaeo,  &
        beta2max_diatoms   , beta2max_sp        , beta2max_phaeo     ,  &
        mu_max_diatoms     , mu_max_sp          , mu_max_phaeo       ,  &
        grow_Tdep_diatoms  , grow_Tdep_sp       , grow_Tdep_phaeo    ,  &
        fr_graze_diatoms   , fr_graze_sp        , fr_graze_phaeo     ,  &
        mort_pre_diatoms   , mort_pre_sp        , mort_pre_phaeo     ,  &
        mort_Tdep_diatoms  , mort_Tdep_sp       , mort_Tdep_phaeo    ,  &
        k_exude_diatoms    , k_exude_sp         , k_exude_phaeo      ,  &
        K_Nit_diatoms      , K_Nit_sp           , K_Nit_phaeo        ,  &
        K_Am_diatoms       , K_Am_sp            , K_Am_phaeo         ,  &
        K_Sil_diatoms      , K_Sil_sp           , K_Sil_phaeo        ,  &
        K_Fe_diatoms       , K_Fe_sp            , K_Fe_phaeo         ,  &
        f_don_protein      , kn_bac_protein     , f_don_Am_protein   ,  &
        f_doc_s            , f_doc_l            , f_exude_s          ,  &
        f_exude_l          , k_bac_s            , k_bac_l            ,  &
        T_max              , fsal               , op_dep_min         ,  &
        fr_graze_s         , fr_graze_e         , fr_mort2min        ,  &
        fr_dFe             , k_nitrif           , t_iron_conv        ,  &
        max_loss           , max_dfe_doc1       , fr_resp_s          ,  &
        y_sk_DMS           , t_sk_conv          , t_sk_ox            ,  &
        algaltype_diatoms  , algaltype_sp       , algaltype_phaeo    ,  &
        nitratetype        , ammoniumtype       , silicatetype       ,  &
        dmspptype          , dmspdtype          , humtype            ,  &
        doctype_s          , doctype_l          , dontype_protein    ,  &
        fedtype_1          , feptype_1          , zaerotype_bc1      ,  &
        zaerotype_bc2      , zaerotype_dust1    , zaerotype_dust2    ,  &
        zaerotype_dust3    , zaerotype_dust4    , ratio_C2N_diatoms  ,  &
        ratio_C2N_sp       , ratio_C2N_phaeo    , ratio_chl2N_diatoms,  & 
        ratio_chl2N_sp     , ratio_chl2N_phaeo  , F_abs_chl_diatoms  ,  &
        F_abs_chl_sp       , F_abs_chl_phaeo    , ratio_C2N_proteins 

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_aero_out=tr_aero)
      call icepack_query_parameters(ktherm_out=ktherm, shortwave_out=shortwave, &
           scale_bgc_out=scale_bgc, skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, &
           dEdd_algae_out=dEdd_algae, solve_zbgc_out=solve_zbgc, phi_snow_out=phi_snow, &
           bgc_flux_type_out=bgc_flux_type, grid_o_out=grid_o, l_sk_out=l_sk, &
           initbio_frac_out=initbio_frac, frazil_scav_out=frazil_scav, &
           algal_vel_out=algal_vel, R_dFe2dust_out=R_dFe2dust, &
           dustFe_sol_out=dustFe_sol, T_max_out=T_max, fsal_out=fsal, &
           op_dep_min_out=op_dep_min, fr_graze_s_out=fr_graze_s, &
           fr_graze_e_out=fr_graze_e, fr_mort2min_out=fr_mort2min, &
           fr_dFe_out=fr_dFe, k_nitrif_out=k_nitrif, t_iron_conv_out=t_iron_conv, &
           max_loss_out=max_loss, max_dfe_doc1_out=max_dfe_doc1, &
           fr_resp_s_out=fr_resp_s, y_sk_DMS_out=y_sk_DMS, t_sk_conv_out=t_sk_conv, &
           t_sk_ox_out=t_sk_ox)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      tr_brine        = .false.  ! brine height differs from ice height
      tr_zaero        = .false.  ! z aerosol tracers
      modal_aero      = .false.  ! use modal aerosol treatment of aerosols
      restore_bgc     = .false.  ! restore bgc if true
      solve_zsal      = .false.  ! update salinity tracer profile from solve_S_dt
      bgc_data_type   = 'default'! source of bgc data
      tr_bgc_PON      = .false.  !---------------------------------------------   
      tr_bgc_Nit      = .false.  ! biogeochemistry (skl or zbgc)
      tr_bgc_C        = .false.  ! if skl_bgc = .true. then skl
      tr_bgc_chl      = .false.  ! if z_tracers = .true. then vertically resolved
      tr_bgc_Sil      = .false.  ! if z_tracers + solve_zbgc = .true. then
      tr_bgc_Am       = .false.  ! vertically resolved with reactions  
      tr_bgc_DMS      = .false.  !------------------------------------------------
      tr_bgc_DON      = .false.  ! 
      tr_bgc_hum      = .false.  !
      tr_bgc_Fe       = .false.  ! 
      tr_bgc_N        = .true.   !

      ! z biology parameters  
      ratio_Si2N_diatoms = 1.8_dbl_kind    ! algal Si to N (mol/mol)                       
      ratio_Si2N_sp      = c0              ! diatoms, small plankton, phaeocystis
      ratio_Si2N_phaeo   = c0
      ratio_S2N_diatoms  = 0.03_dbl_kind   ! algal S  to N (mol/mol)
      ratio_S2N_sp       = 0.03_dbl_kind 
      ratio_S2N_phaeo    = 0.03_dbl_kind
      ratio_Fe2C_diatoms = 0.0033_dbl_kind ! algal Fe to C  (umol/mol)
      ratio_Fe2C_sp      = 0.0033_dbl_kind
      ratio_Fe2C_phaeo   = p1
      ratio_Fe2N_diatoms = 0.023_dbl_kind  ! algal Fe to N  (umol/mol)
      ratio_Fe2N_sp      = 0.023_dbl_kind
      ratio_Fe2N_phaeo   = 0.7_dbl_kind
      ratio_Fe2DON       = 0.023_dbl_kind  ! Fe to N of DON (nmol/umol)
      ratio_Fe2DOC_s     = p1              ! Fe to C of DOC (nmol/umol) saccharids
      ratio_Fe2DOC_l     = 0.033_dbl_kind  ! Fe to C of DOC (nmol/umol) lipids
      tau_min            = 5200.0_dbl_kind ! rapid mobile to stationary exchanges (s)
      tau_max            = 1.73e5_dbl_kind ! long time mobile to stationary exchanges (s)
      chlabs_diatoms     = 0.03_dbl_kind   ! chl absorption (1/m/(mg/m^3))
      chlabs_sp          = 0.01_dbl_kind
      chlabs_phaeo       = 0.05_dbl_kind
      alpha2max_low_diatoms = 0.8_dbl_kind ! light limitation (1/(W/m^2))  
      alpha2max_low_sp      = 0.67_dbl_kind
      alpha2max_low_phaeo   = 0.67_dbl_kind
      beta2max_diatoms   = 0.018_dbl_kind  ! light inhibition (1/(W/m^2))  
      beta2max_sp        = 0.0025_dbl_kind
      beta2max_phaeo     = 0.01_dbl_kind
      mu_max_diatoms     = 1.2_dbl_kind    ! maximum growth rate (1/day) 
      mu_max_sp          = 0.851_dbl_kind
      mu_max_phaeo       = 0.851_dbl_kind
      grow_Tdep_diatoms  = 0.06_dbl_kind ! Temperature dependence of growth (1/C)
      grow_Tdep_sp       = 0.06_dbl_kind
      grow_Tdep_phaeo    = 0.06_dbl_kind
      fr_graze_diatoms   = 0.01_dbl_kind ! Fraction grazed
      fr_graze_sp        = p1
      fr_graze_phaeo     = p1
      mort_pre_diatoms   = 0.007_dbl_kind! Mortality (1/day)
      mort_pre_sp        = 0.007_dbl_kind
      mort_pre_phaeo     = 0.007_dbl_kind
      mort_Tdep_diatoms  = 0.03_dbl_kind ! T dependence of mortality (1/C)
      mort_Tdep_sp       = 0.03_dbl_kind
      mort_Tdep_phaeo    = 0.03_dbl_kind
      k_exude_diatoms    = c0            ! algal exudation (1/d)
      k_exude_sp         = c0
      k_exude_phaeo      = c0
      K_Nit_diatoms      = c1            ! nitrate half saturation (mmol/m^3)
      K_Nit_sp           = c1
      K_Nit_phaeo        = c1
      K_Am_diatoms       = 0.3_dbl_kind  ! ammonium half saturation (mmol/m^3)
      K_Am_sp            = 0.3_dbl_kind
      K_Am_phaeo         = 0.3_dbl_kind
      K_Sil_diatoms      = 4.0_dbl_kind  ! silicate half saturation (mmol/m^3)
      K_Sil_sp           = c0
      K_Sil_phaeo        = c0
      K_Fe_diatoms       = c1            ! iron half saturation (nM)
      K_Fe_sp            = 0.2_dbl_kind
      K_Fe_phaeo         = p1
      f_don_protein      = 0.6_dbl_kind  ! fraction of spilled grazing to proteins           
      kn_bac_protein     = 0.03_dbl_kind ! Bacterial degredation of DON (1/d)                
      f_don_Am_protein   = 0.25_dbl_kind ! fraction of remineralized DON to ammonium         
      f_doc_s            = 0.4_dbl_kind  ! fraction of mortality to DOC 
      f_doc_l            = 0.4_dbl_kind
      f_exude_s          = c1            ! fraction of exudation to DOC
      f_exude_l          = c1
      k_bac_s            = 0.03_dbl_kind ! Bacterial degredation of DOC (1/d)
      k_bac_l            = 0.03_dbl_kind
      algaltype_diatoms  = c0            ! ------------------
      algaltype_sp       = p5            !
      algaltype_phaeo    = p5            !
      nitratetype        = -c1           ! mobility type between
      ammoniumtype       = c1            ! stationary <-->  mobile
      silicatetype       = -c1           !
      dmspptype          = p5            !
      dmspdtype          = -c1           !
      humtype            = c1            !
      doctype_s          = p5            !
      doctype_l          = p5            !
      dontype_protein    = p5            !
      fedtype_1          = p5            !
      feptype_1          = p5            !
      zaerotype_bc1      = c1            !
      zaerotype_bc2      = c1            !
      zaerotype_dust1    = c1            !
      zaerotype_dust2    = c1            !
      zaerotype_dust3    = c1            !
      zaerotype_dust4    = c1            !--------------------
      ratio_C2N_diatoms  = 7.0_dbl_kind  ! algal C to N ratio (mol/mol)
      ratio_C2N_sp       = 7.0_dbl_kind
      ratio_C2N_phaeo    = 7.0_dbl_kind
      ratio_chl2N_diatoms= 2.1_dbl_kind  ! algal chlorophyll to N ratio (mg/mmol)
      ratio_chl2N_sp     = 1.1_dbl_kind
      ratio_chl2N_phaeo  = 0.84_dbl_kind
      F_abs_chl_diatoms  = 2.0_dbl_kind  ! scales absorbed radiation for dEdd
      F_abs_chl_sp       = 4.0_dbl_kind
      F_abs_chl_phaeo    = 5.0_dbl_kind
      ratio_C2N_proteins = 7.0_dbl_kind  ! ratio of C to N in proteins (mol/mol)       

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif 

      print*,'Reading zbgc_nml'
      do while (nml_error > 0)
         read(nu_nml, nml=zbgc_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nu_nml)
      if (nml_error /= 0) then
         print*,'error reading zbgc namelist'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      !-----------------------------------------------------------------
      ! resolve conflicts
      !-----------------------------------------------------------------
      ! zsalinity and brine
      !-----------------------------------------------------------------
      if (solve_zsal .and. TRZS == 0) then
         write(nu_diag,*) 'WARNING: solve_zsal=T but 0 zsalinity tracers'
         write(nu_diag,*) 'WARNING: setting solve_zsal = F'
         solve_zsal = .false.      
      elseif (solve_zsal .and. nblyr < 1)  then
         write(nu_diag,*) 'WARNING: solve_zsal=T but 0 zsalinity tracers'
         write(nu_diag,*) 'WARNING: setting solve_zsal = F'
         solve_zsal = .false.     
      endif 

      if (solve_zsal .and. ((.not. tr_brine) .or. (ktherm /= 1))) then
         write(nu_diag,*) 'WARNING: solve_zsal needs tr_brine=T and ktherm=1'
         write(nu_diag,*) 'WARNING: setting tr_brine=T and ktherm=1'
         tr_brine = .true.
         ktherm = 1
      endif

      if (tr_brine .and. TRBRI == 0 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no brine height compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      elseif (tr_brine .and. nblyr < 1 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no biology layers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      endif 

         write(nu_diag,1010) ' tr_brine                  = ', tr_brine
      if (tr_brine) then
         write(nu_diag,1005) ' phi_snow                  = ', phi_snow
      endif
      if (solve_zsal) then
         write(nu_diag,1010) ' solve_zsal                = ', solve_zsal
         write(nu_diag,1000) ' grid_oS                   = ', grid_oS
         write(nu_diag,1005) ' l_skS                     = ', l_skS
      endif

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      if (.not. tr_brine) then
         if (solve_zbgc) then
            write(nu_diag,*) 'WARNING: tr_brine = F and solve_zbgc = T'
            write(nu_diag,*) 'WARNING: setting solve_zbgc = F'
            solve_zbgc = .false.
         endif
         if (skl_bgc) then
            write(nu_diag,*) 'WARNING: tr_brine = F and skl_bgc = T'
            write(nu_diag,*) 'WARNING: setting skl_bgc = F'
            skl_bgc = .false.
         endif
         if (tr_zaero) then
            write(nu_diag,*) 'WARNING: tr_brine = F and tr_zaero = T'
            write(nu_diag,*) 'WARNING: setting tr_zaero = F'
            tr_zaero = .false.
         endif
      endif

      if ((skl_bgc .AND. solve_zbgc) .or. (skl_bgc .AND. z_tracers)) then
         print*, 'ERROR: skl_bgc and (solve_zbgc or z_tracers) are both true'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (skl_bgc .AND. tr_zaero) then
         write(nu_diag,*) 'WARNING: skl bgc does not use vertical tracers'
         write(nu_diag,*) 'WARNING: setting tr_zaero = F'
         tr_zaero = .false.
      endif

      if (dEdd_algae .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) 'WARNING: dEdd_algae = T but shortwave /= dEdd'
         write(nu_diag,*) 'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (dEdd_algae .AND. (.NOT. tr_bgc_N) .AND. (.NOT. tr_zaero)) then 
         write(nu_diag,*) 'WARNING: need tr_bgc_N or tr_zaero for dEdd_algae'
         write(nu_diag,*) 'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (modal_aero .AND. (.NOT. tr_zaero) .AND. (.NOT. tr_aero)) then
         modal_aero = .false.
      endif
         
      if (modal_aero .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) 'WARNING: modal_aero = T but shortwave /= dEdd'
         write(nu_diag,*) 'WARNING: setting modal_aero = F'
         modal_aero = .false.
      endif
      if (n_algae > icepack_max_algae) then
         print*, 'error:number of algal types exceeds icepack_max_algae'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (n_doc > icepack_max_doc) then
         print*, 'error:number of algal types exceeds icepack_max_doc'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (n_dic > icepack_max_dic) then
         print*, 'error:number of dic types exceeds icepack_max_dic'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (n_don > icepack_max_don) then
         print*, 'error:number of don types exceeds icepack_max_don'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (n_fed > icepack_max_fe) then
         print*, 'error:number of dissolved fe types exceeds icepack_max_fe'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (n_fep > icepack_max_fe) then
         print*, 'error:number of particulate fe types exceeds icepack_max_fe'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if ((TRBGCS == 0 .and. skl_bgc) .or. (TRALG == 0 .and. skl_bgc)) then
         write(nu_diag,*) &
            'WARNING: skl_bgc=T but 0 bgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting skl_bgc = F'
         skl_bgc = .false.
      endif

      if ((TRBGCZ == 0 .and. solve_zbgc) .or. (TRALG == 0 .and. solve_zbgc)) then
         write(nu_diag,*) &
            'WARNING: solve_zbgc=T but 0 zbgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zbgc = F'
         solve_zbgc = .false.
      endif

      if (solve_zbgc .and. .not. z_tracers) z_tracers = .true.
      if (skl_bgc .or. solve_zbgc) then
         tr_bgc_N         = .true.   ! minimum NP biogeochemistry
         tr_bgc_Nit       = .true.
      else
         tr_bgc_N         = .false.
         tr_bgc_C         = .false.
         tr_bgc_chl       = .false.
         tr_bgc_Nit       = .false.
         tr_bgc_Am        = .false.
         tr_bgc_Sil       = .false.
         tr_bgc_hum       = .false.
         tr_bgc_DMS       = .false.
         tr_bgc_PON       = .false.
         tr_bgc_DON       = .false.
         tr_bgc_Fe        = .false.
      endif

      !-----------------------------------------------------------------
      ! z layer aerosols
      !-----------------------------------------------------------------
      if (tr_zaero .and. .not. z_tracers) z_tracers = .true.

      if (n_zaero > icepack_max_aero) then
         print*, 'error:number of z aerosols exceeds icepack_max_aero'
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif         

      if (skl_bgc .and. n_bgc < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of bgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',n_bgc
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (solve_zbgc .and. n_bgc < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of zbgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',n_bgc
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (tr_zaero .and. TRZAERO <  1) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of TRZAERO > 0'
         write (nu_diag,*) 'in order to solve z aerosols:',TRZAERO
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      !-----------------------------------------------------------------
      ! end conflict resolution
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! set Icepack values
      !-----------------------------------------------------------------

      call icepack_init_parameters(ktherm_in=ktherm, shortwave_in=shortwave, &
           scale_bgc_in=scale_bgc, skl_bgc_in=skl_bgc, z_tracers_in=z_tracers, &
           dEdd_algae_in=dEdd_algae, solve_zbgc_in=solve_zbgc, &
           bgc_flux_type_in=bgc_flux_type, grid_o_in=grid_o, l_sk_in=l_sk, &
           initbio_frac_in=initbio_frac, frazil_scav_in=frazil_scav, &
           grid_oS_in=grid_oS, l_skS_in=l_skS, phi_snow_in=phi_snow, &
           algal_vel_in=algal_vel, R_dFe2dust_in=R_dFe2dust, &
           dustFe_sol_in=dustFe_sol, T_max_in=T_max, fsal_in=fsal, &
           op_dep_min_in=op_dep_min, fr_graze_s_in=fr_graze_s, &
           fr_graze_e_in=fr_graze_e, fr_mort2min_in=fr_mort2min, &
           fr_dFe_in=fr_dFe, k_nitrif_in=k_nitrif, t_iron_conv_in=t_iron_conv, &
           max_loss_in=max_loss, max_dfe_doc1_in=max_dfe_doc1, fr_resp_in=fr_resp, &
           fr_resp_s_in=fr_resp_s, y_sk_DMS_in=y_sk_DMS, t_sk_conv_in=t_sk_conv, &
           t_sk_ox_in=t_sk_ox, modal_aero_in=modal_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! initialize zbgc tracer indices
      !----------------------------------------------------------------- 

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

      nt_bgc_S = 0
      if (solve_zsal) then       ! .true. only if tr_brine = .true.
          nt_bgc_S = ntrcr + 1
          ntrcr = ntrcr + nblyr
          do k = 1,nblyr
             trcr_depend(nt_bgc_S + k - 1) = 2 + nt_fbri + ntd
             trcr_base  (nt_bgc_S,1) = c0  ! default: ice area
             trcr_base  (nt_bgc_S,2) = c1 
             trcr_base  (nt_bgc_S,3) = c0  
             n_trcr_strata(nt_bgc_S) = 1
             nt_strata(nt_bgc_S,1) = nt_fbri
             nt_strata(nt_bgc_S,2) = 0
          enddo
      endif 

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      nbtrcr = 0
      nbtrcr_sw = 0

      ! vectors of size icepack_max_algae
      nlt_bgc_N(:) = 0
      nlt_bgc_C(:) = 0
      nlt_bgc_chl(:) = 0
      nt_bgc_N(:) = 0
      nt_bgc_C(:) = 0
      nt_bgc_chl(:) = 0

      ! vectors of size icepack_max_dic
      nlt_bgc_DIC(:) = 0
      nt_bgc_DIC(:) = 0

      ! vectors of size icepack_max_doc
      nlt_bgc_DOC(:) = 0
      nt_bgc_DOC(:) = 0

      ! vectors of size icepack_max_don
      nlt_bgc_DON(:) = 0
      nt_bgc_DON(:) = 0

      ! vectors of size icepack_max_fe 
      nlt_bgc_Fed(:) = 0
      nlt_bgc_Fep(:) = 0
      nt_bgc_Fed(:) = 0
      nt_bgc_Fep(:) = 0

      ! vectors of size icepack_max_aero
      nlt_zaero(:) = 0
      nlt_zaero_sw(:) = 0
      nt_zaero(:) = 0
      nt_zaero_sw(:) = 0

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
      nt_zbgc_frac  = 0

      !-----------------------------------------------------------------
      ! Define array parameters
      !-----------------------------------------------------------------

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

      f_don(1) = f_don_protein
      kn_bac(1) = kn_bac_protein
      f_don_Am(1) = f_don_Am_protein

      f_exude(1) = f_exude_s
      f_exude(2) = f_exude_l
      k_bac(1) = k_bac_s
      k_bac(2) = k_bac_l

      dictype(:) = -c1
      
      algaltype(1) = algaltype_diatoms
      algaltype(2) = algaltype_sp
      algaltype(3) = algaltype_phaeo

      doctype(1) = doctype_s
      doctype(2) = doctype_l
 
      dontype(1) = dontype_protein

      fedtype(1) = fedtype_1
      feptype(1) = feptype_1

      zaerotype(1) = zaerotype_bc1
      zaerotype(2) = zaerotype_bc2
      zaerotype(3) = zaerotype_dust1
      zaerotype(4) = zaerotype_dust2
      zaerotype(5) = zaerotype_dust3
      zaerotype(6) = zaerotype_dust4

      call icepack_init_zbgc ( &
           R_Si2N_in=R_Si2N, &
           R_S2N_in=R_S2N, R_Fe2C_in=R_Fe2C, R_Fe2N_in=R_Fe2N, R_C2N_in=R_C2N, &
           R_chl2N_in=R_chl2N, F_abs_chl_in=F_abs_chl, R_Fe2DON_in=R_Fe2DON, &
           R_C2N_DON_in=R_C2N_DON, &
           R_Fe2DOC_in=R_Fe2DOC, &
           chlabs_in=chlabs, alpha2max_low_in=alpha2max_low, beta2max_in=beta2max, &
           mu_max_in=mu_max, grow_Tdep_in=grow_Tdep, fr_graze_in=fr_graze, &
           mort_pre_in=mort_pre, &
           mort_Tdep_in=mort_Tdep, k_exude_in=k_exude, &
           K_Nit_in=K_Nit, K_Am_in=K_Am, K_sil_in=K_Sil, K_Fe_in=K_Fe, &
           f_don_in=f_don, kn_bac_in=kn_bac, f_don_Am_in=f_don_Am, f_exude_in=f_exude, &
           k_bac_in=k_bac, &
           fr_resp_in=fr_resp, algal_vel_in=algal_vel, R_dFe2dust_in=R_dFe2dust, &
           dustFe_sol_in=dustFe_sol, T_max_in=T_max, fr_mort2min_in=fr_mort2min, &
           fr_dFe_in=fr_dFe, op_dep_min_in=op_dep_min, &
           fr_graze_s_in=fr_graze_s, fr_graze_e_in=fr_graze_e, &
           k_nitrif_in=k_nitrif, t_iron_conv_in=t_iron_conv, &
           max_loss_in=max_loss, max_dfe_doc1_in=max_dfe_doc1, &
           fr_resp_s_in=fr_resp_s, y_sk_DMS_in=y_sk_DMS, &
           t_sk_conv_in=t_sk_conv, t_sk_ox_in=t_sk_ox, fsal_in=fsal)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      if (skl_bgc) then

         nk = 1
         nt_depend = 0

         if (dEdd_algae) then 
           nlt_chl_sw = 1
           nbtrcr_sw = nilyr+nslyr+2  ! only the bottom layer will be nonzero
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
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_N(mm),    nlt_bgc_N(mm), &
                               algaltype(mm),   nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_N(mm)) = mm
         enddo   ! mm
      endif ! tr_bgc_N

      if (tr_bgc_Nit) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Nit,      nlt_bgc_Nit,   &
                               nitratetype,     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Nit) = icepack_max_algae + 1
      endif ! tr_bgc_Nit
         
      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do mm = 1,n_algae      
       !     call init_bgc_trcr(nk,              nt_fbri,       &
       !                        nt_bgc_C(mm),    nlt_bgc_C(mm), &
       !                        algaltype(mm),   nt_depend,     &
       !                        ntrcr,           nbtrcr,        &
       !                        bgc_tracer_type, trcr_depend,   &
       !                        trcr_base,       n_trcr_strata, &
       !                        nt_strata,       bio_index)
       !     bio_index_o(nlt_bgc_C(mm)) = icepack_max_algae + 1 + mm
       !  enddo   ! mm

         do mm = 1, n_doc
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DOC(mm),  nlt_bgc_DOC(mm), &
                               doctype(mm),     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DOC(mm)) = icepack_max_algae + 1 + mm
         enddo   ! mm
         do mm = 1, n_dic
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DIC(mm),  nlt_bgc_DIC(mm), &
                               dictype(mm),     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DIC(mm)) = icepack_max_algae + icepack_max_doc + 1 + mm
         enddo   ! mm
      endif      ! tr_bgc_C

      if (tr_bgc_chl) then
         do mm = 1, n_algae
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_chl(mm),  nlt_bgc_chl(mm), &
                               algaltype(mm),   nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_chl(mm)) = icepack_max_algae + 1 + icepack_max_doc + icepack_max_dic + mm
         enddo   ! mm
      endif      ! tr_bgc_chl

      if (tr_bgc_Am) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Am,       nlt_bgc_Am,    &
                               ammoniumtype,    nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Am) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 2
      endif    
      if (tr_bgc_Sil) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Sil,      nlt_bgc_Sil,   &
                               silicatetype,    nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Sil) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 3
      endif    
      if (tr_bgc_DMS) then   ! all together
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DMSPp,    nlt_bgc_DMSPp, &
                               dmspptype,       nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPp) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 4

            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DMSPd,    nlt_bgc_DMSPd, &
                               dmspdtype,       nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPd) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 5

            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DMS,      nlt_bgc_DMS,   &
                               dmspdtype,       nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMS) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 6
      endif    
      if (tr_bgc_PON) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_PON,      nlt_bgc_PON, &
                               nitratetype,     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_PON) =  2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 7
      endif
      if (tr_bgc_DON) then
         do mm = 1, n_don
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DON(mm),  nlt_bgc_DON(mm), &
                               dontype(mm),     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DON(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_DON
      if (tr_bgc_Fe) then
         do mm = 1, n_fed
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Fed(mm),  nlt_bgc_Fed(mm), &
                               fedtype(mm),     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fed(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic &
                                         + icepack_max_don + 7 + mm
         enddo   ! mm
         do mm = 1, n_fep
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Fep(mm),  nlt_bgc_Fep(mm), &
                               feptype(mm),     nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fep(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic &
                                         + icepack_max_don + icepack_max_fe + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_Fe 
  
      if (tr_bgc_hum) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_hum,      nlt_bgc_hum,   &
                               humtype,         nt_depend,     &
                               ntrcr,           nbtrcr,        &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_hum) =   2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic &
                                         + icepack_max_don + 2*icepack_max_fe + icepack_max_aero 
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
               call init_bgc_trcr(nk,              nt_fbri,       &
                                  nt_zaero(mm),    nlt_zaero(mm), &
                                  zaerotype(mm),   nt_depend,     &
                                  ntrcr,           nbtrcr,        &
                                  bgc_tracer_type, trcr_depend,   &
                                  trcr_base,       n_trcr_strata, &
                                  nt_strata,       bio_index)
               bio_index_o(nlt_zaero(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic &
                                          + icepack_max_don + 2*icepack_max_fe + 7 + mm
            enddo   ! mm
         endif      ! tr_zaero

!echmod keep trcr indices etc here but move zbgc_frac_init, zbgc_init_frac, tau_ret, tau_rel to icepack
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
      ! set values in icepack
      !-----------------------------------------------------------------

      call icepack_init_zbgc( &
           zbgc_init_frac_in=zbgc_init_frac, tau_ret_in=tau_ret, tau_rel_in=tau_rel, &
           zbgc_frac_init_in=zbgc_frac_init, bgc_tracer_type_in=bgc_tracer_type)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_init_tracer_sizes( &
           n_algae_in=n_algae,                                                                   &
           n_DOC_in=n_DOC,             n_DON_in=n_DON,               n_DIC_in=n_DIC,             &
           n_fed_in=n_fed,             n_fep_in=n_fep,               n_zaero_in=n_zaero,         &
           ntrcr_in=ntrcr, ntrcr_o_in=ntrcr_o, nbtrcr_in=nbtrcr, nbtrcr_sw_in=nbtrcr_sw)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_init_tracer_flags( &
           tr_brine_in  =tr_brine, &
           tr_bgc_Nit_in=tr_bgc_Nit, tr_bgc_Am_in =tr_bgc_Am,  tr_bgc_Sil_in=tr_bgc_Sil,   &
           tr_bgc_DMS_in=tr_bgc_DMS, tr_bgc_PON_in=tr_bgc_PON,                             &
           tr_bgc_N_in  =tr_bgc_N,   tr_bgc_C_in  =tr_bgc_C,   tr_bgc_chl_in=tr_bgc_chl,   &
           tr_bgc_DON_in=tr_bgc_DON, tr_bgc_Fe_in =tr_bgc_Fe,  tr_zaero_in  =tr_zaero,     &
           tr_bgc_hum_in=tr_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_init_tracer_indices( &
           nt_fbri_in=nt_fbri,                                                                   &  
           nt_bgc_Nit_in=nt_bgc_Nit,   nt_bgc_Am_in=nt_bgc_Am,       nt_bgc_Sil_in=nt_bgc_Sil,   &
           nt_bgc_DMS_in=nt_bgc_DMS,   nt_bgc_PON_in=nt_bgc_PON,     nt_bgc_S_in=nt_bgc_S,       &
           nt_bgc_N_in=nt_bgc_N,       nt_bgc_chl_in=nt_bgc_chl,     nt_bgc_hum_in=nt_bgc_hum,   &
           nt_bgc_DOC_in=nt_bgc_DOC,   nt_bgc_DON_in=nt_bgc_DON,     nt_bgc_DIC_in=nt_bgc_DIC,   &
           nt_zaero_in=nt_zaero,       nt_bgc_DMSPp_in=nt_bgc_DMSPp, nt_bgc_DMSPd_in=nt_bgc_DMSPd, &
           nt_bgc_Fed_in=nt_bgc_Fed,   nt_bgc_Fep_in=nt_bgc_Fep,     nt_zbgc_frac_in=nt_zbgc_frac, &
           nlt_chl_sw_in=nlt_chl_sw,   nlt_bgc_Sil_in=nlt_bgc_Sil,   nlt_zaero_sw_in=nlt_zaero_sw, &
           nlt_bgc_N_in=nlt_bgc_N,     nlt_bgc_Nit_in=nlt_bgc_Nit,   nlt_bgc_Am_in=nlt_bgc_Am,   &
           nlt_bgc_DMS_in=nlt_bgc_DMS, nlt_bgc_DMSPp_in=nlt_bgc_DMSPp,                           &
           nlt_bgc_DMSPd_in=nlt_bgc_DMSPd,                                                       &
           nlt_bgc_DIC_in=nlt_bgc_DIC, nlt_bgc_DOC_in=nlt_bgc_DOC,   nlt_bgc_PON_in=nlt_bgc_PON, &
           nlt_bgc_DON_in=nlt_bgc_DON, nlt_bgc_Fed_in=nlt_bgc_Fed,   nlt_bgc_Fep_in=nlt_bgc_Fep, &
           nlt_bgc_chl_in=nlt_bgc_chl, nlt_bgc_hum_in=nlt_bgc_hum,   nlt_zaero_in=nlt_zaero,     &
           bio_index_o_in=bio_index_o, bio_index_in=bio_index)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! final consistency checks
      !----------------------------------------------------------------- 
      if (nbtrcr > icepack_max_nbtrcr) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr > icepack_max_nbtrcr'
         write (nu_diag,*) 'nbtrcr, icepack_max_nbtrcr:',nbtrcr, icepack_max_nbtrcr
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif
      if (.NOT. dEdd_algae) nbtrcr_sw = 1

      if (nbtrcr_sw > max_nsw) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr_sw > max_nsw'
         write (nu_diag,*) 'nbtrcr_sw, max_nsw:',nbtrcr_sw, max_nsw
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif

      if (ntrcr > max_ntrcr) then
         write(nu_diag,*) 'max_ntrcr < number of namelist tracers'
         write(nu_diag,*) 'max_ntrcr = ',max_ntrcr,' ntrcr = ',ntrcr
         call icedrv_system_abort(file=__FILE__,line=__LINE__)
      endif                               

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (skl_bgc) then

         write(nu_diag,1010) ' skl_bgc                   = ', skl_bgc
         write(nu_diag,1030) ' bgc_flux_type             = ', bgc_flux_type
         write(nu_diag,1010) ' restore_bgc               = ', restore_bgc
         write(nu_diag,*)    ' bgc_data_type             = ', &
                               trim(bgc_data_type)
         write(nu_diag,1020) ' number of bio tracers     = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc          = ', n_doc
         write(nu_diag,1020) ' number of dic          = ', n_dic
         write(nu_diag,1020) ' number of don          = ', n_don
         write(nu_diag,1020) ' number of fed          = ', n_fed
         write(nu_diag,1020) ' number of fep          = ', n_fep
         write(nu_diag,1010) ' tr_bgc_N               = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_C               = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_chl             = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_Nit             = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_Am              = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_Sil             = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum             = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_DMS             = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON             = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON             = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe              = ', tr_bgc_Fe 
        
      elseif (z_tracers) then

         write(nu_diag,*)    ' bgc_data_type             = ', &
                               trim(bgc_data_type)
         write(nu_diag,1010) ' dEdd_algae                = ', dEdd_algae  
         write(nu_diag,1010) ' modal_aero                = ', modal_aero  
         write(nu_diag,1010) ' scale_bgc                 = ', scale_bgc
         write(nu_diag,1010) ' solve_zbgc                = ', solve_zbgc
         write(nu_diag,1020) ' number of ztracers        = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc             = ', n_doc
         write(nu_diag,1020) ' number of dic             = ', n_dic
         write(nu_diag,1020) ' number of fed             = ', n_fed
         write(nu_diag,1020) ' number of fep             = ', n_fep
         write(nu_diag,1020) ' number of aerosols        = ', n_zaero
         write(nu_diag,1010) ' tr_zaero                  = ', tr_zaero
         write(nu_diag,1010) ' tr_bgc_Nit                = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_N                  = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_Am                 = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_C                  = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_Sil                = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum                = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_chl                = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_DMS                = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON                = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON                = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe                 = ', tr_bgc_Fe 
         write(nu_diag,1000) ' grid_o                    = ', grid_o
         write(nu_diag,1005) ' l_sk                      = ', l_sk
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac
         write(nu_diag,1000) ' frazil_scav               = ', frazil_scav  

      endif  ! skl_bgc or solve_bgc

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

      end subroutine init_zbgc

!=======================================================================

      subroutine init_bgc_trcr(nk,              nt_fbri,       &
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

      character(len=*), parameter :: subname='(init_bgc_trcr)'

      !--------

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

      end subroutine init_bgc_trcr

!=======================================================================

      end module icedrv_init_column

!=======================================================================
