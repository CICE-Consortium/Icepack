!  SVN:$Id: icepack_zbgc.F90 1226 2017-05-22 22:45:03Z tcraig $
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
      use icepack_constants, only: c1,  c2, p5, c0, p1, puny
      use icepack_zbgc_shared, only: zbgc_init_frac
      use icepack_zbgc_shared, only: zbgc_frac_init
      use icepack_zbgc_shared, only: bgc_tracer_type, remap_zbgc
      use icepack_zbgc_shared, only: regrid_stationary
      use icepack_zbgc_shared, only: R_S2N, R_Si2N, R_Fe2C, R_Fe2N, R_Fe2DON, R_Fe2DOC
      use icepack_zbgc_shared, only: chlabs, alpha2max_low, beta2max
      use icepack_zbgc_shared, only: mu_max, grow_Tdep, fr_graze
      use icepack_zbgc_shared, only: mort_pre, mort_Tdep, k_exude
      use icepack_zbgc_shared, only: K_Nit, K_Am, K_Sil, K_Fe
      use icepack_zbgc_shared, only: f_don, kn_bac, f_don_Am
      use icepack_zbgc_shared, only: f_exude, k_bac
      use icepack_zbgc_shared, only: tau_ret, tau_rel

      implicit none 

      private
      public :: add_new_ice_bgc, &
                lateral_melt_bgc, &
                icepack_init_bgc, &
                icepack_init_zbgc, &
                icepack_init_bgc_trcr, &
                icepack_biogeochemistry, &
                icepack_init_OceanConcArray, &
                icepack_init_ocean_conc

!=======================================================================

      contains

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (dt,         nblyr,                &
                                  ncat,       nilyr,      nltrcr,   &
                                  bgrid,      cgrid,      igrid,    &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vsnon1,   &
                                  vi0new,                           &
                                  ntrcr,      trcrn,      nbtrcr,   &
                                  sss,        ocean_bio,  flux_bio, &
                                  hsurp,      l_stop,   &
                                  stop_label, l_conservation_check)

      use icepack_constants, only: c0, c1, puny, depressT
      use icepack_itd, only: column_sum, column_conservation_check
      use icepack_tracers, only: tr_brine, nt_fbri, nt_sice, nt_qice, nt_Tsfc
      use icepack_parameters, only: solve_zsal
      use icepack_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         nblyr   , & ! number of bio layers
         ncat     , & ! number of thickness categories
         nilyr    , & ! number of ice layers
         nltrcr, & ! number of zbgc tracers
         nbtrcr  , & ! number of biology tracers
         ntrcr       ! number of tracers in use

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)

      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! intiial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen           ! volume per unit area of ice          (m)

      real (kind=dbl_kind), intent(in) :: &
         vsnon1          ! category 1 snow volume per unit area (m)

      real (kind=dbl_kind), dimension (:,:), &
         intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         sss              !sea surface salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         vi0_init    , & ! volume of new ice added to cat 1 (intial)
         vi0new          ! volume of new ice added to cat 1

      real (kind=dbl_kind), intent(in) :: &
         hsurp           ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (:), &
         intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         ocean_bio       ! ocean concentration of biological tracer

      logical (kind=log_kind), intent(in) :: &
         l_conservation_check

      logical (kind=log_kind), intent(inout) :: &
         l_stop     
        
      character (char_len), intent(inout) :: stop_label

! local

      integer (kind=int_kind) :: &
         location    , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         n           , & ! ice category index
         k               ! ice layer index

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
         vice_new        ! vicen_init + vsurp

      real (kind=dbl_kind) :: &
         Tmlts       ! melting temperature (oC)

      character (len=char_len) :: &
         fieldid         ! field identifier

      !-----------------------------------------------------------------     
      ! brine
      !-----------------------------------------------------------------
      vbrin(:) = c0
      do n = 1, ncat
         vbrin(n) = vicen_init(n)
         if (tr_brine) vbrin(n) =  trcrn(nt_fbri,n)*vicen_init(n)
      enddo
     
      call column_sum (ncat,  vbrin,  vbri_init)

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
         if (tr_brine .and. vicen(n) > c0) then
            trcrn(nt_fbri,n) = vbrin(n)/vicen(n)
         elseif (tr_brine .and. vicen(n) <= c0) then
            trcrn(nt_fbri,n) = c1
         endif

         if (nltrcr > 0) then 
            location = 1  
            call adjust_tracer_profile(nbtrcr,   dt, ntrcr, &
                                       aicen_init(n),       &
                                       vbrin(n),            &
                                       vice_new,            &
                                       trcrn(:,n),          &
                                       vtmp,                &
                                       vsurp,        sss,   &
                                       nilyr,        nblyr, &
                                       solve_zsal,   bgrid, & 
                                       cgrid,               &
                                       ocean_bio,    igrid, &
                                       location,            &
                                       l_stop,     stop_label)
            if (l_stop) return
         endif       ! nltrcr       
      endif          ! hsurp > 0
      enddo          ! n

      !-----------------------------------------------------------------
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------
       
      if (vi0new > c0) then

         vbri1    = vbrin(1) 
         vbrin(1) = vbrin(1) + vi0new
         if (tr_brine .and. vicen(1) > c0) then
            trcrn(nt_fbri,1) = vbrin(1)/vicen(1)
         elseif (tr_brine .and. vicen(1) <= c0) then
            trcrn(nt_fbri,1) = c1
         endif
       
      ! Diffuse_bio handles concentration changes from ice growth/melt
      ! ice area changes
      ! add salt throughout, location = 0

         if (nltrcr > 0) then 
            location = 0  
            call adjust_tracer_profile(nbtrcr,  dt,     ntrcr,  &
                                       aicen(1),                &
                                       vbrin(1),                &
                                       vicen(1),                &
                                       trcrn(:,1),              &
                                       vbri1,                   &
                                       vi0new,          sss,    &
                                       nilyr,           nblyr,  &
                                       solve_zsal,      bgrid,  &
                                       cgrid,                   &
                                       ocean_bio,       igrid,  &
                                       location,                &
                                       l_stop,     stop_label)
            if (l_stop) return

            if (solve_zsal .and. vsnon1 .le. c0) then
               Tmlts = -trcrn(nt_sice,1)*depressT
               trcrn(nt_Tsfc,1) =  calculate_Tin_from_qin(trcrn(nt_qice,1),Tmlts)
            endif        ! solve_zsal 
         endif           ! nltrcr > 0
      endif              ! vi0new > 0

      if (tr_brine .and. l_conservation_check) then
         call column_sum (ncat,   vbrin,  vbri_final)

         fieldid = 'vbrin, add_new_ice_bgc'
         call column_conservation_check (fieldid,                  &
                                         vbri_init, vbri_final,    &
                                         puny,      l_stop)

         if (l_stop) then
            stop_label = 'add_new_ice_bgc: Column conservation error'
            return
         endif
      endif   ! l_conservation_check

      end subroutine add_new_ice_bgc

!=======================================================================

! When sea ice melts laterally, flux bgc to ocean

      subroutine lateral_melt_bgc (dt,                 &
                                   ncat,     nblyr,    &
                                   rside,    vicen,    &
                                   trcrn,    fzsal,    &
                                   flux_bio, nbltrcr)

      use icepack_tracers, only: nt_fbri, nt_bgc_S, bio_index
      use icepack_parameters, only: solve_zsal, rhosi
      use icepack_constants, only: c1, p001

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nblyr , & ! number of bio layers
         nbltrcr   ! number of biology tracers

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcrn     ! tracer array

      real (kind=dbl_kind), intent(in) :: &
         rside     ! fraction of ice that melts laterally

      real (kind=dbl_kind), intent(inout) :: &
         fzsal     ! salt flux from layer Salinity (kg/m^2/s)
  
      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! layer index
         m     , & !
         n         ! category index

      real (kind=dbl_kind) :: &
         zspace    ! bio grid spacing

      zspace = c1/(real(nblyr,kind=dbl_kind))

      if (solve_zsal) then
         do n = 1, ncat
         do k = 1,nblyr
            fzsal = fzsal + rhosi*trcrn(nt_fbri,n) &
                  * vicen(n)*p001*zspace*trcrn(nt_bgc_S+k-1,n) &
                  * rside/dt
         enddo
         enddo
      endif

      do m = 1, nbltrcr
         do n = 1, ncat
         do k = 1, nblyr+1
            flux_bio(m) = flux_bio(m) + trcrn(nt_fbri,n) &
                        * vicen(n)*zspace*trcrn(bio_index(m)+k-1,n) &
                        * rside/dt
         enddo
         enddo
      enddo

      end subroutine lateral_melt_bgc 

!=======================================================================
!
! Add new ice tracers to the ice bottom and adjust the vertical profile 
!
! author: Nicole Jeffery, LANL

      subroutine adjust_tracer_profile (nbtrcr, dt, ntrcr, &
                                        aicen,      vbrin, &
                                        vicen,      trcrn, &
                                        vtmp,              &
                                        vsurp,      sss,   &
                                        nilyr,      nblyr, &
                                        solve_zsal, bgrid, &
                                        cgrid,      ocean_bio, &
                                        igrid,      location, &
                                        l_stop,     stop_label)

      use icepack_constants, only: c1, c0
      use icepack_tracers, only: nt_sice, nt_bgc_S, bio_index 
      use icepack_parameters, only: min_salin, salt_loss

      integer (kind=int_kind), intent(in) :: &
         location          , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         ntrcr             , & ! number of tracers in use
         nilyr             , & ! number of ice layers
         nbtrcr            , & ! number of biology tracers
         nblyr                 ! number of biology layers

      real (kind=dbl_kind), intent(in) :: &
         dt              ! timestep (s)

      real (kind=dbl_kind), intent(in) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume of ice
         sss     , & ! ocean salinity (ppt)
       ! hsurp   , & ! flags new ice added to each cat
         vsurp   , & ! volume of new ice added to each cat
         vtmp        ! total volume of new and old ice

      real (kind=dbl_kind), dimension (nbtrcr), intent(in) :: &
         ocean_bio

      real (kind=dbl_kind), intent(in) :: &
         vbrin       ! fbri*volume per unit area of ice  (m)
       
      logical (kind=log_kind), intent(in) :: &
         solve_zsal 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid       ! zbio grid

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid       ! zsal grid

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid       ! CICE grid

      real (kind=dbl_kind), dimension (ntrcr), &
         intent(inout) :: &
         trcrn       ! ice tracers
      
      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return
        
      character (char_len), intent(inout) :: stop_label

      ! local variables

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0, &      ! temporary, remapped tracers
         trtmp          ! temporary, remapped tracers
        
      real (kind=dbl_kind) :: &
         hin     , & ! ice height
         hinS_new, & ! brine height
         temp_S         

      integer (kind=int_kind) :: &
         k, m 

      real (kind=dbl_kind), dimension (nblyr+1) ::  &     
         C_stationary      ! stationary bulk concentration*h (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr) ::  &     
         S_stationary      ! stationary bulk concentration*h (ppt*m)

      real(kind=dbl_kind) :: &
         top_conc     , & ! salinity or bgc ocean concentration of frazil
         fluxb        , & ! needed for regrid (set to zero here)
         hbri_old     , & ! previous timestep brine height
         hbri             ! brine height 

      trtmp0(:) = c0
      trtmp(:) = c0
      fluxb = c0

      if (location == 1 .and. vbrin > c0) then  ! add frazil to bottom

         hbri     = vbrin
         hbri_old = vtmp
         if (solve_zsal) then
            top_conc = sss * salt_loss
            do k = 1, nblyr 
               S_stationary(k) = trcrn(nt_bgc_S+k-1)* hbri_old
            enddo
            call regrid_stationary (S_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,                  &
                                    nblyr-1,      top_conc, &
                                    bgrid(2:nblyr+1), fluxb,&
                                    l_stop,       stop_label)
            if (l_stop) return
            do k = 1, nblyr 
               trcrn(nt_bgc_S+k-1) =  S_stationary(k)/hbri
               trtmp0(nt_sice+k-1) = trcrn(nt_bgc_S+k-1)
            enddo
         endif  ! solve_zsal

         do m = 1, nbtrcr
            top_conc = ocean_bio(m)*zbgc_init_frac(m)
            do k = 1, nblyr+1 
               C_stationary(k) = trcrn(bio_index(m) + k-1)* hbri_old
            enddo !k
            call regrid_stationary (C_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,                  &
                                    nblyr,        top_conc, &
                                    igrid,        fluxb,    &
                                    l_stop,       stop_label)
            if (l_stop) return
            do k = 1, nblyr+1 
               trcrn(bio_index(m) + k-1) =  C_stationary(k)/hbri
            enddo !k                  
         enddo !m

         if (solve_zsal) then
            if (aicen > c0) then
               hinS_new  = vbrin/aicen
               hin       = vicen/aicen
            else
               hinS_new  = c0
               hin       = c0
            endif                   ! aicen
            temp_S    = min_salin   ! bio to cice
            call remap_zbgc(ntrcr,           nilyr,    &
                            nt_sice,                   &
                            trtmp0(1:ntrcr), trtmp,    &
                            1,               nblyr,    &
                            hin,             hinS_new, &
                            cgrid(2:nilyr+1),          &
                            bgrid(2:nblyr+1), temp_S,  &
                            l_stop,           stop_label)
            do k = 1, nilyr
               trcrn(nt_sice+k-1) = trtmp(nt_sice+k-1)   
            enddo        ! k
         endif           ! solve_zsal

      elseif (vbrin > c0) then   ! add frazil throughout  location == 0 .and.

         do k = 1, nblyr+1
            if (solve_zsal .and. k < nblyr + 1) then
               trcrn(nt_bgc_S+k-1) = (trcrn(nt_bgc_S+k-1) * vtmp &
                                          + sss*salt_loss * vsurp) / vbrin
               trtmp0(nt_sice+k-1) = trcrn(nt_bgc_S+k-1)
            endif                    ! solve_zsal
            do m = 1, nbtrcr
               trcrn(bio_index(m) + k-1) = (trcrn(bio_index(m) + k-1) * vtmp &
                         + ocean_bio(m)*zbgc_init_frac(m) * vsurp) / vbrin
            enddo
         enddo

         if (solve_zsal) then
            if (aicen > c0) then
               hinS_new  = vbrin/aicen
               hin       = vicen/aicen
            else
               hinS_new  = c0
               hin       = c0
            endif              !aicen
            temp_S    = min_salin   ! bio to cice
            call remap_zbgc(ntrcr,        nilyr,    &
                         nt_sice,                   &
                         trtmp0(1:ntrcr), trtmp,    &
                         1,               nblyr,    &
                         hin,             hinS_new, &
                         cgrid(2:nilyr+1),          &        
                         bgrid(2:nblyr+1),temp_S,   &
                         l_stop,           stop_label)
            do k = 1, nilyr
               trcrn(nt_sice+k-1) = trtmp(nt_sice+k-1)   
            enddo        !k
         endif   ! solve_zsal

      endif     ! location

      end subroutine adjust_tracer_profile

!=======================================================================

      subroutine icepack_init_bgc(dt, ncat, nblyr, nilyr, ntrcr_o, &
         cgrid, igrid, ntrcr, nbtrcr, &
         sicen, trcrn, sss, ocean_bio_all, &
         l_stop, stop_label)

      use icepack_constants, only: c0, c1, c2, p1, p15, p5
      use icepack_tracers, only: nt_fbri, nt_bgc_S, nt_sice, nt_zbgc_frac
      use icepack_tracers, only: bio_index_o,  bio_index  
      use icepack_parameters, only: scale_bgc, ktherm, skl_bgc, solve_zsal

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr_o,& ! number of tracers not including bgc
         ntrcr , & ! number of tracers in use
         nbtrcr    ! number of bio tracers in use
 
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

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      character (len=*), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index 
         n     , & ! category index 
         mm    , & ! bio tracer index
         ki    , & ! loop index
         ks        ! 

      real (kind=dbl_kind), dimension (ntrcr+2) :: & 
         trtmp     ! temporary, remapped tracers   
      
      real (kind=dbl_kind) :: & 
         dvssl , & ! volume of snow surface layer (m)
         dvint     ! volume of snow interior      (m)

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant 
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------   
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

            if (scale_bgc .and. solve_zsal) then ! bulk concentration (mmol or mg per m^3)
               do n = 1,ncat
               do mm = 1,nbtrcr
                  do k = 2, nblyr
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) = &
                          (p5*(trcrn(nt_bgc_S+k-1-ntrcr_o,n)+ trcrn(nt_bgc_S+k-2-ntrcr_o,n)) &
                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  enddo  !k
                  trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
                  trcrn(bio_index(mm)-ntrcr_o,n) = (trcrn(nt_bgc_S-ntrcr_o,n) &
                                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  trcrn(bio_index(mm)+nblyr-ntrcr_o,n) = (trcrn(nt_bgc_S+nblyr-1-ntrcr_o,n) &
                                               / sss*ocean_bio_all(bio_index_o(mm)))
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo ! mm
               enddo ! n 
    
            elseif (scale_bgc .and. ktherm == 2) then
               trtmp(:) = c0
               do n = 1,ncat     
                  call remap_zbgc(nilyr,            nilyr,    &
                                  1,                          &
                                  sicen(:,n),       trtmp,    &
                                  0,                nblyr+1,  &
                                  c1,               c1,       &
                                  cgrid(2:nilyr+1),           &
                                  igrid(1:nblyr+1),           &
                                  sicen(1,n),                 &
                                  l_stop,           stop_label)
                  if (l_stop) return

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

      end subroutine icepack_init_bgc

!=======================================================================

      subroutine icepack_init_zbgc (nblyr, nilyr, nslyr, &
                 n_algae, n_zaero, n_doc, n_dic, n_don, n_fed, n_fep, &
                 trcr_base, trcr_depend, n_trcr_strata, nt_strata, nbtrcr_sw, &
                 tr_brine, nt_fbri, ntrcr, nbtrcr, nt_bgc_Nit, nt_bgc_Am, &
                 nt_bgc_Sil, nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, nt_bgc_N, &
                 nt_bgc_C, nt_bgc_chl, nt_bgc_DOC, nt_bgc_DON, nt_bgc_DIC, & 
                 nt_zaero, nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_Fed, nt_bgc_Fep, &
                 nt_zbgc_frac, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS, &
                 tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, tr_bgc_chl, &
                 tr_bgc_DON, tr_bgc_Fe, tr_zaero, nlt_zaero_sw, nlt_chl_sw, &
                 nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
                 nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
                 nlt_bgc_C, nlt_bgc_chl, nlt_bgc_DIC, nlt_bgc_DOC, &
                 nlt_bgc_PON, nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, &
                 nlt_zaero, &
                 nt_bgc_hum, nlt_bgc_hum, tr_bgc_hum, solve_zsal, &
                 skl_bgc, z_tracers, dEdd_algae, solve_zbgc, &
                 frazil_scav, initbio_frac, bio_index_o, bio_index, ntrcr_o, &
                 max_algae, max_doc, max_dic, max_don, max_fe, &
                 ratio_Si2N_diatoms, ratio_Si2N_sp, ratio_Si2N_phaeo, &
                 ratio_S2N_diatoms, ratio_S2N_sp, ratio_S2N_phaeo, &
                 ratio_Fe2C_diatoms, ratio_Fe2C_sp, ratio_Fe2C_phaeo, &
                 ratio_Fe2N_diatoms, ratio_Fe2N_sp, ratio_Fe2N_phaeo, &
                 ratio_Fe2DON, ratio_Fe2DOC_s,  ratio_Fe2DOC_l, & 
                 chlabs_diatoms, chlabs_sp, chlabs_phaeo, &    
                 alpha2max_low_diatoms, alpha2max_low_sp, alpha2max_low_phaeo, &  
                 beta2max_diatoms, beta2max_sp, beta2max_phaeo, &    
                 mu_max_diatoms, mu_max_sp, mu_max_phaeo, &      
                 grow_Tdep_diatoms, grow_Tdep_sp, grow_Tdep_phaeo, &      
                 fr_graze_diatoms, fr_graze_sp, fr_graze_phaeo, &    
                 mort_pre_diatoms, mort_pre_sp, mort_pre_phaeo, &        
                 mort_Tdep_diatoms, mort_Tdep_sp, mort_Tdep_phaeo, &
                 k_exude_diatoms, k_exude_sp, k_exude_phaeo, &   
                 K_Nit_diatoms, K_Nit_sp, K_Nit_phaeo, &     
                 K_Am_diatoms, K_Am_sp, K_Am_phaeo, &     
                 K_Sil_diatoms, K_Sil_sp, K_Sil_phaeo, &     
                 K_Fe_diatoms, K_Fe_sp, K_Fe_phaeo, &  
                 f_don_protein, kn_bac_protein, & 
                 f_don_Am_protein ,f_doc_s, f_doc_l, &
                 f_exude_s, f_exude_l, k_bac_s,  k_bac_l, &
                 algaltype_diatoms, algaltype_sp, algaltype_phaeo, &
                 doctype_s, doctype_l, dontype_protein, &
                 fedtype_1, feptype_1, zaerotype_bc1, zaerotype_bc2, &
                 zaerotype_dust1, zaerotype_dust2, zaerotype_dust3, &
                 zaerotype_dust4, &
                 ratio_C2N_diatoms, ratio_C2N_sp, ratio_C2N_phaeo, &
                 ratio_chl2N_diatoms, ratio_chl2N_sp, ratio_chl2N_phaeo, &
                 F_abs_chl_diatoms, F_abs_chl_sp, F_abs_chl_phaeo, &
                 ratio_C2N_proteins, &
                 nitratetype, ammoniumtype, dmspptype, dmspdtype, &
                 silicatetype, humtype, tau_min, tau_max)
                    
      use icepack_constants,   only: c1, p5, c0, c2
      use icepack_zbgc_shared, only: algaltype, doctype, dictype, dontype, fedtype, feptype, zaerotype
      use icepack_zbgc_shared, only: R_C2N, R_chl2N, F_abs_chl, R_C2N_DON
      use icepack_tracers,     only: max_aero

      integer (kind=int_kind), intent(in) :: &
         nblyr     , & ! number of bio/brine layers per category 
         nilyr     , & ! number of ice layers per category
         nslyr     , & ! number of snow layers per category
         n_zaero   , & ! number of z aerosols in use 
         n_algae   , & ! number of algae in use 
         n_doc     , & ! number of DOC pools in use
         n_dic     , & ! number of DIC pools in use
         n_don     , & ! number of DON pools in use
         n_fed     , & ! number of Fe  pools in use dissolved Fe
         n_fep     , & ! number of Fe  pools in use particulate Fe
         max_algae , &
         max_doc   , &
         max_dic   , &
         max_don   , &
         max_fe

      integer (kind=int_kind), intent(inout) :: &
         ntrcr_o,     & ! number of non-bio tracers in use
         ntrcr,       & ! number of tracers
         nbtrcr,      & ! number of bgc tracers in use
         nbtrcr_sw      ! size of shorwave tracer vector

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

      logical (kind=log_kind), intent(in) :: &
         tr_brine,       & ! if .true., brine height differs from ice thickness
         tr_bgc_S,       & ! if .true., use zsalinity
         tr_zaero,       & ! if .true., black carbon is tracers  (n_zaero)
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
         tr_bgc_hum,     & ! if .true., humic material as tracer
         solve_zsal,     & ! if true, update salinity profile from solve_S_dt
         z_tracers,      & ! if .true., bgc or aerosol tracers are vertically resolved
         solve_zbgc,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae,     & ! if .true., algal absorption of Shortwave is computed in the 
         skl_bgc           ! if true, solve skeletal biochemistry

       integer (kind=int_kind), intent(out) :: &
         nt_fbri,      & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_bgc_Nit,   & ! nutrients  
         nt_bgc_Am,    & ! 
         nt_bgc_Sil,   & !
         nt_bgc_DMSPp, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd, & ! 
         nt_bgc_DMS,   & ! 
         nt_bgc_PON,   & ! zooplankton and detritus 
         nt_bgc_hum,   & ! humic material 
                         ! bio layer indicess
         nlt_bgc_Nit,  & ! nutrients  
         nlt_bgc_Am,   & ! 
         nlt_bgc_Sil,  & !
         nlt_bgc_DMSPp,& ! trace gases (skeletal layer)
         nlt_bgc_DMSPd,& ! 
         nlt_bgc_DMS,  & ! 
         nlt_bgc_PON,  & ! zooplankton and detritus 
         nlt_bgc_hum,  & ! humic material 
         nlt_chl_sw,   & ! points to total chla in trcrn_sw
         nt_zbgc_frac, & ! fraction of tracer in the mobile phase
         nt_bgc_S        ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl,& ! diatoms, phaeocystis, pico/small  
         nlt_bgc_N ,& ! diatoms, phaeocystis, pico/small   
         nlt_bgc_C ,& ! diatoms, phaeocystis, pico/small   
         nlt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_DOC,   & !  dissolved organic carbon  
         nlt_bgc_DOC     !  dissolved organic carbon

      integer (kind=int_kind), dimension(:), intent(out) :: & 
         nt_bgc_DON,   & !  dissolved organic nitrogen
         nlt_bgc_DON     !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_DIC,    & !  dissolved inorganic carbon
         nlt_bgc_DIC      !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(:), intent(out) :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep,     & !  particulate iron
         nlt_bgc_Fed,    & !  dissolved iron
         nlt_bgc_Fep       !  particulate iron

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_zaero,    & !  black carbon and other aerosols 
         nlt_zaero,   & !  black carbon and other aerosols
         nlt_zaero_sw   
    
      integer (kind=int_kind), dimension(:), intent(out) :: &   
         bio_index_o , & ! nlt  to appropriate value in ocean data array
         bio_index       ! nlt to nt

      real (kind=dbl_kind), intent(in) :: &
         initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

      real (kind=dbl_kind), intent(in) :: &
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
        tau_min           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
        tau_max           , &   ! long time mobile to stationary exchanges (s) = 2 days
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
        K_Nit_sp        , &  
        K_Nit_phaeo      , &  
        K_Am_diatoms     , & ! ammonium half saturation (mmol/m^3)
        K_Am_sp         , &   
        K_Am_phaeo       , &   
        K_Sil_diatoms    , & ! silicate half saturation (mmol/m^3)
        K_Sil_sp        , &   
        K_Sil_phaeo      , &   
        K_Fe_diatoms     , & ! iron half saturation (nM)
        K_Fe_sp         , &   
        K_Fe_phaeo       , &    
        f_don_protein    , & ! fraction of spilled grazing to proteins          
        kn_bac_protein   , & ! Bacterial degredation of DON (1/d)               
        f_don_Am_protein , & ! fraction of remineralized DON to ammonium        
        f_doc_s         , & ! fraction of mortality to DOC 
        f_doc_l         , &   
        f_exude_s        , & ! fraction of exudation to DOC
        f_exude_l        , & 
        k_bac_s         , & ! Bacterial degredation of DOC (1/d)
        k_bac_l         , & 
        algaltype_diatoms  , & ! mobility type
        algaltype_sp       , & !
        algaltype_phaeo    , & !
        nitratetype        , & !
        ammoniumtype       , & !
        silicatetype       , & !
        dmspptype         , & !
        dmspdtype         , & !
        humtype           , & !
        doctype_s         , & !
        doctype_l         , & !
        dontype_protein    , & !
        fedtype_1         , & !
        feptype_1         , & !
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

#if (1 == 0)
      ! local variables

      integer (kind=int_kind) :: &
        k, mm    , & ! loop index  
        ntd      , & ! for tracer dependency calculation
        nk       , & !
        nt_depend

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
      R_C2N(1) = ratio_C2N_proteins
     
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
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
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
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Nit,      nlt_bgc_Nit,   &
                                      nitratetype,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Nit) = max_algae + 1
      endif ! tr_bgc_Nit
         
      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do mm = 1,n_algae      
       !     call icepack_init_bgc_trcr(nk,              nt_fbri,       &
       !                               nt_bgc_C(mm),    nlt_bgc_C(mm), &
       !                               algaltype(mm),   nt_depend,     &
       !                               ntrcr,           nbtrcr,        &
       !                               bgc_tracer_type, trcr_depend,   &
       !                               trcr_base,       n_trcr_strata, &
       !                               nt_strata,       bio_index)
       !     bio_index_o(nlt_bgc_C(mm)) = max_algae + 1 + mm
       !  enddo   ! mm

         do mm = 1, n_doc
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DOC(mm),  nlt_bgc_DOC(mm), &
                                      doctype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DOC(mm)) = max_algae + 1 + mm
         enddo   ! mm
         do mm = 1, n_dic
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DIC(mm),  nlt_bgc_DIC(mm), &
                                      dictype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DIC(mm)) = max_algae + max_doc + 1 + mm
         enddo   ! mm
      endif      ! tr_bgc_C

      if (tr_bgc_chl) then
         do mm = 1, n_algae
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_chl(mm),  nlt_bgc_chl(mm), &
                                      algaltype(mm),   nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_chl(mm)) = max_algae + 1 + max_doc + max_dic + mm
         enddo   ! mm
      endif      ! tr_bgc_chl

      if (tr_bgc_Am) then
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Am,       nlt_bgc_Am,    &
                                      ammoniumtype,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Am) = 2*max_algae + max_doc + max_dic + 2
      endif    
      if (tr_bgc_Sil) then
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Sil,      nlt_bgc_Sil,   &
                                      silicatetype,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Sil) = 2*max_algae + max_doc + max_dic + 3
      endif    
      if (tr_bgc_DMS) then   ! all together
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMSPp,    nlt_bgc_DMSPp, &
                                      dmspptype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPp) = 2*max_algae + max_doc + max_dic + 4

            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMSPd,    nlt_bgc_DMSPd, &
                                      dmspdtype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPd) = 2*max_algae + max_doc + max_dic + 5

            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMS,      nlt_bgc_DMS,   &
                                      dmspdtype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMS) = 2*max_algae + max_doc + max_dic + 6
      endif    
      if (tr_bgc_PON) then
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_PON,      nlt_bgc_PON, &
                                      nitratetype,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_PON) =  2*max_algae + max_doc + max_dic + 7
      endif
      if (tr_bgc_DON) then
         do mm = 1, n_don
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DON(mm),  nlt_bgc_DON(mm), &
                                      dontype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DON(mm)) = 2*max_algae + max_doc + max_dic + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_DON
      if (tr_bgc_Fe) then
         do mm = 1, n_fed
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Fed(mm),  nlt_bgc_Fed(mm), &
                                      fedtype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fed(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + 7 + mm
         enddo   ! mm
         do mm = 1, n_fep
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Fep(mm),  nlt_bgc_Fep(mm), &
                                      feptype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fep(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + max_fe + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_Fe 
  
      if (tr_bgc_hum) then
            call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_hum,      nlt_bgc_hum,   &
                                      humtype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
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
               call icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                         nt_zaero(mm),    nlt_zaero(mm), &
                                         zaerotype(mm),   nt_depend,     &
                                         ntrcr,           nbtrcr,        &
                                         bgc_tracer_type, trcr_depend,   &
                                         trcr_base,       n_trcr_strata, &
                                         nt_strata,       bio_index)
               bio_index_o(nlt_zaero(mm)) = 2*max_algae + max_doc + max_dic &
                                          + max_don + 2*max_fe + 7 + mm
            enddo   ! mm
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

#endif
      end subroutine icepack_init_zbgc
!=======================================================================

      subroutine icepack_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc,          nlt_bgc,       &
                                      bgctype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)

      use icepack_constants, only: c0, c1

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

      subroutine icepack_biogeochemistry(dt, &
                           ntrcr, nbtrcr,  &
                           upNO, upNH, iDi, iki, zfswin, &
                           zsal_tot, darcy_V, grow_net,  &
                           PP_net, hbri,dhbr_bot, dhbr_top, Zoo,&
                           fbio_snoice, fbio_atmice, ocean_bio, &
                           first_ice, fswpenln, bphi, bTiz, ice_bio_net,  &
                           snow_bio_net, fswthrun, Rayleigh_criteria, &
                           sice_rho, fzsal, fzsal_g, &
                           bgrid, igrid, icgrid, cgrid,  &
                           nblyr, nilyr, nslyr, n_algae, n_zaero, ncat, &
                           n_doc, n_dic,  n_don, n_fed, n_fep,  &
                           meltbn, melttn, congeln, snoicen, &
                           sst, sss, fsnow, meltsn, hmix, salinz, &
                           hin_old, flux_bio, flux_bio_atm, &
                           aicen_init, vicen_init, aicen, vicen, vsnon, &
                           aice0, trcrn, vsnon_init, skl_bgc, &
                           max_algae, max_nbtrcr, &
                           l_stop, stop_label)

      use icepack_algae, only: zbio, sklbio
      use icepack_brine, only: preflushing_changes, compute_microS_mushy
      use icepack_brine, only: update_hbrine, compute_microS 
      use icepack_parameters, only: solve_zsal, z_tracers, phi_snow
      use icepack_tracers, only: nt_fbri, tr_brine
      use icepack_tracers, only: nt_bgc_S, nt_qice, nt_sice, nt_zbgc_frac, bio_index 
      use icepack_constants, only: c0, c1, puny
      use icepack_zsalinity, only: zsalinity

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat, &
         nilyr, &
         nslyr, &
         nblyr, &
         ntrcr, &
         nbtrcr, &
         n_algae, n_zaero, &
         n_doc, n_dic,  n_don, n_fed, n_fep, &
         max_algae, max_nbtrcr

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         bgrid         , &  ! biology nondimensional vertical grid points
         igrid         , &  ! biology vertical interface points
         cgrid         , &  ! CICE vertical coordinate   
         icgrid        , &  ! interface grid for CICE (shortwave variable)
         ocean_bio     , &  ! contains all the ocean bgc tracer concentrations
         fbio_snoice   , &  ! fluxes from snow to ice
         fbio_atmice   , &  ! fluxes from atm to ice
         dhbr_top      , &  ! brine top change
         dhbr_bot      , &  ! brine bottom change
         darcy_V       , &  ! darcy velocity positive up (m/s)
         hin_old       , &  ! old ice thickness
         sice_rho      , &  ! avg sea ice density  (kg/m^3) 
         ice_bio_net   , &  ! depth integrated tracer (mmol/m^2) 
         snow_bio_net  , &  ! depth integrated snow tracer (mmol/m^2) 
         flux_bio     ! all bio fluxes to ocean

      logical (kind=log_kind), dimension (:), intent(inout) :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

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
         zsal_tot       , & ! Total ice salinity in per grid cell (g/m^2) 
         fzsal          , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g        , & ! Total gravity drainage flux
         upNO           , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         ! ammonium uptake rate (mmol/m^2/d) times aice

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria    ! .true. means Ra_c was reached  

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fswpenln        ! visible SW entering ice layers (W m-2)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         meltsn      , & ! snow melt in category n (m)
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen     , & ! snow-ice formation in category n (m)
         salinz      , & ! initial salinity  profile (ppt) 
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
         hmix    , & ! mixed layer depth (m)
         fsnow       ! snowfall rate (kg/m^2 s)

      logical (kind=log_kind), intent(in) :: &
         skl_bgc       ! if true, solve skeletal biochemistry

      logical (kind=log_kind), intent(inout) :: &  
         l_stop          ! if true, abort the model

      character (len=*), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k              , & ! vertical index
         n, mm              ! thickness category index

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
         iTin            ! Temperature on the interface grid (oC)

      real (kind=dbl_kind) :: & 
         sloss            ! brine flux contribution from surface runoff (g/m^2)

      ! for bgc sk
      real (kind=dbl_kind) :: & 
         dh_bot_chl  , & ! Chlorophyll may or may not flush
         dh_top_chl  , & ! Chlorophyll may or may not flush
         darcy_V_chl     

      l_stop = .false.

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------
         hin_old(n) = c0
         if (aicen_init(n) > puny) then 
            hin_old(n) = vicen_init(n) &
                                / aicen_init(n)
         else
            first_ice(n) = .true.
            if (tr_brine) trcrn(nt_fbri,n) = c1
            do mm = 1,nbtrcr
               trcrn(nt_zbgc_frac-1+mm,n) = zbgc_frac_init(mm)
            enddo
            if (n == 1) Rayleigh_criteria = .false.
            if (solve_zsal) trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n) = c0
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
               call preflushing_changes  (n,  aicen  (n),   &
                                 vicen   (n), vsnon  (n),   &
                                 meltbn  (n), melttn (n),   &
                                 congeln (n), snoicen(n),   &
                                 hin_old (n), dhice,        & 
                                 trcrn(nt_fbri,n),          &
                                 dhbr_top(n), dhbr_bot(n),  &
                                 hbr_old,     hin,          &
                                 hsn,         first_ice(n), &
                                 l_stop,      stop_label)

               if (l_stop) return

               if (solve_zsal)  then  

                  call compute_microS (n,         nilyr,       nblyr,             &
                                bgrid,            cgrid,       igrid,             &
                                trcrn(1:ntrcr,n), hin_old(n),  hbr_old,           &
                                sss,              sst,         bTiz(:,n),         &
                                iTin,             bphi(:,n),   kavg,              &
                                bphi_o,           phi_snow,    Rayleigh_criteria, &
                                first_ice(n),     bSin,        brine_sal,         &
                                brine_rho,        iphin,       ibrine_rho,        &
                                ibrine_sal,       sice_rho(n), sloss,             &
                                salinz(1:nilyr),  l_stop,      stop_label)

                  if (l_stop) return
               else     

                 ! Requires the average ice permeability = kavg(:)
                 ! and the surface ice porosity = zphi_o(:)
                 ! computed in "compute_microS" or from "thermosaline_vertical"

                  iDi(:,n) = c0

                  call compute_microS_mushy (n,   nilyr,         nblyr,       &
                                   bgrid,         cgrid,         igrid,       &
                                   trcrn(:,n),    hin_old(n),    hbr_old,     &
                                   sss,           sst,           bTiz(:,n),   & 
                                   iTin(:),       bphi(:,n),     kavg,        &
                                   bphi_o,        phi_snow,      bSin(:),     &
                                   brine_sal(:),  brine_rho(:),  iphin(:),    &
                                   ibrine_rho(:), ibrine_sal(:), sice_rho(n), &
                                   iDi(:,n),      l_stop,        stop_label)

               endif ! solve_zsal  

               call update_hbrine (meltbn  (n), melttn(n),   &
                                   meltsn  (n), dt,          &
                                   hin,         hsn,         &
                                   hin_old (n), hbrin,       &

                                   hbr_old,     phi_snow,    &
                                   trcrn(nt_fbri,n),         &
                                   snoicen(n),               &
                                   dhbr_top(n), dhbr_bot(n), &
                                   dh_top_chl,  dh_bot_chl,  & 
                                   kavg,        bphi_o,      &
                                   darcy_V (n), darcy_V_chl, &  
                                   bphi(2,n),   aice0,       &
                                   dh_direct)
               
               hbri = hbri + hbrin * aicen(n)  

               if (solve_zsal) then 

                  call zsalinity (n,             dt,                  &
                                  nilyr,         bgrid,               & 
                                  cgrid,         igrid,               &
                                  trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n), &
                                  trcrn(nt_qice:nt_qice+nilyr-1,n),   &
                                  trcrn(nt_sice:nt_sice+nilyr-1,n),   &
                                  ntrcr,         trcrn(nt_fbri,n),    &
                                  bSin,          bTiz(:,n),           &
                                  bphi(:,n),     iphin,               &
                                  iki(:,n),      hbr_old,             &
                                  hbrin,         hin,                 &
                                  hin_old(n),    iDi(:,n),            &
                                  darcy_V(n),    brine_sal,           & 
                                  brine_rho,     ibrine_sal,          & 
                                  ibrine_rho,    dh_direct,           &
                                  Rayleigh_criteria,                  &
                                  first_ice(n),  sss,                 &
                                  sst,           dhbr_top(n),         &
                                  dhbr_bot(n),                        &
                                  l_stop,        stop_label,          &
                                  fzsal,         fzsal_g,             &
                                  bphi_o,        nblyr,               & 
                                  vicen(n),      aicen_init(n),       &
                                  zsal_tot) 

                  if (l_stop) return

               endif  ! solve_zsal

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
                          bio_index(1:nbtrcr),   aicen_init(n),          &
                          vicen_init(n),         vsnon_init(n),          &
                          vicen(n),              vsnon(n),               &
                          aicen(n),              flux_bio_atm(1:nbtrcr), &
                          n,                     n_algae,                &
                          n_doc,                 n_dic,                  &
                          n_don,                                         &
                          n_fed,                 n_fep,                  &
                          n_zaero,               first_ice(n),           &
                          hin_old(n),            ocean_bio(1:nbtrcr),    &
                          bphi(:,n),             iphin,                  &     
                          iDi(:,n),              sss,                    &
                          fswpenln(:,n),                                 &
                          dhbr_top(n),           dhbr_bot(n),            &
                          dh_top_chl,            dh_bot_chl,             &
                          zfswin(:,n),                                   &
                          hbrin,                 hbr_old,                &
                          darcy_V(n),            darcy_V_chl,            &
                          bgrid,                 cgrid,                  &
                          igrid,                 icgrid,                 &
                          bphi_o,                                        &
                          dhice,                 iTin,                   &
                          Zoo(:,n),                                      &
                          flux_bio(1:nbtrcr),    dh_direct,              &
                          upNO,                  upNH,                   &
                          fbio_snoice,           fbio_atmice,            &
                          PP_net,                ice_bio_net (1:nbtrcr), &
                          snow_bio_net(1:nbtrcr),grow_net,               &
                          l_stop,                stop_label)
            
               if (l_stop) return
     
            elseif (skl_bgc) then

               call sklbio (dt,                      ntrcr,               &
                            nilyr,                                        &
                            nbtrcr,                  n_algae,             &
                            n_zaero,                 n_doc,               &
                            n_dic,                   n_don,               &
                            n_fed,                   n_fep,               &
                            flux_bio (1:nbtrcr),     ocean_bio(1:nbtrcr), &
                            hmix,                    aicen    (n),        &
                            meltbn   (n),            congeln  (n),        &
                            fswthrun (n),            first_ice(n),        &
                            trcrn    (1:ntrcr,n),    hin,                 &
                            PP_net,                  upNO,                &
                            upNH,                    grow_net,            &
                            l_stop,                  stop_label)

               if (l_stop) return

            endif  ! skl_bgc

            first_ice(n) = .false.

         endif             ! aicen > puny
      enddo                ! ncat

      end subroutine icepack_biogeochemistry

!=======================================================================

! basic initialization for ocean_bio_all

      subroutine icepack_init_OceanConcArray(max_nbtrcr, &
          max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
          nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, ocean_bio_all, hum)

      use icepack_zbgc_shared, only: R_CHL2N

      integer (kind=int_kind), intent(in) :: &
         max_algae   , & ! maximum number of algal types 
         max_dic     , & ! maximum number of dissolved inorganic carbon types 
         max_doc     , & ! maximum number of dissolved organic carbon types
         max_don     , & ! maximum number of dissolved organic nitrogen types
         max_fe      , & ! maximum number of iron types
         max_aero    , & ! maximum number of aerosols 
         max_nbtrcr      ! maximum number of bio tracers

      real (kind=dbl_kind), intent(in) :: &
         nit         , & ! ocean nitrate (mmol/m^3)          
         amm         , & ! ammonia/um (mmol/m^3)
         sil         , & ! silicate (mmol/m^3)
         dmsp        , & ! dmsp (mmol/m^3)
         dms         , & ! dms (mmol/m^3)
         hum             ! humic material (mmol/m^3)

      real (kind=dbl_kind), dimension (max_algae), intent(in) :: &
         algalN          ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (max_doc), intent(in) :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (max_don), intent(in) :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_dic), intent(in) :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_fe), intent(in) :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (max_aero), intent(in) :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_nbtrcr), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      ! local variables

      integer (kind=int_kind) :: &
         k, ks           ! tracer indices

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

      end subroutine icepack_init_OceanConcArray

!=======================================================================

!  Initialize ocean concentration

      subroutine icepack_init_ocean_conc (amm, dmsp, dms, algalN, doc, dic, don, &
             fed, fep, hum, nit, sil, zaeros, max_dic, max_don, max_fe, max_aero,&
             CToN, CToN_DON)

      use icepack_zbgc_shared, only: R_C2N, R_C2N_DON

      integer (kind=int_kind), intent(in) :: &
        max_dic, &
        max_don, &
        max_fe, &
        max_aero

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

      integer (kind=int_kind) :: &
        k 

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
            dic(k) = c1
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
 

      end subroutine icepack_init_ocean_conc

!=======================================================================

      end module icepack_zbgc

!=======================================================================
