!=======================================================================

! Diagnostic information output during run for biogeochemistry
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      module icedrv_diagnostics_bgc

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, nu_diag_out
      use icedrv_constants, only: c0, mps_to_cmpdy, c100, p5, c1
      use icedrv_domain_size, only: nx
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icepack_intfc, only: icepack_max_algae, icepack_max_aero, icepack_max_fe
      use icepack_intfc, only: icepack_max_doc, icepack_max_don
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private
      public :: hbrine_diags, bgc_diags, zsal_diags

!=======================================================================

      contains

!=======================================================================

!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Nicole Jeffery, LANL

      subroutine hbrine_diags ()
              
      use icedrv_arrays_column, only: darcy_V
      use icedrv_diagnostics, only: nx_names
      use icedrv_domain_size, only: nilyr
      use icedrv_state, only: aice, aicen, vicen, vice, trcr, trcrn

      ! local variables

      integer (kind=int_kind) :: &
         k, n

      integer (kind=int_kind) :: &
         ktherm

      integer (kind=int_kind) :: &
         nt_sice, &
         nt_fbri

      ! fields at diagnostic points
      real (kind=dbl_kind) :: &
         phinS, phinS1, pdarcy_V, pfbri

      real (kind=dbl_kind), dimension(nilyr) :: &
         pSin, pSin1

      character(len=*), parameter :: subname='(hbrine_diags)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_parameters(ktherm_out=ktherm)
         call icepack_query_tracer_indices(nt_sice_out=nt_sice, nt_fbri_out=nt_fbri)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Dynamic brine height
      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)

         do n = 1, nx
            phinS1   = c0
            phinS    = c0
            pfbri    = trcrn(n,nt_fbri,1)
            pdarcy_V = darcy_V(n,1)
            if (aice (n)  > c0) phinS  = trcr (n,nt_fbri  )*vice (n  )/aice (n  )
            if (aicen(n,1)> c0) phinS1 = trcrn(n,nt_fbri,1)*vicen(n,1)/aicen(n,1)
            do k = 1 , nilyr
               pSin (k) = trcr (n,nt_sice+k-1  )
               pSin1(k) = trcrn(n,nt_sice+k-1,1)
            enddo

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

         write(nu_diag_out+n-1,899) nx_names(n)
        
         write(nu_diag_out+n-1,*) '------ hbrine ------'
         write(nu_diag_out+n-1,900) 'hbrine, (m)        = ',phinS
         write(nu_diag_out+n-1,900) 'fbri, cat1 (m)     = ',pfbri
         write(nu_diag_out+n-1,900) 'hbrine cat1, (m)   = ',phinS1
         write(nu_diag_out+n-1,900) 'darcy_V cat1, (m/s)= ',pdarcy_V
         if (ktherm == 2) then          
            write(nu_diag_out+n-1,*) '                         '
            write(nu_diag_out+n-1,*) '------ Thermosaline Salinity ------'
            write(nu_diag_out+n-1,803) 'Sice1 cat1 S (ppt)'
            write(nu_diag_out+n-1,*) '---------------------------------------------------'
            write(nu_diag_out+n-1,802) (pSin1(k), k = 1,nilyr)
            write(nu_diag_out+n-1,*) '                         '
            write(nu_diag_out+n-1,803) 'Sice bulk S (ppt) '
            write(nu_diag_out+n-1,*) '---------------------------------------------------'
            write(nu_diag_out+n-1,802) (pSin(k), k = 1,nilyr)
            write(nu_diag_out+n-1,*) '                         '
         endif
      enddo ! nx

802   format (f24.17,2x)
803   format (a25,2x)
899   format (43x,a24)
900   format (a25,2x,f24.17)

      end subroutine hbrine_diags

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Nicole Jeffery, LANL

      subroutine bgc_diags ()

      use icedrv_arrays_column, only: ocean_bio, zfswin, fbio_atmice, fbio_snoice
      use icedrv_arrays_column, only: Zoo, grow_net, ice_bio_net, trcrn_sw
      use icedrv_domain_size,   only: ncat, nblyr, n_algae, n_zaero
      use icedrv_domain_size,   only: n_doc, n_don, n_fed, n_fep, nilyr, nslyr
      use icedrv_flux,  only: flux_bio, flux_bio_atm
      use icedrv_state, only: vicen, vice, trcr

      ! local variables

      integer (kind=int_kind) :: &
         k, n, nn, kk, klev

      logical (kind=log_kind) :: &
         skl_bgc, z_tracers, dEdd_algae

      ! fields at diagnostic points
      real (kind=dbl_kind) :: &
         pNit_sk, pAm_sk, pSil_sk, phum_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac, &
         pflux_NO, pflux_Am,  phum_ac, &
         pflux_snow_NO, pflux_snow_Am,  &
         pflux_atm_NO, pflux_atm_Am,  pgrow_net, &
         pflux_hum

      logical (kind=log_kind) :: &
         tr_bgc_DMS, tr_bgc_PON, &
         tr_bgc_N, tr_bgc_C, tr_bgc_DON, tr_zaero, tr_bgc_hum, &
         tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_Fe

      integer (kind=int_kind) :: &
         nt_fbri, nt_sice, nt_bgc_nit, nt_bgc_am, nt_bgc_sil, &
         nt_bgc_hum, nt_bgc_pon, nt_bgc_dmspp, nt_bgc_dmspd, nt_bgc_dms, &
         nlt_bgc_hum, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, nlt_chl_sw, &
         nlt_bgc_DMSPp, nlt_bgc_DMS
      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_n, nlt_bgc_N
      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nt_bgc_doc, nlt_bgc_DOC
      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nt_bgc_don, nlt_bgc_DON 
      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero, nlt_zaero, nlt_zaero_sw
      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nt_bgc_fed, nt_bgc_fep, nlt_bgc_Fed, nlt_bgc_Fep

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         pN_ac, pN_tot, pN_sk, pflux_N
      real (kind=dbl_kind), dimension(icepack_max_doc) :: &
         pDOC_ac, pDOC_sk
      real (kind=dbl_kind), dimension(icepack_max_don) :: &
         pDON_ac, pDON_sk
      real (kind=dbl_kind), dimension(icepack_max_fe ) :: &
         pFed_ac,  pFed_sk, pFep_ac, pFep_sk 
      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
        pflux_zaero, pflux_snow_zaero, pflux_atm_zaero, &
        pflux_atm_zaero_s

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(2) :: &
         pNOs, pAms, pPONs, phums
      real (kind=dbl_kind), dimension(2,icepack_max_algae) :: &
         pNs
      real (kind=dbl_kind), dimension(2,icepack_max_doc) :: &
         pDOCs
      real (kind=dbl_kind), dimension(2,icepack_max_don) :: &
         pDONs
      real (kind=dbl_kind), dimension(2,icepack_max_fe ) :: &
         pFeds, pFeps 
      real (kind=dbl_kind), dimension(2,icepack_max_aero) :: &
         pzaeros
      real (kind=dbl_kind), dimension(nblyr+1) :: &
         pNO, pAm, pPON, pzfswin, pZoo, phum
      real (kind=dbl_kind), dimension(nblyr+1,icepack_max_algae) :: &
         pN
      real (kind=dbl_kind), dimension(nblyr+1,icepack_max_aero) :: &
         pzaero
      real (kind=dbl_kind), dimension(nblyr+1,icepack_max_doc) :: &
         pDOC
      real (kind=dbl_kind), dimension(nblyr+1,icepack_max_don) :: &
         pDON
      real (kind=dbl_kind), dimension(nblyr+1,icepack_max_fe ) :: &
         pFed, pFep 
      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace
      real (kind=dbl_kind), dimension (nslyr+nilyr+2) :: &
         pchlsw
      real (kind=dbl_kind), dimension(nslyr+nilyr+2,icepack_max_aero) :: &
         pzaerosw

      character(len=*), parameter :: subname='(bgc_diags)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(skl_bgc_out=skl_bgc, &
         z_tracers_out=z_tracers, dEdd_algae_out=dEdd_algae)
      call icepack_query_tracer_flags( &
         tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
         tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, &
         tr_bgc_DON_out=tr_bgc_DON, tr_zaero_out=tr_zaero, &
         tr_bgc_hum_out=tr_bgc_hum,tr_bgc_Sil_out=tr_bgc_Sil, &
         tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Am_out=tr_bgc_Am, &
         tr_bgc_Fe_out=tr_bgc_Fe)
      call icepack_query_tracer_indices( &
         nt_fbri_out=nt_fbri, nt_sice_out=nt_sice, nt_zaero_out=nt_zaero, &
         nt_bgc_n_out=nt_bgc_n, nt_bgc_doc_out=nt_bgc_doc, &
         nt_bgc_don_out=nt_bgc_don, nt_bgc_sil_out=nt_bgc_sil, &
         nt_bgc_fed_out=nt_bgc_fed, nt_bgc_fep_out=nt_bgc_fep, &
         nt_bgc_nit_out=nt_bgc_nit, nt_bgc_am_out=nt_bgc_am, &
         nt_bgc_hum_out=nt_bgc_hum, nt_bgc_pon_out=nt_bgc_pon, &
         nt_bgc_dmspp_out=nt_bgc_dmspp, nt_bgc_dmspd_out=nt_bgc_dmspd, &
         nt_bgc_dms_out=nt_bgc_dms, nlt_chl_sw_out=nlt_chl_sw, &
         nlt_bgc_N_out=nlt_bgc_N, nlt_zaero_out=nlt_zaero, &
         nlt_zaero_sw_out=nlt_zaero_sw, &
         nlt_bgc_Fed_out=nlt_bgc_Fed, nlt_bgc_Fep_out=nlt_bgc_Fep, &
         nlt_bgc_hum_out=nlt_bgc_hum, nlt_bgc_Nit_out=nlt_bgc_Nit, &
         nlt_bgc_Am_out=nlt_bgc_Am, nlt_bgc_Sil_out=nlt_bgc_Sil, &
         nlt_bgc_DOC_out=nlt_bgc_DOC, nlt_bgc_DON_out=nlt_bgc_DON, &
         nlt_bgc_DMSPp_out=nlt_bgc_DMSPp, nlt_bgc_DMS_out=nlt_bgc_DMS)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = zspace(1)*p5
      zspace(nblyr+1) = zspace(nblyr+1)*p5
      klev = 1+nilyr+nslyr

      !-----------------------------------------------------------------
      ! biogeochemical state of the ice
      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)

         do n = 1, nx
            pAm_ac   = c0
            pSil_ac  = c0
            phum_ac  = c0
            pDMSP_ac = c0
            pDMS_ac  = c0
            pN_ac(:)  = c0
            pDOC_ac(:)= c0
            pDON_ac(:)= c0
            pFed_ac(:)= c0
            pFep_ac(:)= c0
            pNit_ac = c0
            if (tr_bgc_N) then
               do k = 1,n_algae
                  pN_ac(k)    = ocean_bio(n,nlt_bgc_N(k))
               enddo  !n_algae
            endif    !tr_bgc_N
            if (tr_bgc_C) then
               do k = 1,n_doc
                  pDOC_ac(k)    = ocean_bio(n,nlt_bgc_DOC(k))
               enddo  !n_algae
            endif    !tr_bgc_N
            if (tr_bgc_DON) then
               do k = 1,n_don
                  pDON_ac(k)    = ocean_bio(n,nlt_bgc_DON(k))
               enddo 
            endif
            if (tr_bgc_Fe ) then
               do k = 1,n_fed 
                  pFed_ac (k)   = ocean_bio(n,nlt_bgc_Fed (k))
               enddo 
               do k = 1,n_fep 
                  pFep_ac (k)   = ocean_bio(n,nlt_bgc_Fep (k))
               enddo 
            endif
            if (tr_bgc_Nit) pNit_ac  = ocean_bio(n,nlt_bgc_Nit)
            if (tr_bgc_Am)  pAm_ac   = ocean_bio(n,nlt_bgc_Am)
            if (tr_bgc_Sil) pSil_ac  = ocean_bio(n,nlt_bgc_Sil)
            if (tr_bgc_hum) phum_ac  = ocean_bio(n,nlt_bgc_hum)
            if (tr_bgc_DMS) then
               pDMSP_ac = ocean_bio(n,nlt_bgc_DMSPp)
               pDMS_ac  = ocean_bio(n,nlt_bgc_DMS)
            endif

            ! fluxes in mmol/m^2/d
            ! concentrations are bulk in mmol/m^3
            ! iron is in 10^-3 mmol/m^3  (nM)

            if (skl_bgc) then
               pNit_sk   = c0
               pAm_sk    = c0
               pSil_sk   = c0
               phum_sk   = c0
               pDMSPp_sk = c0
               pDMSPd_sk = c0
               pDMS_sk   = c0
               pN_sk  (:) = c0
               pflux_N(:) = c0
               pDOC_sk(:) = c0
               pDON_sk(:) = c0
               pFed_sk(:) = c0
               pFep_sk(:) = c0
               pgrow_net =  grow_net(n)

               do k = 1,n_algae
                  pN_sk(k)       = trcr    (n, nt_bgc_N(k))
                  pflux_N(k)     = flux_bio(n,nlt_bgc_N(k))*mps_to_cmpdy/c100
               enddo
               if (tr_bgc_C) then
                  do k = 1,n_doc
                     pDOC_sk(k)  = trcr    (n,nt_bgc_DOC(k))
                  enddo
               endif
               if (tr_bgc_DON) then
                  do k = 1,n_don
                     pDON_sk(k)  = trcr    (n,nt_bgc_DON(k))
                  enddo
               endif
               if (tr_bgc_Fe ) then
                  do k = 1,n_fed
                     pFed_sk (k) = trcr    (n,nt_bgc_Fed(k))
                  enddo
                  do k = 1,n_fep
                     pFep_sk (k) = trcr    (n,nt_bgc_Fep(k))
                  enddo
               endif
               if (tr_bgc_Nit) then
                  pNit_sk        = trcr    (n, nt_bgc_Nit)
                  pflux_NO       = flux_bio(n,nlt_bgc_Nit)*mps_to_cmpdy/c100
               endif
               if (tr_bgc_Am) then
                  pAm_sk         = trcr    (n, nt_bgc_Am)
                  pflux_Am       = flux_bio(n,nlt_bgc_Am)*mps_to_cmpdy/c100
               endif
               if (tr_bgc_Sil) then
                  pSil_sk        = trcr    (n, nt_bgc_Sil)
               endif
               if (tr_bgc_hum) then
                  phum_sk        = trcr    (n, nt_bgc_hum)
                  pflux_hum      = flux_bio(n,nlt_bgc_hum)*mps_to_cmpdy/c100
               endif
               if (tr_bgc_DMS) then
                  pDMSPp_sk      = trcr   (n,nt_bgc_DMSPp)
                  pDMSPd_sk      = trcr   (n,nt_bgc_DMSPd)
                  pDMS_sk        = trcr   (n,nt_bgc_DMS)
               endif

            elseif (z_tracers) then   ! zbgc
               pflux_NO             = c0
               pN_tot           (:) = c0
               pflux_Am             = c0
               pflux_hum            = c0
               pflux_atm_Am         = c0
               pflux_snow_Am        = c0
               pflux_N          (:) = c0
               pflux_NO             = c0
               pflux_atm_NO         = c0
               pflux_snow_NO        = c0
               pflux_zaero      (:) = c0
               pflux_atm_zaero_s(:) = c0
               pflux_atm_zaero  (:) = c0
               pflux_snow_zaero (:) = c0

               if (tr_bgc_Nit) then
                  pflux_NO       = flux_bio(n,nlt_bgc_Nit)*mps_to_cmpdy/c100
                  pflux_atm_NO   = fbio_atmice(n,nlt_bgc_Nit)*mps_to_cmpdy/c100
                  pflux_snow_NO  = fbio_snoice(n,nlt_bgc_Nit)*mps_to_cmpdy/c100
               endif
               if (tr_bgc_Am) then
                  pflux_Am       = flux_bio(n,nlt_bgc_Am)*mps_to_cmpdy/c100
                  pflux_atm_Am   = fbio_atmice(n,nlt_bgc_Am)*mps_to_cmpdy/c100
                  pflux_snow_Am  = fbio_snoice(n,nlt_bgc_Am)*mps_to_cmpdy/c100
               endif
               if (tr_bgc_hum) then
                  pflux_hum       = flux_bio(n,nlt_bgc_hum)*mps_to_cmpdy/c100
               endif
               if (tr_bgc_N)  then
                  do k = 1,n_algae
                     pflux_N(k)   = flux_bio(n,nlt_bgc_N(k))*mps_to_cmpdy/c100
                  enddo
               endif
               if (tr_zaero)  then
                  do k = 1,n_zaero
                     pflux_zaero(k)      = flux_bio(n,nlt_zaero(k))*mps_to_cmpdy/c100
                     pflux_atm_zaero_s(k)= flux_bio_atm(n,nlt_zaero(k))*mps_to_cmpdy/c100 !*aice
                     pflux_atm_zaero(k)  = fbio_atmice(n,nlt_zaero(k))*mps_to_cmpdy/c100
                     pflux_snow_zaero(k) = fbio_snoice(n,nlt_zaero(k))*mps_to_cmpdy/c100
                  enddo
               endif
               do k = 1, nblyr+1
                  pzfswin(k) = c0
                  pZoo   (k) = c0
                  do nn = 1,ncat
                     pzfswin(k) = pzfswin(k) + zfswin(n,k,nn)*vicen(n,nn)
                     pZoo   (k) = pZoo   (k) + Zoo   (n,k,nn)*vicen(n,nn)
                  enddo !nn
                  if (vice(n) > c0) then
                     pzfswin(k) = pzfswin(k)/vice(n)
                     pZoo   (k) = pZoo   (k)/vice(n)
                  endif !vice
                  pAm   (k  ) = c0
                  pN    (k,:) = c0
                  pDOC  (k,:) = c0
                  pDON  (k,:) = c0
                  pFed  (k,:) = c0
                  pFep  (k,:) = c0
                  pzaero(k,:) = c0
                  pPON  (k  ) = c0
                  phum  (k  ) = c0
                  pNO   (k  ) = c0
                  if (tr_bgc_Nit) pNO(k) =  trcr(n,nt_bgc_Nit+k-1)
                  if (tr_bgc_Am)  pAm(k) =  trcr(n,nt_bgc_Am+k-1)
                  if (tr_bgc_N) then
                     do nn = 1, n_algae
                        pN    (k,nn) = trcr(n,nt_bgc_N(nn)+k-1)
                     enddo
                  endif
                  if (tr_bgc_C) then
                     do nn = 1, n_doc
                        pDOC  (k,nn) = trcr(n,nt_bgc_DOC(nn)+k-1)
                     enddo
                  endif
                  if (tr_bgc_DON) then
                     do nn = 1, n_don
                        pDON  (k,nn) = trcr(n,nt_bgc_DON(nn)+k-1)
                     enddo
                  endif
                  if (tr_bgc_Fe)  then
                     do nn = 1, n_fed
                        pFed  (k,nn) = trcr(n,nt_bgc_Fed(nn)+k-1)
                     enddo
                     do nn = 1, n_fep
                        pFep  (k,nn) = trcr(n,nt_bgc_Fep(nn)+k-1)
                     enddo
                  endif
                  if (tr_zaero) then
                     do nn = 1, n_zaero
                        pzaero(k,nn) = trcr(n,nt_zaero(nn)+k-1)
                     enddo
                  endif
                  if (tr_bgc_PON) pPON(k) = trcr(n,nt_bgc_PON+k-1)
                  if (tr_bgc_hum) phum(k) = trcr(n,nt_bgc_hum+k-1)
               enddo  ! k
               if (tr_bgc_N) then
                  do nn = 1,n_algae
                     pN_tot(nn) = ice_bio_net(n,nlt_bgc_N(nn))
                  enddo
                  pgrow_net =  grow_net(n)
               endif  ! tr_bgc_N
               do k = 1 , 2  ! snow concentration
                  pAms   (k  ) = c0
                  pNs    (k,:) = c0
                  pDOCs  (k,:) = c0
                  pDONs  (k,:) = c0
                  pFeds  (k,:)= c0
                  pFeps  (k,:)= c0
                  pzaeros(k,:) = c0
                  pPONs  (k  ) = c0
                  phums  (k  ) = c0
                  pNOs   (k  ) = c0
                  if (tr_bgc_Nit) pNOs(k) = trcr(n,nt_bgc_Nit+nblyr+k)
                  if (tr_bgc_Am)  pAms(k) = trcr(n,nt_bgc_Am +nblyr+k)
                  if (tr_bgc_N) then
                     do nn = 1, n_algae
                        pNs    (k,nn) = trcr(n,nt_bgc_N(nn)+nblyr+k)
                     enddo
                  endif
                  if (tr_bgc_C) then
                     do nn = 1, n_doc
                        pDOCs  (k,nn) = trcr(n,nt_bgc_DOC(nn)+nblyr+k)
                     enddo
                  endif
                  if (tr_bgc_DON) then
                     do nn = 1, n_don
                        pDONs  (k,nn) = trcr(n,nt_bgc_DON(nn)+nblyr+k)
                     enddo
                  endif
                  if (tr_bgc_Fe ) then
                     do nn = 1, n_fed
                        pFeds  (k,nn) = trcr(n,nt_bgc_Fed(nn)+nblyr+k)
                     enddo
                     do nn = 1, n_fep
                        pFeps  (k,nn) = trcr(n,nt_bgc_Fep(nn)+nblyr+k)
                     enddo
                  endif
                  if (tr_zaero) then
                     do nn = 1, n_zaero
                        pzaeros(k,nn) = trcr(n,nt_zaero(nn)+nblyr+k)
                     enddo
                  endif
                  if (tr_bgc_PON)pPONs(k) =trcr(n,nt_bgc_PON+nblyr+k)
                  if (tr_bgc_hum)phums(k) =trcr(n,nt_bgc_hum+nblyr+k)
               enddo   ! k
            endif      ! ztracers
            pchlsw(:) = c0
            pzaerosw(:,:) = c0
            if (dEdd_algae) then
               do k = 0, klev
                  if (tr_bgc_N) pchlsw(k+1) = trcrn_sw(n,nlt_chl_sw+k,1)
                  if (tr_zaero) then
                     do nn = 1,n_zaero
                        pzaerosw(k+1,nn) =  trcrn_sw(n,nlt_zaero_sw(nn) + k,1)
                     enddo
                  endif
               enddo
            endif              ! dEdd_algae

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

       if (z_tracers) then
         write(nu_diag_out+n-1,803) 'zfswin PAR  '
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,802) (pzfswin(k), k = 1,nblyr+1)
         write(nu_diag_out+n-1,*) '      '
         write(nu_diag_out+n-1,803) 'Losses: Zoo(mmol/m^3)  '
         write(nu_diag_out+n-1,803) '        Brine Conc.       ',' Brine Conc'
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,802) (pZoo(k), k = 1,nblyr+1)
         write(nu_diag_out+n-1,*) '      '
       endif
       if (tr_bgc_Nit) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '    nitrate conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag_out+n-1,900) 'Ocean conc       = ',pNit_ac
         write(nu_diag_out+n-1,900) 'ice-ocean flux   = ',pflux_NO
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',pNit_sk
         elseif (z_tracers) then
           write(nu_diag_out+n-1,900) 'atm-ice flux     = ',pflux_atm_NO
           write(nu_diag_out+n-1,900) 'snow-ice flux    = ',pflux_snow_NO
           write(nu_diag_out+n-1,*) '             snow + ice conc'
           write(nu_diag_out+n-1,803) '    nitrate'
           write(nu_diag_out+n-1,802) (pNOs(k), k = 1,2)
           write(nu_diag_out+n-1,802) (pNO(k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '    '
         endif
      endif
      if (tr_bgc_PON .and. z_tracers) then
           write(nu_diag_out+n-1,*) '---------------------------------------------------'
           write(nu_diag_out+n-1,*) '    PON snow + ice conc. (mmol/m^3)'
           write(nu_diag_out+n-1,803) '    PON'
           write(nu_diag_out+n-1,802) (pPONs(k), k = 1,2)
           write(nu_diag_out+n-1,802) (pPON(k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) ' '
      endif
      if (tr_bgc_hum) then
           write(nu_diag_out+n-1,*) '---------------------------------------------------'
           write(nu_diag_out+n-1,*) '    hum snow + ice conc. (mmolC/m^3)'
           write(nu_diag_out+n-1,900) 'Ocean conc       = ',phum_ac
           write(nu_diag_out+n-1,900) 'ice-ocean flux   = ',pflux_hum
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',phum_sk
         elseif (z_tracers) then
           write(nu_diag_out+n-1,803) '    hum'
           write(nu_diag_out+n-1,802) (phums(k), k = 1,2)
           write(nu_diag_out+n-1,802) (phum(k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) ' '
         endif
      endif
      if (tr_bgc_Am) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '    ammonium conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag_out+n-1,900) 'Ocean conc       = ',pAm_ac
         write(nu_diag_out+n-1,900) 'ice-ocean flux   = ',pflux_Am
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',pAm_sk
         elseif (z_tracers) then
           write(nu_diag_out+n-1,900) 'atm-ice flux     = ',pflux_atm_Am
           write(nu_diag_out+n-1,900) 'snow-ice flux    = ',pflux_snow_Am
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  ammonium'
           write(nu_diag_out+n-1,802) (pAms(k), k = 1,2)
           write(nu_diag_out+n-1,802) (pAm(k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '       '
         endif
      endif
      if (tr_bgc_N) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,900) 'tot algal growth (1/d) = ',pgrow_net
       do kk = 1,n_algae
         write(nu_diag_out+n-1,*) '  algal conc. (mmol N/m^3) or flux (mmol N/m^2/d)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900) 'Ocean conc           = ',pN_ac(kk)
         write(nu_diag_out+n-1,900) 'ice-ocean flux       = ',pflux_N(kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',pN_sk(kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,900) 'Tot ice (mmolN/m^2) = ',pN_tot(kk)
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  algal N'
           write(nu_diag_out+n-1,802) (pNs(k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pN(k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '         '
         endif
       enddo
      endif
      if (tr_bgc_C) then
       do kk = 1,1 !n_doc
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  DOC conc. (mmol C/m^3)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pDOC_ac(kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pDOC_sk(kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  DOC'
           write(nu_diag_out+n-1,802) (pDOCs(k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pDOC(k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
      endif
      if (tr_bgc_DON) then
       do kk = 1,n_don
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  DON conc. (mmol N/m^3)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pDON_ac(kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pDON_sk(kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  DON'
           write(nu_diag_out+n-1,802) (pDONs(k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pDON(k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
      endif
      if (tr_bgc_Fe ) then
       do kk = 1,n_fed
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) ' dFe  conc. (nM)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pFed_ac (kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pFed_sk (kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  Fed'
           write(nu_diag_out+n-1,802) (pFeds (k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pFed (k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
       do kk = 1,n_fep
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) ' pFe  conc. (nM)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pFep_ac (kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pFep_sk (kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  Fep'
           write(nu_diag_out+n-1,802) (pFeps (k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pFep (k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
      endif
      if (tr_bgc_DMS) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '    DMS (mmol/m^3)      '
         write(nu_diag_out+n-1,900) 'Ocean DMSP    = ',pDMSP_ac
         write(nu_diag_out+n-1,900) 'Ocean DMS     = ',pDMS_ac
         if (skl_bgc) then
          write(nu_diag_out+n-1,900) 'Ice DMSPp    = ',pDMSPp_sk
          write(nu_diag_out+n-1,900) 'Ice DMSPd    = ',pDMSPd_sk
          write(nu_diag_out+n-1,900) 'Ice DMS      = ',pDMS_sk
         endif
      endif
      if (tr_zaero .and. z_tracers) then
       do kk = 1,n_zaero
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  aerosol conc. (kg/m^3) or flux (kg/m^2/d)'
         write(nu_diag_out+n-1,1020) '  type: ',kk
         write(nu_diag_out+n-1,900) 'Atm source flux     = ',pflux_atm_zaero_s(kk)
         write(nu_diag_out+n-1,900) 'ice-ocean flux*aice = ',pflux_zaero(kk)
         write(nu_diag_out+n-1,900) 'atm-ice flux*aice   = ',pflux_atm_zaero(kk)
         write(nu_diag_out+n-1,900) 'snow-ice flux*aice  = ',pflux_snow_zaero(kk)
         write(nu_diag_out+n-1,*) '             snow + ice conc.'
         write(nu_diag_out+n-1,803) ' aerosol'
         write(nu_diag_out+n-1,802) (pzaeros(k,kk), k = 1,2)
         write(nu_diag_out+n-1,802) (pzaero(k,kk), k = 1,nblyr+1)
         write(nu_diag_out+n-1,*) '            '
       enddo
      endif
      if (dEdd_algae) then
          if (tr_zaero) then
          do kk = 1,n_zaero
          write(nu_diag_out+n-1,*) '---------------------------------------------------'
          write(nu_diag_out+n-1,*) '  Cat 1 aerosol conc. (kg/m^3) on delta-Eddington grid  '
          write(nu_diag_out+n-1,802) (pzaerosw(k,kk), k = 1,klev +1)
          enddo
          endif
         if (tr_bgc_N) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  Cat 1 chl (mg/m^3) on delta-Eddington grid  '
         write(nu_diag_out+n-1,802) (pchlsw(k), k = 1,klev +1)
         endif
      endif
      enddo                  ! nx

802   format (f24.17,2x)
803   format (a25,2x)
900   format (a25,2x,f24.17)
1020  format (a30,2x,i6)    ! integer

      end subroutine bgc_diags

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Nicole Jeffery, LANL

      subroutine zsal_diags ()

      use icedrv_arrays_column, only: fzsal, fzsal_g, sice_rho, bTiz
      use icedrv_arrays_column, only: iDi, bphi, dhbr_top, dhbr_bot, darcy_V
      use icedrv_domain_size, only: nblyr, ncat, nilyr
      use icedrv_state, only: aicen, aice, vice, trcr, trcrn, vicen, vsno

      ! local variables

      integer (kind=int_kind) :: &
         k, n, nn

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(nx) :: &
         phinS, phinS1, &
         phbrn,pdh_top1,pdh_bot1, psice_rho, pfzsal, &
         pfzsal_g, pdarcy_V1

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(nblyr+2) :: &
         pphin, pphin1
      real (kind=dbl_kind), dimension(nblyr) :: &
         pSin, pSice, pSin1
      real (kind=dbl_kind), dimension(nblyr+1) :: &
         pbTiz, piDin

      real (kind=dbl_kind) :: &
         rhosi, rhow, rhos

      logical (kind=log_kind) :: tr_brine
      integer (kind=int_kind) :: nt_fbri, nt_bgc_S, nt_sice

      character(len=*), parameter :: subname='(zsal_diags)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_parameters(rhosi_out=rhosi, rhow_out=rhow, rhos_out=rhos)
         call icepack_query_tracer_flags(tr_brine_out=tr_brine)
         call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, &
             nt_bgc_S_out=nt_bgc_S, nt_sice_out=nt_sice)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! salinity and microstructure of the ice
      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)

         do n = 1, nx
            pfzsal    = fzsal(n)
            pfzsal_g  = fzsal_g(n)
            phinS     = c0
            phinS1    = c0
            phbrn     = c0
            psice_rho = c0
            pdh_top1  = c0
            pdh_bot1  = c0
            pdarcy_V1 = c0
            do nn = 1, ncat
               psice_rho = psice_rho(n) + sice_rho(n,nn)*aicen(n,nn)
            enddo
            if (aice(n) > c0) psice_rho = psice_rho/aice(n)
            if (tr_brine .and. aice(n) > c0) &
               phinS = trcr(n,nt_fbri)*vice(n)/aice(n)
            if (aicen(n,1) > c0) then
               if (tr_brine) phinS1 = trcrn(n,nt_fbri,1) &
                                    * vicen(n,1)/aicen(n,1)
               pdh_top1  = dhbr_top(n,1)
               pdh_bot1  = dhbr_bot(n,1)
               pdarcy_V1 = darcy_V(n,1)
            endif
            if (tr_brine .AND. aice(n) > c0) &
               phbrn = (c1 - rhosi/rhow)*vice(n)/aice(n) &
                           - rhos /rhow *vsno(n)/aice(n)
            do k = 1, nblyr+1
               pbTiz(k) = c0
               piDin(k) = c0
               do nn = 1,ncat
                  pbTiz(k) = pbTiz(k) + bTiz(n,k,nn)*vicen(n,nn)
                  piDin(k) = piDin(k) +  iDi(n,k,nn)*vicen(n,nn)
               enddo
               if (vice(n) > c0) then
                  pbTiz(k) = pbTiz(k)/vice(n)
                  piDin(k) = piDin(k)/vice(n)
               endif
            enddo                 ! k
            do k = 1, nblyr+2
               pphin (k) = c0
               pphin1(k) = c0
               if (aicen(n,1) > c0) pphin1(k) = bphi(n,k,1)
               do nn = 1,ncat
                  pphin(k) = pphin(k) + bphi(n,k,nn)*vicen(n,nn)
               enddo
               if (vice(n) > c0) pphin(k) = pphin(k)/vice(n)
            enddo
            do k = 1,nblyr
               pSin (k) = c0
               pSin1(k) = c0
               pSin (k) = trcr(n,nt_bgc_S+k-1)
               if (aicen(n,1) > c0) pSin1(k) = trcrn(n,nt_bgc_S+k-1,1)
            enddo
            do k = 1,nilyr
               pSice(k) = trcr(n,nt_sice+k-1)
            enddo
         enddo                  ! nx

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) '      Brine height       '
        write(nu_diag_out+n-1,900) 'hbrin                   = ',phinS
        write(nu_diag_out+n-1,900) 'hbrin cat 1             = ',phinS1
        write(nu_diag_out+n-1,900) 'Freeboard               = ',phbrn
        write(nu_diag_out+n-1,900) 'dhbrin cat 1 top        = ',pdh_top1
        write(nu_diag_out+n-1,900) 'dhbrin cat 1 bottom     = ',pdh_bot1
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) '     zSalinity         '
        write(nu_diag_out+n-1,900) 'Avg density (kg/m^3)   = ',psice_rho
        write(nu_diag_out+n-1,900) 'Salt flux (kg/m^2/s)   = ',pfzsal
        write(nu_diag_out+n-1,900) 'Grav. Drain. Salt flux = ',pfzsal_g
        write(nu_diag_out+n-1,900) 'Darcy V cat 1 (m/s)    = ',pdarcy_V1
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) ' Top down bgc Layer Model'
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'bTiz ice temp'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pbTiz(k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'iDi diffusivity'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (piDin(k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'bphi porosity'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pphin(k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'phi1 porosity'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pphin1(k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'zsal cat 1'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pSin1(k), k = 1,nblyr)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'zsal Avg S'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pSin(k), k = 1,nblyr)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'Sice Ice S'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pSice(k), k = 1,nilyr)
        write(nu_diag_out+n-1,*) '                         '

802   format (f24.17,2x)
803   format (a25,2x)
900   format (a25,2x,f24.17)

      end subroutine zsal_diags

!=======================================================================

      end module icedrv_diagnostics_bgc

!=======================================================================
