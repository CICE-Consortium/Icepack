!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      module icepack_drv_diagnostics_bgc

      use icepack_drv_kinds
      use icepack_drv_constants, only: c0, nu_diag, nu_diag_out
      use icepack_drv_calendar, only: diagfreq, istep1, istep
      use icepack_drv_domain_size, only: nx

      implicit none
      private
      public :: hbrine_diags, bgc_diags, zsal_diags
      save

!=======================================================================

      contains

!=======================================================================

!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Nicole Jeffery, LANL

      subroutine hbrine_diags (dt)
              
      use icepack_drv_arrays_column, only: darcy_V
      use icepack_drv_diagnostics, only: nx_names
      use icepack_drv_domain_size, only: ncat, nltrcr, nilyr
      use icepack_drv_state, only: aice, aicen, vicen, vice, trcr, trcrn
      use icepack_drv_tracers, only: nt_sice, nt_fbri
      use icepack_drv_parameters, only: ktherm

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, k, n

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(nx) :: &
         phinS, phinS1, pdarcy_V, pfbri

      real (kind=dbl_kind), dimension(nx,nilyr) :: &
         pSin, pSin1

      !-----------------------------------------------------------------
      ! Dynamic brine height
      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)

         do n = 1, nx
               i = n
               phinS1(n) = c0             
               phinS(n) = c0
               pfbri(n) = trcrn(i,nt_fbri,1) 
               pdarcy_V(n) = darcy_V(i,1)
               if (aice(i) > c0) &
                       phinS(n) = trcr(i,nt_fbri)*vice(i)/aice(i)
               if (aicen(i,1)> c0) &
                       phinS1(n) = trcrn(i,nt_fbri,1)*vicen(i,1)/&
                                                aicen(i,1)                    
               do k = 1,nilyr
                  pSin1(n,k) = trcrn(i,nt_sice+k-1,1)
                  pSin(n,k)  = trcr(i,nt_sice+k-1)
               enddo
        ! enddo                  ! nx

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

         write(nu_diag_out+n-1,899) nx_names(n)
        
         write(nu_diag_out+n-1,*) '------ hbrine ------'
         write(nu_diag_out+n-1,900) 'hbrine, (m)        = ',phinS(n)
         write(nu_diag_out+n-1,900) 'fbri, cat1 (m)     = ',pfbri(n)
         write(nu_diag_out+n-1,900) 'hbrine cat1, (m)   = ',phinS1(n)
         write(nu_diag_out+n-1,900) 'darcy_V cat1, (m/s)= ',pdarcy_V(n)
         if (ktherm == 2) then          
            write(nu_diag_out+n-1,*) '                         '
            write(nu_diag_out+n-1,*) '------ Thermosaline Salinity ------'
            write(nu_diag_out+n-1,803) 'Sice1 cat1 S (ppt)'
            write(nu_diag_out+n-1,*) '---------------------------------------------------'
            write(nu_diag_out+n-1,802) (pSin1(n,k), k = 1,nilyr)              
            write(nu_diag_out+n-1,*) '                         '
            write(nu_diag_out+n-1,803) 'Sice bulk S (ppt) '
            write(nu_diag_out+n-1,*) '---------------------------------------------------'
            write(nu_diag_out+n-1,802) (pSin(n,k), k = 1,nilyr)              
            write(nu_diag_out+n-1,*) '                         '
         endif
       enddo ! n

802 format (f24.17,2x)
803 format (a25,2x)
899    format (43x,a24)
900 format (a25,2x,f24.17)

      end subroutine hbrine_diags

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Nicole Jeffery, LANL

      subroutine bgc_diags (dt)

      use icepack_drv_arrays_column, only: ocean_bio, zfswin, fbio_atmice, fbio_snoice
      use icepack_drv_arrays_column, only: Zoo, grow_net, ice_bio_net, trcrn_sw
      use icepack_drv_parameters,    only: max_algae, max_aero, max_dic, max_doc, max_don
      use icepack_drv_parameters,    only: max_fe, dEdd_algae, skl_bgc, z_tracers
      use icepack_drv_constants,     only: c0, mps_to_cmpdy, c100, p5, c1, secday
      use icepack_drv_domain_size,   only: ncat, nltrcr, nblyr, n_algae, n_zaero
      use icepack_drv_domain_size,   only: n_dic, n_doc, n_don, n_fed, n_fep, nilyr, nslyr
      use icepack_drv_flux,  only: flux_bio, flux_bio_atm
      use icepack_drv_state, only: aice, vicen, vice, trcr
      use icepack_drv_tracers, only: tr_bgc_DMS, tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C
      use icepack_drv_tracers, only: tr_bgc_DON, tr_zaero, tr_bgc_hum
      use icepack_drv_tracers, only: tr_bgc_Fe

      use icepack_drv_tracers, only: tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil
      use icepack_drv_tracers, only: nt_fbri, nt_sice, nt_zaero
      use icepack_drv_tracers, only: nt_bgc_n, nt_bgc_doc, nt_bgc_don, nt_bgc_fed, nt_bgc_fep
      use icepack_drv_tracers, only: nt_bgc_nit, nt_bgc_am, nt_bgc_sil, nt_bgc_hum, nt_bgc_pon
      use icepack_drv_tracers, only: nt_bgc_dmspp, nt_bgc_dmspd, nt_bgc_dms
      use icepack_drv_tracers, only: nlt_bgc_N, nlt_zaero, nlt_chl_sw, nlt_zaero_sw
      use icepack_drv_tracers, only: nlt_bgc_Fed, nlt_bgc_Fep  
      use icepack_drv_tracers, only: nlt_bgc_hum
      use icepack_drv_tracers, only: nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil
      use icepack_drv_tracers, only: nlt_bgc_DOC, nlt_bgc_DON, nlt_bgc_DMSPp, nlt_bgc_DMS

      !use icepack_drv_tracers, only: nt_bgc_N, nt_bgc_C, nt_bgc_chl, nt_bgc_Am
      !use icepack_drv_tracers, only: nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil
      !use icepack_drv_tracers, only: nt_bgc_PON, nt_bgc_DON, nt_bgc_DIC, nt_bgc_DOC, nt_zaero, nt_bgc_Fed
      !use icepack_drv_tracers, only: nt_bgc_Fep
      !use icepack_drv_tracers, only: nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, nlt_bgc_C, nlt_bgc_chl
      !use icepack_drv_tracers, only: nlt_bgc_DIC, nlt_bgc_DOC, nlt_bgc_PON
      !use icepack_drv_tracers, only: nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, nlt_zaero
      !use icepack_drv_tracers, only: nt_bgc_hum,  nlt_bgc_hum, nlt_chl_sw
      !use icepack_drv_tracers, only: nlt_zaero_sw

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, k, n, nn, kk, klev
      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(nx) :: &
         pNit_sk, pAm_sk, pSil_sk, phum_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac, &
         pflux_NO, pflux_Am,  phum_ac, &
         pflux_snow_NO, pflux_snow_Am,  &
         pflux_atm_NO, pflux_atm_Am,  pgrow_net, &
         pflux_hum

      real (kind=dbl_kind), dimension(nx,max_algae) :: &
         pN_ac, pN_tot, pN_sk, pflux_N
      real (kind=dbl_kind), dimension(nx,max_doc) :: &
         pDOC_ac, pDOC_sk
      real (kind=dbl_kind), dimension(nx,max_don) :: &
         pDON_ac, pDON_sk
      real (kind=dbl_kind), dimension(nx,max_fe ) :: &
         pFed_ac,  pFed_sk, pFep_ac, pFep_sk 
      real (kind=dbl_kind), dimension(nx,max_aero) :: &
        pflux_zaero, pflux_snow_zaero, pflux_atm_zaero, &
        pflux_atm_zaero_s

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(nx,2) :: &
         pNOs, pAms, pPONs, phums
      real (kind=dbl_kind), dimension(nx,2,max_algae) :: &
         pNs
      real (kind=dbl_kind), dimension(nx,2,max_doc) :: &
         pDOCs
      real (kind=dbl_kind), dimension(nx,2,max_don) :: &
         pDONs
      real (kind=dbl_kind), dimension(nx,2,max_fe ) :: &
         pFeds, pFeps 
      real (kind=dbl_kind), dimension(nx,2,max_aero) :: &
         pzaeros
      real (kind=dbl_kind), dimension(nx,nblyr+1) :: &
         pNO, pAm, pPON, pzfswin, pZoo, phum
      real (kind=dbl_kind), dimension(nx,nblyr+1,max_algae) :: &
         pN
      real (kind=dbl_kind), dimension(nx,nblyr+1,max_aero) :: &
         pzaero
      real (kind=dbl_kind), dimension(nx,nblyr+1,max_doc) :: &
         pDOC
      real (kind=dbl_kind), dimension(nx,nblyr+1,max_don) :: &
         pDON
      real (kind=dbl_kind), dimension(nx,nblyr+1,max_fe ) :: &
         pFed, pFep 
      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace
      real (kind=dbl_kind), dimension (nx,nslyr+nilyr+2) :: & 
         pchlsw
      real (kind=dbl_kind), dimension(nx,nslyr+nilyr+2,max_aero) :: &
         pzaerosw

      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = zspace(1)*p5
      zspace(nblyr+1) = zspace(nblyr+1)*p5    

      klev = 1+nilyr+nslyr
      !-----------------------------------------------------------------
      ! biogeochemical state of the ice
      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)

         do n = 1, nx
               i = n
               pAm_ac(n)   = c0
               pSil_ac(n)  = c0
               phum_ac(n)  = c0
               pDMSP_ac(n) = c0
               pDMS_ac(n)  = c0
               pN_ac(n,:)  = c0
               pDOC_ac(n,:)= c0
               pDON_ac(n,:)= c0
               pFed_ac(n,:)= c0
               pFep_ac(n,:)= c0
               pNit_ac(n) = c0
               if (tr_bgc_N) then
                 do k = 1,n_algae
                    pN_ac(n,k)    = ocean_bio(i,nlt_bgc_N(k)) 
                 enddo  !n_algae
               endif    !tr_bgc_N
               if (tr_bgc_C) then
                 do k = 1,n_doc
                    pDOC_ac(n,k)    = ocean_bio(i,nlt_bgc_DOC(k)) 
                 enddo  !n_algae
               endif    !tr_bgc_N
               if (tr_bgc_DON) then
                 do k = 1,n_don
                    pDON_ac(n,k)    = ocean_bio(i,nlt_bgc_DON(k)) 
                 enddo 
               endif
               if (tr_bgc_Fe ) then
                 do k = 1,n_fed 
                    pFed_ac (n,k)   = ocean_bio(i,nlt_bgc_Fed (k)) 
                 enddo 
                 do k = 1,n_fep 
                    pFep_ac (n,k)   = ocean_bio(i,nlt_bgc_Fep (k)) 
                 enddo 
               endif
               if (tr_bgc_Nit) &
               pNit_ac(n)  = ocean_bio(i,nlt_bgc_Nit)  ! nit(i)
               if (tr_bgc_Am) &
               pAm_ac(n)   = ocean_bio(i,nlt_bgc_Am)   ! amm(i)
               if (tr_bgc_Sil) &
               pSil_ac(n)  = ocean_bio(i,nlt_bgc_Sil)  ! sil(i)
               if (tr_bgc_hum) &
               phum_ac(n)  = ocean_bio(i,nlt_bgc_hum)  ! hum(i)
               if (tr_bgc_DMS) then
               pDMSP_ac(n) = ocean_bio(i,nlt_bgc_DMSPp)! dmsp(i)
               pDMS_ac(n)  = ocean_bio(i,nlt_bgc_DMS)  ! dms(i)
               endif

               ! fluxes in mmol/m^2/d
               ! concentrations are bulk in mmol/m^3
               ! iron is in 10^-3 mmol/m^3  (nM)

               if (skl_bgc) then
                  pNit_sk(n)   = c0
                  pAm_sk(n)    = c0
                  pSil_sk(n)   = c0
                  phum_sk(n)   = c0
                  pDMSPp_sk(n) = c0
                  pDMSPd_sk(n) = c0
                  pDMS_sk(n)   = c0
                  pN_sk(n,:)   = c0
                  pflux_N(n,:) = c0
                  pDOC_sk(n,:) = c0
                  pDON_sk(n,:) = c0
                  pFed_sk(n,:) = c0
                  pFep_sk(n,:) = c0

                  do k = 1,n_algae
                    pN_sk(n,k)       = trcr    (i,nt_bgc_N(k))
                    pflux_N(n,k)     = flux_bio(i,nlt_bgc_N(k))*mps_to_cmpdy/c100
                  enddo
                  if (tr_bgc_C) then
                     do k = 1,n_doc
                        pDOC_sk(n,k) = trcr    (i,nt_bgc_DOC(k))
                     enddo
                  endif
                  if (tr_bgc_DON) then
                     do k = 1,n_don
                        pDON_sk(n,k) = trcr    (i,nt_bgc_DON(k))
                     enddo
                  endif
                  if (tr_bgc_Fe ) then
                     do k = 1,n_fed
                        pFed_sk (n,k)= trcr    (i,nt_bgc_Fed(k))
                     enddo
                     do k = 1,n_fep
                        pFep_sk (n,k)= trcr    (i,nt_bgc_Fep(k))
                     enddo
                  endif
                  if (tr_bgc_Nit) then
                     pNit_sk(n)  = trcr    (i, nt_bgc_Nit)
                     pflux_NO(n) = flux_bio(i,nlt_bgc_Nit)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_Am) then
                     pAm_sk(n)   = trcr    (i, nt_bgc_Am)
                     pflux_Am(n) = flux_bio(i,nlt_bgc_Am)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_Sil) then
                     pSil_sk(n)  = trcr    (i, nt_bgc_Sil)
                  endif
                  if (tr_bgc_hum) then
                     phum_sk(n)  = trcr    (i, nt_bgc_hum)
                     pflux_hum(n)= flux_bio(i,nlt_bgc_hum)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_DMS) then
                     pDMSPp_sk(n) = trcr   (i,nt_bgc_DMSPp)
                     pDMSPd_sk(n) = trcr   (i,nt_bgc_DMSPd)
                     pDMS_sk  (n) = trcr   (i,nt_bgc_DMS)
                  endif

               elseif (z_tracers) then   ! zbgc
                  pflux_NO(n) = c0
                  pN_tot(n,:) = c0
                  pflux_Am(n) = c0
                  pflux_hum(n) = c0
                  pflux_atm_Am(n) = c0
                  pflux_snow_Am(n) = c0
                  pflux_N(n,:) = c0
                  pflux_NO(n) = c0
                  pflux_atm_NO(n) = c0
                  pflux_snow_NO(n) = c0
                  pflux_zaero(n,:) = c0
                  pflux_atm_zaero_s(n,:) = c0
                  pflux_atm_zaero(n,:) = c0
                  pflux_snow_zaero(n,:) = c0
                  if (tr_bgc_Nit) then
                    pflux_NO(n)       = flux_bio(i,nlt_bgc_Nit)*mps_to_cmpdy/c100
                    pflux_atm_NO(n)   = fbio_atmice(i,nlt_bgc_Nit)*mps_to_cmpdy/c100
                    pflux_snow_NO(n)  = fbio_snoice(i,nlt_bgc_Nit)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_Am) then
                    pflux_Am(n)       = flux_bio(i,nlt_bgc_Am)*mps_to_cmpdy/c100
                    pflux_atm_Am(n)   = fbio_atmice(i,nlt_bgc_Am)*mps_to_cmpdy/c100
                    pflux_snow_Am(n)  = fbio_snoice(i,nlt_bgc_Am)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_hum) then
                    pflux_hum(n)       = flux_bio(i,nlt_bgc_hum)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_N)  then
                     do k = 1,n_algae
                     pflux_N(n,k)     = flux_bio(i,nlt_bgc_N(k))*mps_to_cmpdy/c100
                     enddo
                  endif
                  if (tr_zaero)  then
                     do k = 1,n_zaero
                     pflux_zaero(n,k)      = flux_bio(i,nlt_zaero(k))*mps_to_cmpdy/c100
                     pflux_atm_zaero_s(n,k)= flux_bio_atm(i,nlt_zaero(k))*mps_to_cmpdy/c100 !*aice
                     pflux_atm_zaero(n,k)  = fbio_atmice(i,nlt_zaero(k))*mps_to_cmpdy/c100
                     pflux_snow_zaero(n,k) = fbio_snoice(i,nlt_zaero(k))*mps_to_cmpdy/c100
                     enddo
                  endif
                  do k = 1, nblyr+1
                   pzfswin(n,k) = c0
                   pZoo(n,k) = c0
                   do nn = 1,ncat
                     pzfswin(n,k) = pzfswin(n,k) + zfswin(i,k,nn)*vicen(i,nn)
                     pZoo(n,k) = pZoo(n,k) + Zoo(i,k,nn)*vicen(i,nn)
                   enddo !nn
                   if (vice(i) > c0) then
                     pzfswin(n,k) = pzfswin(n,k)/vice(i)
                     pZoo(n,k)    = pZoo(n,k)/vice(i)
                   endif !vice
                   pAm(n,k) = c0
                   pN(n,k,:) = c0
                   pDOC(n,k,:) = c0
                   pDON(n,k,:) = c0
                   pFed(n,k,:) = c0
                   pFep(n,k,:) = c0
                   pzaero(n,k,:) = c0
                   pPON(n,k) = c0
                   phum(n,k) = c0
                   pNO(n,k) = c0
                   if (tr_bgc_Nit) pNO(n,k) =  trcr(i,nt_bgc_Nit+k-1)
                   if (tr_bgc_Am) pAm(n,k) =  trcr(i,nt_bgc_Am+k-1)
                   if (tr_bgc_N) then
                     do nn = 1, n_algae
                        pN(n,k,nn)   =  trcr(i,nt_bgc_N(nn)+k-1)
                     enddo
                   endif
                   if (tr_bgc_C) then
                     do nn = 1, n_doc
                        pDOC(n,k,nn)   =  trcr(i,nt_bgc_DOC(nn)+k-1)
                     enddo
                   endif
                   if (tr_bgc_DON) then
                     do nn = 1, n_don
                        pDON(n,k,nn)   =  trcr(i,nt_bgc_DON(nn)+k-1)
                     enddo
                   endif
                   if (tr_bgc_Fe)  then
                     do nn = 1, n_fed
                        pFed(n,k,nn)   =  trcr(i,nt_bgc_Fed(nn)+k-1)
                     enddo
                     do nn = 1, n_fep
                        pFep(n,k,nn)   =  trcr(i,nt_bgc_Fep(nn)+k-1)
                     enddo
                   endif
                   if (tr_zaero) then
                     do nn = 1, n_zaero
                        pzaero(n,k,nn)   =  trcr(i,nt_zaero(nn)+k-1)
                     enddo
                   endif
                   if (tr_bgc_PON) pPON(n,k) =  trcr(i,nt_bgc_PON+k-1)
                   if (tr_bgc_hum) phum(n,k) =  trcr(i,nt_bgc_hum+k-1)
                 enddo  !k
                 if (tr_bgc_N) then
                   do nn = 1,n_algae
                      pN_tot(n,nn) = ice_bio_net(i,nlt_bgc_N(nn))
                   enddo
                   pgrow_net(n) =  grow_net(i)
                 endif !tr_bgc_N
                 do k = 1,2  !snow concentration
                   pAms(n,k) = c0
                   pNs(n,k,:) = c0
                   pDOCs(n,k,:) = c0
                   pDONs(n,k,:) = c0
                   pFeds (n,k,:)= c0
                   pFeps (n,k,:)= c0
                   pzaeros(n,k,:) = c0
                   pPONs(n,k) = c0
                   phums(n,k) = c0
                   pNOs(n,k) = c0
                   if (tr_bgc_Nit) pNOs(n,k)  = trcr(i,nt_bgc_Nit+nblyr+k)
                   if (tr_bgc_Am) pAms(n,k) = trcr(i,nt_bgc_Am+nblyr+k)
                   if (tr_bgc_N) then
                     do nn = 1, n_algae
                       pNs(n,k,nn) =  trcr(i,nt_bgc_N(nn)+nblyr+k)
                     enddo
                   endif
                   if (tr_bgc_C) then
                     do nn = 1, n_doc
                       pDOCs(n,k,nn) =  trcr(i,nt_bgc_DOC(nn)+nblyr+k)
                     enddo
                   endif
                   if (tr_bgc_DON) then
                     do nn = 1, n_don
                       pDONs(n,k,nn) =  trcr(i,nt_bgc_DON(nn)+nblyr+k)
                     enddo
                   endif
                   if (tr_bgc_Fe ) then
                     do nn = 1, n_fed
                       pFeds(n,k,nn) =  trcr(i,nt_bgc_Fed(nn)+nblyr+k)
                     enddo
                     do nn = 1, n_fep
                       pFeps(n,k,nn) =  trcr(i,nt_bgc_Fep(nn)+nblyr+k)
                     enddo
                   endif
                   if (tr_zaero) then
                     do nn = 1, n_zaero
                       pzaeros(n,k,nn) =  trcr(i,nt_zaero(nn)+nblyr+k)
                     enddo
                   endif
                   if (tr_bgc_PON)pPONs(n,k) =trcr(i,nt_bgc_PON+nblyr+k)
                   if (tr_bgc_hum)phums(n,k) =trcr(i,nt_bgc_hum+nblyr+k)
                 enddo   !k
               endif
               pchlsw(n,:) = c0
               pzaerosw(n,:,:) = c0
               if (dEdd_algae) then
                 do k = 0, klev
                    if (tr_bgc_N) pchlsw(n,k+1) = trcrn_sw(i,nlt_chl_sw+k,1)
                    if (tr_zaero) then
                    do nn = 1,n_zaero
                       pzaerosw(n,k+1,nn) =  trcrn_sw(i,nlt_zaero_sw(nn) + k,1)
                    enddo
                    endif
                  enddo
               endif              ! dEdd_algae

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

       if (z_tracers) then
         write(nu_diag_out+n-1,803) 'zfswin PAR  '
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,802) (pzfswin(n,k), k = 1,nblyr+1)
         write(nu_diag_out+n-1,*) '      '
         write(nu_diag_out+n-1,803) 'Losses: Zoo(mmol/m^3)  '
         write(nu_diag_out+n-1,803) '        Brine Conc.       ',' Brine Conc'
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,802) (pZoo(n,k), k = 1,nblyr+1)
         write(nu_diag_out+n-1,*) '      '
       endif
       if (tr_bgc_Nit) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '    nitrate conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag_out+n-1,900) 'Ocean conc       = ',pNit_ac(n)
         write(nu_diag_out+n-1,900) 'ice-ocean flux   = ',pflux_NO(n)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',pNit_sk(n)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,900) 'atm-ice flux     = ',pflux_atm_NO(n)
           write(nu_diag_out+n-1,900) 'snow-ice flux    = ',pflux_snow_NO(n)
           write(nu_diag_out+n-1,*) '             snow + ice conc'
           write(nu_diag_out+n-1,803) '    nitrate'
           write(nu_diag_out+n-1,802) (pNOs(n,k), k = 1,2)
           write(nu_diag_out+n-1,802) (pNO(n,k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '    '
         endif
      endif
      if (tr_bgc_PON .and. z_tracers) then
           write(nu_diag_out+n-1,*) '---------------------------------------------------'
           write(nu_diag_out+n-1,*) '    PON snow + ice conc. (mmol/m^3)'
           write(nu_diag_out+n-1,803) '    PON'
           write(nu_diag_out+n-1,802) (pPONs(n,k), k = 1,2)
           write(nu_diag_out+n-1,802) (pPON(n,k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) ' '
      endif
      if (tr_bgc_hum) then
           write(nu_diag_out+n-1,*) '---------------------------------------------------'
           write(nu_diag_out+n-1,*) '    hum snow + ice conc. (mmolC/m^3)'
           write(nu_diag_out+n-1,900) 'Ocean conc       = ',phum_ac(n)
           write(nu_diag_out+n-1,900) 'ice-ocean flux   = ',pflux_hum(n)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',phum_sk(n)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,803) '    hum'
           write(nu_diag_out+n-1,802) (phums(n,k), k = 1,2)
           write(nu_diag_out+n-1,802) (phum(n,k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) ' '
         endif
      endif
      if (tr_bgc_Am) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '    ammonium conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag_out+n-1,900) 'Ocean conc       = ',pAm_ac(n)
         write(nu_diag_out+n-1,900) 'ice-ocean flux   = ',pflux_Am(n)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',pAm_sk(n)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,900) 'atm-ice flux     = ',pflux_atm_Am(n)
           write(nu_diag_out+n-1,900) 'snow-ice flux    = ',pflux_snow_Am(n)
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  ammonium'
           write(nu_diag_out+n-1,802) (pAms(n,k), k = 1,2)
           write(nu_diag_out+n-1,802) (pAm(n,k), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '       '
         endif
      endif
      if (tr_bgc_N) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,900) 'tot algal growth (1/d) = ',pgrow_net(n)
       do kk = 1,n_algae
         write(nu_diag_out+n-1,*) '  algal conc. (mmol N/m^3) or flux (mmol N/m^2/d)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900) 'Ocean conc           = ',pN_ac(n,kk)
         write(nu_diag_out+n-1,900) 'ice-ocean flux       = ',pflux_N(n,kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900) 'Bulk ice conc.   = ',pN_sk(n,kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,900) 'Tot ice (mmolN/m^2) = ',pN_tot(n,kk)
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  algal N'
           write(nu_diag_out+n-1,802) (pNs(n,k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pN(n,k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '         '
         endif
       enddo
      endif
      if (tr_bgc_C) then
       do kk = 1,1 !n_doc
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  DOC conc. (mmol C/m^3)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pDOC_ac(n,kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pDOC_sk(n,kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  DOC'
           write(nu_diag_out+n-1,802) (pDOCs(n,k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pDOC(n,k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
      endif
      if (tr_bgc_DON) then
       do kk = 1,n_don
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  DON conc. (mmol N/m^3)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pDON_ac(n,kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pDON_sk(n,kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  DON'
           write(nu_diag_out+n-1,802) (pDONs(n,k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pDON(n,k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
      endif
      if (tr_bgc_Fe ) then
       do kk = 1,n_fed
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) ' dFe  conc. (nM)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pFed_ac (n,kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pFed_sk (n,kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  Fed'
           write(nu_diag_out+n-1,802) (pFeds (n,k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pFed (n,k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
       do kk = 1,n_fep
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) ' pFe  conc. (nM)'
         write(nu_diag_out+n-1,1020) '  type:', kk
         write(nu_diag_out+n-1,900)  'Ocean conc       = ',pFep_ac (n,kk)
         if (skl_bgc) then
           write(nu_diag_out+n-1,900)'Bulk ice conc.   = ',pFep_sk (n,kk)
         elseif (z_tracers) then
           write(nu_diag_out+n-1,*) '             snow + ice conc.'
           write(nu_diag_out+n-1,803) '  Fep'
           write(nu_diag_out+n-1,802) (pFeps (n,k,kk), k = 1,2)
           write(nu_diag_out+n-1,802) (pFep (n,k,kk), k = 1,nblyr+1)
           write(nu_diag_out+n-1,*) '      '
         endif
       enddo
      endif
      if (tr_bgc_DMS) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '    DMS (mmol/m^3)      '
         write(nu_diag_out+n-1,900) 'Ocean DMSP    = ',pDMSP_ac(n)
         write(nu_diag_out+n-1,900) 'Ocean DMS     = ',pDMS_ac(n)
         if (skl_bgc) then
          write(nu_diag_out+n-1,900) 'Ice DMSPp    = ',pDMSPp_sk(n)
          write(nu_diag_out+n-1,900) 'Ice DMSPd    = ',pDMSPd_sk(n)
          write(nu_diag_out+n-1,900) 'Ice DMS      = ',pDMS_sk(n)
         endif
      endif
      if (tr_zaero .and. z_tracers) then
       do kk = 1,n_zaero
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  aerosol conc. (kg/m^3) or flux (kg/m^2/d)'
         write(nu_diag_out+n-1,1020) '  type: ',kk
         write(nu_diag_out+n-1,900) 'Atm source flux     = ',pflux_atm_zaero_s(n,kk)
         write(nu_diag_out+n-1,900) 'ice-ocean flux*aice = ',pflux_zaero(n,kk)
         write(nu_diag_out+n-1,900) 'atm-ice flux*aice   = ',pflux_atm_zaero(n,kk)
         write(nu_diag_out+n-1,900) 'snow-ice flux*aice  = ',pflux_snow_zaero(n,kk)
         write(nu_diag_out+n-1,*) '             snow + ice conc.'
         write(nu_diag_out+n-1,803) ' aerosol'
         write(nu_diag_out+n-1,802) (pzaeros(n,k,kk), k = 1,2)
         write(nu_diag_out+n-1,802) (pzaero(n,k,kk), k = 1,nblyr+1)
         write(nu_diag_out+n-1,*) '            '
       enddo
      endif
      if (dEdd_algae) then
          if (tr_zaero) then
          do kk = 1,n_zaero
          write(nu_diag_out+n-1,*) '---------------------------------------------------'
          write(nu_diag_out+n-1,*) '  Cat 1 aerosol conc. (kg/m^3) on delta-Eddington grid  '
          write(nu_diag_out+n-1,802) (pzaerosw(n,k,kk), k = 1,klev +1)
          enddo
          endif
         if (tr_bgc_N) then
         write(nu_diag_out+n-1,*) '---------------------------------------------------'
         write(nu_diag_out+n-1,*) '  Cat 1 chl (mg/m^3) on delta-Eddington grid  '
         write(nu_diag_out+n-1,802) (pchlsw(n,k), k = 1,klev +1)
         endif
      endif
      enddo                  ! nx

802 format (f24.17,2x)
803 format (a25,2x)
900 format (a25,2x,f24.17)
1020 format (a30,2x,i6)    ! integer

      end subroutine bgc_diags

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Nicole Jeffery, LANL

      subroutine zsal_diags (dt)

      use icepack_drv_arrays_column, only: fzsal, fzsal_g, sice_rho, bTiz
      use icepack_drv_arrays_column, only: iDi, bphi, dhbr_top, dhbr_bot, darcy_V
      use icepack_drv_constants, only: rhos, rhoi, rhow, c1
      use icepack_drv_domain_size, only: nblyr, ncat, nilyr
      use icepack_drv_state, only: aicen, aice, vice, trcr, trcrn, vicen, vsno
      use icepack_drv_parameters, only: rhosi
      use icepack_drv_tracers, only: tr_brine, nt_fbri, nt_bgc_S, nt_sice

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, k, n, nn

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(nx) :: &
         phinS, phinS1,&
         phbrn,pdh_top1,pdh_bot1, psice_rho, pfzsal, &
         pfzsal_g, pdarcy_V1

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(nx,nblyr+2) :: &
         pphin, pgrid, pphin1
      real (kind=dbl_kind), dimension(nx,nblyr) :: &
         pSin, pSice, pSin1

      real (kind=dbl_kind), dimension(nx,nblyr+1) :: &
         pbTiz, piDin

      !-----------------------------------------------------------------
      ! salinity and microstructure  of the ice
      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)

         do n = 1, nx
               i = n

               pfzsal(n) = fzsal(i)
               pfzsal_g(n) = fzsal_g(i)
               phinS(n) = c0
               phinS1(n) = c0
               phbrn(n) = c0
               psice_rho(n) = c0
               pdh_top1(n) = c0
               pdh_bot1(n) = c0
               pdarcy_V1(n) = c0
               do nn = 1,ncat
                  psice_rho(n) = psice_rho(n) + sice_rho(i,nn)*aicen(i,nn)
               enddo
               if (aice(i) > c0) &
                  psice_rho(n) = psice_rho(n)/aice(i)
               if (tr_brine .and. aice(i) > c0) &
                  phinS(n) = trcr(i,nt_fbri)*vice(i)/aice(i)

               if (aicen(i,1)> c0) then
                  if (tr_brine) phinS1(n) = trcrn(i,nt_fbri,1) &
                                          * vicen(i,1)/aicen(i,1)
                  pdh_top1(n) = dhbr_top(i,1)
                  pdh_bot1(n) = dhbr_bot(i,1)
                  pdarcy_V1(n) = darcy_V(i,1)
               endif
               if (tr_brine .AND. aice(i) > c0) &
                  phbrn(n) = (c1 - rhosi/rhow)*vice(i)/aice(i) &
                                 - rhos/rhow  *vsno(i)/aice(i)
               do k = 1, nblyr+1
                  pbTiz(n,k) = c0
                  piDin(n,k) = c0
                  do nn = 1,ncat
                     pbTiz(n,k) = pbTiz(n,k) + bTiz(i,k,nn)*vicen(i,nn)
                     piDin(n,k) = piDin(n,k) +  iDi(i,k,nn)*vicen(i,nn)
                  enddo
                  if (vice(i) > c0) then
                     pbTiz(n,k) = pbTiz(n,k)/vice(i)
                     piDin(n,k) = piDin(n,k)/vice(i)
                  endif
               enddo                 ! k
               do k = 1, nblyr+2
                  pphin(n,k) = c0
                  pphin1(n,k) = c0
                  if (aicen(i,1) > c0) pphin1(n,k) = bphi(i,k,1)
                  do nn = 1,ncat
                     pphin(n,k) = pphin(n,k) + bphi(i,k,nn)*vicen(i,nn)
                  enddo
                  if (vice(i) > c0) then
                     pphin(n,k) = pphin(n,k)/vice(i)
                  endif
               enddo
               do k = 1,nblyr
                  pSin(n,k) = c0
                  pSin1(n,k) = c0
                  pSin(n,k)= trcr(i,nt_bgc_S+k-1)
                  if (aicen(i,1) > c0) pSin1(n,k) = trcrn(i,nt_bgc_S+k-1,1)
               enddo
               do k = 1,nilyr
                  pSice(n,k) = trcr(i,nt_sice+k-1)
               enddo

         enddo                  ! nx

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) '      Brine height       '
        write(nu_diag_out+n-1,900) 'hbrin                   = ',phinS(n)
        write(nu_diag_out+n-1,900) 'hbrin cat 1             = ',phinS1(n)
        write(nu_diag_out+n-1,900) 'Freeboard               = ',phbrn(n)
        write(nu_diag_out+n-1,900) 'dhbrin cat 1 top        = ',pdh_top1(n)
        write(nu_diag_out+n-1,900) 'dhbrin cat 1 bottom     = ',pdh_bot1(n)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) '     zSalinity         '
        write(nu_diag_out+n-1,900) 'Avg density (kg/m^3)   = ',psice_rho(n)
        write(nu_diag_out+n-1,900) 'Salt flux (kg/m^2/s)   = ',pfzsal(n)
        write(nu_diag_out+n-1,900) 'Grav. Drain. Salt flux = ',pfzsal_g(n)
        write(nu_diag_out+n-1,900) 'Darcy V cat 1 (m/s)    = ',pdarcy_V1(n)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) ' Top down bgc Layer Model'
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'bTiz ice temp'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pbTiz(n,k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'iDi diffusivity'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (piDin(n,k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'bphi porosity'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pphin(n,k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'phi1 porosity'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pphin1(n,k), k = 1,nblyr+1)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'zsal cat 1'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pSin1(n,k), k = 1,nblyr)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'zsal Avg S'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pSin(n,k), k = 1,nblyr)
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,803) 'Sice Ice S'
        write(nu_diag_out+n-1,*) '---------------------------------------------------'
        write(nu_diag_out+n-1,802) (pSice(n,k), k = 1,nilyr)
        write(nu_diag_out+n-1,*) '                         '

802 format (f24.17,2x)
803 format (a25,2x)
900 format (a25,2x,f24.17)

      end subroutine zsal_diags

!=======================================================================

      end module icepack_drv_diagnostics_bgc

!=======================================================================
