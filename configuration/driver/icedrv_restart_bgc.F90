!=========================================================================
!
! Restart routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
      module icedrv_restart_bgc

      use icedrv_kinds
      use icedrv_constants
      use icedrv_domain_size, only: ncat, nblyr, nx
      use icepack_intfc, only: icepack_max_algae, icepack_max_doc, icepack_max_don
      use icepack_intfc, only: icepack_max_dic, icepack_max_aero, icepack_max_fe
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icedrv_system, only: icedrv_system_abort

      implicit none

      private
      public ::  write_restart_bgc, &
                 read_restart_bgc

!=======================================================================

      contains

!=======================================================================
!
! Dumps all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_bgc()

      use icedrv_arrays_column, only: Rayleigh_criteria, Rayleigh_real
      use icedrv_domain_size, only: n_algae, n_doc, n_dic
      use icedrv_domain_size, only: n_don, n_zaero, n_fed, n_fep
      use icedrv_flux, only: sss, nit, amm, sil, dmsp, dms, algalN
      use icedrv_flux, only: doc, don, dic, fed, fep, zaeros, hum
      use icedrv_state, only: trcrn
      use icedrv_restart, only:  write_restart_field

      ! local variables

      integer (kind=int_kind) :: &
       i, k             , & ! horizontal, vertical indices
       mm                   ! n_algae

      logical (kind=log_kind) :: skl_bgc, solve_zsal, z_tracers

      integer (kind=int_kind) :: nbtrcr
      logical (kind=log_kind) :: tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_hum
      logical (kind=log_kind) :: tr_bgc_DMS, tr_bgc_PON, tr_bgc_N, tr_bgc_C
      logical (kind=log_kind) :: tr_bgc_DON, tr_bgc_Fe,  tr_zaero , tr_bgc_chl
      integer (kind=int_kind) :: nt_bgc_S, nt_bgc_Nit, nt_bgc_AM, nt_bgc_Sil
      integer (kind=int_kind) :: nt_bgc_hum, nt_bgc_PON
      integer (kind=int_kind) :: nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_DMS
      integer (kind=int_kind) :: nt_zbgc_frac
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_N
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_chl
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_C
      integer (kind=int_kind), dimension(icepack_max_doc)   :: nt_bgc_DOC
      integer (kind=int_kind), dimension(icepack_max_don)   :: nt_bgc_DON
      integer (kind=int_kind), dimension(icepack_max_dic)   :: nt_bgc_DIC
      integer (kind=int_kind), dimension(icepack_max_aero)  :: nt_zaero
      integer (kind=int_kind), dimension(icepack_max_fe)    :: nt_bgc_Fed
      integer (kind=int_kind), dimension(icepack_max_fe)    :: nt_bgc_Fep

      character(len=*), parameter :: subname='(write_restart_bgc)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(skl_bgc_out=skl_bgc)
      call icepack_query_parameters(solve_zsal_out=solve_zsal)
      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(tr_bgc_Nit_out=tr_bgc_Nit, & 
         tr_bgc_Am_out=tr_bgc_Am, &
         tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_hum_out=tr_bgc_hum, &
         tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
         tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, &
         tr_bgc_DON_out=tr_bgc_DON, tr_bgc_Fe_out=tr_bgc_Fe,  &
         tr_zaero_out=tr_zaero, tr_bgc_chl_out=tr_bgc_chl)
      call icepack_query_tracer_indices( nt_bgc_S_out=nt_bgc_S, &
         nt_bgc_N_out=nt_bgc_N, nt_bgc_AM_out=nt_bgc_AM, &
         nt_bgc_chl_out=nt_bgc_chl, nt_bgc_C_out=nt_bgc_C, &
         nt_bgc_DOC_out=nt_bgc_DOC, &
         nt_bgc_DIC_out=nt_bgc_DIC, nt_bgc_Nit_out=nt_bgc_Nit, &
         nt_bgc_Sil_out=nt_bgc_Sil, nt_bgc_hum_out=nt_bgc_hum, &
         nt_bgc_DMSPp_out=nt_bgc_DMSPp, nt_bgc_DMSPd_out=nt_bgc_DMSPd, &
         nt_bgc_DMS_out=nt_bgc_DMS, nt_bgc_PON_out=nt_bgc_PON, &
         nt_bgc_DON_out=nt_bgc_DON, nt_bgc_Fed_out=nt_bgc_Fed, &
         nt_bgc_Fep_out=nt_bgc_Fep, nt_zbgc_frac_out=nt_zbgc_frac, &
         nt_zaero_out=nt_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Salinity and extras
      !-----------------------------------------------------------------
      if (solve_zsal) then

      do k = 1, nblyr
        call write_restart_field(nu_dump,trcrn(:,nt_bgc_S+k-1,:),ncat) 
      enddo
    
      call write_restart_field(nu_dump,sss,1) 

         do i = 1, nx
            if (Rayleigh_criteria(i)) then
                Rayleigh_real    (i) = c1
            elseif (.NOT. Rayleigh_criteria(i)) then
                Rayleigh_real    (i) = c0
            endif
         enddo

      call write_restart_field(nu_dump,Rayleigh_real,1) 

      endif ! solve_zsal

      !-----------------------------------------------------------------
      ! Skeletal layer BGC
      !-----------------------------------------------------------------

      if (skl_bgc) then
         do k = 1, n_algae
           call write_restart_field(nu_dump,trcrn(:,nt_bgc_N(k),:),ncat) 
           if (tr_bgc_chl) &
           call write_restart_field(nu_dump,trcrn(:,nt_bgc_chl(k),:),ncat) 
         enddo
         if (tr_bgc_C)  then
!           do k = 1, n_algae
!             call write_restart_field(nu_dump,trcrn(:,nt_bgc_C(k),:),ncat) 
!           enddo
            do k = 1, n_doc
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DOC(k),:),ncat)
            enddo
            do k = 1, n_dic
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DIC(k),:),ncat)
            enddo
         endif
         if (tr_bgc_Nit) &
             call write_restart_field(nu_dump,trcrn(:,nt_bgc_Nit,:),ncat)
         if (tr_bgc_Am) &
             call write_restart_field(nu_dump,trcrn(:,nt_bgc_Am,:),ncat)
         if (tr_bgc_Sil) &
             call write_restart_field(nu_dump,trcrn(:,nt_bgc_Sil,:),ncat)
         if (tr_bgc_hum) &
             call write_restart_field(nu_dump,trcrn(:,nt_bgc_hum,:),ncat)
         if (tr_bgc_DMS) then
            call write_restart_field(nu_dump,trcrn(:,nt_bgc_DMSPp,:),ncat)
            call write_restart_field(nu_dump,trcrn(:,nt_bgc_DMSPd,:),ncat)
            call write_restart_field(nu_dump,trcrn(:,nt_bgc_DMS,:),ncat)
         endif
         if (tr_bgc_PON) &
             call write_restart_field(nu_dump,trcrn(:,nt_bgc_PON,:),ncat)
         if (tr_bgc_DON)  then
            do k = 1, n_don
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DON(k),:),ncat)
            enddo
         endif
         if (tr_bgc_Fe )  then
            do k = 1, n_fed 
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Fed(k),:),ncat)
            enddo
            do k = 1, n_fep 
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Fep(k),:),ncat)
            enddo
         endif

      elseif (z_tracers) then

      !-----------------------------------------------------------------
      ! Z layer BGC
      !-----------------------------------------------------------------

         if (tr_bgc_Nit) then
            write(nu_diag,*) 'z bgc restart: min/max Nitrate'
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Nit+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_N) then
            do mm = 1, n_algae
               do k = 1, nblyr+3
                  call write_restart_field(nu_dump,trcrn(:,nt_bgc_N(mm)+k-1,:),ncat)
               enddo
               if (tr_bgc_chl) then
                  do k = 1, nblyr+3
                     call write_restart_field(nu_dump,trcrn(:,nt_bgc_chl(mm)+k-1,:),ncat)
                  enddo
               endif
            enddo   ! n_algae
         endif      ! tr_bgc_N
         if (tr_bgc_C) then
! algal C is not yet distinct from algal N
!           do mm = 1, n_algae
!           do k = 1, nblyr+3
!              call write_restart_field(nu_dump,trcrn(:,nt_bgc_C(mm)+k-1,:),ncat)
!           enddo
!           enddo
            do mm = 1, n_doc
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DOC(mm)+k-1,:),ncat)
            enddo
            enddo
            do mm = 1, n_dic
            do k = 1, nblyr+3
              call write_restart_field(nu_dump,trcrn(:,nt_bgc_DIC(mm)+k-1,:),ncat)
            enddo
            enddo
         endif  !tr_bgc_C
         if (tr_bgc_Am) then
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Am+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_Sil) then
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Sil+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_hum) then
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_hum+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_DMS) then
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DMSPp+k-1,:),ncat)
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DMSPd+k-1,:),ncat)
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DMS+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_PON) then
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_PON+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_DON) then
            do mm = 1, n_don
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_DON(mm)+k-1,:),ncat)
            enddo
            enddo
         endif
         if (tr_bgc_Fe ) then
            do mm = 1, n_fed
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Fed(mm)+k-1,:),ncat)
            enddo
            enddo
            do mm = 1, n_fep
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_bgc_Fep(mm)+k-1,:),ncat)
            enddo
            enddo
         endif
         if (tr_zaero) then
            do mm = 1, n_zaero
            do k = 1, nblyr+3
               call write_restart_field(nu_dump,trcrn(:,nt_zaero(mm)+k-1,:),ncat)
            enddo
            enddo
         endif
         if (nbtrcr > 0) then
            do mm = 1, nbtrcr
               call write_restart_field(nu_dump,trcrn(:,nt_zbgc_frac+mm-1,:),ncat)
            enddo
         endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

         if (tr_bgc_N) then
            do k = 1, n_algae
               call write_restart_field(nu_dump,algalN(:,k),1)
            enddo  ! k
         endif
         if (tr_bgc_C) then
            do k = 1, n_doc
               call write_restart_field(nu_dump,doc(:,k),1)
            enddo  ! k
            do k = 1, n_dic
               call write_restart_field(nu_dump,dic(:,k),1)
            enddo  ! k
         endif      
         if (tr_bgc_Nit) &
             call write_restart_field(nu_dump,nit,1)
         if (tr_bgc_Am) &
             call write_restart_field(nu_dump,amm,1)
         if (tr_bgc_Sil) &
             call write_restart_field(nu_dump,sil,1)
         if (tr_bgc_hum) &
             call write_restart_field(nu_dump,hum,1)
         if (tr_bgc_DMS) then
             call write_restart_field(nu_dump,dmsp,1)
             call write_restart_field(nu_dump,dms,1)
         endif
         if (tr_bgc_DON) then
            do k = 1, n_don
               call write_restart_field(nu_dump,don(:,k),1)
            enddo  ! k
         endif
         if (tr_bgc_Fe ) then
            do k = 1, n_fed
               call write_restart_field(nu_dump,fed(:,k),1)
            enddo  ! k
            do k = 1, n_fep
               call write_restart_field(nu_dump,fep(:,k),1)
            enddo  ! k
         endif
         if (tr_zaero) then
            do k = 1, n_zaero
               call write_restart_field(nu_dump,zaeros(:,k),1)
            enddo  ! k
         endif

      endif  ! skl_bgc or z_tracers

      end subroutine write_restart_bgc

!=======================================================================
!
! Reads all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_bgc()

      use icedrv_arrays_column, only: Rayleigh_real, Rayleigh_criteria
      use icedrv_domain_size, only: n_algae, n_doc, n_dic
      use icedrv_domain_size, only: n_don, n_zaero, n_fed, n_fep
      use icedrv_flux, only: sss, nit, amm, sil, dmsp, dms, algalN
      use icedrv_flux, only: doc, don, dic, fed, fep, zaeros, hum
      use icedrv_state, only: trcrn
      use icedrv_restart, only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, k, & ! indices
         mm      ! n_algae

      logical (kind=log_kind) :: skl_bgc, solve_zsal, z_tracers

      integer (kind=int_kind) :: nbtrcr
      logical (kind=log_kind) :: tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_hum
      logical (kind=log_kind) :: tr_bgc_DMS, tr_bgc_PON, tr_bgc_N, tr_bgc_C
      logical (kind=log_kind) :: tr_bgc_DON, tr_bgc_Fe,  tr_zaero , tr_bgc_chl
      integer (kind=int_kind) :: nt_bgc_S, nt_bgc_Nit, nt_bgc_AM, nt_bgc_Sil
      integer (kind=int_kind) :: nt_bgc_hum, nt_bgc_PON
      integer (kind=int_kind) :: nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_DMS
      integer (kind=int_kind) :: nt_zbgc_frac
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_N
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_chl
      integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_C
      integer (kind=int_kind), dimension(icepack_max_doc)   :: nt_bgc_DOC
      integer (kind=int_kind), dimension(icepack_max_don)   :: nt_bgc_DON
      integer (kind=int_kind), dimension(icepack_max_dic)   :: nt_bgc_DIC
      integer (kind=int_kind), dimension(icepack_max_aero)  :: nt_zaero
      integer (kind=int_kind), dimension(icepack_max_fe)    :: nt_bgc_Fed
      integer (kind=int_kind), dimension(icepack_max_fe)    :: nt_bgc_Fep

      character(len=*), parameter :: subname='(read_restart_bgc)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(skl_bgc_out=skl_bgc)
      call icepack_query_parameters(solve_zsal_out=solve_zsal)
      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr)


      call icepack_query_tracer_flags(tr_bgc_Nit_out=tr_bgc_Nit, & 
         tr_bgc_Am_out=tr_bgc_Am, &
         tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_hum_out=tr_bgc_hum, &
         tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
         tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, &
         tr_bgc_DON_out=tr_bgc_DON, tr_bgc_Fe_out=tr_bgc_Fe,  &
         tr_zaero_out=tr_zaero, tr_bgc_chl_out=tr_bgc_chl)
      call icepack_query_tracer_indices( nt_bgc_S_out=nt_bgc_S, &
         nt_bgc_N_out=nt_bgc_N, nt_bgc_AM_out=nt_bgc_AM, &
         nt_bgc_chl_out=nt_bgc_chl, nt_bgc_C_out=nt_bgc_C, &
         nt_bgc_DOC_out=nt_bgc_DOC, &
         nt_bgc_DIC_out=nt_bgc_DIC, nt_bgc_Nit_out=nt_bgc_Nit, &
         nt_bgc_Sil_out=nt_bgc_Sil, nt_bgc_hum_out=nt_bgc_hum, &
         nt_bgc_DMSPp_out=nt_bgc_DMSPp, nt_bgc_DMSPd_out=nt_bgc_DMSPd, &
         nt_bgc_DMS_out=nt_bgc_DMS, nt_bgc_PON_out=nt_bgc_PON, &
         nt_bgc_DON_out=nt_bgc_DON, nt_bgc_Fed_out=nt_bgc_Fed, &
         nt_bgc_Fep_out=nt_bgc_Fep, nt_zbgc_frac_out=nt_zbgc_frac, &
         nt_zaero_out=nt_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Salinity and extras
      !-----------------------------------------------------------------

      if (solve_zsal) then

         write(nu_diag,*)'zSalinity restart'
         do k = 1, nblyr
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_S+k-1,:),ncat)
         enddo
 
         write(nu_diag,*) 'min/max sea surface salinity'
         call read_restart_field(nu_restart,sss,1)
         write(nu_diag,*) 'min/max Rayleigh'
         call read_restart_field(nu_restart,Rayleigh_real,1)

         do i = 1, nx  
            if (Rayleigh_real     (i) .GE. c1) then
                Rayleigh_criteria (i) = .true.
            elseif (Rayleigh_real (i) < c1) then
                Rayleigh_criteria (i) = .false.
            endif
         enddo

      endif ! solve_zsal

      !-----------------------------------------------------------------
      ! Skeletal Layer BGC
      !-----------------------------------------------------------------

      if (skl_bgc) then
         write(nu_diag,*) 'skl bgc restart'

         write (nu_diag,*) 'min/max algal N, chl'
         do k = 1, n_algae
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_N(k),:),ncat)
            if (tr_bgc_chl) &
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_chl(k),:),ncat)
         enddo   ! k
         if (tr_bgc_C) then
! algal C is not yet distinct from algal N
!           write (nu_diag,*) 'min/max algal C, DOC, DIC'
!           do k = 1, n_algae
!              call read_restart_field(nu_restart,trcrn(:,nt_bgc_C(k),:),ncat)
!           enddo
            write (nu_diag,*) 'min/max DOC, DIC'
            do k = 1, n_doc
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_DOC(k),:),ncat)
            enddo
            do k = 1, n_dic
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_DIC(k),:),ncat)
            enddo
         endif
         write (nu_diag,*) 'min/max Nit, Am, Sil, hum'
         if (tr_bgc_Nit) &
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_Nit,:),ncat)
         if (tr_bgc_Am) &
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_Am,:),ncat)
         if (tr_bgc_Sil) &
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_Sil,:),ncat)
         if (tr_bgc_hum) &
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_hum,:),ncat)
         if (tr_bgc_DMS) then
            write (nu_diag,*) 'min/max pDMSP, dDMSP, DMS'
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_DMSPp,:),ncat)
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_DMSPd,:),ncat)
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_DMS,:),ncat)
         endif
         write (nu_diag,*) 'min/max PON'
         if (tr_bgc_PON) &
            call read_restart_field(nu_restart,trcrn(:,nt_bgc_PON,:),ncat)
         if (tr_bgc_DON) then
            write (nu_diag,*) 'min/max DON'
            do k = 1, n_don
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_DON(k),:),ncat)
            enddo
         endif
         if (tr_bgc_Fe) then
            write (nu_diag,*) 'min/max dFe, pFe'
            do k = 1, n_fed 
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_Fed (k),:),ncat)
            enddo
            do k = 1, n_fep 
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_Fep (k),:),ncat)
            enddo
         endif

      elseif (z_tracers) then

      !-----------------------------------------------------------------
      ! Z Layer BGC
      !-----------------------------------------------------------------

         if (tr_bgc_Nit) then
            write(nu_diag,*) 'z bgc restart: min/max Nitrate'
            do k=1, nblyr+3
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_Nit+k-1,:),ncat)
            enddo
         endif   ! Nit
         if (tr_bgc_N) then
            do mm = 1, n_algae
               write(nu_diag,*) ' min/max Algal N', n_algae
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_bgc_N(mm)+k-1,:),ncat)
               enddo
               if (tr_bgc_chl) then
                  write(nu_diag,*) ' min/max Algal chla'
                  do k=1, nblyr+3
                     call read_restart_field(nu_restart,trcrn(:,nt_bgc_chl(mm)+k-1,:),ncat)
                  enddo
               endif   ! tr_bgc_chl
            enddo      ! n_algae
         endif         ! tr_bgc_N
         if (tr_bgc_C) then
! algal C is not yet distinct from algal N
!           do mm = 1, n_algae
!              write(nu_diag,*) ' min/max algalC'
!              do k=1, nblyr+3
!                 call read_restart_field(nu_restart,trcrn(:,nt_bgc_C(mm)+k-1,:),ncat)
!              enddo
!           enddo  !mm
            do mm = 1, n_doc
               write(nu_diag,*) ' min/max DOC'
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_bgc_DOC(mm)+k-1,:),ncat)
               enddo
            enddo  !mm
            do mm = 1, n_dic
               write(nu_diag,*) ' min/max DIC'
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_bgc_DIC(mm)+k-1,:),ncat)
               enddo
            enddo  !mm
         endif  ! tr_bgc_C
         if (tr_bgc_Am) then
            write(nu_diag,*) ' min/max ammonium'
            do k=1, nblyr+3
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_Am+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_Sil) then
            write(nu_diag,*) ' min/max silicate'
            do k=1, nblyr+3
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_Sil+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_hum) then
            write(nu_diag,*) ' min/max humic material'
            do k=1, nblyr+3
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_hum+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_DMS) then
            write (nu_diag,*) 'min/max pDMSP, dDMSP, DMS'
            do k=1, nblyr+3
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_DMSPp+k-1,:),ncat)
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_DMSPd+k-1,:),ncat)
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_DMS+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_PON) then
            write(nu_diag,*) ' min/max PON'
            do k=1, nblyr+3
               call read_restart_field(nu_restart,trcrn(:,nt_bgc_PON+k-1,:),ncat)
            enddo
         endif
         if (tr_bgc_DON) then
            do mm = 1, n_don
               write(nu_diag,*) ' min/max DON'
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_bgc_DON(mm)+k-1,:),ncat)
               enddo
            enddo  !mm
         endif
         if (tr_bgc_Fe) then
            do mm = 1, n_fed
               write(nu_diag,*) ' min/max dFe '
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_bgc_Fed (mm)+k-1,:),ncat)
               enddo
            enddo  !mm
            do mm = 1, n_fep
               write(nu_diag,*) ' min/max pFe '
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_bgc_Fep (mm)+k-1,:),ncat)
               enddo
            enddo  !mm
         endif
         if (tr_zaero) then
            do mm = 1, n_zaero
               write(nu_diag,*) ' min/max z aerosols'
               do k=1, nblyr+3
                  call read_restart_field(nu_restart,trcrn(:,nt_zaero(mm)+k-1,:),ncat)
               enddo
            enddo  !mm
         endif
         if (nbtrcr > 0) then
            write (nu_diag,*) 'min/max zbgc_frac'
            do mm = 1, nbtrcr
               call read_restart_field(nu_restart,trcrn(:,nt_zbgc_frac+mm-1,:),ncat)
            enddo
         endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

         write(nu_diag,*) 'mixed layer ocean bgc restart'
         if (tr_bgc_N) then
            write (nu_diag,*) 'min/max algalN'
            do k = 1, n_algae
               call read_restart_field(nu_restart,algalN(:,k),1)
            enddo  ! k
         endif
         if (tr_bgc_C) then
            write (nu_diag,*) 'min/max DOC, DIC'
            do k = 1, n_doc
               call read_restart_field(nu_restart,doc(:,k),1)
            enddo  ! k
            do k = 1, n_dic
               call read_restart_field(nu_restart,dic(:,k),1)
            enddo  ! k
         endif  !tr_bgc_C

         write (nu_diag,*) 'min/max Nit, Am, Sil, hum, DMSP, DMS'
         if (tr_bgc_Nit) &
            call read_restart_field(nu_restart,nit,1)
         if (tr_bgc_Am) &
            call read_restart_field(nu_restart,amm ,1)
         if (tr_bgc_Sil) &
            call read_restart_field(nu_restart,sil,1)
         if (tr_bgc_hum) &
            call read_restart_field(nu_restart,hum,1)
         if (tr_bgc_DMS) then
            call read_restart_field(nu_restart,dmsp,1)
            call read_restart_field(nu_restart,dms,1)
         endif
         if (tr_bgc_DON) then
            write (nu_diag,*) 'min/max DON'
            do k = 1, n_don
               call read_restart_field(nu_restart,don(:,k),1)
            enddo  ! k
         endif
         if (tr_bgc_Fe ) then
            write (nu_diag,*) 'min/max dFe, pFe'
            do k = 1, n_fed
               call read_restart_field(nu_restart,fed(:,k),1)
            enddo  ! k
            do k = 1, n_fep
               call read_restart_field(nu_restart,fep(:,k),1)
            enddo  ! k
         endif
         if (tr_zaero) then
            write (nu_diag,*) 'min/max zaeros'
            do k = 1, n_zaero
               call read_restart_field(nu_restart,zaeros(:,k),1)
            enddo  ! k
         endif

      endif  ! skl_bgc or z_tracers
   
      end subroutine read_restart_bgc

!=======================================================================

      end module icedrv_restart_bgc

!=======================================================================
