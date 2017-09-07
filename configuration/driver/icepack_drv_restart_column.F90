!=========================================================================
!
! Restart routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
      module icepack_drv_restart_column

      use icepack_kinds_mod
      use icepack_drv_constants
      use icepack_drv_domain_size, only: ncat, nilyr, nslyr, nblyr, nx
      use icepack_drv_restart, only: read_restart_field, write_restart_field

      implicit none
      save

      private
      public ::  write_restart_bgc,       read_restart_bgc,  &
                 write_restart_hbrine,    read_restart_hbrine

!=======================================================================

      contains

!=======================================================================

      subroutine read_restart_hbrine()

! Reads all values needed for hbrine
! author Elizabeth C. Hunke, LANL

      use icepack_drv_arrays_column, only: first_ice_real, first_ice
!      use icepack_drv_fileunits, only: nu_diag, nu_restart_hbrine
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_fbri
      use icepack_drv_restart, only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, n ! horizontal indices

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      write(nu_diag,*) 'brine restart'

!      call read_restart_field(nu_restart_hbrine,0,trcrn(:,nt_fbri,:),'ruf8', &
!                              'fbrn',ncat,diag,field_loc_center,field_type_scalar)
!      call read_restart_field(nu_restart_hbrine,0,first_ice_real(:,:),'ruf8', &
!                              'first_ice',ncat,diag,field_loc_center,field_type_scalar)

         do i = 1, nx
            do n = 1, ncat
               if (first_ice_real(i,n) >= p5) then
                   first_ice     (i,n) = .true.
               else
                   first_ice     (i,n) = .false.
               endif
            enddo ! ncat
         enddo    ! i 

      end subroutine read_restart_hbrine

!=======================================================================

      subroutine write_restart_hbrine()

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use icepack_drv_arrays_column, only: first_ice, first_ice_real
!      use icepack_drv_fileunits, only: nu_diag, nu_dump_hbrine
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_fbri
      use icepack_drv_restart, only: write_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, n ! horizontal indices

      logical (kind=log_kind) :: diag

      diag = .true.

        do i = 1, nx  
           do n = 1, ncat
              if (first_ice     (i,n)) then
                  first_ice_real(i,n) = c1
              else
                  first_ice_real(i,n) = c0
              endif
           enddo ! n
        enddo    ! i

!      call write_restart_field(nu_dump_hbrine,0,trcrn(:,nt_fbri,:),'ruf8', &
!                               'fbrn',ncat,diag)
!      call write_restart_field(nu_dump_hbrine,0,first_ice_real(:,:),'ruf8', &
!                               'first_ice',ncat,diag)

      end subroutine write_restart_hbrine

!=======================================================================
!
! Dumps all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_bgc()

      use icepack_drv_arrays_column, only: Rayleigh_criteria, Rayleigh_real
      use icepack_drv_domain_size, only: ncat, n_algae, n_doc, n_dic, &
          n_don, n_zaero, n_fed, n_fep
!      use icepack_drv_fileunits, only: nu_diag, nu_dump_bgc
      use icepack_drv_flux, only: sss, nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use icepack_drv_state, only: trcrn
      use icepack_drv_restart, only:  write_restart_field
      use icepack_intfc_tracers, only: nt_bgc_S, nt_bgc_Am, &
          nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_C, nt_bgc_chl, &
          nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
          nt_bgc_PON, nt_bgc_DON, nt_bgc_DOC, nt_bgc_DIC, &
          nt_bgc_N, nt_zaero, nt_bgc_Fed, nt_bgc_Fep, &
          nt_zbgc_frac, nbtrcr,  &
          nt_bgc_Fep, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil,&
          tr_bgc_DMS, tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, &
          tr_bgc_DON, tr_bgc_Fe,  tr_zaero , tr_bgc_chl, &
          nt_bgc_hum, tr_bgc_hum
      use icepack_intfc_shared, only: skl_bgc, solve_zsal

      ! local variables

      integer (kind=int_kind) :: &
       i, k             , & ! horizontal, vertical indices
       mm                   ! n_algae

      real (kind=dbl_kind) :: cszn ! counter for history averaging

      logical (kind=log_kind) :: diag

      character (len=3) :: nchar, ncharb

      integer (kind=int_kind) :: &
         ipoint

      diag = .true.

      !-----------------------------------------------------------------
      ! Salinity and extras
      !-----------------------------------------------------------------
      if (solve_zsal) then 

      do k = 1,nblyr
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_S+k-1,:),'ruf8', &
!                   'zSalinity'//trim(nchar),ncat,diag)
      enddo
    
!      call write_restart_field(nu_dump_bgc,0,sss,'ruf8','sss',1,diag)

         do i = 1, nx
            if (Rayleigh_criteria(i)) then
                Rayleigh_real    (i) = c1
            elseif (.NOT. Rayleigh_criteria(i)) then
                Rayleigh_real    (i) = c0
            endif
         enddo

!      call write_restart_field(nu_dump_bgc,0,Rayleigh_real,'ruf8','Rayleigh',1,diag)

      endif ! solve_zsal

      !-----------------------------------------------------------------
      ! Skeletal layer BGC
      !-----------------------------------------------------------------

      if (skl_bgc) then
         do k = 1, n_algae
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_N(k),:), &
!                                  'ruf8','bgc_N'//trim(nchar),ncat,diag)
!            if (tr_bgc_chl) &
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_chl(k),:), &
!                                  'ruf8','bgc_chl'//trim(nchar),ncat,diag)
         enddo
        if (tr_bgc_C)  then
          ! do k = 1, n_algae
          !  write(nchar,'(i3.3)') k
          !  call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_C(k),:), &
          !                        'ruf8','bgc_C'//trim(nchar),ncat,diag)
          ! enddo
           do k = 1, n_doc
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DOC(k),:), &
!                                  'ruf8','bgc_DOC'//trim(nchar),ncat,diag)
           enddo
           do k = 1, n_dic
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DIC(k),:), &
!                                  'ruf8','bgc_DIC'//trim(nchar),ncat,diag)
           enddo
         endif
!         if (tr_bgc_Nit) &
!         call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Nit,:), &
!                                  'ruf8','bgc_Nit',ncat,diag)
!         if (tr_bgc_Am) &
!         call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Am,:), &
!                                  'ruf8','bgc_Am',ncat,diag)
!         if (tr_bgc_Sil) &
!         call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Sil,:), &
!                                  'ruf8','bgc_Sil',ncat,diag)
!         if (tr_bgc_hum) &
!         call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_hum,:), &
!                                  'ruf8','bgc_hum',ncat,diag)
!         if (tr_bgc_DMS) then
!           call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DMSPp,:), &
!                                  'ruf8','bgc_DMSPp',ncat,diag)
!           call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DMSPd,:), &
!                                  'ruf8','bgc_DMSPd',ncat,diag)
!           call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DMS,:), &
!                                  'ruf8','bgc_DMS',ncat,diag)
!         endif
!         if (tr_bgc_PON) &
!         call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_PON,:), &
!                                  'ruf8','bgc_PON',ncat,diag)
      
        if (tr_bgc_DON)  then
           do k = 1, n_don
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DON(k),:), &
!                                  'ruf8','bgc_DON'//trim(nchar),ncat,diag)
           enddo
         endif
        if (tr_bgc_Fe )  then
           do k = 1, n_fed 
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Fed (k),:), &
!                                  'ruf8','bgc_Fed'//trim(nchar),ncat,diag)
           enddo
           do k = 1, n_fep 
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Fep (k),:), &
!                                  'ruf8','bgc_Fep'//trim(nchar),ncat,diag)
           enddo
         endif

      else 

      !-----------------------------------------------------------------
      ! Z layer BGC
      !-----------------------------------------------------------------

         if (tr_bgc_Nit) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Nit+k-1,:),'ruf8', &
!                                 'bgc_Nit'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_N) then
         do mm = 1,n_algae
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0, &
!                                  trcrn(:,nt_bgc_N(mm)+k-1,:),'ruf8', &
!                                 'bgc_N'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         if (tr_bgc_chl) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_bgc_chl(mm)+k-1,:),'ruf8', &
!                                 'bgc_chl'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         endif
         enddo   !n_algae
         endif   ! tr_bgc_N
         if (tr_bgc_C) then
        ! do mm = 1,n_algae
        ! write(ncharb, '(i3.3)') mm
        ! do k = 1,nblyr+3
        !    write(nchar,'(i3.3)') k
        !    call write_restart_field(nu_dump_bgc,0,  &
        !                          trcrn(:,nt_bgc_C(mm)+k-1,:),'ruf8', &
        !                         'bgc_C'//trim(ncharb)//trim(nchar),ncat,diag)
        ! enddo
        ! enddo
         do mm = 1,n_doc
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_bgc_DOC(mm)+k-1,:),'ruf8', &
!                                 'bgc_DOC'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         do mm = 1,n_dic
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_bgc_DIC(mm)+k-1,:),'ruf8', &
!                                 'bgc_DIC'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif  !tr_bgc_C
         if (tr_bgc_Am) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Am+k-1,:),'ruf8', &
!                                 'bgc_Am'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_Sil) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_Sil+k-1,:),'ruf8', &
!                                 'bgc_Sil'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_hum) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_hum+k-1,:),'ruf8', &
!                                 'bgc_hum'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_DMS) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DMSPp+k-1,:),'ruf8', &
!                                 'bgc_DMSPp'//trim(nchar),ncat,diag)
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DMSPd+k-1,:),'ruf8', &
!                                 'bgc_DMSPd'//trim(nchar),ncat,diag)
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_DMS+k-1,:),'ruf8', &
!                                 'bgc_DMS'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_PON) then
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,trcrn(:,nt_bgc_PON+k-1,:),'ruf8', &
!                                 'bgc_PON'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_DON) then
         do mm = 1,n_don
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_bgc_DON(mm)+k-1,:),'ruf8', &
!                                 'bgc_DON'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         if (tr_bgc_Fe ) then
         do mm = 1,n_fed
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_bgc_Fed(mm)+k-1,:),'ruf8', &
!                                 'bgc_Fed'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         do mm = 1,n_fep
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_bgc_Fep(mm)+k-1,:),'ruf8', &
!                                 'bgc_Fep'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         if (tr_zaero) then
         do mm = 1,n_zaero
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
!            write(nchar,'(i3.3)') k
!            call write_restart_field(nu_dump_bgc,0,  &
!                                  trcrn(:,nt_zaero(mm)+k-1,:),'ruf8', &
!                                 'zaero'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         do mm = 1,nbtrcr
          write(nchar,'(i3.3)') mm
!          call write_restart_field(nu_dump_bgc,0,   &
!                                trcrn(:,nt_zbgc_frac+mm-1,:),'ruf8', &
!                                'zbgc_frac'//trim(nchar),ncat,diag)
         enddo
      endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      if (tr_bgc_N) then
      do k = 1,n_algae
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,algalN(:,k),'ruf8','algalN'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_C) then
      do k = 1,n_doc
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,doc(:,k),'ruf8','doc'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_dic
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,dic(:,k),'ruf8','dic'//trim(nchar),1,diag)
      enddo  !k
      endif      
!      if (tr_bgc_Nit) &
!      call write_restart_field(nu_dump_bgc,0,nit,   'ruf8','nit',   1,diag)
!      if (tr_bgc_Am) &
!      call write_restart_field(nu_dump_bgc,0,amm,   'ruf8','amm',   1,diag)
!      if (tr_bgc_Sil) &
!      call write_restart_field(nu_dump_bgc,0,sil,   'ruf8','sil',   1,diag)
!      if (tr_bgc_hum) &
!      call write_restart_field(nu_dump_bgc,0,hum,   'ruf8','hum',   1,diag)
!      if (tr_bgc_DMS) then
!        call write_restart_field(nu_dump_bgc,0,dmsp,  'ruf8','dmsp',  1,diag)
!        call write_restart_field(nu_dump_bgc,0,dms,   'ruf8','dms',   1,diag)
!      endif
      if (tr_bgc_DON) then
      do k = 1,n_don
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,don(:,k),'ruf8','don'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_Fe ) then
      do k = 1,n_fed
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,fed(:,k),'ruf8','fed'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_fep
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,fep(:,k),'ruf8','fep'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_zaero) then
      do k = 1,n_zaero
!          write(nchar,'(i3.3)') k
!          call write_restart_field(nu_dump_bgc,0,zaeros(:,k),'ruf8','zaeros'//trim(nchar),1,diag)
      enddo  !k
      endif

      end subroutine write_restart_bgc

!=======================================================================
!
! Reads all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_bgc()

      use icepack_drv_arrays_column, only: Rayleigh_real, Rayleigh_criteria
      use icepack_drv_domain_size, only: ncat, n_algae, n_doc, n_dic,&
          n_don, n_zaero, n_fed, n_fep
!      use icepack_drv_fileunits, only: nu_diag, nu_restart_bgc
      use icepack_drv_flux, only: sss, nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use icepack_drv_state, only: trcrn
      use icepack_drv_restart, only: read_restart_field
      use icepack_intfc_tracers, only: nt_bgc_S, nt_bgc_Am, &
          nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_C, nt_bgc_chl, &
          nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
          nt_bgc_PON, nt_bgc_DON, nt_bgc_DOC, nt_bgc_DIC, &
          nt_bgc_N, nt_zaero, nt_bgc_Fed, nt_bgc_Fep, &
          nt_zbgc_frac, nbtrcr,  &
          nt_bgc_Fep, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil,&
          tr_bgc_DMS, tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, &
          tr_bgc_DON, tr_bgc_Fe,  tr_zaero , tr_bgc_chl, &
          nt_bgc_hum, tr_bgc_hum
      use icepack_intfc_shared, only: skl_bgc

      ! local variables

      integer (kind=int_kind) :: &
         i, k, & ! indices
         mm      ! n_algae

      logical (kind=log_kind) :: diag

      character (len=3) :: nchar, ncharb

      diag = .true.

      !-----------------------------------------------------------------
      ! Salinity and extras
      !-----------------------------------------------------------------

!      if (restart_zsal) then 

      write(nu_diag,*)'zSalinity restart'
      do k = 1,nblyr
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_S+k-1,:),'ruf8', &
!              'zSalinity'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
 
      write(nu_diag,*) 'sea surface salinity'
!      call read_restart_field(nu_restart_bgc,0,sss,'ruf8','sss',1,diag)
!      call read_restart_field(nu_restart_bgc,0,Rayleigh_real,'ruf8','Rayleigh',1,diag)

         do i = 1, nx  
            if (Rayleigh_real     (i) .GE. c1) then
                Rayleigh_criteria (i) = .true.
            elseif (Rayleigh_real (i) < c1) then
                Rayleigh_criteria (i) = .false.
            endif
         enddo
!      endif ! restart_zsal

      !-----------------------------------------------------------------
      ! Skeletal Layer BGC
      !-----------------------------------------------------------------
!      if (restart_bgc) then

      if (skl_bgc) then
       write(nu_diag,*) 'skl bgc restart'

       do k = 1, n_algae
!          write(nchar,'(i3.3)') k
!          call read_restart_field(nu_restart_bgc,0, &
!                 trcrn(:,nt_bgc_N(k),:), &
!                 'ruf8','bgc_N'//trim(nchar),ncat,diag)
!          if (tr_bgc_chl) &
!          call read_restart_field(nu_restart_bgc,0,  &
!                 trcrn(:,nt_bgc_chl(k),:), &
!                 'ruf8','bgc_chl'//trim(nchar),ncat,diag)
       enddo   !k
       if (tr_bgc_C) then
         ! do k = 1, n_algae
         !     write(nchar,'(i3.3)') k
         !     call read_restart_field(nu_restart_bgc,0,  &
         !        trcrn(:,nt_bgc_C(k),:), &
         !        'ruf8','bgc_C'//trim(nchar),ncat,diag)
         ! enddo
          do k = 1, n_doc
!              write(nchar,'(i3.3)') k
!              call read_restart_field(nu_restart_bgc,0,  &
!                 trcrn(:,nt_bgc_DOC(k),:), &
!                 'ruf8','bgc_DOC'//trim(nchar),ncat,diag)
          enddo
          do k = 1, n_dic
!              write(nchar,'(i3.3)') k
!              call read_restart_field(nu_restart_bgc,0,  &
!                 trcrn(:,nt_bgc_DIC(k),:), &
!                 'ruf8','bgc_DIC'//trim(nchar),ncat,diag)
          enddo
       endif
!       if (tr_bgc_Nit) &
!       call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_Nit,:), &
!           'ruf8','bgc_Nit',ncat,diag)
!       if (tr_bgc_Am) &
!       call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_Am,:), &
!           'ruf8','bgc_Am',ncat,diag)
!       if (tr_bgc_Sil) &
!       call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_Sil,:), &
!           'ruf8','bgc_Sil',ncat,diag)
!       if (tr_bgc_hum) &
!       call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_hum,:), &
!           'ruf8','bgc_hum',ncat,diag)
!       if(tr_bgc_DMS) then
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_DMSPp,:), &
!           'ruf8','bgc_DMSPp',ncat,diag)
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_DMSPd,:), &
!           'ruf8','bgc_DMSPd',ncat,diag)
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_DMS,:), &
!           'ruf8','bgc_DMS',ncat,diag)
!       endif
!       if (tr_bgc_PON) &
!       call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_PON,:), &
!           'ruf8','bgc_PON',ncat,diag)
       if (tr_bgc_DON) then
          do k = 1, n_don
!              write(nchar,'(i3.3)') k
!              call read_restart_field(nu_restart_bgc,0,  &
!                 trcrn(:,nt_bgc_DON(k),:), &
!                 'ruf8','bgc_DON'//trim(nchar),ncat,diag)
          enddo
       endif
       if (tr_bgc_Fe) then
          do k = 1, n_fed 
!              write(nchar,'(i3.3)') k
!              call read_restart_field(nu_restart_bgc,0,  &
!                 trcrn(:,nt_bgc_Fed (k),:), &
!                 'ruf8','bgc_Fed'//trim(nchar),ncat,diag)
          enddo
          do k = 1, n_fep 
!              write(nchar,'(i3.3)') k
!              call read_restart_field(nu_restart_bgc,0,  &
!                 trcrn(:,nt_bgc_Fep (k),:), &
!                 'ruf8','bgc_Fep'//trim(nchar),ncat,diag)
          enddo
       endif

!      else

      !-----------------------------------------------------------------
      ! Z Layer BGC
      !-----------------------------------------------------------------

      if (tr_bgc_Nit) then
      write(nu_diag,*) 'z bgc restart: min/max Nitrate'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_Nit+k-1,:),'ruf8', &
!              'bgc_Nit'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif   !Nit
      if (tr_bgc_N) then
      do mm = 1,n_algae
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max Algal N'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0, &
!               trcrn(:,nt_bgc_N(mm)+k-1,:),'ruf8', &
!              'bgc_N'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      if (tr_bgc_chl) then
      write(nu_diag,*) ' min/max Algal chla'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_bgc_chl(mm)+k-1,:),'ruf8', &
!              'bgc_chl'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif   ! tr_bgc_chl
      enddo   !n_algae
      endif   ! tr_bgc_N
      if (tr_bgc_C) then
     ! do mm = 1,n_algae
     ! write(ncharb,'(i3.3)') mm
     ! write(nu_diag,*) ' min/max Algal C'
     ! do k=1,nblyr+3
     !    write(nchar,'(i3.3)') k
     !    call read_restart_field(nu_restart_bgc,0,  &
     !          trcrn(:,nt_bgc_C(mm)+k-1,:),'ruf8', &
     !         'bgc_C'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
     ! enddo
     ! enddo  !mm
      do mm = 1,n_doc
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max DOC'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_bgc_DOC(mm)+k-1,:),'ruf8', &
!              'bgc_DOC'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      do mm = 1,n_dic
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max DIC'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_bgc_DIC(mm)+k-1,:),'ruf8', &
!              'bgc_DIC'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif  ! tr_bgc_C
      if (tr_bgc_Am) then
      write(nu_diag,*) ' min/max ammonium'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_Am+k-1,:),'ruf8', &
!              'bgc_Am'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_Sil) then
      write(nu_diag,*) ' min/max silicate'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_Sil+k-1,:),'ruf8', &
!              'bgc_Sil'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_hum) then
      write(nu_diag,*) ' min/max humic material'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_hum+k-1,:),'ruf8', &
!              'bgc_hum'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_DMS) then
      do k=1,nblyr+3
!      write(nu_diag,*) ' min/max DMSPp'
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_DMSPp+k-1,:),'ruf8', &
!              'bgc_DMSPp'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
!      write(nu_diag,*) ' min/max DMSPd'
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_DMSPd+k-1,:),'ruf8', &
!              'bgc_DMSPd'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
!      write(nu_diag,*) ' min/max DMS'
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_DMS+k-1,:),'ruf8', &
!              'bgc_DMS'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_PON) then
      write(nu_diag,*) ' min/max PON'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,trcrn(:,nt_bgc_PON+k-1,:),'ruf8', &
!              'bgc_PON'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_DON) then
      do mm = 1,n_don
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max DON'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_bgc_DON(mm)+k-1,:),'ruf8', &
!              'bgc_DON'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      if (tr_bgc_Fe) then
      do mm = 1,n_fed
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max dFe '
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_bgc_Fed (mm)+k-1,:),'ruf8', &
!              'bgc_Fed'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      do mm = 1,n_fep
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max pFe '
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_bgc_Fep (mm)+k-1,:),'ruf8', &
!              'bgc_Fep'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      if (tr_zaero) then
      do mm = 1,n_zaero
      write(ncharb,'(i3.3)') mm
      write(nu_diag,*) ' min/max z aerosols'
      do k=1,nblyr+3
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,  &
!               trcrn(:,nt_zaero(mm)+k-1,:),'ruf8', &
!              'zaero'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      do mm = 1,nbtrcr
!          write(nchar,'(i3.3)') mm
!          call read_restart_field(nu_restart_bgc,0, &
!               trcrn(:,nt_zbgc_frac+mm-1,:),'ruf8', &
!               'zbgc_frac'//trim(nchar),ncat,diag)
      enddo
!      endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      write(nu_diag,*) 'mixed layer ocean bgc restart'
      if (tr_bgc_N) then
      do k = 1,n_algae
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,algalN(:,k),'ruf8','algalN'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_C) then
      do k = 1,n_doc
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,doc(:,k),'ruf8','doc'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_dic
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,dic(:,k),'ruf8','dic'//trim(nchar),1,diag)
      enddo  !k
      endif  !tr_bgc_C

!      if (tr_bgc_Nit) &
!      call read_restart_field(nu_restart_bgc,0,nit   ,'ruf8','nit'   ,&
!                              1,diag,field_loc_center,field_type_scalar)
!      if (tr_bgc_Am) &
!      call read_restart_field(nu_restart_bgc,0,amm   ,'ruf8','amm'   ,&
!                              1,diag,field_loc_center,field_type_scalar)
!      if (tr_bgc_Sil) &
!      call read_restart_field(nu_restart_bgc,0,sil   ,'ruf8','sil'   ,&
!                              1,diag,field_loc_center,field_type_scalar)
!      if (tr_bgc_hum) &
!      call read_restart_field(nu_restart_bgc,0,hum   ,'ruf8','hum'   ,&
!                              1,diag,field_loc_center,field_type_scalar)
!      if (tr_bgc_DMS) then
!        call read_restart_field(nu_restart_bgc,0,dmsp  ,'ruf8','dmsp'  ,1,diag)
!        call read_restart_field(nu_restart_bgc,0,dms   ,'ruf8','dms'   ,1,diag)
!      endif
      if (tr_bgc_DON) then
      do k = 1,n_don
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,don(:,k),'ruf8','don'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_Fe ) then
      do k = 1,n_fed
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,fed(:,k),'ruf8','fed'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_fep
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,fep(:,k),'ruf8','fep'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_zaero) then
      do k = 1,n_zaero
!         write(nchar,'(i3.3)') k
!         call read_restart_field(nu_restart_bgc,0,zaeros(:,k),'ruf8','zaeros'//trim(nchar),1,diag)
      enddo  !k
      endif
      endif  ! restart_bgc
   
      end subroutine read_restart_bgc

!=======================================================================

      end module icepack_drv_restart_column

!=======================================================================
