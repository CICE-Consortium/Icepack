!=======================================================================

! Read and write ice model restart files
!
! authors Elizabeth C. Hunke, LANL

      module icepack_drv_restart

      use icepack_kinds_mod
      use icepack_drv_constants, only: nu_rst_pointer, nu_diag, nu_restart, nu_dump
!      use icepack_drv_restart
      use icepack_drv_restart_shared, only: &
          restart, restart_dir, restart_file, lenstr

      implicit none
      private
      public :: dumpfile, restartfile, &
                read_restart_field, write_restart_field, final_restart
      save

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for a restart
! author Elizabeth C. Hunke, LANL

      subroutine dumpfile

      use icepack_drv_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, year_init
      use icepack_intfc_shared, only: oceanmixed_ice
      use icepack_drv_constants, only: c0, c1, nu_diag, nu_rst_pointer, nu_dump
      use icepack_drv_domain_size, only: nilyr, nslyr, ncat, nx
      use icepack_drv_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          sst, frzmlt, coszen
      !use icepack_drv_read_write, only: ice_open, ice_write
      use icepack_drv_state, only: aicen, vicen, vsnon, trcrn, uvel, vvel
      use icepack_intfc_tracers, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno

      ! local variables

      integer (kind=int_kind) :: &
          i, k, n, &              ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind), dimension (nx) :: &
         work1

      character (len=3) :: nchar

      ! construct path/file
      iyear = nyr + year_init - 1
      imonth = month
      iday = mday
      
      write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
          restart_dir(1:lenstr(restart_dir)), &
          restart_file(1:lenstr(restart_file)),'.', &
          iyear,'-',month,'-',mday,'-',sec
      
      !cn need to enable this call, it writes the binary restart file
      !cn corresponds to open(nu_dump,file=filename,form='unformatted')
      !            call ice_open(nu_dump,filename,0)
      open(nu_dump,file=filename,form='unformatted')
      write(nu_dump) istep1,time,time_forc
      write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      
      diag = .true.
      
      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to binary files.  All other
      ! tracers are written to their own dump/restart binary files.
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump,0,aicen(:,:),'ruf8','aicen',ncat,diag)
      call write_restart_field(nu_dump,0,vicen(:,:),'ruf8','vicen',ncat,diag)
      call write_restart_field(nu_dump,0,vsnon(:,:),'ruf8','vsnon',ncat,diag)
!cn this is surface temperature
      call write_restart_field(nu_dump,0,trcrn(:,nt_Tsfc,:),'ruf8','Tsfcn',ncat,diag)
      !write(*,*) trcrn(:,nt_Tsfc,:) !cn
!cn this is ice salinity in each nilyr
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,nt_sice+k-1,:),'ruf8', &
                                 'sice'//trim(nchar),ncat,diag)
      enddo
!cn this is ice enthalpy in each nilyr
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,nt_qice+k-1,:),'ruf8', &
                                 'qice'//trim(nchar),ncat,diag)
      enddo
!cn this is snow enthalpy in each nslyr
      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,nt_qsno+k-1,:),'ruf8', &
                                 'qsno'//trim(nchar),ncat,diag)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
#ifdef CCSMCOUPLED
      call write_restart_field(nu_dump,0,coszen,'ruf8','coszen',1,diag)
#endif
      call write_restart_field(nu_dump,0,scale_factor,'ruf8','scale_factor',1,diag)
      call write_restart_field(nu_dump,0,swvdr,'ruf8','swvdr',1,diag)
      call write_restart_field(nu_dump,0,swvdf,'ruf8','swvdf',1,diag)
      call write_restart_field(nu_dump,0,swidr,'ruf8','swidr',1,diag)
      call write_restart_field(nu_dump,0,swidf,'ruf8','swidf',1,diag)

!cn probably write all the tracers right here????


      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------
      if (oceanmixed_ice) then
!         call write_restart_field(nu_dump,0,sst,'ruf8','sst',1,diag)
!         call write_restart_field(nu_dump,0,frzmlt,'ruf8','frzmlt',1,diag)
      endif

      end subroutine dumpfile

!=======================================================================

! Restarts from a dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile (ice_ic)

      use icepack_drv_calendar, only: istep0, istep1, time, time_forc, calendar, npt
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc_shared, only: oceanmixed_ice
      use icepack_drv_constants, only: c0, p5, nu_diag, nu_rst_pointer, nu_restart
      use icepack_drv_domain_size, only: nilyr, nslyr, ncat, &
          max_ntrcr, nx
      use icepack_drv_flux, only: swvdr, swvdf, swidr, swidf, &
          sst, frzmlt, coszen, scale_factor
      use icepack_drv_init, only: tmask
!      use icepack_drv_read_write, only: ice_open, ice_read, ice_read_global
      use icepack_drv_state, only: trcr_depend, aice, vice, vsno, trcr, &
          aice0, aicen, vicen, vsnon, trcrn, aice_init, uvel, vvel, &
          trcr_base, nt_strata, n_trcr_strata
      use icepack_intfc_tracers, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno

      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, k, n              ! counting indices

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

      real (kind=dbl_kind), dimension (nx) :: &
         work1

      character (len=3) :: nchar

      if (present(ice_ic)) then 
         filename = trim(ice_ic)
      else
!cn need to do something here ....
!cn probably make sure there is a default for ice_ic up stream and require it as an arg here
!cn it is probably iced
        stop 'no ice_ic present'
      endif

      write(nu_diag,*) 'Using restart dump=', trim(filename)
      !            call ice_open(nu_restart,trim(filename),0)
      open(nu_restart,file=filename,form='unformatted')
      read (nu_restart) istep0,time,time_forc
      write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc

      istep1 = istep0

      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      write(nu_diag,*) ' min/max area, vol ice, vol snow, Tsfc'

      call read_restart_field(nu_restart,0,aicen,'ruf8', &
              'aicen',ncat,diag)
      call read_restart_field(nu_restart,0,vicen,'ruf8', &
              'vicen',ncat,diag)
      call read_restart_field(nu_restart,0,vsnon,'ruf8', &
              'vsnon',ncat,diag)
      call read_restart_field(nu_restart,0,trcrn(:,nt_Tsfc,:),'ruf8', &
              'Tsfcn',ncat,diag)

      write(nu_diag,*) 'min/max sice for each layer'
      do k=1,nilyr
        write(nchar,'(i3.3)') k
        call read_restart_field(nu_restart,0,trcrn(:,nt_sice+k-1,:),'ruf8', &
            'sice'//trim(nchar),ncat,diag)
      enddo
      
      write(nu_diag,*) 'min/max qice for each layer'
      do k=1,nilyr
        write(nchar,'(i3.3)') k
        call read_restart_field(nu_restart,0,trcrn(:,nt_qice+k-1,:),'ruf8', &
            'qice'//trim(nchar),ncat,diag)
      enddo
      
      write(nu_diag,*) 'min/max qsno for each layer'
      do k=1,nslyr
        write(nchar,'(i3.3)') k
        call read_restart_field(nu_restart,0,trcrn(:,nt_qsno+k-1,:),'ruf8', &
            'qsno'//trim(nchar),ncat,diag)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      write(nu_diag,*) 'radiation fields'

#ifdef CCSMCOUPLED
      call read_restart_field(nu_restart,0,coszen,'ruf8', &
           'coszen',1,diag)
#endif
      call read_restart_field(nu_restart,0,scale_factor,'ruf8',&
          'scale_factor',1,diag)
      call read_restart_field(nu_restart,0,swvdr,'ruf8', &
           'swvdr',1,diag)
      call read_restart_field(nu_restart,0,swvdf,'ruf8', &
           'swvdf',1,diag)
      call read_restart_field(nu_restart,0,swidr,'ruf8', &
           'swidr',1,diag)
      call read_restart_field(nu_restart,0,swidf,'ruf8', &
           'swidf',1,diag)

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------

      if (oceanmixed_ice) then

        write(nu_diag,*) 'min/max sst, frzmlt'

        call read_restart_field(nu_restart,0,sst,'ruf8', &
            'sst',1,diag)
        call read_restart_field(nu_restart,0,frzmlt,'ruf8', &
            'frzmlt',1,diag)
      endif

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      do i = 1, nx
         if (tmask(i)) &
         call icepack_aggregate (ncat,               &
                                aicen(i,:),  &
                                trcrn(i,:,:),&
                                vicen(i,:),  &
                                vsnon(i,:),  &
                                aice (i),  &
                                trcr (i,:),  &
                                vice (i),  &
                                vsno (i),  &
                                aice0(i),  &
                                max_ntrcr,          &
                                trcr_depend,        &
                                trcr_base,          &
                                n_trcr_strata,      &
                                nt_strata)

         aice_init(i) = aice(i)
      enddo

      end subroutine restartfile

!=======================================================================

! Reads a single restart field
! author David A Bailey, NCAR

      subroutine read_restart_field(nu,nrec,work,atype,vname,ndim3, &
                                    diag)

      use icepack_drv_domain_size, only: nx
!      use icepack_drv_read_write, only: ice_read, ice_read_nc

      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number (not used for netcdf)
           ndim3         , & ! third dimension
           nrec              ! record number (0 for sequential access)

      integer (kind=int_kind) :: i

      real (kind=dbl_kind), dimension(nx,ndim3), &
           intent(inout) :: &
           work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
           atype             ! format for output array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: vname

      ! local variables

      integer (kind=int_kind) :: &
        n,     &      ! number of dimensions for variable
        varid, &      ! variable id
        status        ! status variable from netCDF routine

      real (kind=dbl_kind), dimension(nx) :: &
           work2              ! input array (real, 8-byte)

         write(nu_diag,*) 'vname ',trim(vname)
            do n=1,ndim3
!cn               read(nu) (((work_g4(i,j,k),i=1,nx_global),j=1,ny_global),&
!cn                                                         k=1,nblyr+2)
              read(nu) (work2(i),i=1,nx)
!              call ice_read(nu,nrec,work2,atype,diag)
               work(:,n) = work2(:)
            enddo

      end subroutine read_restart_field
      
!=======================================================================

! Writes a single restart field.
! author David A Bailey, NCAR

      subroutine write_restart_field(nu,nrec,work,atype,vname,ndim3,diag)

      use icepack_drv_domain_size, only: nx
!      use icepack_drv_read_write, only: ice_write, ice_write_nc

      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number
           ndim3         , & ! third dimension
           nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx,ndim3), &
           intent(in) :: &
           work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
           atype             ! format for output array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      character (len=*), intent(in)  :: vname

      ! local variables

      integer (kind=int_kind) :: &
        n,     &      ! dimension counter
        varid, &      ! variable id
        status        ! status variable from netCDF routine
      
      integer (kind=int_kind) :: i

      real (kind=dbl_kind), dimension(nx) :: &
           work2              ! input array (real, 8-byte)

         do n=1,ndim3
            work2(:) = work(:,n)
!cn need to enable this ice_write....
!cn probably with something like the following copied from cice with atype=ruf8
!cn if this is right, then we need to do a lot of clean up here...
!cn            write(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
            write(nu) (work2(i),i=1,nx)
!               call ice_write(nu,nrec,work2,atype,diag)
         enddo

      end subroutine write_restart_field

!=======================================================================

! Finalize the restart file.
! author David A Bailey, NCAR

      subroutine final_restart()

      use icepack_drv_calendar, only: istep1, time, time_forc

      integer (kind=int_kind) :: status

         close(nu_dump)

         write(nu_diag,*) 'Restart read/written ',istep1,time,time_forc

      end subroutine final_restart

!=======================================================================

      end module icepack_drv_restart

!=======================================================================
