!=======================================================================

! Read and write ice model restart files
!
! authors Elizabeth C. Hunke, LANL

      module icepack_drv_restart

      use icepack_kinds_mod
      use icepack_drv_constants, only: nu_diag, nu_restart, nu_dump
      use icepack_drv_restart_shared, only: &
          restart, restart_dir, restart_file, lenstr

      implicit none
      private :: write_restart_pond_topo, read_restart_pond_topo, &
          write_restart_age,       read_restart_age, &
          write_restart_FY,        read_restart_FY, &  
          write_restart_lvl,       read_restart_lvl, & 
          write_restart_pond_cesm, read_restart_pond_cesm, & 
          write_restart_pond_lvl,  read_restart_pond_lvl, &
          write_restart_aero,      read_restart_aero

      public :: dumpfile, restartfile, &
                read_restart_field, write_restart_field, final_restart, &
                write_restart_field_cn, read_restart_field_cn, &
                write_restart_hbrine, read_restart_hbrine
      save
!cn future stuff for writing the number of tracers in file
#if 0
      integer (kind=int_kind), parameter :: f_ntrcr = 7 ! number of tracers !cn
      integer (kind=int_kind), dimension(f_ntrcr) :: f_trcr !list of tracers in file !cn
#endif

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
      use icepack_drv_constants, only: nu_diag, nu_dump
      use icepack_drv_domain_size, only: nilyr, nslyr, ncat, nx
      use icepack_drv_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          sst, frzmlt, coszen
      use icepack_drv_state, only: aicen, vicen, vsnon, trcrn, uvel, vvel
      use icepack_intfc_tracers, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno,&
          tr_iage, tr_FY, tr_lvl, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_aero

      ! local variables

      integer (kind=int_kind) :: &
          i, k, n, &              ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      character (len=3) :: nchar

      ! construct path/file
      iyear = nyr + year_init - 1
      imonth = month
      iday = mday
      
      write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
          restart_dir(1:lenstr(restart_dir)), &
          restart_file(1:lenstr(restart_file)),'.', &
          iyear,'-',month,'-',mday,'-',sec
      
!cn future stuff for writing the number of tracers in file
#if 0
      f_trcr = 0
      ! ice age tracer    
      if (tr_iage)  f_trcr(1) = 1
      ! first-year area tracer
      if (tr_FY)   f_trcr(2) = 1 
      ! level ice tracer
      if (tr_lvl)   f_trcr(3) = 1 
      ! CESM melt ponds
      if (tr_pond_cesm)   f_trcr(4) = 1 
      ! level-ice melt ponds
      if (tr_pond_lvl)  f_trcr(5) = 1
      ! topographic melt ponds
      if (tr_pond_topo)   f_trcr(6) = 1
      ! ice aerosol
      if (tr_aero)   f_trcr(7) = 1
#endif

      open(nu_dump,file=filename,form='unformatted')
      !cn future stuff for writing the number of tracers in file
      !cn write(nu_dump) (f_trcr(i),i=1,f_ntrcr)
      write(nu_dump) istep1,time,time_forc
      write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      
      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to binary files.  All other
      ! tracers are written to their own dump/restart binary files.
      !-----------------------------------------------------------------

      call write_restart_field_cn(nu_dump,aicen(:,:),ncat)
      call write_restart_field_cn(nu_dump,vicen(:,:),ncat)
      call write_restart_field_cn(nu_dump,vsnon(:,:),ncat)
!this is surface temperature
      call write_restart_field_cn(nu_dump,trcrn(:,nt_Tsfc,:),ncat)
!this is ice salinity in each nilyr
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field_cn(nu_dump,trcrn(:,nt_sice+k-1,:),ncat)
      enddo
!this is ice enthalpy in each nilyr
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field_cn(nu_dump,trcrn(:,nt_qice+k-1,:),ncat)
      enddo
!this is snow enthalpy in each nslyr
      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call write_restart_field_cn(nu_dump,trcrn(:,nt_qsno+k-1,:),ncat)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
#ifdef CCSMCOUPLED
      call write_restart_field_cn(nu_dump,coszen,1)
#endif
      call write_restart_field_cn(nu_dump,scale_factor,1)
      call write_restart_field_cn(nu_dump,swvdr,1)
      call write_restart_field_cn(nu_dump,swvdf,1)
      call write_restart_field_cn(nu_dump,swidr,1)
      call write_restart_field_cn(nu_dump,swidf,1)

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------
      if (oceanmixed_ice) then
!         call write_restart_field(nu_dump,0,sst,'ruf8','sst',1,diag)
!         call write_restart_field(nu_dump,0,frzmlt,'ruf8','frzmlt',1,diag)
      endif

      ! tracers
      ! ice age tracer    
      if (tr_iage) then 
        call write_restart_age()
      endif
      ! first-year area tracer
      if (tr_FY) then
        call write_restart_FY()
      endif
      ! level ice tracer
      if (tr_lvl) then
        call write_restart_lvl()
      endif
      ! CESM melt ponds
      if (tr_pond_cesm) then
        call write_restart_pond_cesm()
      endif
      ! level-ice melt ponds
      if (tr_pond_lvl) then
        call write_restart_pond_lvl()
      endif
      ! topographic melt ponds
      if (tr_pond_topo) then
        call write_restart_pond_topo()
      endif

      if (tr_aero) then ! ice aerosol
        call write_restart_aero()
      endif

      end subroutine dumpfile

!=======================================================================

! Restarts from a dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile (ice_ic)

      use icepack_drv_calendar, only: istep0, istep1, time, time_forc, calendar, npt
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc_shared, only: oceanmixed_ice
      use icepack_drv_constants, only: c0, p5, nu_diag, nu_restart
      use icepack_drv_domain_size, only: nilyr, nslyr, ncat, &
          max_ntrcr, nx
      use icepack_drv_flux, only: swvdr, swvdf, swidr, swidf, &
          sst, frzmlt, coszen, scale_factor
      use icepack_drv_init, only: tmask
      use icepack_drv_state, only: trcr_depend, aice, vice, vsno, trcr, &
          aice0, aicen, vicen, vsnon, trcrn, aice_init, uvel, vvel, &
          trcr_base, nt_strata, n_trcr_strata
      use icepack_intfc_tracers, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno,&
          tr_iage, tr_FY, tr_lvl, tr_pond_cesm, tr_pond_lvl, &
          tr_pond_topo, tr_aero, tr_brine

      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, k, n              ! counting indices

      character(len=char_len_long) :: &
         filename, filename0

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
      open(nu_restart,file=filename,form='unformatted')
      !cn read(nu_restart) (f_trcr(i),i=1,f_ntrcr)
      read (nu_restart) istep0,time,time_forc
      write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc

      istep1 = istep0
 
      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      write(nu_diag,*) ' min/max area, vol ice, vol snow, Tsfc'

      call read_restart_field_cn(nu_restart,aicen,ncat)
      call read_restart_field_cn(nu_restart,vicen,ncat)
      call read_restart_field_cn(nu_restart,vsnon,ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_Tsfc,:),ncat)

      write(nu_diag,*) 'min/max sice for each layer'
      do k=1,nilyr
        write(nchar,'(i3.3)') k
        call read_restart_field_cn(nu_restart,trcrn(:,nt_sice+k-1,:),ncat)
      enddo
      
      write(nu_diag,*) 'min/max qice for each layer'
      do k=1,nilyr
        write(nchar,'(i3.3)') k
        call read_restart_field_cn(nu_restart,trcrn(:,nt_qice+k-1,:),ncat)
      enddo
      
      write(nu_diag,*) 'min/max qsno for each layer'
      do k=1,nslyr
        write(nchar,'(i3.3)') k
        call read_restart_field_cn(nu_restart,trcrn(:,nt_qsno+k-1,:),ncat)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      write(nu_diag,*) 'radiation fields'

#ifdef CCSMCOUPLED
      call read_restart_field_cn(nu_restart,0,coszen,'ruf8', &
           'coszen',1,diag)
#endif
      call read_restart_field_cn(nu_restart,scale_factor,1)
      call read_restart_field_cn(nu_restart,swvdr,1)
      call read_restart_field_cn(nu_restart,swvdf,1)
      call read_restart_field_cn(nu_restart,swidr,1)
      call read_restart_field_cn(nu_restart,swidf,1)

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------

      if (oceanmixed_ice) then

        !write(nu_diag,*) 'min/max sst, frzmlt'

        !call read_restart_field(nu_restart,0,sst,'ruf8', &
        !    'sst',1,diag)
        !call read_restart_field(nu_restart,0,frzmlt,'ruf8', &
        !    'frzmlt',1,diag)
      endif

      ! tracers
      ! ice age tracer    
      if (tr_iage) then 
        call read_restart_age()
      endif
      ! first-year area tracer
      if (tr_FY) then
        call read_restart_FY()
      endif
      ! level ice tracer
      if (tr_lvl) then
        call read_restart_lvl()
      endif
      ! CESM melt ponds
      if (tr_pond_cesm) then
        call read_restart_pond_cesm()
      endif
      ! level-ice melt ponds
      if (tr_pond_lvl) then
        call read_restart_pond_lvl()
      endif
      ! topographic melt ponds
      if (tr_pond_topo) then
        call read_restart_pond_topo()
      endif
      ! ice aerosol
      if (tr_aero) then
        call read_restart_aero() 
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

!cn this gets called again upon returning...
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

      subroutine read_restart_field_cn(nu,work,ndim3)

      use icepack_drv_domain_size, only: nx

      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number (not used for netcdf)
           ndim3             ! third dimension

      real (kind=dbl_kind), dimension(nx,ndim3), &
           intent(inout) :: &
           work              ! input array (real, 8-byte)

      ! local variables

      integer (kind=int_kind) :: &
        n,     &      ! number of dimensions for variable
        i

      real (kind=dbl_kind), dimension(nx) :: &
           work2              ! input array (real, 8-byte)

      do n=1,ndim3
        read(nu) (work2(i),i=1,nx)
        work(:,n) = work2(:)
      enddo
      
    end subroutine read_restart_field_cn
      
!=======================================================================

! Reads a single restart field
! author David A Bailey, NCAR

      subroutine read_restart_field(nu,nrec,work,atype,vname,ndim3, diag)

      use icepack_drv_domain_size, only: nx

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
        read(nu) (work2(i),i=1,nx)
        work(:,n) = work2(:)
      enddo
      
    end subroutine read_restart_field
      
!=======================================================================

! Writes a single restart field.
! author David A Bailey, NCAR

      subroutine write_restart_field_cn(nu,work,ndim3)

      use icepack_drv_domain_size, only: nx

      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number
           ndim3         

      real (kind=dbl_kind), dimension(nx,ndim3), &
           intent(in) :: &
           work              ! input array (real, 8-byte)

      ! local variables

      integer (kind=int_kind) :: &
        n,     &      ! dimension counter
        i
      
      real (kind=dbl_kind), dimension(nx) :: &
          work2              ! input array (real, 8-byte)
      
      do n=1,ndim3
        work2(:) = work(:,n)
        write(nu) (work2(i),i=1,nx)
      enddo
      
    end subroutine write_restart_field_cn
     
!=======================================================================

! Writes a single restart field.
! author David A Bailey, NCAR

      subroutine write_restart_field(nu,nrec,work,atype,vname,ndim3,diag)

      use icepack_drv_domain_size, only: nx

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
        write(nu) (work2(i),i=1,nx)
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

! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_topo()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_apnd, nt_hpnd, nt_ipnd
      use icepack_drv_domain_size, only: ncat

      call write_restart_field_cn(nu_dump,trcrn(:,nt_apnd,:),ncat)
      call write_restart_field_cn(nu_dump,trcrn(:,nt_hpnd,:),ncat)
      call write_restart_field_cn(nu_dump,trcrn(:,nt_ipnd,:),ncat)

      end subroutine write_restart_pond_topo

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_topo()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_apnd, nt_hpnd, nt_ipnd
      use icepack_drv_domain_size, only: ncat

      write(nu_diag,*) 'min/max topo ponds'

      call read_restart_field_cn(nu_restart,trcrn(:,nt_apnd,:),ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_hpnd,:),ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_ipnd,:),ncat)

      end subroutine read_restart_pond_topo

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_age()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_iage
      use icepack_drv_domain_size, only: ncat

      call write_restart_field_cn(nu_dump,trcrn(:,nt_iage,:),ncat)

      end subroutine write_restart_age

!=======================================================================

! Reads all values needed for an ice age restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_age()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_iage
      use icepack_drv_domain_size, only: ncat

      write(nu_diag,*) 'min/max age (s)'

      call read_restart_field_cn(nu_restart,trcrn(:,nt_iage,:),ncat)

      end subroutine read_restart_age

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_FY()

      use icepack_drv_flux, only: frz_onset
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_FY
      use icepack_drv_domain_size, only: ncat

      call write_restart_field_cn(nu_dump,trcrn(:,nt_FY,:),ncat)
      call write_restart_field_cn(nu_dump,frz_onset,1)

      end subroutine write_restart_FY

!=======================================================================

! Reads all values needed for an ice FY restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_FY()

      use icepack_drv_flux, only: frz_onset
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_FY
      use icepack_drv_domain_size, only: ncat

      write(nu_diag,*) 'min/max first-year ice area'

      call read_restart_field_cn(nu_restart,trcrn(:,nt_FY,:),ncat)

      write(nu_diag,*) 'min/max frz_onset'

      call read_restart_field_cn(nu_restart,frz_onset,1)

      end subroutine read_restart_FY

!=======================================================================

! Dumps all values needed for restarting
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_lvl()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_alvl, nt_vlvl
      use icepack_drv_domain_size, only: ncat

      call write_restart_field_cn(nu_dump,trcrn(:,nt_alvl,:),ncat)
      call write_restart_field_cn(nu_dump,trcrn(:,nt_vlvl,:),ncat)

      end subroutine write_restart_lvl

!=======================================================================

! Reads all values needed for an ice lvl restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_lvl()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_alvl, nt_vlvl
      use icepack_drv_domain_size, only: ncat

      write(nu_diag,*) 'min/max level ice area, volume'

      call read_restart_field_cn(nu_restart,trcrn(:,nt_alvl,:),ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_vlvl,:),ncat)

      end subroutine read_restart_lvl

!=======================================================================!

! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_cesm()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_apnd, nt_hpnd
      use icepack_drv_domain_size, only: ncat

      call write_restart_field_cn(nu_dump,trcrn(:,nt_apnd,:),ncat)
      call write_restart_field_cn(nu_dump,trcrn(:,nt_hpnd,:),ncat)

      end subroutine write_restart_pond_cesm

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_cesm()

      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_apnd, nt_hpnd
      use icepack_drv_domain_size, only: ncat

      write(nu_diag,*) 'min/max cesm ponds'

      call read_restart_field_cn(nu_restart,trcrn(:,nt_apnd,:),ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_hpnd,:),ncat)

      end subroutine read_restart_pond_cesm

!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_pond_lvl()

      use icepack_drv_arrays_column, only: dhsn, ffracn
      use icepack_drv_flux, only: fsnow
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_apnd, nt_hpnd, nt_ipnd
      use icepack_drv_domain_size, only: ncat

      call write_restart_field_cn(nu_dump,trcrn(:,nt_apnd,:),ncat)
      call write_restart_field_cn(nu_dump,trcrn(:,nt_hpnd,:),ncat)
      call write_restart_field_cn(nu_dump,trcrn(:,nt_ipnd,:),ncat)
      call write_restart_field_cn(nu_dump,fsnow(:),1)
      call write_restart_field_cn(nu_dump,dhsn(:,:),ncat)
      call write_restart_field_cn(nu_dump,ffracn(:,:),ncat)

      end subroutine write_restart_pond_lvl

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL

      subroutine read_restart_pond_lvl()

      use icepack_drv_arrays_column, only: dhsn, ffracn
      use icepack_drv_flux, only: fsnow
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_apnd, nt_hpnd, nt_ipnd
      use icepack_drv_domain_size, only: ncat

      write(nu_diag,*) 'min/max level-ice ponds'

      call read_restart_field_cn(nu_restart,trcrn(:,nt_apnd,:),ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_hpnd,:),ncat)
      call read_restart_field_cn(nu_restart,trcrn(:,nt_ipnd,:),ncat)
      call read_restart_field_cn(nu_restart,fsnow(:),1)
      call read_restart_field_cn(nu_restart,dhsn(:,:),ncat)
      call read_restart_field_cn(nu_restart,ffracn(:,:),ncat)

      end subroutine read_restart_pond_lvl

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine write_restart_aero()

      use icepack_drv_domain_size, only: n_aero
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_aero
      use icepack_drv_domain_size, only: ncat

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices

      character (len=3)       :: nchar
    
      write(nu_diag,*) 'write_restart_aero (aerosols)'

      do k = 1, n_aero
       write(nchar,'(i3.3)') k
       call write_restart_field_cn(nu_dump, &
            trcrn(:,nt_aero  +(k-1)*4,:), &
            ncat)
       call write_restart_field_cn(nu_dump, &
            trcrn(:,nt_aero+1+(k-1)*4,:), &
            ncat)
       call write_restart_field_cn(nu_dump, &
            trcrn(:,nt_aero+2+(k-1)*4,:), &
            ncat)
       call write_restart_field_cn(nu_dump, &
            trcrn(:,nt_aero+3+(k-1)*4,:), &
            ncat)
      enddo

      end subroutine write_restart_aero

!=======================================================================

! Reads all values needed for an ice aerosol restart
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine read_restart_aero()

      use icepack_drv_domain_size, only: n_aero
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_aero
      use icepack_drv_domain_size, only: ncat

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices

      character (len=3)       :: nchar

      !-----------------------------------------------------------------

      write(nu_diag,*) 'read_restart_aero (aerosols)'

      do k = 1, n_aero
       write(nchar,'(i3.3)') k
       call read_restart_field_cn(nu_restart, trcrn(:,nt_aero  +(k-1)*4,:), ncat)
       call read_restart_field_cn(nu_restart, trcrn(:,nt_aero+1+(k-1)*4,:), ncat)
       call read_restart_field_cn(nu_restart, trcrn(:,nt_aero+2+(k-1)*4,:), ncat)
       call read_restart_field_cn(nu_restart, trcrn(:,nt_aero+3+(k-1)*4,:), ncat)
      enddo

      end subroutine read_restart_aero

!=======================================================================

      subroutine write_restart_hbrine()

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use icepack_drv_arrays_column, only: first_ice, first_ice_real
!      use icepack_drv_fileunits, only: nu_diag, nu_dump_hbrine
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_fbri
      !cn use icepack_drv_restart, only: write_restart_field
      use icepack_drv_constants, only: c1, c0
      use icepack_drv_domain_size, only: ncat, nx

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

      subroutine read_restart_hbrine()

! Reads all values needed for hbrine
! author Elizabeth C. Hunke, LANL

      use icepack_drv_arrays_column, only: first_ice_real, first_ice
!      use icepack_drv_fileunits, only: nu_diag, nu_restart_hbrine
      use icepack_drv_state, only: trcrn
      use icepack_intfc_tracers, only: nt_fbri
      use icepack_drv_constants, only: p5
      use icepack_drv_domain_size, only: ncat, nx

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

      end module icepack_drv_restart

!=======================================================================
