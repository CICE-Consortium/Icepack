!=======================================================================

! Read and write ice model restart files
!
! authors Elizabeth C. Hunke, LANL

      module icedrv_restart

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, nu_restart, nu_dump
      use icedrv_constants, only: c0, c1, p5
      use icedrv_restart_shared, only: restart, restart_dir, restart_file, lenstr
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private :: write_restart_pond_topo, read_restart_pond_topo, &
          write_restart_age,       read_restart_age, &
          write_restart_FY,        read_restart_FY, &  
          write_restart_lvl,       read_restart_lvl, & 
          write_restart_pond_cesm, read_restart_pond_cesm, & 
          write_restart_pond_lvl,  read_restart_pond_lvl, &
          write_restart_fsd,       read_restart_fsd, &
          write_restart_aero,      read_restart_aero

      public :: dumpfile, restartfile, final_restart, &
                write_restart_field, read_restart_field

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for a restart
! author Elizabeth C. Hunke, LANL

      subroutine dumpfile

      use icedrv_calendar, only: sec, month, mday, nyr, istep1
      use icedrv_calendar, only: time, time_forc, year_init
      use icedrv_domain_size, only: nilyr, nslyr, ncat
      use icedrv_forcing, only: oceanmixed_ice
      use icedrv_flux, only: scale_factor, swvdr, swvdf, swidr, swidf
      use icedrv_flux, only: sst, frzmlt
      use icedrv_state, only: aicen, vicen, vsnon, trcrn

      ! local variables

      integer (kind=int_kind) :: &
         k, &              ! counting indices
         iyear

      integer (kind=int_kind) :: &
         nt_Tsfc, nt_sice, nt_qice, nt_qsno

      logical (kind=log_kind) :: &
         tr_iage, tr_FY, tr_lvl, tr_aero, tr_brine, &
         tr_pond_topo, tr_pond_cesm, tr_pond_lvl, tr_fsd
!         solve_zsal, skl_bgc, z_tracers

      character(len=char_len_long) :: filename
      character(len=*), parameter :: subname='(dumpfile)'

      ! construct path/file
      iyear = nyr + year_init - 1
      
      write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
         restart_dir(1:lenstr(restart_dir)), &
         restart_file(1:lenstr(restart_file)),'.', &
         iyear,'-',month,'-',mday,'-',sec
      
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
           nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
           tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_brine_out=tr_brine, &
           tr_pond_topo_out=tr_pond_topo, tr_pond_cesm_out=tr_pond_cesm, &
           tr_pond_lvl_out=tr_pond_lvl,tr_fsd_out=tr_fsd)
!      call icepack_query_parameters(solve_zsal_out=solve_zsal, &
!         skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      open(nu_dump,file=filename,form='unformatted')
      write(nu_dump) istep1,time,time_forc
      write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      
      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to binary files.  All other
      ! tracers are written to their own dump/restart binary files.
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump,aicen(:,:),ncat)
      call write_restart_field(nu_dump,vicen(:,:),ncat)
      call write_restart_field(nu_dump,vsnon(:,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_Tsfc,:),ncat)
      do k=1,nilyr
         call write_restart_field(nu_dump,trcrn(:,nt_sice+k-1,:),ncat)
      enddo
      do k=1,nilyr
         call write_restart_field(nu_dump,trcrn(:,nt_qice+k-1,:),ncat)
      enddo
      do k=1,nslyr
         call write_restart_field(nu_dump,trcrn(:,nt_qsno+k-1,:),ncat)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,scale_factor,1)
      call write_restart_field(nu_dump,swvdr,1)
      call write_restart_field(nu_dump,swvdf,1)
      call write_restart_field(nu_dump,swidr,1)
      call write_restart_field(nu_dump,swidf,1)

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------
      if (oceanmixed_ice) then
         call write_restart_field(nu_dump,sst,1)
         call write_restart_field(nu_dump,frzmlt,1)
      endif

      ! tracers
      if (tr_iage)      call write_restart_age()       ! ice age tracer
      if (tr_FY)        call write_restart_FY()        ! first-year area tracer
      if (tr_lvl)       call write_restart_lvl()       ! level ice tracer
      if (tr_pond_cesm) call write_restart_pond_cesm() ! CESM melt ponds
      if (tr_pond_lvl)  call write_restart_pond_lvl()  ! level-ice melt ponds
      if (tr_pond_topo) call write_restart_pond_topo() ! topographic melt ponds
      if (tr_aero)      call write_restart_aero()      ! ice aerosol
      if (tr_brine)     call write_restart_hbrine()    ! brine height
      if (tr_fsd)       call write_restart_fsd()       ! floe size distribution
! called in icedrv_RunMod.F90 to prevent circular dependencies
!      if (solve_zsal .or. skl_bgc .or. z_tracers) &
!                        call write_restart_bgc         ! biogeochemistry

      end subroutine dumpfile

!=======================================================================

! Restarts from a dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile (ice_ic)

      use icedrv_calendar, only: istep0, istep1, time, time_forc
      use icepack_intfc, only: icepack_aggregate
      use icedrv_domain_size, only: nilyr, nslyr, ncat
      use icedrv_domain_size, only: max_ntrcr, nx
      use icedrv_flux, only: swvdr, swvdf, swidr, swidf
      use icedrv_flux, only: sst, frzmlt, scale_factor
      use icedrv_forcing, only: oceanmixed_ice
      use icedrv_init, only: tmask
      use icedrv_state, only: trcr_depend, aice, vice, vsno, trcr
      use icedrv_state, only: aice0, aicen, vicen, vsnon, trcrn, aice_init
      use icedrv_state, only: trcr_base, nt_strata, n_trcr_strata

      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, k              ! counting indices

      integer (kind=int_kind) :: &
         nt_Tsfc, nt_sice, nt_qice, nt_qsno

      logical (kind=log_kind) :: &
         tr_iage, tr_FY, tr_lvl, tr_aero, tr_brine, &
         tr_pond_topo, tr_pond_cesm, tr_pond_lvl, tr_fsd

      character(len=char_len_long) :: filename
      character(len=*), parameter :: subname='(restartfile)'

      if (present(ice_ic)) then 
         filename = trim(ice_ic)
      else
         call icedrv_system_abort(string=subname//'no ice_ic present', &
                                  file=__FILE__,line=__LINE__)
      endif

      write(nu_diag,*) 'Using restart dump=', trim(filename)
      open(nu_restart,file=filename,form='unformatted')
      read (nu_restart) istep0,time,time_forc
      write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
           nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
           tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_brine_out=tr_brine, &
           tr_pond_topo_out=tr_pond_topo, tr_pond_cesm_out=tr_pond_cesm, &
           tr_pond_lvl_out=tr_pond_lvl,tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      istep1 = istep0

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      write(nu_diag,*) 'min/max area, vol ice, vol snow, Tsfc'

      call read_restart_field(nu_restart,aicen,ncat)
      call read_restart_field(nu_restart,vicen,ncat)
      call read_restart_field(nu_restart,vsnon,ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_Tsfc,:),ncat)

      write(nu_diag,*) 'min/max sice for each layer'
      do k=1,nilyr
         call read_restart_field(nu_restart,trcrn(:,nt_sice+k-1,:),ncat)
      enddo
      
      write(nu_diag,*) 'min/max qice for each layer'
      do k=1,nilyr
         call read_restart_field(nu_restart,trcrn(:,nt_qice+k-1,:),ncat)
      enddo
      
      write(nu_diag,*) 'min/max qsno for each layer'
      do k=1,nslyr
         call read_restart_field(nu_restart,trcrn(:,nt_qsno+k-1,:),ncat)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      write(nu_diag,*) 'min/max radiation fields'

      call read_restart_field(nu_restart,scale_factor,1)
      call read_restart_field(nu_restart,swvdr,1)
      call read_restart_field(nu_restart,swvdf,1)
      call read_restart_field(nu_restart,swidr,1)
      call read_restart_field(nu_restart,swidf,1)

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------

      if (oceanmixed_ice) then
         write(nu_diag,*) 'min/max sst, frzmlt'
         call read_restart_field(nu_restart,sst,1)
         call read_restart_field(nu_restart,frzmlt,1)
      endif

      ! tracers
      if (tr_iage)      call read_restart_age()       ! ice age tracer
      if (tr_FY)        call read_restart_FY()        ! first-year area tracer
      if (tr_lvl)       call read_restart_lvl()       ! level ice tracer
      if (tr_pond_cesm) call read_restart_pond_cesm() ! CESM melt ponds
      if (tr_pond_lvl)  call read_restart_pond_lvl()  ! level-ice melt ponds
      if (tr_pond_topo) call read_restart_pond_topo() ! topographic melt ponds
      if (tr_aero)      call read_restart_aero()      ! ice aerosol
      if (tr_brine)     call read_restart_hbrine      ! brine height
      if (tr_fsd)       call read_restart_fsd()       ! floe size distribution
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
         call icepack_aggregate (ncat=ncat,          &
                                 aicen=aicen(i,:),   &
                                 trcrn=trcrn(i,:,:), &
                                 vicen=vicen(i,:),   &
                                 vsnon=vsnon(i,:),   &
                                 aice=aice (i),      &
                                 trcr=trcr (i,:),    &
                                 vice=vice (i),      &
                                 vsno=vsno (i),      &
                                 aice0=aice0(i),     &
                                 ntrcr=max_ntrcr,    &
                                 trcr_depend=trcr_depend, &
                                 trcr_base=trcr_base,     &
                                 n_trcr_strata=n_trcr_strata, &
                                 nt_strata=nt_strata)

         aice_init(i) = aice(i)
      enddo
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      end subroutine restartfile

!=======================================================================

! Reads a single restart field
! author Chris Newman, LANL

      subroutine read_restart_field(nu,work,ndim)

      use icedrv_domain_size, only: nx

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number (not used for netcdf)
         ndim             ! number of dimensions

      real (kind=dbl_kind), dimension(nx,ndim), intent(inout) :: &
         work              ! input array (real, 8-byte)

      ! local variables

      integer (kind=int_kind) :: &
         n, i               ! loop indices

      real (kind=dbl_kind), dimension(nx) :: &
         work2              ! input array (real, 8-byte)

      real (kind=dbl_kind) :: &
        minw, maxw          ! diagnostics

      character(len=*), parameter :: subname='(read_restart_field)'

      do n = 1, ndim
         read(nu) (work2(i), i=1,nx)
         work(:,n) = work2(:)
      enddo

      minw = minval(work)
      maxw = maxval(work)
      write(nu_diag,*) minw, maxw
      
      end subroutine read_restart_field
      
!=======================================================================

! Writes a single restart field.
! author Chris Newman, LANL

      subroutine write_restart_field(nu,work,ndim)

      use icedrv_domain_size, only: nx

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         ndim              ! number of dimensions

      real (kind=dbl_kind), dimension(nx,ndim), intent(in) :: &
         work              ! input array (real, 8-byte)

      ! local variables

      integer (kind=int_kind) :: &
         n, i              ! loop indices
      
      real (kind=dbl_kind), dimension(nx) :: &
         work2             ! input array (real, 8-byte)
      
      character(len=*), parameter :: subname='(write_restart_field)'

      do n = 1, ndim
        work2(:) = work(:,n)
        write(nu) (work2(i), i=1,nx)
      enddo
      
      end subroutine write_restart_field

!=======================================================================

! Finalize the restart file.
! author David A. Bailey, NCAR

      subroutine final_restart()

      use icedrv_calendar, only: istep1, time, time_forc

      character(len=*), parameter :: subname='(final_restart)'

      close(nu_dump)
      write(nu_diag,*) 'Restart read/written ',istep1,time,time_forc

      end subroutine final_restart

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_topo()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd 
      character(len=*), parameter :: subname='(write_restart_pond_topo)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_apnd,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_hpnd,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_ipnd,:),ncat)

      end subroutine write_restart_pond_topo

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_topo()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd 
      character(len=*), parameter :: subname='(read_restart_pond_topo)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'min/max topo ponds'

      call read_restart_field(nu_restart,trcrn(:,nt_apnd,:),ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_hpnd,:),ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_ipnd,:),ncat)

      end subroutine read_restart_pond_topo

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_age()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_iage
      character(len=*), parameter :: subname='(write_restart_age)'

      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_iage,:),ncat)

      end subroutine write_restart_age

!=======================================================================

! Reads all values needed for an ice age restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_age()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_iage
      character(len=*), parameter :: subname='(read_restart_age)'

      write(nu_diag,*) 'min/max age (s)'

      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call read_restart_field(nu_restart,trcrn(:,nt_iage,:),ncat)

      end subroutine read_restart_age

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_fsd()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat, nfsd
      integer (kind=int_kind) :: nt_fsd, k
      character(len=*), parameter :: subname='(write_restart_fsd)'

      call icepack_query_tracer_indices(nt_fsd_out=nt_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      do k = 1, nfsd
          call write_restart_field(nu_dump,trcrn(:,nt_fsd+k-1,:),ncat)
      end do

      end subroutine write_restart_fsd

!=======================================================================

! Reads all values needed for a FSD restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_fsd()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat, nfsd
      integer (kind=int_kind) :: nt_fsd, k
      character(len=*), parameter :: subname='(read_restart_fsd)'

      write(nu_diag,*) 'min/max fsd'

      call icepack_query_tracer_indices(nt_fsd_out=nt_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      do k =1, nfsd
          call read_restart_field(nu_restart,trcrn(:,nt_fsd+k-1,:),ncat)
      end do

      end subroutine read_restart_fsd

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_FY()

      use icedrv_flux, only: frz_onset
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_FY
      character(len=*), parameter :: subname='(write_restart_FY)'

      call icepack_query_tracer_indices(nt_FY_out=nt_FY)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_FY,:),ncat)
      call write_restart_field(nu_dump,frz_onset,1)

      end subroutine write_restart_FY

!=======================================================================

! Reads all values needed for an ice FY restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_FY()

      use icedrv_flux, only: frz_onset
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_FY
      character(len=*), parameter :: subname='(read_restart_FY)'

      call icepack_query_tracer_indices(nt_FY_out=nt_FY)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'min/max first-year ice area'

      call read_restart_field(nu_restart,trcrn(:,nt_FY,:),ncat)

      write(nu_diag,*) 'min/max frz_onset'

      call read_restart_field(nu_restart,frz_onset,1)

      end subroutine read_restart_FY

!=======================================================================

! Dumps all values needed for restarting
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_lvl()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_alvl, nt_vlvl
      character(len=*), parameter :: subname='(write_restart_lvl)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_alvl,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_vlvl,:),ncat)

      end subroutine write_restart_lvl

!=======================================================================

! Reads all values needed for an ice lvl restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_lvl()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_alvl, nt_vlvl
      character(len=*), parameter :: subname='(read_restart_lvl)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'min/max level ice area, volume'

      call read_restart_field(nu_restart,trcrn(:,nt_alvl,:),ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_vlvl,:),ncat)

      end subroutine read_restart_lvl

!=======================================================================!

! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_cesm()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_apnd, nt_hpnd
      character(len=*), parameter :: subname='(write_restart_pond_cesm)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_apnd,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_hpnd,:),ncat)

      end subroutine write_restart_pond_cesm

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_cesm()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_apnd, nt_hpnd
      character(len=*), parameter :: subname='(read_restart_pond_cesm)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'min/max cesm ponds'

      call read_restart_field(nu_restart,trcrn(:,nt_apnd,:),ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_hpnd,:),ncat)

      end subroutine read_restart_pond_cesm

!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_pond_lvl()

      use icedrv_arrays_column, only: dhsn, ffracn
      use icedrv_flux, only: fsnow
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*), parameter :: subname='(write_restart_pond_lvl)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_apnd,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_hpnd,:),ncat)
      call write_restart_field(nu_dump,trcrn(:,nt_ipnd,:),ncat)
      call write_restart_field(nu_dump,fsnow(:),1)
      call write_restart_field(nu_dump,dhsn(:,:),ncat)
      call write_restart_field(nu_dump,ffracn(:,:),ncat)

      end subroutine write_restart_pond_lvl

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL

      subroutine read_restart_pond_lvl()

      use icedrv_arrays_column, only: dhsn, ffracn
      use icedrv_flux, only: fsnow
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*), parameter :: subname='(write_restart_pond_lvl)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'min/max level-ice ponds'

      call read_restart_field(nu_restart,trcrn(:,nt_apnd,:),ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_hpnd,:),ncat)
      call read_restart_field(nu_restart,trcrn(:,nt_ipnd,:),ncat)
      call read_restart_field(nu_restart,fsnow(:),1)
      call read_restart_field(nu_restart,dhsn(:,:),ncat)
      call read_restart_field(nu_restart,ffracn(:,:),ncat)

      end subroutine read_restart_pond_lvl

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth Hunke, LANL
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine write_restart_aero()

      use icedrv_domain_size, only: n_aero
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices
      integer (kind=int_kind) :: nt_aero
      character(len=*), parameter :: subname='(write_restart_aero)'
    
      call icepack_query_tracer_indices(nt_aero_out=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'write_restart_aero (aerosols)'

      do k = 1, n_aero
         call write_restart_field(nu_dump, trcrn(:,nt_aero  +(k-1)*4,:), ncat)
         call write_restart_field(nu_dump, trcrn(:,nt_aero+1+(k-1)*4,:), ncat)
         call write_restart_field(nu_dump, trcrn(:,nt_aero+2+(k-1)*4,:), ncat)
         call write_restart_field(nu_dump, trcrn(:,nt_aero+3+(k-1)*4,:), ncat)
      enddo

      end subroutine write_restart_aero

!=======================================================================

! Reads all values needed for an ice aerosol restart
!
! authors Elizabeth Hunke, LANL
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine read_restart_aero()

      use icedrv_domain_size, only: n_aero
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices
      integer (kind=int_kind) :: nt_aero
      character(len=*), parameter :: subname='(read_restart_aero)'

      !-----------------------------------------------------------------

      call icepack_query_tracer_indices(nt_aero_out=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'read_restart_aero (aerosols)'

      do k = 1, n_aero
         call read_restart_field(nu_restart, trcrn(:,nt_aero  +(k-1)*4,:), ncat)
         call read_restart_field(nu_restart, trcrn(:,nt_aero+1+(k-1)*4,:), ncat)
         call read_restart_field(nu_restart, trcrn(:,nt_aero+2+(k-1)*4,:), ncat)
         call read_restart_field(nu_restart, trcrn(:,nt_aero+3+(k-1)*4,:), ncat)
      enddo

      end subroutine read_restart_aero

!=======================================================================

      subroutine write_restart_hbrine()

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use icedrv_arrays_column, only: first_ice, first_ice_real
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat, nx

      ! local variables

      integer (kind=int_kind) :: &
         i, n ! horizontal indices
      integer (kind=int_kind) :: nt_fbri
      character(len=*), parameter :: subname='(write_restart_hbrine)'

      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

        do i = 1, nx  
           do n = 1, ncat
              if (first_ice     (i,n)) then
                  first_ice_real(i,n) = c1
              else
                  first_ice_real(i,n) = c0
              endif
           enddo ! n
        enddo    ! i

        call write_restart_field(nu_dump,trcrn(:,nt_fbri,:),ncat)
        call write_restart_field(nu_dump,first_ice_real(:,:),ncat)

      end subroutine write_restart_hbrine

!=======================================================================

      subroutine read_restart_hbrine()

! Reads all values needed for hbrine
! author Elizabeth C. Hunke, LANL

      use icedrv_arrays_column, only: first_ice_real, first_ice
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat, nx

      ! local variables

      integer (kind=int_kind) :: &
         i, n ! horizontal indices
      integer (kind=int_kind) :: nt_fbri
      character(len=*), parameter :: subname='(read_restart_hbrine)'

      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'read brine restart'

      call read_restart_field(nu_restart,trcrn(:,nt_fbri,:),ncat)
      call read_restart_field(nu_restart,first_ice_real(:,:),ncat)

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

      end module icedrv_restart

!=======================================================================
