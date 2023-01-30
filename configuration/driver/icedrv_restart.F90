!=======================================================================

! Read and write ice model restart files
!
! authors Elizabeth C. Hunke, LANL

      module icedrv_restart

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, nu_restart, nu_dump
      use icedrv_constants, only: c0, c1, p5
      use icedrv_restart_shared, only: restart, restart_dir, restart_file, lenstr
      use icedrv_restart_shared, only: restart_format
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters
      use icedrv_system, only: icedrv_system_abort
      use icepack_intfc, only: icepack_warnings_add, warnstr
#ifdef USE_NETCDF
      use netcdf
#endif

      implicit none
      private :: &
          write_restart_age,       read_restart_age, &
          write_restart_FY,        read_restart_FY, &
          write_restart_lvl,       read_restart_lvl, &
          write_restart_pond_lvl,  read_restart_pond_lvl, &
          write_restart_pond_topo, read_restart_pond_topo, &
          write_restart_snow,      read_restart_snow, &
          write_restart_fsd,       read_restart_fsd, &
          write_restart_iso,       read_restart_iso, &
          write_restart_aero,      read_restart_aero, &
          define_rest_field
      
      integer (kind=int_kind), private :: &
         ncid     ! ID for NetCDF file
         
      public :: dumpfile, restartfile, final_restart, &
                write_restart_field, read_restart_field

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read restart data files ..
!     Which can either be in unformatted Fortran format
!     (restart_format = 'bin') or NetCDF (restart_format = 'nc')      
!=======================================================================

! Dumps all values needed for a restart
! authors Elizabeth C. Hunke, LANL
!         David Clemens-Sewall, NCAR

      subroutine dumpfile

      use icedrv_calendar, only: sec, month, mday, nyr, istep1
      use icedrv_calendar, only: time, time_forc, year_init
      use icedrv_domain_size, only: nilyr, nslyr, ncat, n_aero, nfsd, nx, n_iso
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
         tr_iage, tr_FY, tr_lvl, tr_iso, tr_aero, tr_brine, &
         tr_pond_topo, tr_pond_lvl, tr_snow, tr_fsd
!         solve_zsal, skl_bgc, z_tracers

      character(len=char_len_long) :: filename
      character(len=*), parameter :: subname='(dumpfile)'
      
      ! local variables if we're writing in NetCDF
      integer (kind=int_kind) :: &
      dimid_ni,   & ! netCDF identifiers
      dimid_ncat, & !
      dimid_aero, &
      dimid_floe, &
      iflag,      & ! netCDF creation flag
      status        ! status variable from netCDF routine

      integer (kind=int_kind), allocatable :: dims(:),dims1(:)

      character (len=3) :: nchar


      ! get year
      iyear = nyr + year_init - 1
      
      ! query which tracers are active and their indices
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
         nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_iso_out=tr_iso, &
         tr_brine_out=tr_brine, &
         tr_pond_topo_out=tr_pond_topo, &
         tr_pond_lvl_out=tr_pond_lvl,tr_snow_out=tr_snow,tr_fsd_out=tr_fsd)
!      call icepack_query_parameters(solve_zsal_out=solve_zsal, &
!         skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)

      if (restart_format == 'bin') then
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
            restart_dir(1:lenstr(restart_dir)), &
            restart_file(1:lenstr(restart_file)),'.', &
            iyear,'-',month,'-',mday,'-',sec

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
         if (tr_pond_lvl)  call write_restart_pond_lvl()  ! level-ice melt ponds
         if (tr_pond_topo) call write_restart_pond_topo() ! topographic melt ponds
         if (tr_snow)      call write_restart_snow()      ! snow metamorphosis tracers
         if (tr_iso)       call write_restart_iso()       ! ice isotopes
         if (tr_aero)      call write_restart_aero()      ! ice aerosols
         if (tr_brine)     call write_restart_hbrine()    ! brine height
         if (tr_fsd)       call write_restart_fsd()       ! floe size distribution
   ! called in icedrv_RunMod.F90 to prevent circular dependencies
   !      if (solve_zsal .or. skl_bgc .or. z_tracers) &
   !                        call write_restart_bgc         ! biogeochemistry
      else if (restart_format == 'nc') then
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
         restart_dir(1:lenstr(restart_dir)), &
         restart_file(1:lenstr(restart_file)),'.', &
         iyear,'-',month,'-',mday,'-',sec
         filename = trim(filename) // '.nc'

         ! Open NetCDF restart file in define mode
         iflag = nf90_clobber
         status = nf90_create(trim(filename), iflag, ncid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname, &
            file=__FILE__,line= __LINE__)
         ! Set time attributes for restart file
         status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
         status = nf90_put_att(ncid,nf90_global,'time',time)
         status = nf90_put_att(ncid,nf90_global,'time_forc',time_forc)
         status = nf90_put_att(ncid,nf90_global,'nyr',iyear)
         status = nf90_put_att(ncid,nf90_global,'month',month)
         status = nf90_put_att(ncid,nf90_global,'mday',mday)
         status = nf90_put_att(ncid,nf90_global,'sec',sec)
         ! Set dimensions for restart file
         status = nf90_def_dim(ncid,'ni',nx,dimid_ni)
         status = nf90_def_dim(ncid,'ncat',ncat,dimid_ncat)
         !status = nf90_def_dim(ncid,'naero',n_aero,dimid_aero)
         !status = nf90_def_dim(ncid,'nfloe',nfsd,dimid_floe) 

         !-----------------------------------------------------------------
         ! state variables field creation
         !-----------------------------------------------------------------
         
         allocate(dims(2))
         allocate(dims1(1))
         dims(1) = dimid_ni
         dims(2) = dimid_ncat
         dims1(1) = dimid_ni 
         call define_rest_field(ncid,'aicen',dims)
         call define_rest_field(ncid,'vicen',dims)
         call define_rest_field(ncid,'vsnon',dims)
         call define_rest_field(ncid,'Tsfcn',dims)
         do k=1,nilyr
            write(nchar,'(i3.3)') k
            call define_rest_field(ncid,'sice'//trim(nchar),dims)
            call define_rest_field(ncid,'qice'//trim(nchar),dims)
         enddo
         do k=1,nslyr
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'qsno'//trim(nchar),dims)
         enddo

         !-----------------------------------------------------------------
         ! radiation fields
         !-----------------------------------------------------------------
         call define_rest_field(ncid,'scale_factor',dims1)
         call define_rest_field(ncid,'swvdr',dims1)
         call define_rest_field(ncid,'swvdf',dims1)
         call define_rest_field(ncid,'swidr',dims1)
         call define_rest_field(ncid,'swidf',dims1)
         
         !-----------------------------------------------------------------
         ! for mixed layer model
         !-----------------------------------------------------------------
         if (oceanmixed_ice) then
            call define_rest_field(ncid,'sst',dims1)
            call define_rest_field(ncid,'frzmlt',dims1)
         endif
         
         ! tracers
         if (tr_iage)      call define_rest_field(ncid,'iage',dims)  ! ice age tracer
         if (tr_FY) then                                             ! first-year area tracer
            call define_rest_field(ncid,'FY',dims)
            call define_rest_field(ncid,'frz_onset',dims1)
         endif
         if (tr_lvl) then                                            ! level ice tracer
            call define_rest_field(ncid,'alvl',dims)
            call define_rest_field(ncid,'vlvl',dims)
         endif
         if (tr_pond_lvl) then                                       ! level-ice melt ponds
            call define_rest_field(ncid,'apnd',dims)
            call define_rest_field(ncid,'hpnd',dims)
            call define_rest_field(ncid,'ipnd',dims)
            call define_rest_field(ncid,'dhsn',dims)
            call define_rest_field(ncid,'ffracn',dims)
            call define_rest_field(ncid,'fsnow',dims1)
         endif
         if (tr_pond_topo) then                                      ! topographic melt ponds
            call define_rest_field(ncid,'apnd',dims)
            call define_rest_field(ncid,'hpnd',dims)
            call define_rest_field(ncid,'ipnd',dims)
         endif
         if (tr_snow) then                                           ! snow metamorphosis tracers
            do k=1,nslyr
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'smice'//trim(nchar),dims)    
               call define_rest_field(ncid,'smliq'//trim(nchar),dims)
               call define_rest_field(ncid,'rhos'//trim(nchar),dims)
               call define_rest_field(ncid,'rsnw'//trim(nchar),dims)
            enddo
         endif
         if (tr_iso) then                                            ! ice isotopes
            do k=1,n_iso
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'isosno'//trim(nchar),dims) 
               call define_rest_field(ncid,'isoice'//trim(nchar),dims)
            enddo
         endif
         if (tr_aero) then                                           ! ice aerosols
            do k = 1, n_aero
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'aerosnossl'//nchar,dims)
               call define_rest_field(ncid,'aerosnoint'//nchar,dims)
               call define_rest_field(ncid,'aeroicessl'//nchar,dims)
               call define_rest_field(ncid,'aeroiceint'//nchar,dims)
            enddo
         endif
         if (tr_brine) then                                          ! brine height
            call define_rest_field(ncid,'fbri',dims)
            call define_rest_field(ncid,'first_ice',dims)
         endif
         if (tr_fsd) then                                            ! floe size distribution
            do k = 1, nfsd
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'fsd'//nchar,dims)
            enddo 
         endif
         status = nf90_enddef(ncid)

         ! Populate fields in NetCDF restart file

         status = nf90_close(ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname, &
            file=__FILE__,line= __LINE__)
      else
         write(warnstr,*) subname, 'Restart format must be either "bin" or "nc", no restart file written'
         call icepack_warnings_add(warnstr)
         call icepack_warnings_flush(nu_diag)
      endif

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
         tr_iage, tr_FY, tr_lvl, tr_iso, tr_aero, tr_brine, &
         tr_pond_topo, tr_pond_lvl, tr_snow, tr_fsd

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
           tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_iso_out=tr_iso, &
           tr_brine_out=tr_brine, &
           tr_pond_topo_out=tr_pond_topo, &
           tr_pond_lvl_out=tr_pond_lvl,tr_snow_out=tr_snow,tr_fsd_out=tr_fsd)
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
      if (tr_pond_lvl)  call read_restart_pond_lvl()  ! level-ice melt ponds
      if (tr_pond_topo) call read_restart_pond_topo() ! topographic melt ponds
      if (tr_snow)      call read_restart_snow()      ! snow metamorphosis tracers
      if (tr_iso)       call read_restart_iso()       ! ice isotopes
      if (tr_aero)      call read_restart_aero()      ! ice aerosols
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
        minw, maxw, sumw    ! diagnostics

      character(len=*), parameter :: subname='(read_restart_field)'

      do n = 1, ndim
         read(nu) (work2(i), i=1,nx)
         work(:,n) = work2(:)
      enddo

      minw = minval(work)
      maxw = maxval(work)
      sumw = sum(work)
      write(nu_diag,*) subname, minw, maxw, sumw

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

      real (kind=dbl_kind) :: &
        minw, maxw, sumw    ! diagnostics

      character(len=*), parameter :: subname='(write_restart_field)'

      do n = 1, ndim
        work2(:) = work(:,n)
        write(nu) (work2(i), i=1,nx)
      enddo

      minw = minval(work)
      maxw = maxval(work)
      sumw = sum(work)
      write(nu_diag,*) subname, minw, maxw, sumw

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

! Dumps values needed to restart snow redistribution/metamorphism tracers
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_snow()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: nslyr, ncat

      integer (kind=int_kind) :: nt_smice, nt_smliq, nt_rhos, nt_rsnw, k
      character(len=*), parameter :: subname='(write_restart_snow)'

      call icepack_query_tracer_indices(nt_smice_out=nt_smice, &
           nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      do k=1,nslyr
         call write_restart_field(nu_dump,trcrn(:,nt_smice+k-1,:),ncat)
         call write_restart_field(nu_dump,trcrn(:,nt_smliq+k-1,:),ncat)
         call write_restart_field(nu_dump,trcrn(:,nt_rhos +k-1,:),ncat)
         call write_restart_field(nu_dump,trcrn(:,nt_rsnw +k-1,:),ncat)
      enddo

      end subroutine write_restart_snow

!=======================================================================

! Reads all values needed to restart snow redistribution/metamorphism
!
! authors Elizabeth C. Hunke, LANL

      subroutine read_restart_snow()

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: nslyr, ncat

      integer (kind=int_kind) :: nt_smice, nt_smliq, nt_rhos, nt_rsnw, k
      character(len=*), parameter :: subname='(read_restart_snow)'

      call icepack_query_tracer_indices(nt_smice_out=nt_smice, &
           nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'min/max snow metamorphosis tracers'

      do k=1,nslyr
         call read_restart_field(nu_restart,trcrn(:,nt_smice+k-1,:),ncat)
         call read_restart_field(nu_restart,trcrn(:,nt_smliq+k-1,:),ncat)
         call read_restart_field(nu_restart,trcrn(:,nt_rhos +k-1,:),ncat)
         call read_restart_field(nu_restart,trcrn(:,nt_rsnw +k-1,:),ncat)
      enddo

      end subroutine read_restart_snow

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

! Dumps all values needed for restarting
!
! authors Elizabeth Hunke, LANL
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine write_restart_iso()

      use icedrv_domain_size, only: n_iso
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices
      integer (kind=int_kind) :: nt_isosno, nt_isoice
      character(len=*), parameter :: subname='(write_restart_iso)'

      call icepack_query_tracer_indices(nt_isosno_out=nt_isosno)
      call icepack_query_tracer_indices(nt_isoice_out=nt_isoice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'write_restart_iso (isotopes)'

      do k = 1, n_iso
         call write_restart_field(nu_dump, trcrn(:,nt_isosno+(k-1),:), ncat)
         call write_restart_field(nu_dump, trcrn(:,nt_isoice+(k-1),:), ncat)
      enddo

      end subroutine write_restart_iso

!=======================================================================

! Reads all values needed for an ice isotope restart
!
! authors Elizabeth Hunke, LANL
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine read_restart_iso()

      use icedrv_domain_size, only: n_iso
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices
      integer (kind=int_kind) :: nt_isosno, nt_isoice
      character(len=*), parameter :: subname='(read_restart_iso)'

      !-----------------------------------------------------------------

      call icepack_query_tracer_indices(nt_isosno_out=nt_isosno)
      call icepack_query_tracer_indices(nt_isoice_out=nt_isoice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(nu_diag,*) 'read_restart_iso (isotopes)'

      do k = 1, n_iso
         call read_restart_field(nu_restart, trcrn(:,nt_isosno+(k-1),:), ncat)
         call read_restart_field(nu_restart, trcrn(:,nt_isoice+(k-1),:), ncat)
      enddo

      end subroutine read_restart_iso

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
      subroutine define_rest_field(ncid, vname, dims)

! Defines a field in NetCDF restart file
! author Chris Riedel, NCAR

         character (len=*)      , intent(in)  :: vname
         integer (kind=int_kind), intent(in)  :: dims(:)
         integer (kind=int_kind), intent(in)  :: ncid
   
         integer (kind=int_kind) :: varid
   
         integer (kind=int_kind) :: &
           status        ! status variable from netCDF routine
   
         status = nf90_def_var(ncid,trim(vname),nf90_double,dims,varid)
   
         end subroutine define_rest_field
!=======================================================================

      end module icedrv_restart

!=======================================================================
