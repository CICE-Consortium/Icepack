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
          write_restart_aero,      read_restart_aero

      character (len=3), private :: nchar !

      logical (kind=log_kind), private :: diag ! netCDF diagnostics flag

#ifdef USE_NETCDF
      private :: &
         define_rest_field,         write_restart_field_nc2D, &
         write_restart_field_nc1D,  ice_write_nc2D, &
         ice_write_nc1D,            read_restart_field_net2D, &
         read_restart_field_net1D,  ice_read_nc2D, &
         ice_read_nc1D

      integer (kind=int_kind), private :: &
         ncid     ! ID for NetCDF file
#endif

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
      use icedrv_flux, only: sst, frzmlt, frz_onset, fsnow
      use icedrv_state, only: aicen, vicen, vsnon, trcrn
      use icedrv_arrays_column, only: dhsn, ffracn
      use icedrv_arrays_column, only: first_ice, first_ice_real

      ! local variables

      integer (kind=int_kind) :: &
         k, &              ! counting indices
         iyear

      integer (kind=int_kind) :: &
         nt_Tsfc, nt_sice, nt_qice, nt_qsno

      logical (kind=log_kind) :: &
         tr_iage, tr_FY, tr_lvl, tr_iso, tr_aero, tr_brine, &
         tr_pond_topo, tr_pond_lvl, tr_snow, tr_fsd
!         skl_bgc, z_tracers

      integer (kind=int_kind) :: dims(2)

      character(len=char_len_long) :: filename
      character(len=*), parameter :: subname='(dumpfile)'

#ifdef USE_NETCDF
      ! local variables if we're writing in NetCDF
      integer (kind=int_kind) :: &
      dimid_ni,   & ! netCDF identifiers
      dimid_ncat, & !
      iflag,      & ! netCDF creation flag
      status        ! status variable from netCDF routine
#endif

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
!      call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__,line= __LINE__)

      dims(:) = 0
      if (restart_format == 'bin') then
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
            restart_dir(1:lenstr(restart_dir)), &
            restart_file(1:lenstr(restart_file)),'.', &
            iyear,'-',month,'-',mday,'-',sec

         open(nu_dump,file=filename,form='unformatted')
         write(nu_dump) istep1,time,time_forc
      else if (restart_format == 'nc') then
#ifdef USE_NETCDF
         ! Change this if you want diagnostic output
         diag = .false.

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

         ! Get dimension IDs
         dims(1) = dimid_ni
         dims(2) = dimid_ncat
#else
   call icedrv_system_abort(string=subname//' ERROR: restart_format = "nc" requires USE_NETCDF',file=__FILE__,line=__LINE__)
#endif
      endif

      write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to binary files.  All other
      ! tracers are written to their own dump/restart binary files.
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump,aicen(:,:),ncat,'aicen',dims)
      call write_restart_field(nu_dump,vicen(:,:),ncat,'vicen',dims)
      call write_restart_field(nu_dump,vsnon(:,:),ncat,'vsnon',dims)
      call write_restart_field(nu_dump,trcrn(:,nt_Tsfc,:),ncat,'Tsfcn',dims)
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,trcrn(:,nt_sice+k-1,:),ncat,&
            'sice'//trim(nchar),dims)
      enddo
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,trcrn(:,nt_qice+k-1,:),ncat,&
            'qice'//trim(nchar),dims)
      enddo
      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,trcrn(:,nt_qsno+k-1,:),ncat,&
            'qsno'//trim(nchar),dims)
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,scale_factor,1,'scale_factor',dims)
      call write_restart_field(nu_dump,swvdr,1,'swvdr',dims)
      call write_restart_field(nu_dump,swvdf,1,'swvdf',dims)
      call write_restart_field(nu_dump,swidr,1,'swidr',dims)
      call write_restart_field(nu_dump,swidf,1,'swidf',dims)

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------
      if (oceanmixed_ice) then
         call write_restart_field(nu_dump,sst,1,'sst',dims)
         call write_restart_field(nu_dump,frzmlt,1,'frzmlt',dims)
      endif

      ! tracers
      if (tr_iage)      call write_restart_age(dims)       ! ice age tracer
      if (tr_FY)        call write_restart_FY(dims)        ! first-year area tracer
      if (tr_lvl)       call write_restart_lvl(dims)       ! level ice tracer
      if (tr_pond_lvl)  call write_restart_pond_lvl(dims)  ! level-ice melt ponds
      if (tr_pond_topo) call write_restart_pond_topo(dims) ! topographic melt ponds
      if (tr_snow)      call write_restart_snow(dims)      ! snow metamorphosis tracers
      if (tr_iso)       call write_restart_iso(dims)       ! ice isotopes
      if (tr_aero)      call write_restart_aero(dims)      ! ice aerosols
      if (tr_brine)     call write_restart_hbrine(dims)    ! brine height
      if (tr_fsd)       call write_restart_fsd(dims)       ! floe size distribution
   ! called in icedrv_RunMod.F90 to prevent circular dependencies
   !      if (skl_bgc .or. z_tracers) &
   !                        call write_restart_bgc         ! biogeochemistry

      ! Close netcdf file
      if (restart_format == 'nc') then
#ifdef USE_NETCDF
         status = nf90_close(ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname, &
            file=__FILE__,line= __LINE__)
#endif
      endif

      end subroutine dumpfile

!=======================================================================

! Restarts from a dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile (ice_ic)

      use icedrv_calendar, only: istep0, istep1, time, time_forc
      use icepack_intfc, only: icepack_aggregate
      use icedrv_domain_size, only: nilyr, nslyr, ncat, nfsd, nblyr
      use icedrv_domain_size, only: max_ntrcr, nx, n_iso, n_aero
      use icedrv_flux, only: swvdr, swvdf, swidr, swidf, Tf
      use icedrv_flux, only: sst, frzmlt, scale_factor
      use icedrv_flux, only: frz_onset,fsnow
      use icedrv_forcing, only: oceanmixed_ice
      use icedrv_init, only: tmask
      use icedrv_state, only: trcr_depend, aice, vice, vsno, trcr
      use icedrv_state, only: aice0, aicen, vicen, vsnon, trcrn, aice_init
      use icedrv_state, only: trcr_base, nt_strata, n_trcr_strata
      use icedrv_restart_shared, only: restart_format
      use icedrv_arrays_column, only: dhsn, ffracn, hin_max
      use icedrv_arrays_column, only: first_ice, first_ice_real
      use icepack_tracers, only: ntrcr, nbtrcr

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

      ! local variable for reading from a netcdf file
      integer (kind=int_kind) :: &
         status

      ! Query tracers
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

      if (present(ice_ic)) then
         filename = trim(ice_ic)
      else
         call icedrv_system_abort(string=subname//'no ice_ic present', &
                                  file=__FILE__,line=__LINE__)
      endif

      write(nu_diag,*) 'Using restart dump=', trim(filename)

      if (restart_format == 'bin') then
         open(nu_restart,file=filename,form='unformatted')
         read (nu_restart) istep0,time,time_forc
      else if (restart_format == 'nc') then
#ifdef USE_NETCDF
         ! set this to .true. for netcdf diagnostic output
         diag = .false.
         ! Open restart files
         status = nf90_open(trim(filename), nf90_nowrite, ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//'Couldnt open netcdf file', &
                                    file=__FILE__,line=__LINE__)

         status = nf90_get_att(ncid, nf90_global, 'istep1', istep0)
         status = nf90_get_att(ncid, nf90_global, 'time', time)
         status = nf90_get_att(ncid, nf90_global, 'time_forc', time_forc)
#else
         call icedrv_system_abort(string=subname//' ERROR: restart_format = "nc" requires USE_NETCDF',file=__FILE__,line=__LINE__)
#endif
      else
         call icedrv_system_abort(string=subname//'unrecognized restart format', &
                                 file=__FILE__,line=__LINE__)
      endif

      write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc

      istep1 = istep0

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      write(nu_diag,*) 'min/max area, vol ice, vol snow, Tsfc'

      call read_restart_field(nu_restart,aicen,ncat,'aicen')
      call read_restart_field(nu_restart,vicen,ncat,'vicen')
      call read_restart_field(nu_restart,vsnon,ncat,'vsnon')
      call read_restart_field(nu_restart,trcrn(:,nt_Tsfc,:),ncat,'Tsfcn')

      write(nu_diag,*) 'min/max sice for each layer'
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,trcrn(:,nt_sice+k-1,:),ncat, &
         'sice'//trim(nchar))
      enddo

      write(nu_diag,*) 'min/max qice for each layer'
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,trcrn(:,nt_qice+k-1,:),ncat, &
         'qice'//trim(nchar))
      enddo

      write(nu_diag,*) 'min/max qsno for each layer'
      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,trcrn(:,nt_qsno+k-1,:),ncat, &
         'qsno'//trim(nchar))
      enddo

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      write(nu_diag,*) 'min/max radiation fields'

      call read_restart_field(nu_restart,scale_factor,1,'scale_factor')
      call read_restart_field(nu_restart,swvdr,1,'swvdr')
      call read_restart_field(nu_restart,swvdf,1,'swvdf')
      call read_restart_field(nu_restart,swidr,1,'swidr')
      call read_restart_field(nu_restart,swidf,1,'swidf')

      !-----------------------------------------------------------------
      ! for mixed layer model
      !-----------------------------------------------------------------

      if (oceanmixed_ice) then
         write(nu_diag,*) 'min/max sst, frzmlt'
         call read_restart_field(nu_restart,sst,1,'sst')
         call read_restart_field(nu_restart,frzmlt,1,'frzmlt')
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
                                 nt_strata=nt_strata, &
                                 Tf = Tf(i))

         aice_init(i) = aice(i)
      enddo
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

      end subroutine restartfile

!=======================================================================

! Reads a single restart field
! author Chris Newman, LANL

      subroutine read_restart_field(nu,work,ndim,name)

      use icedrv_domain_size, only: nx

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number (not used for netcdf)
         ndim             ! number of dimensions

      real (kind=dbl_kind), dimension(nx,ndim), intent(inout) :: &
         work              ! input array (real, 8-byte)

      character (len=*), intent(in), optional :: &
         name

      ! local variables

      integer (kind=int_kind) :: &
         n, i               ! loop indices

      real (kind=dbl_kind), dimension(nx) :: &
         work2              ! input array (real, 8-byte)

      real (kind=dbl_kind) :: &
        minw, maxw, sumw    ! diagnostics

      character(len=*), parameter :: subname='(read_restart_field)'

      if (restart_format == 'bin') then
         do n = 1, ndim
            read(nu) (work2(i), i=1,nx)
            work(:,n) = work2(:)
         enddo
      else if (restart_format == 'nc') then
#ifdef USE_NETCDF
         if (ndim == 1) then
            call read_restart_field_net1D(ncid,1,work,name,1,diag)
         else
            call read_restart_field_net2D(ncid,1,work,name,ndim,diag)
         endif
#else
         call icedrv_system_abort(string=subname//' ERROR: restart_format = "nc" requires USE_NETCDF',file=__FILE__,line=__LINE__)
#endif
      else
         call icedrv_system_abort(string=subname//'unrecognized restart format', &
                                  file=__FILE__,line=__LINE__)
      endif

      minw = minval(work)
      maxw = maxval(work)
      sumw = sum(work)
      write(nu_diag,*) subname, minw, maxw, sumw

      end subroutine read_restart_field

!=======================================================================

! Writes a single restart field.
! author Chris Newman, LANL

      subroutine write_restart_field(nu,work,ndim,vname,dims)

      use icedrv_domain_size, only: nx

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         ndim              ! number of dimensions

      real (kind=dbl_kind), dimension(nx,ndim), intent(in) :: &
         work              ! input array (real, 8-byte)

      character (len=*), intent(in), optional :: &
         vname             ! variable name for netcdf output

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      ! local variables

      integer (kind=int_kind) :: &
         n, i              ! loop indices

      real (kind=dbl_kind), dimension(nx) :: &
         work2             ! input array (real, 8-byte)

      real (kind=dbl_kind) :: &
        minw, maxw, sumw    ! diagnostics

      character(len=*), parameter :: subname='(write_restart_field)'

#ifdef USE_NETCDF
      ! local variables if we're writing in NetCDF
      integer (kind=int_kind) :: &
         status, &          ! Status variable from netCDF routine
         varid              ! Variable ID
#endif

      if (restart_format == 'bin') then
         do n = 1, ndim
         work2(:) = work(:,n)
         write(nu) (work2(i), i=1,nx)
         enddo
      else if (restart_format == 'nc') then
#ifdef USE_NETCDF
         ! Check whether the variable already exists in the netcdf file
         status = nf90_inq_varid(ncid,trim(vname),varid)
         if (status == nf90_enotvar) then ! if not add it to the file
            status = nf90_redef(ncid) ! Put file into define mode
            if (ndim == 1) then
               ! Handle 1D cases
               status = nf90_def_var(ncid,trim(vname),nf90_double,dims(1:1),varid)
            else
               status = nf90_def_var(ncid,trim(vname),nf90_double,dims,varid)
            endif
         endif
         status = nf90_enddef(ncid) ! put file into data mode
         if (ndim == 1) then
            call write_restart_field_nc1D(ncid,0,work,'ruf8',vname,1,diag)
         else
            call write_restart_field_nc2D(ncid,0,work,'ruf8',vname,ndim,diag)
         endif
#else
   call icedrv_system_abort(string=subname//' ERROR: restart_format = "nc" requires USE_NETCDF',file=__FILE__,line=__LINE__)
#endif
      else
         write (nu_diag,*) 'WARNING: Restart format must be either "bin" or "nc", no restart file written'
      endif

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

      subroutine write_restart_pond_topo(dims)

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*), parameter :: subname='(write_restart_pond_topo)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_apnd,:),ncat,'apnd',dims)
      call write_restart_field(nu_dump,trcrn(:,nt_hpnd,:),ncat,'hpnd',dims)
      call write_restart_field(nu_dump,trcrn(:,nt_ipnd,:),ncat,'ipnd',dims)

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

      call read_restart_field(nu_restart,trcrn(:,nt_apnd,:),ncat,'apnd')
      call read_restart_field(nu_restart,trcrn(:,nt_hpnd,:),ncat,'hpnd')
      call read_restart_field(nu_restart,trcrn(:,nt_ipnd,:),ncat,'ipnd')

      end subroutine read_restart_pond_topo

!=======================================================================

! Dumps values needed to restart snow redistribution/metamorphism tracers
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_snow(dims)

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: nslyr, ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_smice, nt_smliq, nt_rhos, nt_rsnw, k
      character(len=*), parameter :: subname='(write_restart_snow)'

      call icepack_query_tracer_indices(nt_smice_out=nt_smice, &
           nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,trcrn(:,nt_smice+k-1,:),ncat,'smice'//trim(nchar),dims)
         call write_restart_field(nu_dump,trcrn(:,nt_smliq+k-1,:),ncat,'smliq'//trim(nchar),dims)
         call write_restart_field(nu_dump,trcrn(:,nt_rhos +k-1,:),ncat,'rhos'//trim(nchar),dims)
         call write_restart_field(nu_dump,trcrn(:,nt_rsnw +k-1,:),ncat,'rsnw'//trim(nchar),dims)
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
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,trcrn(:,nt_smice+k-1,:),ncat,'smice'//trim(nchar))
         call read_restart_field(nu_restart,trcrn(:,nt_smliq+k-1,:),ncat,'smliq'//trim(nchar))
         call read_restart_field(nu_restart,trcrn(:,nt_rhos +k-1,:),ncat,'rhos'//trim(nchar))
         call read_restart_field(nu_restart,trcrn(:,nt_rsnw +k-1,:),ncat,'rsnw'//trim(nchar))
      enddo

      end subroutine read_restart_snow

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_age(dims)

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_iage
      character(len=*), parameter :: subname='(write_restart_age)'

      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_iage,:),ncat,'iage',dims)

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

      call read_restart_field(nu_restart,trcrn(:,nt_iage,:),ncat,'iage')

      end subroutine read_restart_age

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_fsd(dims)

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat, nfsd

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_fsd, k
      character(len=*), parameter :: subname='(write_restart_fsd)'

      call icepack_query_tracer_indices(nt_fsd_out=nt_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      do k = 1, nfsd
         write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump,trcrn(:,nt_fsd+k-1,:),ncat,&
            'fsd'//trim(nchar),dims)
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
         write(nchar,'(i3.3)') k
          call read_restart_field(nu_restart,trcrn(:,nt_fsd+k-1,:),ncat,'fsd'//trim(nchar))
      end do

      end subroutine read_restart_fsd

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_FY(dims)

      use icedrv_flux, only: frz_onset
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_FY
      character(len=*), parameter :: subname='(write_restart_FY)'

      call icepack_query_tracer_indices(nt_FY_out=nt_FY)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_FY,:),ncat,'FY',dims)
      call write_restart_field(nu_dump,frz_onset,1,'frz_onset',dims)

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

      call read_restart_field(nu_restart,trcrn(:,nt_FY,:),ncat,'FY')

      write(nu_diag,*) 'min/max frz_onset'

      call read_restart_field(nu_restart,frz_onset,1,'frz_onset')

      end subroutine read_restart_FY

!=======================================================================

! Dumps all values needed for restarting
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_lvl(dims)

      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_alvl, nt_vlvl
      character(len=*), parameter :: subname='(write_restart_lvl)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_alvl,:),ncat,'alvl',dims)
      call write_restart_field(nu_dump,trcrn(:,nt_vlvl,:),ncat,'vlvl',dims)

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

      call read_restart_field(nu_restart,trcrn(:,nt_alvl,:),ncat,'alvl')
      call read_restart_field(nu_restart,trcrn(:,nt_vlvl,:),ncat,'vlvl')

      end subroutine read_restart_lvl

!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_pond_lvl(dims)

      use icedrv_arrays_column, only: dhsn, ffracn
      use icedrv_flux, only: fsnow
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*), parameter :: subname='(write_restart_pond_lvl)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call write_restart_field(nu_dump,trcrn(:,nt_apnd,:),ncat,'apnd',dims)
      call write_restart_field(nu_dump,trcrn(:,nt_hpnd,:),ncat,'hpnd',dims)
      call write_restart_field(nu_dump,trcrn(:,nt_ipnd,:),ncat,'ipnd',dims)
      call write_restart_field(nu_dump,fsnow(:),1,'fsnow',dims)
      call write_restart_field(nu_dump,dhsn(:,:),ncat,'dhsn',dims)
      call write_restart_field(nu_dump,ffracn(:,:),ncat,'ffracn',dims)

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

      call read_restart_field(nu_restart,trcrn(:,nt_apnd,:),ncat,'apnd')
      call read_restart_field(nu_restart,trcrn(:,nt_hpnd,:),ncat,'hpnd')
      call read_restart_field(nu_restart,trcrn(:,nt_ipnd,:),ncat,'ipnd')
      call read_restart_field(nu_restart,fsnow(:),1,'fsnow')
      call read_restart_field(nu_restart,dhsn(:,:),ncat,'dhsn')
      call read_restart_field(nu_restart,ffracn(:,:),ncat,'ffracn')

      end subroutine read_restart_pond_lvl

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth Hunke, LANL
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine write_restart_aero(dims)

      use icedrv_domain_size, only: n_aero
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

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
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump, trcrn(:,nt_aero  +(k-1)*4,:),ncat,&
            'aerosnossl'//trim(nchar),dims)
         call write_restart_field(nu_dump, trcrn(:,nt_aero+1+(k-1)*4,:),ncat,&
            'aerosnoint'//trim(nchar),dims)
         call write_restart_field(nu_dump, trcrn(:,nt_aero+2+(k-1)*4,:),ncat,&
            'aeroicessl'//trim(nchar),dims)
         call write_restart_field(nu_dump, trcrn(:,nt_aero+3+(k-1)*4,:),ncat,&
            'aeroiceint'//trim(nchar),dims)
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
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart, trcrn(:,nt_aero  +(k-1)*4,:),ncat,'aerosnossl'//trim(nchar))
         call read_restart_field(nu_restart, trcrn(:,nt_aero+1+(k-1)*4,:),ncat,'aerosnoint'//trim(nchar))
         call read_restart_field(nu_restart, trcrn(:,nt_aero+2+(k-1)*4,:),ncat,'aeroicessl'//trim(nchar))
         call read_restart_field(nu_restart, trcrn(:,nt_aero+3+(k-1)*4,:),ncat,'aeroiceint'//trim(nchar))
      enddo

      end subroutine read_restart_aero

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth Hunke, LANL
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine write_restart_iso(dims)

      use icedrv_domain_size, only: n_iso
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

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
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump, trcrn(:,nt_isosno+(k-1),:),ncat,&
            'isosno'//trim(nchar),dims)
         call write_restart_field(nu_dump, trcrn(:,nt_isoice+(k-1),:),ncat,&
            'isoice'//trim(nchar),dims)
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
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart, trcrn(:,nt_isosno+(k-1),:),ncat,'isosno'//trim(nchar))
         call read_restart_field(nu_restart, trcrn(:,nt_isoice+(k-1),:),ncat,'isoice'//trim(nchar))
      enddo

      end subroutine read_restart_iso

!=======================================================================

      subroutine write_restart_hbrine(dims)

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use icedrv_arrays_column, only: first_ice, first_ice_real
      use icedrv_state, only: trcrn
      use icedrv_domain_size, only: ncat, nx

      integer (kind=int_kind), intent(in), optional :: &
         dims(:)           ! netcdf dimension IDs

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

        call write_restart_field(nu_dump,trcrn(:,nt_fbri,:),ncat,'fbri',dims)
        call write_restart_field(nu_dump,first_ice_real(:,:),ncat,'first_ice',dims)

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

      call read_restart_field(nu_restart,trcrn(:,nt_fbri,:),ncat,'fbri')
      call read_restart_field(nu_restart,first_ice_real(:,:),ncat,'first_ice')

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
#ifdef USE_NETCDF
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
      subroutine write_restart_field_nc2D(ncid,nrec,work,atype,vname,ndim3,diag)

! Write a 2D field (nx, ncat) in NetCDF restart
! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx
      !use ice_fileunits, only: nu_diag
      !use ice_read_write, only: ice_write, ice_write_nc

      integer (kind=int_kind), intent(in) :: &
            ncid            , & ! unit number
            ndim3         , & ! third dimension
            nrec              ! record number (0 for sequential access)
      real (kind=dbl_kind), intent(in), dimension(nx,ncat) :: &
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


      status = nf90_inq_varid(ncid,trim(vname),varid)
      write(nu_diag,*) 'Writing out ',trim(vname)
      call ice_write_nc2D(ncid, 1, varid, work, diag)

      end subroutine write_restart_field_nc2D
!=======================================================================
      subroutine write_restart_field_nc1D(ncid,nrec,work,atype,vname,ndim3,diag)

! Write a 1D field (nx) in NetCDF restart
! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
            ncid            , & ! unit number
            ndim3         , & ! third dimension
            nrec              ! record number (0 for sequential access)
      real (kind=dbl_kind), intent(in), dimension(nx) :: &
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


      status = nf90_inq_varid(ncid,trim(vname),varid)
      if (status /= 0) then
         write(nu_diag,*) 'Writing out ',trim(vname)
         write(nu_diag,*) 'Erros Status ',status
         call icedrv_system_abort(string='Write out restart', &
            file=__FILE__,line= __LINE__)
      else
         write(nu_diag,*) 'Writing out ',trim(vname)
      endif
      call ice_write_nc1D(ncid, 1, varid, work, diag)

      end subroutine write_restart_field_nc1D
!=======================================================================
      subroutine ice_write_nc2D(fid,  nrec,  varid, work,  diag)

! Write a 2D field (nx, ncat) in NetCDF restart
! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
            fid           , & ! file id
            varid         , & ! variable id
            nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
            diag              ! if true, write diagnostic output

      real (kind=dbl_kind), intent(in), dimension(nx,ncat) :: &
            work              ! output array (real, 8-byte)

      ! local variables

      integer (kind=int_kind) :: &
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension
      integer (kind=int_kind) :: start2(2),count2(2)

      real (kind=dbl_kind) :: &
         amin, amax, asum   ! min, max values and sum of input array

      character (char_len) :: &
         dimname            ! dimension name

      integer (kind=int_kind) :: ny,numDims

      ny = ncat
      start2(1) = 1
      count2(1) = nx
      start2(2) = 1
      count2(2) = ncat

      status = nf90_put_var( fid, varid, work, &
               start=start2, &
               count=count2)

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (diag) then
         amin = minval(work)
         amax = maxval(work)
         !asum = sum   (work)
         write(nu_diag,*) ' min, max =', amin, amax
      endif

      !deallocate(work)

      end subroutine ice_write_nc2D
!=======================================================================
      subroutine ice_write_nc1D(fid,  nrec,  varid, work,  diag)

! Write a 1D field (nx) in NetCDF restart
! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
            fid           , & ! file id
            varid         , & ! variable id
            nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
            diag              ! if true, write diagnostic output

      real (kind=dbl_kind), intent(in), dimension(nx) :: &
            work              ! output array (real, 8-byte)

      ! local variables

      integer (kind=int_kind) :: &
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension
      integer (kind=int_kind) :: start2(1),count2(1)

      real (kind=dbl_kind) :: &
         amin, amax, asum   ! min, max values and sum of input array

      character (char_len) :: &
         dimname            ! dimension name

      integer (kind=int_kind) :: ny,numDims

      ny = ncat
      start2(1) = 1
      count2(1) = nx

      status = nf90_put_var( fid, varid, work, &
               start=start2, &
               count=count2)

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (diag) then
         amin = minval(work)
         amax = maxval(work)
         !asum = sum   (work)
         write(nu_diag,*) ' min, max =', amin, amax
      endif

      !deallocate(work)

      end subroutine ice_write_nc1D
!=======================================================================
      subroutine read_restart_field_net2D(nu,nrec,work,vname,ndim3, &
                     diag)

! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
      nu            , & ! unit number (not used for netcdf)
      ndim3         , & ! third dimension
      nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), intent(inout), dimension(nx,ncat) :: &
      work              ! input array (real, 8-byte)

      logical (kind=log_kind), intent(in) :: &
      diag              ! if true, write diagnostic output

      character (len=*), intent(in)  :: vname

      ! local variables

      integer (kind=int_kind) :: &
      n,     &      ! number of dimensions for variable
      varid, &      ! variable id
      status        ! status variable from netCDF routine

      call ice_read_nc2D(ncid, 1, vname, work, diag)

      end subroutine read_restart_field_net2D
!=======================================================================
      subroutine read_restart_field_net1D(nu,nrec,work,vname,ndim3, &
                  diag)

! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
      nu            , & ! unit number (not used for netcdf)
      ndim3         , & ! third dimension
      nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), intent(inout), dimension(nx) :: &
      work              ! input array (real, 8-byte)

      logical (kind=log_kind), intent(in) :: &
      diag              ! if true, write diagnostic output

      character (len=*), intent(in)  :: vname

      ! local variables

      integer (kind=int_kind) :: &
      n,     &      ! number of dimensions for variable
      varid, &      ! variable id
      status        ! status variable from netCDF routine

      call ice_read_nc1D(ncid, 1, vname, work, diag)

      end subroutine read_restart_field_net1D
!=======================================================================
      subroutine ice_read_nc2D(fid,  nrec,  varname, work,  diag)

! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
            fid           , & ! file id
            nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
            diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: &
            varname           ! field name in netcdf file

      real (kind=dbl_kind), intent(out), dimension(nx,ncat) :: &
            work              ! output array (real, 8-byte)

      integer (kind=int_kind) :: &
         varid          , & ! variable id
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax, asum   ! min, max values and sum of input array

      character (char_len) :: &
         dimname            ! dimension name

      integer (kind=int_kind) :: start2(2),count2(2)

      integer (kind=int_kind) :: ny,numDims

      ny = ncat
      start2(1) = 1
      count2(1) = nx
      start2(2) = 1
      count2(2) = ncat

      status = nf90_inq_varid(fid, trim(varname), varid)

      if (status /= nf90_noerr) then
         call icedrv_system_abort(string='ice_read_nc2D:Var not found '//trim(varname), &
                                    file=__FILE__,line=__LINE__)
      endif

      status = nf90_get_var( fid, varid, work, &
               start=start2, &
               count=count2 )


      if (diag) then
         amin = minval(work)
         amax = maxval(work)
         asum = sum   (work)
         write(nu_diag,*) ' min, max, sum =', amin, amax, asum
      endif

      end subroutine ice_read_nc2D
!=======================================================================
      subroutine ice_read_nc1D(fid,  nrec,  varname, work,  diag)

! author Chris Riedel, NCAR

      use icedrv_domain_size, only: ncat,nx

      integer (kind=int_kind), intent(in) :: &
            fid           , & ! file id
            nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
            diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: &
            varname           ! field name in netcdf file

      real (kind=dbl_kind), intent(out), dimension(nx) :: &
            work              ! output array (real, 8-byte)

      integer (kind=int_kind) :: &
         varid          , & ! variable id
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax, asum   ! min, max values and sum of input array

      character (char_len) :: &
         dimname            ! dimension name

      integer (kind=int_kind) :: start2(1),count2(1)

      integer (kind=int_kind) :: ny,numDims

      ny = ncat
      start2(1) = 1
      count2(1) = nx

      status = nf90_inq_varid(fid, trim(varname), varid)

      if (status /= nf90_noerr) then
         call icedrv_system_abort(string='ice_read_nc2D:Var not found'//trim(varname), &
                                    file=__FILE__,line=__LINE__)
      endif

      status = nf90_get_var( fid, varid, work, &
               start=start2, &
               count=count2 )

      if (diag) then
         amin = minval(work)
         amax = maxval(work)
         asum = sum   (work)
         write(nu_diag,*) ' min, max, sum =', amin, amax, asum
      endif

      end subroutine ice_read_nc1D
#endif
!=======================================================================
      end module icedrv_restart

!=======================================================================
