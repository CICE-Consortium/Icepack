!=======================================================================

! Diagnostic information output during run
!
! authors: T. Craig

      module icedrv_history

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, nu_diag_out
      use icedrv_domain_size, only: nx, ncat, nfsd
      use icedrv_diagnostics, only: nx_names
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private
      public :: history_write, &
                history_close

      ! history output file info

      logical (kind=log_kind), public :: history_cdf      ! flag to turn on cdf history files

      character (len=char_len_long) :: hist_file  ! hist file name

      integer (kind=int_kind) :: ncid             ! cdf file id
      integer (kind=int_kind) :: nxid, ncatid, ntrcrid, nfsdid, timid     ! cdf dim ids
      integer (kind=int_kind) :: timcnt           ! time counter

!=======================================================================

      contains

!=======================================================================

! Writes history information

      subroutine history_write()

      use icedrv_calendar, only: idate0, days_per_year, use_leap_years
      use icedrv_calendar, only: time, time0, secday, istep1, idate, sec
      use icedrv_state, only: aice, vice, vsno, uvel, vvel, divu, shear, strength
      use icedrv_state, only: trcr, trcrn
      use icedrv_state, only: aicen, vicen, vsnon
      use icedrv_flux, only: evap, fsnow, frain, frazil
      use icedrv_flux, only: fswabs, flw, flwout, fsens, fsurf, flat
      use icedrv_flux, only: Tair, Qa, fsw, fcondtop
      use icedrv_flux, only: meltt, meltb, meltl, snoice
      use icedrv_flux, only: dsnow, congel, sst, sss, Tf, fhocn
      use icedrv_arrays_column, only: d_afsd_newi, d_afsd_latg, d_afsd_latm, d_afsd_wave, d_afsd_weld
#ifdef USE_NETCDF
      use netcdf
#endif

      ! local variables

      logical (kind=log_kind), save :: &
         first_call = .true.             ! first call flag

      integer (kind=int_kind) :: &
         n, &                            ! counters
         ntrcr, &                        ! tracer count from icepack
         dimid1(1), dimid2(2), dimid3(3), dimid4(4), & ! cdf dimids
         start1(1), start2(2), start3(3), start4(4), & ! cdf start/count arrays
         count1(1), count2(2), count3(3), count4(4), & ! cdf start/count arrays
         varid, &                        ! cdf varid
         status, &                       ! cdf status flag
         iflag                           ! history file attributes

      character (len=8) :: &
         cdate                           ! date string

      real (kind=dbl_kind) :: &
         value                           ! temporary
      real (kind=dbl_kind),allocatable :: &
         value1(:), value2(:,:), value3(:,:,:), value4(:,:,:,:)  ! temporary

      integer (kind=dbl_kind), parameter :: num_2d = 32
      character(len=16), parameter :: fld_2d(num_2d) = &
         (/ 'aice            ', 'vice            ', 'vsno            ', &
            'uvel            ', 'vvel            ', 'divu            ', &
            'shear           ', 'strength        ',                     &
            'evap            ', 'fsnow           ', 'frazil          ', &
            'fswabs          ', 'flw             ', 'flwout          ', &
            'fsens           ', 'fsurf           ', 'flat            ', &
            'frain           ', 'Tair            ', 'Qa              ', &
            'fsw             ', 'fcondtop        ', 'meltt           ', &
            'meltb           ', 'meltl           ', 'snoice          ', &
            'dsnow           ', 'congel          ', 'sst             ', &
            'sss             ', 'Tf              ', 'fhocn           '    /)

      integer (kind=dbl_kind), parameter :: num_3d_ncat = 3
      character(len=16), parameter :: fld_3d_ncat(num_3d_ncat) = &
         (/ 'aicen           ', 'vicen           ', 'vsnon           ' /)

      logical (kind=log_kind) :: &
         tr_fsd                          ! flag for tracing fsd

      integer (kind=dbl_kind), parameter :: num_3d_nfsd = 5
      character(len=16), parameter :: fld_3d_nfsd(num_3d_nfsd) = &
         (/ 'd_afsd_newi     ', 'd_afsd_latg     ', 'd_afsd_latm     ', &
            'd_afsd_wave     ', 'd_afsd_weld     ' /)

      integer (kind=dbl_kind), parameter :: num_3d_ntrcr = 1
      character(len=16), parameter :: fld_3d_ntrcr(num_3d_ntrcr) = &
         (/ 'trcr            ' /)

      integer (kind=dbl_kind), parameter :: num_4d_ncat_ntrcr = 1
      character(len=16), parameter :: fld_4d_ncat_ntrcr(num_4d_ncat_ntrcr) = &
         (/ 'trcrn           ' /)

      character (len=char_len_long) :: tmpstr

      character(len=*), parameter :: subname='(history_write)'

#ifdef USE_NETCDF
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      if (first_call) then
         timcnt = 0
         write(hist_file,'(a,i8.8,a)') './history/icepack.h.',idate,'.nc'
         iflag = nf90_clobber
         status = nf90_create(trim(hist_file),iflag,ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: nf90_create '//trim(hist_file))

         ! nx columns dimension
         status = nf90_def_dim(ncid,'ni',nx,nxid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_dim ni')
         status = nf90_def_var(ncid,'ni',NF90_INT,nxid,varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var ni')
         do n = 1,nx
            write(tmpstr,'(a,i3.3)') 'column_name_',n
            status = nf90_put_att(ncid,varid,trim(tmpstr),trim(nx_names(n)))
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_att columns names')
         enddo

         ! ncat category dimension
         status = nf90_def_dim(ncid,'ncat',ncat,ncatid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_dim ncat')
         status = nf90_def_var(ncid,'ncat',NF90_INT,ncatid,varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var ncat')

         ! ntrcr dimension
         status = nf90_def_dim(ncid,'ntrcr',ntrcr,ntrcrid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_dim ntrcr')
         status = nf90_def_var(ncid,'ntrcr',NF90_INT,ntrcrid,varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var ntrcr')

         if (tr_fsd) then
            ! nfsd category dimension
            status = nf90_def_dim(ncid,'nfsd',nfsd,nfsdid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_dim nfsd')
            status = nf90_def_var(ncid,'nfsd',NF90_INT,nfsdid,varid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var nfsd')
         endif

         ! time dimension
         status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,timid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_dim time')
         status = nf90_def_var(ncid,'time',NF90_DOUBLE,timid,varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var time')
         status = nf90_put_att(ncid,varid,'long_name','model time')
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_att time long_name')
         write(cdate,'(i8.8)') idate0
         write(tmpstr,'(a,a,a,a,a,a,a,a)') 'days since ', &
            cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
         status = nf90_put_att(ncid,varid,'units',trim(tmpstr))
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_att time units')
         if (days_per_year == 360) then
            status = nf90_put_att(ncid,varid,'calendar','360_day')
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_att calendar 360_day')
         elseif (days_per_year == 365 .and. .not.use_leap_years ) then
            status = nf90_put_att(ncid,varid,'calendar','NoLeap')
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_att calendar noleap')
         elseif (use_leap_years) then
            status = nf90_put_att(ncid,varid,'calendar','Gregorian')
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_att calendar gregorian')
         else
            call icedrv_system_abort(string=subname//' ERROR: invalid calendar settings')
         endif
         status = nf90_def_var(ncid,'timestep',NF90_INT,timid,varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var timestep')
         status = nf90_def_var(ncid,'date',NF90_DOUBLE,timid,varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: def_var date')

         ! 2d fields

         dimid2(1) = nxid
         dimid2(2) = timid

         do n = 1,num_2d
            status = nf90_def_var(ncid,trim(fld_2d(n)),NF90_DOUBLE,dimid2,varid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR in def_var '//trim(fld_2d(n)))
         enddo

         ! 3d ncat fields

         dimid3(1) = nxid
         dimid3(2) = ncatid
         dimid3(3) = timid

         do n = 1,num_3d_ncat
            status = nf90_def_var(ncid,trim(fld_3d_ncat(n)),NF90_DOUBLE,dimid3,varid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR in def_var '//trim(fld_3d_ncat(n)))
         enddo

         if (tr_fsd) then
            ! 3d nfsd fields

            dimid3(1) = nxid
            dimid3(2) = nfsdid
            dimid3(3) = timid

            do n = 1,num_3d_nfsd
               status = nf90_def_var(ncid,trim(fld_3d_nfsd(n)),NF90_DOUBLE,dimid3,varid)
               if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR in def_var '//trim(fld_3d_nfsd(n)))
            enddo
         endif

         ! 3d ntrcr fields

         dimid3(1) = nxid
         dimid3(2) = ntrcrid
         dimid3(3) = timid

         do n = 1,num_3d_ntrcr
            status = nf90_def_var(ncid,trim(fld_3d_ntrcr(n)),NF90_DOUBLE,dimid3,varid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR in def_var '//trim(fld_3d_ntrcr(n)))
         enddo

         ! 4d ncat ntrcr fields

         dimid4(1) = nxid
         dimid4(2) = ntrcrid
         dimid4(3) = ncatid
         dimid4(4) = timid

         do n = 1,num_4d_ncat_ntrcr
            status = nf90_def_var(ncid,trim(fld_4d_ncat_ntrcr(n)),NF90_DOUBLE,dimid4,varid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR in def_var '//trim(fld_4d_ncat_ntrcr(n)))
         enddo

         ! enddef

         status = nf90_enddef(ncid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR in nf90_enddef')

         ! static dimension variables

         status = nf90_inq_varid(ncid,'ni',varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'ni')
         status = nf90_put_var(ncid,varid,(/(n,n=1,nx)/))
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'ni')

         status = nf90_inq_varid(ncid,'ncat',varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'ncat')
         status = nf90_put_var(ncid,varid,(/(n,n=1,ncat)/))
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'ncat')

         status = nf90_inq_varid(ncid,'ntrcr',varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'ntrcr')
         status = nf90_put_var(ncid,varid,(/(n,n=1,ntrcr)/))
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'ntrcr')

         if (tr_fsd) then
            status = nf90_inq_varid(ncid,'nfsd',varid)
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'nfsd')
            status = nf90_put_var(ncid,varid,(/(n,n=1,nfsd)/))
            if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'nfsd')
         endif

      endif

      first_call = .false.

      ! Time

      timcnt = timcnt + 1

      status = nf90_inq_varid(ncid,'time',varid)
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'time')
      value = (time-time0)/secday
      status = nf90_put_var(ncid,varid,value,start=(/timcnt/))
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'time')

      status = nf90_inq_varid(ncid,'timestep',varid)
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'timestep')
      status = nf90_put_var(ncid,varid,istep1,start=(/timcnt/))
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'timestep')

      status = nf90_inq_varid(ncid,'date',varid)
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//'date')
      value = real(idate,kind=dbl_kind) + real(sec,kind=dbl_kind)/(secday)
      status = nf90_put_var(ncid,varid,value,start=(/timcnt/))
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//'date')

      ! 2d fields

      start2(1) = 1
      count2(1) = nx
      start2(2) = timcnt
      count2(2) = 1

      do n = 1,num_2d
         allocate(value2(count2(1),1))

         value2 = -9999._dbl_kind
         if (trim(fld_2d(n)) == 'aice')     value2(1:count2(1),1) = aice(1:count2(1))
         if (trim(fld_2d(n)) == 'vice')     value2(1:count2(1),1) = vice(1:count2(1))
         if (trim(fld_2d(n)) == 'vsno')     value2(1:count2(1),1) = vsno(1:count2(1))
         if (trim(fld_2d(n)) == 'uvel')     value2(1:count2(1),1) = uvel(1:count2(1))
         if (trim(fld_2d(n)) == 'vvel')     value2(1:count2(1),1) = vvel(1:count2(1))
         if (trim(fld_2d(n)) == 'divu')     value2(1:count2(1),1) = divu(1:count2(1))
         if (trim(fld_2d(n)) == 'shear')    value2(1:count2(1),1) = shear(1:count2(1))
         if (trim(fld_2d(n)) == 'strength') value2(1:count2(1),1) = strength(1:count2(1))
         if (trim(fld_2d(n)) == 'evap')     value2(1:count2(1),1) = evap(1:count2(1))
         if (trim(fld_2d(n)) == 'fsnow')    value2(1:count2(1),1) = fsnow(1:count2(1))
         if (trim(fld_2d(n)) == 'frazil')   value2(1:count2(1),1) = frazil(1:count2(1))
         if (trim(fld_2d(n)) == 'fswabs')   value2(1:count2(1),1) = fswabs(1:count2(1))
         if (trim(fld_2d(n)) == 'flw')      value2(1:count2(1),1) = flw(1:count2(1))
         if (trim(fld_2d(n)) == 'flwout')   value2(1:count2(1),1) = flwout(1:count2(1))
         if (trim(fld_2d(n)) == 'fsens')    value2(1:count2(1),1) = fsens(1:count2(1))
         if (trim(fld_2d(n)) == 'fsurf')    value2(1:count2(1),1) = fsurf(1:count2(1))
         if (trim(fld_2d(n)) == 'flat')     value2(1:count2(1),1) = flat(1:count2(1))
         if (trim(fld_2d(n)) == 'frain')    value2(1:count2(1),1) = frain(1:count2(1))
         if (trim(fld_2d(n)) == 'Tair')     value2(1:count2(1),1) = Tair(1:count2(1))
         if (trim(fld_2d(n)) == 'Qa')       value2(1:count2(1),1) = Qa(1:count2(1))
         if (trim(fld_2d(n)) == 'fsw')      value2(1:count2(1),1) = fsw(1:count2(1))
         if (trim(fld_2d(n)) == 'fcondtop') value2(1:count2(1),1) = fcondtop(1:count2(1))
         if (trim(fld_2d(n)) == 'meltt')    value2(1:count2(1),1) = meltt(1:count2(1))
         if (trim(fld_2d(n)) == 'meltb')    value2(1:count2(1),1) = meltb(1:count2(1))
         if (trim(fld_2d(n)) == 'meltl')    value2(1:count2(1),1) = meltl(1:count2(1))
         if (trim(fld_2d(n)) == 'snoice')   value2(1:count2(1),1) = snoice(1:count2(1))
         if (trim(fld_2d(n)) == 'dsnow')    value2(1:count2(1),1) = dsnow(1:count2(1))
         if (trim(fld_2d(n)) == 'congel')   value2(1:count2(1),1) = congel(1:count2(1))
         if (trim(fld_2d(n)) == 'sst')      value2(1:count2(1),1) = sst(1:count2(1))
         if (trim(fld_2d(n)) == 'sss')      value2(1:count2(1),1) = sss(1:count2(1))
         if (trim(fld_2d(n)) == 'Tf')       value2(1:count2(1),1) = Tf(1:count2(1))
         if (trim(fld_2d(n)) == 'fhocn')    value2(1:count2(1),1) = fhocn(1:count2(1))

         status = nf90_inq_varid(ncid,trim(fld_2d(n)),varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//trim(fld_2d(n)))
         status = nf90_put_var(ncid,varid,value2,start=start2,count=count2)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//trim(fld_2d(n)))

         deallocate(value2)
     enddo

      ! 3d ncat fields

      start3(1) = 1
      count3(1) = nx
      start3(2) = 1
      count3(2) = ncat
      start3(3) = timcnt
      count3(3) = 1

      do n = 1,num_3d_ncat
         allocate(value3(count3(1),count3(2),1))

         value3 = -9999._dbl_kind
         if (trim(fld_3d_ncat(n)) == 'aicen') value3(1:count3(1),1:count3(2),1) = aicen(1:count3(1),1:count3(2))
         if (trim(fld_3d_ncat(n)) == 'vicen') value3(1:count3(1),1:count3(2),1) = vicen(1:count3(1),1:count3(2))
         if (trim(fld_3d_ncat(n)) == 'vsnon') value3(1:count3(1),1:count3(2),1) = vsnon(1:count3(1),1:count3(2))

         status = nf90_inq_varid(ncid,trim(fld_3d_ncat(n)),varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//trim(fld_3d_ncat(n)))
         status = nf90_put_var(ncid,varid,value3,start=start3,count=count3)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//trim(fld_3d_ncat(n)))

         deallocate(value3)
     enddo

     if (tr_fsd) then
        ! 3d nfsd fields

        start3(1) = 1
        count3(1) = nx
        start3(2) = 1
        count3(2) = nfsd
        start3(3) = timcnt
        count3(3) = 1

        do n = 1,num_3d_nfsd
           allocate(value3(count3(1),count3(2),1))

           value3 = -9999._dbl_kind
           if (trim(fld_3d_nfsd(n)) == 'd_afsd_newi') value3(1:count3(1),1:count3(2),1) = d_afsd_newi(1:count3(1),1:count3(2))
           if (trim(fld_3d_nfsd(n)) == 'd_afsd_latg') value3(1:count3(1),1:count3(2),1) = d_afsd_latg(1:count3(1),1:count3(2))
           if (trim(fld_3d_nfsd(n)) == 'd_afsd_latm') value3(1:count3(1),1:count3(2),1) = d_afsd_latm(1:count3(1),1:count3(2))
           if (trim(fld_3d_nfsd(n)) == 'd_afsd_wave') value3(1:count3(1),1:count3(2),1) = d_afsd_wave(1:count3(1),1:count3(2))
           if (trim(fld_3d_nfsd(n)) == 'd_afsd_weld') value3(1:count3(1),1:count3(2),1) = d_afsd_weld(1:count3(1),1:count3(2))

           status = nf90_inq_varid(ncid,trim(fld_3d_nfsd(n)),varid)
           if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//trim(fld_3d_nfsd(n)))
           status = nf90_put_var(ncid,varid,value3,start=start3,count=count3)
           if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//trim(fld_3d_nfsd(n)))

           deallocate(value3)
        enddo
     endif

      ! 3d ntrcr fields

      start3(1) = 1
      count3(1) = nx
      start3(2) = 1
      count3(2) = ntrcr
      start3(3) = timcnt
      count3(3) = 1

      do n = 1,num_3d_ntrcr
         allocate(value3(count3(1),count3(2),1))

         value3 = -9999._dbl_kind
         if (trim(fld_3d_ntrcr(n)) == 'trcr') value3(1:count3(1),1:count3(2),1) = trcr(1:count3(1),1:count3(2))

         status = nf90_inq_varid(ncid,trim(fld_3d_ntrcr(n)),varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//trim(fld_3d_ntrcr(n)))
         status = nf90_put_var(ncid,varid,value3,start=start3,count=count3)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//trim(fld_3d_ntrcr(n)))

         deallocate(value3)
     enddo

      ! 4d ncat ntrcr fields

      start4(1) = 1
      count4(1) = nx
      start4(2) = 1
      count4(2) = ntrcr
      start4(3) = 1
      count4(3) = ncat
      start4(4) = timcnt
      count4(4) = 1

      do n = 1,num_4d_ncat_ntrcr
         allocate(value4(count4(1),count4(2),count4(3),1))

         value4 = -9999._dbl_kind
         if (trim(fld_4d_ncat_ntrcr(n)) == 'trcrn') value4(1:count4(1),1:count4(2),1:count4(3),1) = trcrn(1:count4(1),1:count4(2),1:count4(3))

         status = nf90_inq_varid(ncid,trim(fld_4d_ncat_ntrcr(n)),varid)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: inq_var '//trim(fld_4d_ncat_ntrcr(n)))
         status = nf90_put_var(ncid,varid,value4,start=start4,count=count4)
         if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: put_var '//trim(fld_4d_ncat_ntrcr(n)))

         deallocate(value4)
     enddo

#else
      call icedrv_system_abort(string=subname//' ERROR: history requires USE_NETCDF',file=__FILE__,line=__LINE__)
#endif

      end subroutine history_write

!=======================================================================

! Close history file

      subroutine history_close()

#ifdef USE_NETCDF
      use netcdf
#endif

      ! local variables

      integer (kind=int_kind) :: &
         status                          ! cdf status flag

      character(len=*), parameter :: subname='(history_close)'

#ifdef USE_NETCDF
      status = nf90_close(ncid)
      if (status /= nf90_noerr) call icedrv_system_abort(string=subname//' ERROR: nf90_close')
#else
      call icedrv_system_abort(string=subname//' ERROR: history requires USE_NETCDF',file=__FILE__,line=__LINE__)
#endif

      end subroutine history_close

!=======================================================================

      end module icedrv_history

!=======================================================================
