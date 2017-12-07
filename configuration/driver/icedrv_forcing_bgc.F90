!  SVN:$Id: ice_forcing.F90 973 2015-04-15 21:07:21Z akt $
!=======================================================================
!
! Reads and interpolates forcing data for biogeochemistry
!
! authors:  Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL
!
      module icedrv_forcing_bgc

      use icedrv_kinds
      use icedrv_domain_size, only: nx
      use icedrv_tracers, only: bio_index_o
      use icedrv_calendar, only: dt, istep, sec, mday, month, daymo
      use icedrv_constants, only: nu_diag

      implicit none
      private
      public :: get_forcing_bgc, faero_default, init_bgc_data, init_forcing_bgc 
      !cn, get_atm_bgc, fzaero_data, &
                !cn , faero_data, faero_optics

      real (kind=dbl_kind), dimension(365) :: & ! hardwired for now
          sil_data, nit_data

      integer (kind=int_kind) :: &
         bgcrecnum = 0   ! old record number (save between steps)

!=======================================================================

      contains

!=======================================================================

      subroutine init_forcing_bgc

        use icedrv_constants, only: nu_forcing
        use icedrv_arrays_column, only: nit_data_type, sil_data_type
        use icedrv_forcing, only: data_dir 

        integer (kind=int_kind) :: &
            ntime, &
            i

        real(kind=dbl_kind), dimension(365) :: &
            sil, &
            nit

        character (char_len_long) filename

        character(len=*), parameter :: subname='(init_forcing_bgc)'
        
        if (trim(nit_data_type) == 'ISPOL' .or. &
            trim(sil_data_type) == 'ISPOL' .or. &
            trim(nit_data_type) == 'NICE' .or. &
            trim(sil_data_type) == 'NICE') then 
          
           if (trim(nit_data_type) == 'ISPOL' .or. &
               trim(sil_data_type) == 'ISPOL') &
               filename = trim(data_dir)//'/ISPOL_2004/nutrients_daily_ISPOL_WOA_field3.txt'
           if (trim(nit_data_type) == 'NICE' .or. &
               trim(sil_data_type) == 'NICE') &
               filename = trim(data_dir)//'/NICE_2015/nutrients_daily_ISPOL_WOA_field3.txt'

          write (nu_diag,*) 'Reading ',filename

          ntime = 365 !daily

          open (nu_forcing, file=filename, form='formatted')
          read (nu_forcing,*) sil
          read (nu_forcing,*) nit
          close(nu_forcing)

          do i = 1, ntime
            sil_data(i) = sil(i)
            nit_data(i) = nit(i)
          end do

          !write(*,*)sil_data
          !write(*,*)nit_data
          
        end if
      !stop  
      end subroutine init_forcing_bgc
        
!=======================================================================
!
! Read and interpolate annual climatologies of silicate and nitrate.
! Restore model quantities to data if desired.
!
! author: Elizabeth C. Hunke, LANL

      subroutine get_forcing_bgc

      use icedrv_arrays_column, only: ocean_bio_all
      use icedrv_calendar, only:  yday
      use icedrv_constants, only: secday
      use icedrv_flux, only: sss, sil, nit
      use icedrv_forcing, only: interp_coeff
      use icedrv_arrays_column, only: nit_data_type, sil_data_type, bgc_data_dir
      use icedrv_tracers, only: max_algae, max_doc, max_dic
      use icedrv_tracers, only: tr_bgc_Sil, tr_bgc_Nit

      integer (kind=int_kind) :: &
          i, j, k,iblk, & ! horizontal indices !cn remove
          ixm,ixp, ixx, & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth    , & ! middle day of month
          recnum      , & ! record number
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          sec_day     , & !  fix time to noon
          ks              ! bgc tracer index (bio_index_o)

      character (char_len_long) :: & 
         met_file,   &    ! netcdf filename
         fieldname        ! field name in netcdf file

      real (kind=dbl_kind) :: &
          sec1hr,&              ! number of seconds in 1 hour
          c1intp, c2intp

      logical (kind=log_kind) :: readm, read1

      character (char_len_long) :: &        ! input data file names
         nit_file   , & ! nitrate input file
         sil_file       ! silicate input file

      character(len=*), parameter :: subname='(get_forcing_bgc)'

      !cn write(*,*) nit_data_type


      if (.not. trim(nit_data_type)=='ISPOL' .AND. &
          .not. trim(sil_data_type)=='ISPOL' .AND. &
          .not. trim(nit_data_type)=='INICE' .AND. &
          .not. trim(sil_data_type)=='NICE') then


      elseif (trim(nit_data_type) == 'ISPOL' .or. trim(sil_data_type) == 'ISPOL') then 
      
        dataloc = 2                          ! data located at end of interval
        sec1hr = secday                      ! seconds in day
        maxrec = 365                         ! 
        
        ! current record number
        recnum = int(yday)   
        
        ! Compute record numbers for surrounding data (2 on each side)
        ixm = mod(recnum+maxrec-2,maxrec) + 1
        ixx = mod(recnum-1,       maxrec) + 1
        
        recslot = 2
        ixp = -99
        call interp_coeff (recnum, recslot, sec1hr, dataloc, c1intp, c2intp)
        read1 = .false.
        if (istep==1 .or. bgcrecnum .ne. recnum) read1 = .true.
                 
        if (tr_bgc_Sil) then
   
          sil(:) =  c1intp * sil_data(ixm) &
              + c2intp * sil_data(ixx)
        endif

        if (tr_bgc_Nit) then
                   
          nit(:) =  c1intp * nit_data(ixm) &
              + c2intp * nit_data(ixx)
        endif

        do i = 1, nx
          ks = 2*max_algae + max_doc + 3 + max_dic
          ocean_bio_all(i,ks) = sil(i)                       !Sil  
          ks = max_algae + 1
          ocean_bio_all(i,ks) = nit(i)                       !nit
          ks =  2*max_algae + max_doc + 7 + max_dic
          ocean_bio_all(i,ks) = nit(i)                       !PON    
        enddo

      endif !ISPOL
      end subroutine get_forcing_bgc

!=======================================================================

! constant values for atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_default
        
        use icedrv_flux, only: faero_atm
        !use icepack_constants, only: nspint
        !use icepack_intfc_shared, only: max_aero
        !use icepack_intfc_tracers, only: tr_aero
        character(len=*), parameter :: subname='(faero_default)'
        
        faero_atm(:,1) = 1.e-12_dbl_kind ! kg/m^2 s
        faero_atm(:,2) = 1.e-13_dbl_kind
        faero_atm(:,3) = 1.e-14_dbl_kind 
        faero_atm(:,4) = 1.e-14_dbl_kind 
        faero_atm(:,5) = 1.e-14_dbl_kind 
        faero_atm(:,6) = 1.e-14_dbl_kind 
        
      end subroutine faero_default
      
!=======================================================================

! Initialize ocean iron from file
!
! authors: Nicole Jeffery, LANL
      subroutine init_bgc_data (fed1,fep1)
      !cn use ice_read_write, only: ice_open_nc, ice_read_nc, ice_close_nc
      use icedrv_constants, only: c0, p1 !, nu_forcing
      use icedrv_tracers, only: max_fe
      use icedrv_arrays_column, only: fe_data_type, bgc_data_dir

#ifdef ncdf
      use netcdf
#endif
           
      !real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(out) :: &
      real (kind=dbl_kind), dimension(nx, max_fe), intent(out) :: &
           fed1, &  ! first dissolved iron pool (nM)
           fep1    ! first particulate iron pool (nM)

      ! local parameters

      integer (kind=int_kind) :: &
         fid              , & ! file id for netCDF file 
         nbits

      logical (kind=log_kind) :: diag

      character (char_len_long) :: & 
         iron_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      character(len=*), parameter :: subname='(init_bgc_data)'

      nbits = 64              ! double precision data

#if 0
    !-------------------------------------------------------------------
    ! Annual average data from Tagliabue, 2012 (top 50 m average
    ! poisson grid filled on gx1v6
    !-------------------------------------------------------------------

      if (trim(fe_data_type) == 'clim') then
       	diag = .true.   ! write diagnostic information 
        iron_file = trim(bgc_data_dir)//'/dFe_50m_annual_Tagliabue_gx1.nc'

        write (nu_diag,*) ' '
        write (nu_diag,*) 'Dissolved iron ocean concentrations from:'
        write (nu_diag,*) trim(iron_file)
        !cn call ice_open_nc(iron_file,fid)
        open(fid,file=iron_file,form='unformatted')

        !cn fieldname='dFe'
        !cn Not sure what this comment means, do we not need an implied do ver nx, max_fe within the read?
        ! Currently only first fed  value is read
        !cn call ice_read_nc(fid,1,fieldname,fed1,diag)
        read(fid,fed1(:,:))

        where ( fed1(:,:) > 1.e20) fed1(:,:) = p1  

        !cn call ice_close_nc(fid)  
        close(fid)

       	diag = .true.   ! write diagnostic information 
        iron_file = trim(bgc_data_dir)//'/pFe_bathy_gx1.nc'

        write (nu_diag,*) ' '
        write (nu_diag,*) 'Particulate iron ocean concentrations from:'
        write (nu_diag,*) trim(iron_file)
        !cn call ice_open_nc(iron_file,fid)
        open(fid, file = iron_file,form='unformatted')

        !cn fieldname='pFe'
        !cn Not sure what this comment means, do we not need an implied do within the read?
        ! Currently only first fep value is read
        !cn call ice_read_nc(fid,1,fieldname,fep1,diag) 
        read(fid,fep1(:,:)) 
        where ( fep1(:,:) > 1.e20) fep1(:,:) = p1  

        close(fid)  

      endif
#endif
      end subroutine init_bgc_data

!=======================================================================

      end module icedrv_forcing_bgc

!=======================================================================
