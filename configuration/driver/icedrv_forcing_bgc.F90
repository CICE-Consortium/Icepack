!=======================================================================
!
! Reads and interpolates forcing data for biogeochemistry
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL
!
      module icedrv_forcing_bgc

      use icedrv_kinds
      use icedrv_domain_size, only: nx
      use icedrv_calendar, only: secday
      use icedrv_constants, only: nu_forcing, nu_diag
      use icepack_intfc, only: icepack_max_algae, icepack_max_doc
      use icepack_intfc, only: icepack_max_dic
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icedrv_system, only: icedrv_system_abort

      implicit none
      private
      public :: get_forcing_bgc, faero_default, init_forcing_bgc 

      real (kind=dbl_kind), dimension(365) :: & ! hardwired for now
         sil_data, nit_data

!=======================================================================

      contains

!=======================================================================

      subroutine init_forcing_bgc

        use icedrv_forcing, only: data_dir, bgc_data_type, bgc_data_file

        integer (kind=int_kind) :: &
           ntime, &
           i

        real (kind=dbl_kind), dimension(365) :: &
           sil, &
           nit

        character (char_len_long) filename

        character(len=*), parameter :: subname='(init_forcing_bgc)'
        
        if (trim(bgc_data_type) == 'ISPOL' .or. &
            trim(bgc_data_type) == 'NICE') then
          
           if (trim(bgc_data_type) == 'ISPOL') &
           filename = trim(data_dir)//'/ISPOL_2004/'//trim(bgc_data_file)
           if (trim(bgc_data_type) == 'NICE') &
           filename = trim(data_dir)//'/NICE_2015/'//trim(bgc_data_file)

          write (nu_diag,*) 'Reading ',filename

          ntime = 365 ! daily

          open (nu_forcing, file=filename, form='formatted')
          read (nu_forcing,*) sil
          read (nu_forcing,*) nit
          close(nu_forcing)

          do i = 1, ntime
             sil_data(i) = sil(i)
             nit_data(i) = nit(i)
          end do

        end if

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
      use icedrv_flux, only: sil, nit
      use icedrv_forcing, only: interp_coeff, bgc_data_type

      integer (kind=int_kind) :: &
         i,            & ! horizontal indices
         ixm,ixx,      & ! record numbers for neighboring months
         maxrec      , & ! maximum record number
         recslot     , & ! spline slot for current record
         recnum      , & ! record number
         dataloc     , & ! = 1 for data located in middle of time interval
                         ! = 2 for date located at end of time interval
         ks              ! bgc tracer index (bio_index_o)

      real (kind=dbl_kind) :: &
          c1intp, c2intp

      logical (kind=log_kind) :: tr_bgc_Sil, tr_bgc_Nit

      character(len=*), parameter :: subname='(get_forcing_bgc)'

      call icepack_query_tracer_flags(tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_Nit_out=tr_bgc_Nit)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      if (trim(bgc_data_type) == 'ISPOL' .or. &
          trim(bgc_data_type) == 'NICE') then

        dataloc = 2                          ! data located at end of interval
        maxrec = 365                         ! 
        
        ! current record number
        recnum = int(yday)   
        
        ! Compute record numbers for surrounding data (2 on each side)
        ixm = mod(recnum+maxrec-2,maxrec) + 1
        ixx = mod(recnum-1,       maxrec) + 1
        
        recslot = 2
        call interp_coeff (recnum, recslot, secday, dataloc, c1intp, c2intp)
                 
        if (tr_bgc_Sil) then
           sil(:) =  c1intp * sil_data(ixm) + c2intp * sil_data(ixx)
        endif

        if (tr_bgc_Nit) then
           nit(:) =  c1intp * nit_data(ixm) + c2intp * nit_data(ixx)
        endif

        do i = 1, nx
           ks = 2*icepack_max_algae + icepack_max_doc + 3 + icepack_max_dic
           ocean_bio_all(i,ks) = sil(i)                       ! Sil
           ks =   icepack_max_algae + 1
           ocean_bio_all(i,ks) = nit(i)                       ! Nit
           ks = 2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
           ocean_bio_all(i,ks) = nit(i)                       ! PON
        enddo

      endif

      end subroutine get_forcing_bgc

!=======================================================================

! constant values for atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_default
        
      use icedrv_flux, only: faero_atm
      character(len=*), parameter :: subname='(faero_default)'
        
      faero_atm(:,1) = 1.e-12_dbl_kind ! kg/m^2 s
      faero_atm(:,2) = 1.e-13_dbl_kind
      faero_atm(:,3) = 1.e-14_dbl_kind
      faero_atm(:,4) = 1.e-14_dbl_kind
      faero_atm(:,5) = 1.e-14_dbl_kind
      faero_atm(:,6) = 1.e-14_dbl_kind
        
      end subroutine faero_default
      
!=======================================================================

      end module icedrv_forcing_bgc

!=======================================================================
