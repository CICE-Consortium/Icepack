!=======================================================================
! Copyright 1998-2026, Triad National Security, LLC
! All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others
! acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! The full license and distribution policy are available from
! https://github.com/CICE-Consortium
!
!=======================================================================

! Main driver routine for Icepack, the column package for CICE.
! Initializes and steps through the model.
!
! author Elizabeth C. Hunke, LANL
!
      program icedrv

      use icedrv_InitMod
      use icedrv_RunMod
      use icedrv_constants, only: ice_stdout, nu_diag, nu_diag_out
      use icedrv_domain_size, only: nx
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icedrv_system, only: icedrv_system_abort, icedrv_system_flush

      implicit none

      integer n
      logical openflag
      character(len=*), parameter :: subname='(icedrv)'

      !-----------------------------------------------------------------
      ! Initialize Icepack
      !-----------------------------------------------------------------

      call icedrv_initialize

      !-----------------------------------------------------------------
      ! Run Icepack
      !-----------------------------------------------------------------

      call icedrv_run

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(ice_stdout, *) "ICEPACK COMPLETED SUCCESSFULLY "

      inquire(unit=ice_stdout,opened=openflag)
      if (openflag) then
         call icedrv_system_flush(ice_stdout)
         close (ice_stdout)
      endif

      inquire(unit=nu_diag,opened=openflag)
      if (openflag) then
         call icedrv_system_flush(nu_diag)
         close (nu_diag)
      endif

      do n = 1, nx
         inquire(unit=nu_diag_out+n-1,opened=openflag)
         if (openflag) then
            call icedrv_system_flush(nu_diag_out+n-1)
            close (nu_diag_out+n-1)
         endif
      enddo

      end program icedrv

!=======================================================================
