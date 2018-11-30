!=======================================================================
! Copyright (c) 2018, Triad National Security, LLC 
! All rights reserved.
!                
! Copyright 2018. Triad National Security, LLC. This software was 
! produced under U.S. Government contract DE-AC52-06NA25396 for Los 
! Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S. Department of Energy. The U.S.  
! Government has rights to use, reproduce, and distribute this software.  
! NEITHER THE GOVERNMENT NOR TRIAD NATIONAL SECURITY, LLC MAKES ANY  
! WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
! THIS SOFTWARE. If software is modified to produce derivative works, 
! such modified software should be clearly marked, so as not to confuse 
! it with the version available from LANL.
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
      use icedrv_constants, only: ice_stdout, nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icedrv_system, only: icedrv_system_abort

      implicit none

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

      end program icedrv

!=======================================================================
