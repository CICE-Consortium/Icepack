!=======================================================================
! Copyright (c) 2017, Los Alamos National Security, LLC 
! All rights reserved.
!                
! Copyright 2017. Los Alamos National Security, LLC. This software was 
! produced under U.S. Government contract DE-AC52-06NA25396 for Los 
! Alamos National Laboratory (LANL), which is operated by Los Alamos 
! National Security, LLC for the U.S. Department of Energy. The U.S.  
! Government has rights to use, reproduce, and distribute this software.  
! NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY  
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
      program icepackdriver

      use Icepack_drv_InitMod
      use Icepack_drv_RunMod
      use icepack_drv_constants, only: ice_stdout

      implicit none

      !-----------------------------------------------------------------
      ! Initialize Icepack
      !-----------------------------------------------------------------

      call icepack_initialize

      !-----------------------------------------------------------------
      ! Run Icepack
      !-----------------------------------------------------------------

      call icepack_run

      write(ice_stdout, *) "ICEPACK COMPLETED SUCCESSFULLY "

      end program icepackdriver

!=======================================================================
