!  SVN:$Id: icepack_age.F90 1226 2017-05-22 22:45:03Z tcraig $
!=======================================================================
!
! authors Elizabeth Hunke

      module icepack_age

      use icepack_kinds
      use icepack_warnings, only: warnstr, add_warning
      use icepack_warnings, only: set_warning_abort, icepack_aborted

      implicit none

      private
      public :: increment_age

!=======================================================================

      contains

!=======================================================================

!  Increase ice age tracer by timestep length.

      subroutine increment_age (dt, iage)

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), &
         intent(inout) :: &
         iage

      character(len=*),parameter :: subname='(increment_age)'

      iage = iage + dt 

      end subroutine increment_age

!=======================================================================

      end module icepack_age

!=======================================================================
