!=======================================================================

      module icedrv_restart_shared

      use icedrv_kinds
      implicit none
      private
      public :: lenstr

      logical (kind=log_kind), public :: &
         restart          ! if true, initialize using restart file instead of defaults

      character (len=char_len_long), public :: &
         restart_file  , & ! output file for restart dump
         restart_dir       ! directory name for restart dump

!=======================================================================

      contains

!=======================================================================

! Compute length of string by finding first non-blank
! character from the right.

      integer function lenstr(label)

      character*(*) label

      ! local variables

      integer (kind=int_kind) :: &
         length, & ! length of character string
         n         ! loop index

      character(len=*), parameter :: subname='(lenstr)'

      length = len(label)
      do n=length, 1, -1
         if( label(n:n) /= ' ' ) exit
      enddo
      lenstr = n

      end function lenstr

!=======================================================================

      end module icedrv_restart_shared

!=======================================================================
