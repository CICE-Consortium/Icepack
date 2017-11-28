module icepack_warnings

  use icepack_kinds

  implicit none

  private
  save

  ! warning messages
  character(len=char_len_long), dimension(:), allocatable :: warnings
  integer :: nWarnings = 0

  ! abort flag, accessed via set_warning_abort and icepack_aborted
  logical :: warning_abort = .false.

  public :: &
       icepack_clear_warnings, & 
       icepack_get_warnings, & 
       icepack_print_warnings, &
       icepack_flush_warnings, &
       icepack_aborted, &
       set_warning_abort, &
       add_warning, &
       reset_warnings, &
       get_number_warnings, &
       get_warning

!=======================================================================

contains

!=======================================================================

      logical function icepack_aborted()

        icepack_aborted = warning_abort

      end function icepack_aborted

!=======================================================================

      subroutine set_warning_abort(abortflag)

        logical, intent(in) :: abortflag

        warning_abort = abortflag

      end subroutine set_warning_abort

!=======================================================================

      subroutine icepack_clear_warnings()

        call reset_warnings()

      end subroutine icepack_clear_warnings

!=======================================================================
      
      subroutine icepack_get_warnings(warningsOut)

        character(len=char_len_long), dimension(:), allocatable, intent(out) :: &
             warningsOut

        integer :: &
             iWarning

        if (allocated(warningsOut)) deallocate(warningsOut)
        allocate(warningsOut(get_number_warnings()))

        do iWarning = 1, get_number_warnings()
           warningsOut(iWarning) = trim(get_warning(iWarning))
        enddo

      end subroutine icepack_get_warnings

!=======================================================================

      subroutine icepack_print_warnings(nu_diag)

        integer, intent(in) :: nu_diag

        integer :: &
             iWarning

        do iWarning = 1, get_number_warnings()
           write(nu_diag,*) trim(get_warning(iWarning))
        enddo

      end subroutine icepack_print_warnings

!=======================================================================

      subroutine icepack_flush_warnings(nu_diag)

        integer, intent(in) :: nu_diag

        integer :: &
             iWarning

        call icepack_print_warnings(nu_diag)
        call icepack_clear_warnings()

      end subroutine icepack_flush_warnings

!=======================================================================

  subroutine add_warning(warning)

    character(len=*), intent(in) :: &
         warning ! warning to add to array of warnings

    ! number of array elements to increase size of warnings array if that array has run out of space
    integer, parameter :: &
         nWarningsBuffer = 10

    ! temporary array to store previous warnings while warning array is increased in size
    character(len=char_len_long), dimension(:), allocatable :: &
         warningsTmp

    integer :: &
         nWarningsArray, & ! size of warnings array at start
         iWarning ! warning index

    ! check if warnings array is not allocated
    if (.not. allocated(warnings)) then

       ! allocate warning array with number of buffer elements
       allocate(warnings(nWarningsBuffer))

       ! set initial number of nWarnings
       nWarnings = 0

    ! already allocated
    else

       ! find the size of the warnings array at the start
       nWarningsArray = size(warnings)
       
       ! check to see if need more space in warnings array
       if (nWarnings + 1 > nWarningsArray) then
       
          ! allocate the temporary warning storage
          allocate(warningsTmp(nWarningsArray))

          ! copy the warnings to temporary storage
          do iWarning = 1, nWarningsArray
             warningsTmp(iWarning) = trim(warnings(iWarning))
          enddo ! iWarning

          ! increase the size of the warning array by the buffer size
          deallocate(warnings)
          allocate(warnings(nWarningsArray + nWarningsBuffer))

          ! copy back the temporary stored warnings
          do iWarning = 1, nWarningsArray
             warnings(iWarning) = trim(warningsTmp(iWarning))
          enddo ! iWarning

          ! deallocate the temporary storage
          deallocate(warningsTmp)

       endif
          
    endif

    ! increase warning number
    nWarnings = nWarnings + 1

    ! add the new warning
    warnings(nWarnings) = trim(warning)

  end subroutine add_warning

!=======================================================================

  subroutine reset_warnings()

    nWarnings = 0

  end subroutine reset_warnings

!=======================================================================

  function get_number_warnings() result(nWarningsOut)

    integer :: nWarningsOut

    nWarningsOut = nWarnings

  end function get_number_warnings

!=======================================================================

  function get_warning(iWarning) result(warning)

    integer, intent(in) :: iWarning

    character(len=char_len_long) :: warning

    if (iWarning <= nWarnings) then
       warning = warnings(iWarning)
    else
       warning = ""
    endif

  end function get_warning

!=======================================================================

end module icepack_warnings
