
module icepack_warnings

! Provides a logging and abort package for Icepack.
! Icepack has no idea about MPI, OpenMP, or IO.
! Store error message and provide methods for the driver
! to write these messages to a Fortran unit number.
! Needs to be thread safe.  This could be called within
! a threaded or non-threaded region or both.  Need to make
! sure multiple threads are not adding to the warnings
! buffer at the same time.  Also need to make sure warnings
! buffers are not added at the same time messages are
! cleared by a different thread.  Use multiple critical
! regions using the same ID to allow threads to block
! each other during multiple operations.

      use icepack_kinds
      implicit none

      private

      ! warning messages
      character(len=char_len_long), dimension(:), allocatable :: warnings
      integer :: nWarnings = 0
      integer :: nWarningsBuffer = 10 ! incremental number of messages

      ! abort flag, accessed via icepack_warnings_setabort and icepack_warnings_aborted
      logical :: warning_abort = .false.

      ! public string for all subroutines to use
      character(len=char_len_long), public :: warnstr

      public :: &
        icepack_warnings_clear,    &
        icepack_warnings_print,    &
        icepack_warnings_flush,    &
        icepack_warnings_aborted,  &
        icepack_warnings_add,      &
        icepack_warnings_setabort, &
        icepack_warnings_getall

      private :: &
        icepack_warnings_getone

! variables are shared by default
! have warnstr be private
!$OMP THREADPRIVATE(warnstr)

!=======================================================================

contains

!=======================================================================
!autodocument_start icepack_warnings_aborted
! turn on the abort flag in the icepack warnings package
! pass in an optional error message

      logical function icepack_warnings_aborted(instring)

        character(len=*),intent(in), optional :: instring

!autodocument_end

        character(len=*),parameter :: subname='(icepack_warnings_aborted)'

        icepack_warnings_aborted = warning_abort
        if (warning_abort .and. present(instring)) then
           call icepack_warnings_add(subname//' ... '//trim(instring))
        endif

      end function icepack_warnings_aborted

!=======================================================================

      subroutine icepack_warnings_setabort(abortflag,file,line)

        logical, intent(in) :: abortflag
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        character(len=*),parameter :: subname='(icepack_warnings_setabort)'

        ! try to capture just the first setabort call

        if (abortflag) then
          write(warnstr,*) subname,abortflag
          call icepack_warnings_add(warnstr)
          if (present(file)) then
             write(warnstr,*) trim(warnstr)//' :file '//trim(file)
             call icepack_warnings_add(warnstr)
          endif
          if (present(line)) then
             write(warnstr,*) trim(warnstr)//' :line ',line
             call icepack_warnings_add(warnstr)
          endif
        endif

        warning_abort = abortflag

      end subroutine icepack_warnings_setabort

!=======================================================================
!autodocument_start icepack_warnings_clear
! clear all warning messages from the icepack warning buffer

      subroutine icepack_warnings_clear()

!autodocument_end

        character(len=*),parameter :: subname='(icepack_warnings_clear)'

        nWarnings = 0

      end subroutine icepack_warnings_clear

!=======================================================================
!autodocument_start icepack_warnings_clear
! return an array of all the current warning messages

      subroutine icepack_warnings_getall(warningsOut)

        character(len=char_len_long), dimension(:), allocatable, intent(out) :: &
             warningsOut

!autodocument_end

        integer :: iWarning
        character(len=*),parameter :: subname='(icepack_warnings_getall)'

        if (allocated(warningsOut)) deallocate(warningsOut)
        allocate(warningsOut(nWarnings))

        do iWarning = 1, nWarnings
           warningsOut(iWarning) = trim(icepack_warnings_getone(iWarning))
        enddo

      end subroutine icepack_warnings_getall

!=======================================================================
!autodocument_start icepack_warnings_print
! print all warning messages from the icepack warning buffer

      subroutine icepack_warnings_print(iounit)

        integer, intent(in) :: iounit

!autodocument_end

        integer :: iWarning
        character(len=*),parameter :: subname='(icepack_warnings_print)'

        do iWarning = 1, nWarnings
          write(iounit,*) trim(icepack_warnings_getone(iWarning))
        enddo

      end subroutine icepack_warnings_print

!=======================================================================
!autodocument_start icepack_warnings_flush
! print and clear all warning messages from the icepack warning buffer

      subroutine icepack_warnings_flush(iounit)

        integer, intent(in) :: iounit

!autodocument_end

        character(len=*),parameter :: subname='(icepack_warnings_flush)'

!$OMP CRITICAL (omp_warnings)
        if (nWarnings > 0) then
          call icepack_warnings_print(iounit)
        endif
        call icepack_warnings_clear()
!$OMP END CRITICAL (omp_warnings)

      end subroutine icepack_warnings_flush

!=======================================================================

      subroutine icepack_warnings_add(warning)

        character(len=*), intent(in) :: warning ! warning to add to array of warnings

        ! local

        character(len=char_len_long), dimension(:), allocatable :: warningsTmp
        integer :: &
             nWarningsArray, & ! size of warnings array at start
             iWarning ! warning index
        character(len=*),parameter :: subname='(icepack_warnings_add)'

!$OMP CRITICAL (omp_warnings)
        ! check if warnings array is not allocated
        if (.not. allocated(warnings)) then

           ! allocate warning array with number of buffer elements
           allocate(warnings(nWarningsBuffer))

           ! set initial number of nWarnings
           nWarnings = 0

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

              ! increase nWarningsBuffer for next reallocation
              nWarningsBuffer = nWarningsBuffer * 2
           endif

        endif

        ! increase warning number
        nWarnings = nWarnings + 1

        ! add the new warning
        warnings(nWarnings) = trim(warning)
!$OMP END CRITICAL (omp_warnings)

      end subroutine icepack_warnings_add

!=======================================================================

      function icepack_warnings_getone(iWarning) result(warning)

        integer, intent(in) :: iWarning

        character(len=char_len_long) :: warning

        character(len=*),parameter :: subname='(icepack_warnings_getone)'

        if (iWarning <= nWarnings) then
           warning = warnings(iWarning)
        else
           warning = ""
        endif

      end function icepack_warnings_getone

!=======================================================================

end module icepack_warnings
