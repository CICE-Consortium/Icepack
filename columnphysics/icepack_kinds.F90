!=======================================================================

! Defines variable precision for all common data types
! Code originally based on kinds_mod.F in POP
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
! 2006: ECH converted to free source form (F90)

      module icepack_kinds

!=======================================================================

      implicit none
      public

      integer, parameter :: char_len  = 80, &
                            char_len_long  = 256, &
                            log_kind  = kind(.true.), &
                            int_kind  = selected_int_kind(6), &
                            int8_kind = selected_int_kind(13), &
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13), &
                            r16_kind  = selected_real_kind(33,4931)

!=======================================================================

      end module icepack_kinds

!=======================================================================
