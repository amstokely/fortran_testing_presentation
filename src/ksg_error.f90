module ksg_error
   use iso_fortran_env, only : real64
   implicit none

   type :: error_t
      integer :: code = 0
      character(len = :), allocatable :: message
   end type error_t

contains

   function make_error_t(code, message) result(err)
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      type(error_t) :: err

      err%code = code
      err%message = message
   end function make_error_t

   subroutine validate_ksg_parameters(Mx, My, k, err)
      real(real64), intent(in) :: Mx(:)
      real(real64), intent(in) :: My(:)
      integer, intent(in) :: k
      type(error_t), intent(out) :: err

      if (size(Mx) /= size(My)) then
         err = make_error_t(-1, "Input arrays Mx and My must have the same size.")
         return
      end if
      if (k < 1) then
         err = make_error_t(-1, "Parameter k must be a positive integer greater than zero.")
         return
      end if
   end subroutine validate_ksg_parameters


   subroutine handle_ksg_status(error_code, err)
      integer, intent(inout) :: error_code
      type(error_t), intent(out) :: err

      if (error_code /= 0) then
         err = make_error_t(error_code, "An error occurred in the ksg_counts implementation.")
      end if

   end subroutine handle_ksg_status

end module ksg_error