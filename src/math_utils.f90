module math_utils
   use iso_fortran_env, only : real64
   use ksg_error
   implicit none
   private

   public digamma_t

   interface digamma_t
      module procedure make_digamma
   end interface

   real(real64), parameter :: euler_gamma = 0.5772156649015328606_real64

   !==================================================================
   ! Digamma lookup object
   !==================================================================
   type, private :: digamma_t
      private
      real(real64), allocatable :: table(:)
   contains
      procedure :: value => get_value
   end type

contains

   !------------------------------------------------------------------
   ! Factory
   !------------------------------------------------------------------
   function make_digamma(n) result(obj)
      integer, intent(in) :: n
      type(digamma_t) :: obj

      integer :: i

      allocate(obj%table(n+1))
      obj%table = 0.0_real64

      ! psi(1) unused
      obj%table(2) = euler_gamma

      do i = 2, n+1
         obj%table(i) = obj%table(i-1) + 1.0_real64 / real(i-1, real64)
      end do
   end function make_digamma


   !------------------------------------------------------------------
   ! Accessor
   !------------------------------------------------------------------
   function get_value(self, k, err) result(val)
      class(digamma_t), intent(in) :: self
      integer, intent(in) :: k
      type(error_t), intent(out) :: err
      real(real64) :: val

      if (k == 0) then
         err = make_error_t(-1, "digamma is undefined at zero")
         val = 0.0_real64
         return
      end if

      if (k+1 > size(self%table)) then
         err = make_error_t(-1, "Requested digamma value exceeds precomputed table size")
         val = 0.0_real64
         return
      end if

      val = self%table(k+1)
   end function get_value

end module math_utils
