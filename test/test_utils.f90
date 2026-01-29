module test_utils
   use iso_fortran_env, only: real64
   implicit none
contains

   subroutine generate_independent_arrays(a, b, n)
      integer, intent(in) :: n
      real(real64), intent(out) :: a(n), b(n)

      integer :: i, j
      integer, allocatable :: seed(:)
      integer :: seed_size

      call random_seed(size = seed_size)
      allocate(seed(seed_size))

      call system_clock(count = i)
      seed = i + 37 * [(j, j = 1, seed_size)]

      call random_seed(put = seed)
      deallocate(seed)

      do i = 1, n
         call random_number(a(i))
         call random_number(b(i))
      end do
   end subroutine generate_independent_arrays



   subroutine correlated_gaussian_arrays(x, y, n_point, rho)
      integer, intent(in) :: n_point
      real(real64), intent(in) :: rho
      real(real64), intent(out) :: x(n_point), y(n_point)

      real(real64) :: u1(n_point), u2(n_point)
      real(real64) :: r(n_point), theta(n_point)
      real(real64), parameter :: pi  = 3.1415926535897932384626433832795_real64
      real(real64) :: s

      s = sqrt(1.0_real64 - rho*rho)

      call random_number(u1)
      call random_number(u2)

      where (u1 <= 1.0e-12_real64)
         u1 = 1.0e-12_real64
      end where

      r     = sqrt(-2.0_real64 * log(u1))
      theta = 2.0_real64 * pi * u2

      x = r * cos(theta)
      y = r * sin(theta)

      y = rho * x + s * y
   end subroutine correlated_gaussian_arrays

   function pearson_corr(x, y, n) result(r)
      use iso_fortran_env, only: real64
      implicit none

      integer, intent(in) :: n
      real(real64), intent(in) :: x(n), y(n)
      real(real64) :: r

      real(real64) :: mx, my
      real(real64) :: sx, sy, sxy
      integer :: i

      ! Means
      mx = sum(x) / real(n, real64)
      my = sum(y) / real(n, real64)

      ! Accumulate centered sums
      sx  = 0.0_real64
      sy  = 0.0_real64
      sxy = 0.0_real64

      do i = 1, n
         sx  = sx  + (x(i) - mx)**2
         sy  = sy  + (y(i) - my)**2
         sxy = sxy + (x(i) - mx)*(y(i) - my)
      end do

      r = sxy / sqrt(sx * sy)
   end function pearson_corr

   subroutine sin_cos_arrays(x, y, n)
      use iso_fortran_env, only: real64
      implicit none

      integer, intent(in) :: n
      real(real64), intent(out) :: x(n), y(n)

      real(real64) :: theta(n)
      real(real64), parameter :: pi = 3.1415926535897932384626433832795_real64

      ! Uniform angles in [0, 2π)
      call random_number(theta)
      theta = 2.0_real64 * pi * theta

      x = sin(theta)
      y = cos(theta)
   end subroutine sin_cos_arrays



end module test_utils
