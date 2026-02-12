module test_utils
   use iso_c_binding
   use ksg_error
   implicit none

   interface

      subroutine generate_independent_arrays(a, b, n) &
            bind(C, name = "c_generate_independent_arrays")
         import :: c_double, c_int
         real(c_double) :: a(*), b(*)
         integer(c_int), value :: n
      end subroutine


      subroutine correlated_gaussian_arrays(x, y, n_point, rho) &
            bind(C, name = "c_correlated_gaussian_arrays")
         import :: c_double, c_int
         real(c_double) :: x(*), y(*)
         integer(c_int), value :: n_point
         real(c_double), value :: rho
      end subroutine

      function pearson_corr(x, y, n) result(r) &
            bind(C, name = "c_pearson_corr")
         import :: c_double, c_int
         real(c_double) :: x(*), y(*)
         integer(c_int), value :: n
         real(c_double) :: r
      end function pearson_corr


      subroutine sin_cos_arrays(x, y, n) &
            bind(C, name = "c_sin_cos_arrays")
         import :: c_double, c_int
         real(c_double) :: x(*), y(*)
         integer(c_int), value :: n
      end subroutine

      subroutine cpp_throws_ksg_counts_i(Mx, My, n_points, k, mx_counts, my_counts, err) &
            bind(C, name = "c_cpp_throws_ksg_counts")
         import :: c_double, c_int
         real(c_double) :: Mx(*), My(*)
         integer(c_int), value :: k
         integer(c_int), value :: n_points
         integer(c_int) :: mx_counts(*), my_counts(*)
         integer(c_int), intent(out) :: err
      end subroutine
   end interface

contains

   subroutine cpp_throws_ksg_counts(Mx, My, k, mx_counts, my_counts, err)
      use iso_fortran_env, only : real64
      use ksg_error
      implicit none
      real(real64), intent(in) :: Mx(:)
      real(real64), intent(in) :: My(:)
      integer, intent(in) :: k
      integer, intent(out) :: mx_counts(size(Mx)), my_counts(size(Mx))
      type(error_t), intent(out) :: err

      integer(c_int) :: c_mx_counts(size(Mx)), c_my_counts(size(Mx))
      integer(c_int) :: c_err

      call cpp_throws_ksg_counts_i(&
            Mx, &
            My, &
            size(Mx), &
            k, &
            c_mx_counts, c_my_counts, c_err)
      call handle_ksg_status(c_err, err)
      if (err%code /= 0) then
         return
      end if
      mx_counts = c_mx_counts
      my_counts = c_my_counts
   end subroutine cpp_throws_ksg_counts

end module
