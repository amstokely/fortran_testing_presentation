module test_utils
   use iso_c_binding
   implicit none

   interface

      subroutine generate_independent_arrays(a, b, n) &
            bind(C, name="c_generate_independent_arrays")
         import :: c_double, c_int
         real(c_double) :: a(*), b(*)
         integer(c_int), value :: n
      end subroutine


      subroutine correlated_gaussian_arrays(x, y, n_point, rho) &
            bind(C, name="c_correlated_gaussian_arrays")
         import :: c_double, c_int
         real(c_double) :: x(*), y(*)
         integer(c_int), value :: n_point
         real(c_double), value :: rho
      end subroutine

      function pearson_corr(x, y, n) result(r) &
            bind(C, name="c_pearson_corr")
         import :: c_double, c_int
         real(c_double) :: x(*), y(*)
         integer(c_int), value :: n
         real(c_double) :: r
      end function pearson_corr


   subroutine sin_cos_arrays(x, y, n) &
            bind(C, name="c_sin_cos_arrays")
         import :: c_double, c_int
         real(c_double) :: x(*), y(*)
         integer(c_int), value :: n
      end subroutine

   end interface

end module
