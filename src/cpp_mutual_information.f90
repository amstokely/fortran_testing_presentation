module cpp_mutual_information
   use iso_c_binding
   implicit none
   private
   public :: cpp_ksg_counts

   interface
      subroutine cpp_ksg_counts_i(Mx, My, n_points, k, mx_counts, my_counts) &
            bind(C, name = "c_cpp_ksg_counts")

         import :: c_double, c_int

         real(c_double), intent(in) :: Mx(*)
         real(c_double), intent(in) :: My(*)
         integer(c_int), value :: n_points
         integer(c_int), value :: k
         integer(c_int), intent(out) :: mx_counts(*)
         integer(c_int), intent(out) :: my_counts(*)
      end subroutine
   end interface

contains

   subroutine cpp_ksg_counts(Mx, My, k, mx_counts, my_counts)
      use iso_fortran_env, only : real64

      real(real64), intent(in) :: Mx(:)
      real(real64), intent(in) :: My(:)
      integer, intent(in) :: k
      integer, intent(out) :: mx_counts(size(Mx)), my_counts(size(Mx))

      integer(c_int) :: c_mx_counts(size(Mx)), c_my_counts(size(Mx))

      call cpp_ksg_counts_i(&
            Mx, &
            My, &
            size(Mx), &
            k, &
            c_mx_counts, c_my_counts)
      mx_counts = c_mx_counts
      my_counts = c_my_counts
   end subroutine

end module
