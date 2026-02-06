module cpp_mutual_information
   use iso_c_binding
   implicit none
   private
   public :: cpp_ksg_counts

   interface
      subroutine cpp_ksg_counts_i(Mx, My, n_points, k, nx, ny) &
            bind(C, name = "c_cpp_ksg_counts")

         import :: c_double, c_int

         real(c_double), intent(in) :: Mx(*)
         real(c_double), intent(in) :: My(*)
         integer(c_int), value :: n_points
         integer(c_int), value :: k
         integer(c_int), intent(out) :: nx(*)
         integer(c_int), intent(out) :: ny(*)
      end subroutine
   end interface

contains

   subroutine cpp_ksg_counts(Mx, My, k, nx, ny)
      use iso_fortran_env, only : real64

      real(real64), intent(in) :: Mx(:)
      real(real64), intent(in) :: My(:)
      integer, intent(in) :: k
      integer, intent(out) :: nx(size(Mx)), ny(size(Mx))

      integer :: ref_idx
      integer(c_int) :: c_nx(size(Mx)), c_ny(size(Mx))

      call cpp_ksg_counts_i(&
            Mx, &
            My, &
            size(Mx), &
            k, &
            c_nx, c_ny)
      nx = c_nx
      ny = c_ny
   end subroutine

end module
