module cpp_mutual_information
   use iso_c_binding
   implicit none
   private
   public :: cpp_ksg_count

   interface
      subroutine cpp_ksg_count_i(J, n_points, n_dims, ref_idx, k, counts) &
            bind(C, name="c_ksg_count")

         import :: c_double, c_int

         real(c_double), intent(in)  :: J(*)      ! flat row-major buffer
         integer(c_int), value       :: n_points
         integer(c_int), value       :: n_dims
         integer(c_int), value       :: ref_idx   ! 0-based
         integer(c_int), value       :: k
         integer(c_int), intent(out) :: counts(2)
      end subroutine
   end interface

contains

   subroutine cpp_ksg_count(J, ref_idx, k, counts)
      use iso_fortran_env, only: real64
      real(real64), intent(in)  :: J(:, :)   ! column-major Fortran
      integer, intent(in)       :: ref_idx   ! 1-based Fortran
      integer, intent(in)       :: k
      integer, intent(out)      :: counts(:)

      integer(c_int) :: c_counts(2)
      integer(c_int) :: n_points, n_dims
      real(c_double), allocatable :: J_rowmajor(:)

      n_points = size(J, 1)
      n_dims   = size(J, 2)

      allocate(J_rowmajor(n_points * n_dims))
      J_rowmajor = reshape(transpose(J), [n_points * n_dims])

      call cpp_ksg_count_i( J_rowmajor, &
            n_points,   &
            n_dims,     &
            ref_idx - 1, &  ! convert to 0-based
            k,          &
            c_counts )

      counts = c_counts
      deallocate(J_rowmajor)
   end subroutine
end module
