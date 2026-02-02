module mutual_information
   use iso_fortran_env, only: real64
   use cpp_mutual_information
   implicit none

   abstract interface
      subroutine ksg_count_i(Mx, My, ref_idx, k, nx, ny)
         import :: real64
         real(real64), intent(in) :: Mx(:), My(:)
         integer, intent(in)      :: ref_idx
         integer, intent(in)      :: k
         integer, intent(out)     :: nx, ny
      end subroutine
   end interface

   abstract interface
      subroutine ksg_counts_i(Mx, My, k, nx, ny)
         import :: real64
         real(real64), intent(in) :: Mx(:), My(:)
         integer, intent(in)      :: k
         integer, intent(out)     :: nx(size(Mx)), ny(size(My))
      end subroutine
   end interface

contains

   ! --------------------------------------------------------------------------
   subroutine max_norm_from_point(Mx, My, xref_x, xref_y, dists)
      real(real64), intent(in)  :: Mx(:), My(:)
      real(real64), intent(in)  :: xref_x, xref_y
      real(real64), intent(out) :: dists(size(Mx))

      integer :: i
      do i = 1, size(Mx)
         dists(i) = max( abs(Mx(i) - xref_x), &
               abs(My(i) - xref_y) )
      end do
   end subroutine max_norm_from_point

   ! --------------------------------------------------------------------------
   subroutine k_argsort(X, k, idxs)
      real(real64), intent(in)  :: X(:)
      integer,      intent(in)  :: k
      integer,      intent(out) :: idxs(k)

      real(real64) :: best_vals(k)
      integer      :: best_idxs(k)

      integer :: i, pos, n
      real(real64) :: xval

      n = size(X)
      best_vals = huge(0.0_real64)
      best_idxs = -1

      do i = 1, n
         xval = X(i)
         if (xval < best_vals(k)) then
            pos = k
            do while (pos > 1 .and. xval < best_vals(pos-1))
               best_vals(pos) = best_vals(pos-1)
               best_idxs(pos) = best_idxs(pos-1)
               pos = pos - 1
            end do
            best_vals(pos) = xval
            best_idxs(pos) = i
         end if
      end do

      idxs = best_idxs
   end subroutine k_argsort

   ! --------------------------------------------------------------------------
   subroutine max_neighbor_distance(X, xref, idxs, max_dist)
      real(real64), intent(in) :: X(:)
      real(real64), intent(in) :: xref
      integer, intent(in)      :: idxs(:)
      real(real64), intent(out):: max_dist

      integer :: i
      max_dist = -huge(0.0_real64)

      do i = 1, size(idxs)
         max_dist = max(max_dist, abs(X(idxs(i)) - xref))
      end do
   end subroutine max_neighbor_distance

   ! --------------------------------------------------------------------------
   subroutine count_neighbors_within_radius(X, xref, radius, count)
      real(real64), intent(in) :: X(:)
      real(real64), intent(in) :: xref
      real(real64), intent(in) :: radius
      integer, intent(out)     :: count

      integer :: i
      count = 0

      do i = 1, size(X)
         if (abs(X(i) - xref) <= radius) count = count + 1
      end do

      count = count - 1
   end subroutine count_neighbors_within_radius

   ! --------------------------------------------------------------------------
   subroutine nf90_ksg_count(Mx, My, ref_idx, k, nx, ny)
      real(real64), intent(in) :: Mx(:), My(:)
      integer, intent(in)      :: ref_idx
      integer, intent(in)      :: k
      integer, intent(out)     :: nx
      integer, intent(out)     :: ny

      integer :: n_points
      real(real64), allocatable :: dists(:)
      integer, allocatable :: neighbor_idxs(:)
      real(real64) :: max_dist_x, max_dist_y

      n_points = size(Mx)
      allocate(dists(n_points))
      allocate(neighbor_idxs(k + 1))

      call max_norm_from_point(Mx, My, Mx(ref_idx), My(ref_idx), dists)
      call k_argsort(dists, k + 1, neighbor_idxs)

      call max_neighbor_distance(Mx, Mx(ref_idx), neighbor_idxs(2:k+1), max_dist_x)
      call max_neighbor_distance(My, My(ref_idx), neighbor_idxs(2:k+1), max_dist_y)

      call count_neighbors_within_radius(Mx, Mx(ref_idx), max_dist_x, nx)
      call count_neighbors_within_radius(My, My(ref_idx), max_dist_y, ny)

      deallocate(dists, neighbor_idxs)
   end subroutine nf90_ksg_count

   subroutine nf90_ksg_counts(Mx, My, k, nx, ny)
      real(real64), intent(in) :: Mx(:), My(:)
      integer, intent(in)      :: k
      integer, intent(out)     :: nx(size(Mx))
      integer, intent(out)     :: ny(size(My))

      integer :: n_points
      integer :: i

      n_points = size(Mx)

      do i = 1, n_points
         call nf90_ksg_count(Mx, My, i, k, nx(i), ny(i))
      end do
   end subroutine nf90_ksg_counts

   ! --------------------------------------------------------------------------
   subroutine calc_mutual_information(Mx, My, k, mi, ksg_counts_strategy)
      real(real64), intent(in)  :: Mx(:), My(:)
      integer, intent(in)       :: k
      real(real64), intent(out) :: mi
      procedure(ksg_counts_i), optional :: ksg_counts_strategy

      integer :: n_points, i
      integer :: nx(size(Mx)), ny(size(My))
      real(real64), allocatable :: psi(:)
      real(real64), parameter :: gamma = -0.5772156649015328606_real64
      real(real64) :: avg_psi_ksg_sum
      procedure(ksg_counts_i), pointer :: ksg_counts

      if (present(ksg_counts_strategy)) then
         ksg_counts => ksg_counts_strategy
      else
         ksg_counts => nf90_ksg_counts
      end if

      n_points = size(Mx)
      allocate(psi(n_points + 1))

      psi(2) = gamma
      do i = 2, n_points + 1
         psi(i) = psi(i - 1) + 1.0_real64 / real(i - 1, real64)
      end do

      avg_psi_ksg_sum = 0.0_real64
      call ksg_counts(Mx, My, k, nx, ny)
      do i = 1, n_points
         avg_psi_ksg_sum = avg_psi_ksg_sum + &
               (psi(nx(i) + 1) + psi(ny(i) + 1))
      end do

      avg_psi_ksg_sum = avg_psi_ksg_sum / real(n_points, real64)

      mi = psi(k+1) + psi(n_points + 1) - avg_psi_ksg_sum - 1.0_real64 / real(k, real64)

      deallocate(psi)
   end subroutine calc_mutual_information

   ! --------------------------------------------------------------------------
   subroutine calc_generalized_correlation(Mx, My, k, gc, ksg_counts_strategy)
      real(real64), intent(in)  :: Mx(:), My(:)
      integer, intent(in)       :: k
      real(real64), intent(out) :: gc
      procedure(ksg_counts_i), optional :: ksg_counts_strategy

      real(real64) :: mi
      procedure(ksg_counts_i), pointer :: ksg_counts

      if (present(ksg_counts_strategy)) then
         ksg_counts => ksg_counts_strategy
      else
         ksg_counts => nf90_ksg_counts
      end if

      call calc_mutual_information(Mx, My, k, mi, ksg_counts)

      if (mi <= 0.0_real64) then
         gc = 0.0_real64
      else
         gc = sqrt(1.0_real64 - exp(-2.0_real64 * mi))
      end if
   end subroutine calc_generalized_correlation

end module mutual_information

