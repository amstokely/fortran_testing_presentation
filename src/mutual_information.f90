module mutual_information
   use iso_fortran_env, only: real64
   use c_mutual_information
   implicit none
   abstract interface
      subroutine ksg_count_i(J, ref_idx, k, counts)
         import :: real64
         real(real64), intent(in) :: J(:, :)
         integer, intent(in)      :: ref_idx
         integer, intent(in)      :: k
         integer, intent(out)     :: counts(:)
      end subroutine ksg_count_i
   end interface
contains

   subroutine max_norm_from_point(X, xref, dists)
      real(real64), intent(in)  :: X(:, :)
      real(real64), intent(in)  :: xref(:)
      real(real64), intent(out) :: dists(size(X, 1))

      integer :: i, j
      integer :: n_points, n_dims
      real(real64) :: max_diff

      n_points = size(X, 1)
      n_dims   = size(X, 2)

      do i = 1, n_points
         max_diff = 0.0_real64
         do j = 1, n_dims
            max_diff = max(max_diff, abs(X(i, j) - xref(j)))
         end do
         dists(i) = max_diff
      end do
   end subroutine max_norm_from_point


   subroutine k_argsort(X, k, idxs)
      use iso_fortran_env, only: real64
      implicit none

      real(real64), intent(in)  :: X(:)
      integer,      intent(in)  :: k
      integer,      intent(out) :: idxs(k)

      real(real64) :: best_vals(k)
      integer      :: best_idxs(k)

      integer :: i, j, n, pos
      real(real64) :: xval

      n = size(X)

      ! Initialize with +inf
      best_vals = huge(0.0_real64)
      best_idxs = -1

      do i = 1, n
         xval = X(i)

         ! Only consider if better than current worst
         if (xval < best_vals(k)) then

            ! Find insertion position (small k → linear is fastest)
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
   end subroutine


   subroutine max_neighbor_distance(X, xref, idxs, max_dist)
      real(real64), intent(in) :: X(:)
      real(real64), intent(in) :: xref
      integer, intent(in)      :: idxs(:)
      real(real64), intent(out):: max_dist

      integer :: i
      real(real64) :: dist

      max_dist = -huge(0.0_real64)
      do i = 1, size(idxs)
         dist = abs(X(idxs(i)) - xref)
         if (dist > max_dist) then
            max_dist = dist
         end if
      end do
   end subroutine max_neighbor_distance


   subroutine count_neighbors_within_radius(X, xref, radius, count)
      real(real64), intent(in) :: X(:)
      real(real64), intent(in) :: xref
      real(real64), intent(in) :: radius
      integer, intent(out)     :: count

      integer :: i
      count = 0

      do i = 1, size(X)
         if (abs(X(i) - xref) <= radius) then
            count = count + 1
         end if
      end do

      count = count - 1
   end subroutine count_neighbors_within_radius


   subroutine ksg_count(J, ref_idx, k, counts)
      real(real64), intent(in) :: J(:, :)
      integer, intent(in)      :: ref_idx
      integer, intent(in)      :: k
      integer, intent(out)     :: counts(:)

      integer :: n_points
      real(real64), allocatable :: dists(:)
      integer, allocatable :: neighbor_idxs(:)
      real(real64) :: max_dist_x, max_dist_y

      n_points = size(J, 1)
      allocate(dists(n_points))
      allocate(neighbor_idxs(k + 1))

      call max_norm_from_point(J, J(ref_idx, :), dists)
      call k_argsort(dists, k + 1, neighbor_idxs)
      call max_neighbor_distance(J(:, 1), J(ref_idx, 1), neighbor_idxs(2:k + 1), max_dist_x)
      call max_neighbor_distance(J(:, 2), J(ref_idx, 2), neighbor_idxs(2:k + 1), max_dist_y)
      call count_neighbors_within_radius(J(:, 1), J(ref_idx, 1), max_dist_x, counts(1))
      call count_neighbors_within_radius(J(:, 2), J(ref_idx, 2), max_dist_y, counts(2))

      deallocate(dists)
      deallocate(neighbor_idxs)
   end subroutine ksg_count


   subroutine calc_mutual_information(J, k, mi, ksg_count_proc_arg)
      real(real64), intent(in)  :: J(:, :)
      integer, intent(in)       :: k
      real(real64), intent(out) :: mi
      procedure(ksg_count_i), optional :: ksg_count_proc_arg

      integer :: n_points, i
      integer :: counts(2)
      real(real64), allocatable :: psi(:)
      real(real64), parameter :: gamma = -0.5772156649015328606_real64
      real(real64) :: inversePsiIndex, avg_psi_ksg_sum
      procedure(ksg_count_i), pointer :: ksg_count_proc

      if (present(ksg_count_proc_arg)) then
         ksg_count_proc => ksg_count_proc_arg
      else
         ksg_count_proc => ksg_count
      end if

      n_points = size(J, 1)
      allocate(psi(n_points + 1))

      psi(2) = gamma
      do i = 2, n_points + 1
         inversePsiIndex = 1.0_real64 / real(i - 1, real64)
         psi(i) = psi(i - 1) + inversePsiIndex
      end do

      avg_psi_ksg_sum = 0.0_real64
      do i = 1, n_points
         call ksg_count_proc(J, i, k, counts)
         avg_psi_ksg_sum = avg_psi_ksg_sum + &
               (psi(counts(1) + 1) + psi(counts(2) + 1))
      end do

      avg_psi_ksg_sum = avg_psi_ksg_sum / real(n_points, real64)

      mi = psi(k+1) + psi(n_points + 1) - avg_psi_ksg_sum - 1.0_real64 / real(k, real64)

      deallocate(psi)
   end subroutine calc_mutual_information

   subroutine calc_generalized_correlation(J, k, gc, ksg_count_proc_arg)
      real(real64), intent(in)  :: J(:, :)
      integer, intent(in)       :: k
      real(real64), intent(out) :: gc
      procedure(ksg_count_i), optional :: ksg_count_proc_arg

      real(real64) :: mi
      integer :: n_points
      procedure(ksg_count_i), pointer :: ksg_count_proc

      if (present(ksg_count_proc_arg)) then
         ksg_count_proc => ksg_count_proc_arg
      else
         ksg_count_proc => ksg_count
      end if

      n_points = size(J, 1)
      call calc_mutual_information(J, k, mi, ksg_count_proc)
      if (mi <= 0.0_real64) then
         gc = 0.0_real64
         return
      end if
      gc = sqrt(1.0_real64 - exp(-2.0_real64 * mi))
   end subroutine calc_generalized_correlation

end module mutual_information
