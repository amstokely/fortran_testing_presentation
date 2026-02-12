---
marp: true
theme: gaia
paginate: true
---

```fortran
 @TestCase
   type, extends(TestCase) :: ksg_count_fixture
      real(real64) :: Mx(8)
      real(real64) :: My(8)
   contains
      procedure :: setup
   end type ksg_count_fixture

contains

   subroutine setup(this)
      class(ksg_count_fixture), intent(inout) :: this

      this%Mx = [1.0_real64, 4.0_real64, &
            6.0_real64, 3.0_real64, &
            5.0_real64, 2.0_real64, &
            7.0_real64, 8.0_real64]

      this%My = [5.0_real64, 6.0_real64, &
            4.0_real64, 8.0_real64, &
            1.0_real64, 7.0_real64, &
            2.0_real64, 3.0_real64]
   end subroutine setup
```

---

```fortran
 @Test
   subroutine test_max_norm_from_point(this)
      class(ksg_count_fixture), intent(inout) :: this

      real(real64), parameter :: expected_dists(8) = &
            [0.0_real64, 3.0_real64, 5.0_real64, 3.0_real64, &
                  4.0_real64, 2.0_real64, 6.0_real64, 7.0_real64]

      real(real64) :: dists(8)

      call max_norm_from_point(this%Mx, this%My, &
            this%Mx(1), this%My(1), dists)

      call assertEqual( &
            expected_dists, dists, &
            tolerance=1.0e-12_real64)
   end subroutine
```

---

```fortran
  @Test
   subroutine test_k_argsort(this)
      class(ksg_count_fixture), intent(inout) :: this

      real(real64), parameter :: X(8) = &
            [0.0_real64, 3.0_real64, 5.0_real64, 3.0_real64, &
                  4.0_real64, 2.0_real64, 6.0_real64, 7.0_real64]

      integer, parameter :: expected_idxs(4) = [1, 6, 2, 4]
      integer :: idxs(4)

      call k_argsort(X, 4, idxs)

      call assertEqual( &
            expected_idxs, idxs)
   end subroutine
```

---

```fortran
subroutine k_argsort(X, k, idxs)
   real(real64), intent(in)  :: X(:)
   integer,      intent(in)  :: k
   integer,      intent(out) :: idxs(k)

   integer :: i, j, min_idx
   real(real64) :: min_val
   real(real64) :: X_copy(size(X))

   X_copy = X
   do i = 1, k
      min_val = huge(0.0_real64)
      min_idx = -1
      do j = 1, size(X)
         if (X_copy(j) < min_val) then
            min_val = X_copy(j)
            min_idx = j
         end if
      end do
      X_copy(min_idx) = huge(0.0_real64)
      idxs(i) = min_idx
   end do
end subroutine k_argsort
```

---

```fortran
subroutine k_argsort(X, k, idxs)
   real(real64), intent(in) :: X(:)
   integer, intent(in) :: k
   integer, intent(out) :: idxs(k)

   real(real64) :: best_vals(k)
   integer :: best_idxs(k)

   integer :: i, pos, n
   real(real64) :: xval

   n = size(X)
   best_vals = huge(0.0_real64)
   best_idxs = -1

   do i = 1, n
      xval = X(i)
      if (xval < best_vals(k)) then
         pos = k
         do while (pos > 1 .and. xval < best_vals(pos - 1))
            best_vals(pos) = best_vals(pos - 1)
            best_idxs(pos) = best_idxs(pos - 1)
            pos = pos - 1
         end do
         best_vals(pos) = xval
         best_idxs(pos) = i
      end if
   end do

   idxs = best_idxs
end subroutine k_argsort
```

---

```fortran
  @Test
   subroutine test_max_neighbor_distance(this)
      class(ksg_count_fixture), intent(inout) :: this

      integer, parameter :: neighbors(3) = [6, 2, 4]
      real(real64) :: max_dist

      call max_neighbor_distance(this%My, 5.0_real64, neighbors, max_dist)

      call assertEqual( &
            3.0_real64, max_dist, &
            tolerance=1.0e-12_real64)
   end subroutine
```

---

```fortran
   subroutine max_neighbor_distance(X, xref, idxs, max_dist)
      real(real64), intent(in) :: X(:)
      real(real64), intent(in) :: xref
      integer, intent(in) :: idxs(:)
      real(real64), intent(out) :: max_dist

      integer :: i
      max_dist = -huge(0.0_real64)

      do i = 1, size(idxs)
         max_dist = max(max_dist, abs(X(idxs(i)) - xref))
      end do
   end subroutine max_neighbor_distance
```

---


```fortran
  @Test
   subroutine test_count_neighbors_within_radius(this)
      class(ksg_count_fixture), intent(inout) :: this
      integer :: count

      call count_neighbors_within_radius(this%My, 5.0_real64, 3.0_real64, count)

      call assertEqual( &
            6, count)
   end subroutine
```

---

```fortran
   subroutine count_neighbors_within_radius(X, xref, radius, count)
      real(real64), intent(in) :: X(:)
      real(real64), intent(in) :: xref
      real(real64), intent(in) :: radius
      integer, intent(out) :: count

      integer :: i
      count = 0

      do i = 1, size(X)
         if (abs(X(i) - xref) <= radius) count = count + 1
      end do

      count = count - 1
   end subroutine count_neighbors_within_radius
```

---

```fortran
 @Test
   subroutine test_complete_ksg_count(this)
      class(ksg_count_fixture), intent(inout) :: this

      integer :: f90_mx_counts, f90_my_counts

      call f90_ksg_count(this%Mx, this%My, 1, 3, f90_mx_counts, f90_my_counts)

      call assertEqual( &
            [3, 6], [f90_mx_counts, f90_my_counts])
   end subroutine
```

---

```fortran
  subroutine f90_ksg_count(Mx, My, ref_idx, k, mx_counts, my_counts)
      real(real64), intent(in) :: Mx(:), My(:)
      integer, intent(in) :: ref_idx
      integer, intent(in) :: k
      integer, intent(out) :: mx_counts
      integer, intent(out) :: my_counts

      integer :: n_points
      real(real64), allocatable :: dists(:)
      integer, allocatable :: neighbor_idxs(:)
      real(real64) :: max_dist_x, max_dist_y

      n_points = size(Mx)
      allocate(dists(n_points))
      allocate(neighbor_idxs(k + 1))

      call max_norm_from_point(Mx, My, Mx(ref_idx), My(ref_idx), dists)
      call k_argsort(dists, k + 1, neighbor_idxs)

      call max_neighbor_distance(Mx, Mx(ref_idx), neighbor_idxs(2:k + 1), max_dist_x)
      call max_neighbor_distance(My, My(ref_idx), neighbor_idxs(2:k + 1), max_dist_y)

      call count_neighbors_within_radius(Mx, Mx(ref_idx), max_dist_x, mx_counts)
      call count_neighbors_within_radius(My, My(ref_idx), max_dist_y, my_counts)

      deallocate(dists, neighbor_idxs)
   end subroutine f90_ksg_count
```

---

```fortran
  @Test
   subroutine test_ksg_counts_complete(this)
      class(ksg_count_fixture), intent(inout) :: this

      integer :: mx_counts(8), my_counts(8)

      call f90_ksg_counts(this%Mx, this%My, 3, mx_counts, my_counts)

      call assertEqual( &
            [3, 6], [mx_counts(1), my_counts(1)])
   end subroutine
```

---

```fortran
  subroutine f90_ksg_counts(Mx, My, k, mx_counts, my_counts)
      real(real64), intent(in) :: Mx(:), My(:)
      integer, intent(in) :: k
      integer, intent(out) :: mx_counts(size(Mx))
      integer, intent(out) :: my_counts(size(My))

      integer :: n_points
      integer :: i

      n_points = size(Mx)

      do i = 1, n_points
         call f90_ksg_count(Mx, My, i, k, mx_counts(i), my_counts(i))
      end do
   end subroutine f90_ksg_counts
```

---

```fortran
 @testCase
   type, extends(TestCase) :: mutual_information_fixture
      integer :: n_points = 1024
      real(real64), allocatable :: Mx(:), My(:)
   contains
      procedure :: setup
      procedure :: teardown
   end type mutual_information_fixture

contains
   
   subroutine setup(this)
      class(mutual_information_fixture), intent(inout) :: this
      allocate(this%Mx(this%n_points))
      allocate(this%My(this%n_points))
   end subroutine

   subroutine teardown(this)
      class(mutual_information_fixture), intent(inout) :: this
      deallocate(this%Mx)
      deallocate(this%My)
   end subroutine
```

---

```fortran
   @test
   subroutine test_mutual_information_independent(this)
      class(mutual_information_fixture), intent(inout) :: this
      real(real64) :: mi

      call generate_independent_arrays(this%Mx, this%My, this%n_points)

      call calc_mutual_information(&
            this%Mx, this%My, this%n_points / 2, mi)

      @assertEqual(0.0_real64, mi, tolerance = 0.01_real64)
   end subroutine
```

---

```fortran
  @test
   subroutine test_mutual_information_analytic_gaussian(this)
      class(mutual_information_fixture), intent(inout) :: this
      integer, parameter :: n_big = 4096
      integer :: k
      real(real64) :: Mx(n_big), My(n_big)
      real(real64) :: mi

      call correlated_gaussian_arrays(Mx, My, n_big, 0.9_real64)

      k = int(n_big * 0.04_real64 + 0.5_real64)

      call calc_mutual_information(Mx, My, k, mi)

      @assertEqual(0.830366_real64, mi, tolerance = 0.05_real64)
   end subroutine
```

---

```fortran
   @test
   subroutine test_analytic_gaussian_ratio(this)
      class(mutual_information_fixture), intent(inout) :: this
      real(real64) :: mi_02, mi_04, mi_06
      integer :: k_02, k_04, k_06

      call correlated_gaussian_arrays(this%Mx, this%My, this%n_points, 0.9_real64)

      k_02 = int(this%n_points * 0.02_real64 + 0.5_real64)
      k_04 = int(this%n_points * 0.04_real64 + 0.5_real64)
      k_06 = int(this%n_points * 0.06_real64 + 0.5_real64)

      call calc_mutual_information(this%Mx, this%My, k_02, mi_02)
      call calc_mutual_information(this%Mx, this%My, k_04, mi_04)
      call calc_mutual_information(this%Mx, this%My, k_06, mi_06)

      @assertTrue(mi_02 > mi_04)
      @assertTrue(mi_04 > mi_06)
   end subroutine
```

---

```fortran
   subroutine calc_mutual_information(Mx, My, k, mi)
      real(real64), intent(in) :: Mx(:), My(:)
      integer, intent(in) :: k
      real(real64), intent(out) :: mi

      integer :: n_points, i
      integer :: mx_counts(size(Mx)), my_counts(size(My))
      real(real64), allocatable :: psi(:)
      real(real64), parameter :: gamma = -0.5772156649015328606_real64
      real(real64) :: avg_psi_ksg_sum

      n_points = size(Mx)
      allocate(psi(n_points + 1))

      psi = 0.0_real64
      psi(2) = gamma
      do i = 2, n_points + 1
         psi(i) = psi(i - 1) + 1.0_real64 / real(i - 1, real64)
      end do

      avg_psi_ksg_sum = 0.0_real64
      call f90_ksg_counts(Mx, My, k, mx_counts, my_counts)
      do i = 1, n_points
         avg_psi_ksg_sum = avg_psi_ksg_sum + &
               (psi(mx_counts(i) + 1) + psi(my_counts(i) + 1))
      end do

      avg_psi_ksg_sum = avg_psi_ksg_sum / real(n_points, real64)

      mi = psi(k + 1) + psi(n_points + 1) - avg_psi_ksg_sum - 1.0_real64 / real(k, real64)

      deallocate(psi)
   end subroutine calc_mutual_information
```

---

```fortran
  @test
   subroutine test_generalized_correlation_independent(this)
      class(mutual_information_fixture), intent(inout) :: this
      real(real64) :: gc

      call generate_independent_arrays(this%Mx, this%My, this%n_points)

      call calc_generalized_correlation(&
            this%Mx, this%My, this%n_points / 2, gc)

      @assertEqual(0.0_real64, gc, tolerance = 0.15_real64)
   end subroutine
      
   @test
   subroutine test_generalized_correlation_non_negative(this)
      class(mutual_information_fixture), intent(inout) :: this
      real(real64) :: gc, pearson_r

      call correlated_gaussian_arrays(this%Mx, this%My, this%n_points, -1.0_real64)

      call calc_generalized_correlation(&
            this%Mx, this%My, 4, gc)

      pearson_r = pearson_corr(this%Mx, this%My, this%n_points)

      @assertEqual(0.0_real64, gc + pearson_r, tolerance = 0.1_real64)
   end subroutine

   ! ---------------------------------------------------------------------------

   @test
   subroutine test_generalized_correlation_non_linear(this)
      class(mutual_information_fixture), intent(inout) :: this
      real(real64) :: gc

      call sin_cos_arrays(this%Mx, this%My, this%n_points)

      call calc_generalized_correlation(&
            this%Mx, this%My, 4, gc)

      @assertEqual(1.0_real64, gc, tolerance = 0.1_real64)
   end subroutine

```