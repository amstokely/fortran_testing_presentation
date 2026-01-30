---
marp: true
theme: gaia
paginate: true
---

<!-- _class: lead -->

# My Presentation Title
## A short subtitle

<!-- Speaker note: introduce yourself -->

---

## Agenda

- Problem
- Approach
- Demo
- Takeaways

---

## The Problem

Explain the problem clearly in 2–3 sentences.

---

## Key Idea

- The time complexity is $I(X;Y)=\psi(N) +\psi(k) - \langle\psi(n_x) + \psi(n_y)\rangle-\frac{1}{k}$.


---

## Example Code

```fortran
 @Test
   subroutine test_k_argsort(this)
      class(ksg_count_fixture), intent(inout) :: this
      real(real64), parameter :: X(8) = [0.0_real64, 3.0_real64, &
            5.0_real64, 3.0_real64, &
            4.0_real64, 2.0_real64, &
            6.0_real64, 7.0_real64]
      integer, parameter :: expected_idxs(4) = [1, 6, 2, 4]
      integer :: idxs(4)

      call k_argsort(X, 4, idxs)
      @assertEqual(expected_idxs, idxs)
   end subroutine
```

---

```fortran
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
```