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
subroutine test_max_norm_from_point(this)
    class(ksg_count_fixture), intent(inout) :: this
    real(real64), parameter :: expected_dists(8) = &
          [0.0_real64, 3.0_real64, 5.0_real64, 3.0_real64, &
                4.0_real64, 2.0_real64, 6.0_real64, 7.0_real64]
    real(real64) :: dists(8)

    call max_norm_from_point(this%Mx, this%My, &
          this%Mx(1), this%My(1), dists)
    @assertEqual(expected_dists, dists, tolerance = 1.0e-12_real64)
end subroutine

```

---

```fortran
@Test
subroutine test_k_argsort(this)
    class(ksg_count_fixture), intent(inout) :: this
    real(real64), parameter :: X(8) = [&
          0.0_real64, 3.0_real64, &
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
@Test
subroutine test_max_neighbor_distance(this)
    class(ksg_count_fixture), intent(inout) :: this
    integer, parameter :: neighbors(3) = [6, 2, 4]
    real(real64) :: max_dist

    call max_neighbor_distance(this%My, 5.0_real64, neighbors, max_dist)
    @assertEqual(3.0_real64, max_dist, tolerance=1.0e-12_real64)
end subroutine
```

---

```fortran 
@Test
subroutine test_count_neighbors_within_radius(this)
    class(ksg_count_fixture), intent(inout) :: this
    integer :: count

    call count_neighbors_within_radius(this%My, 5.0_real64, 3.0_real64, count)
    @assertEqual(6, count)
end subroutine
```

---

```fortran
@Test
subroutine test_complete_ksg_count(this)
    class(ksg_count_fixture), intent(inout) :: this

    integer :: nf90_nx, nf90_ny
    integer :: cpp_nx, cpp_ny

    call nf90_ksg_count(this%Mx, this%My, 1, 3, nf90_nx, nf90_ny)
    call cpp_ksg_count(this%Mx, this%My, 1, 3, cpp_nx, cpp_ny)

    @assertEqual([3, 6], [nf90_nx, nf90_ny])
    @assertEqual([3, 6], [cpp_nx, cpp_ny])
end subroutine
```

---

```fortran
 @Test
   subroutine test_complete_ksg_counts(this)
      class(ksg_count_fixture), intent(inout) :: this

      integer :: nf90_nx(8), nf90_ny(8)
      integer :: cpp_nx(8), cpp_ny(8)
      integer :: cuda_nx(8), cuda_ny(8)

      call nf90_ksg_counts(this%Mx, this%My, 3, nf90_nx, nf90_ny)
      call cpp_ksg_counts(this%Mx, this%My, 3, cpp_nx, cpp_ny)
      call cuda_ksg_counts(this%Mx, this%My, 3, cuda_nx, cuda_ny)

      @assertEqual([3, 6], [nf90_nx(1), nf90_ny(1)])
      @assertEqual([3, 6], [cpp_nx(1), cpp_ny(1)])
      @assertEqual([3, 6], [cuda_nx(1), cuda_ny(1)])
   end subroutine
```
