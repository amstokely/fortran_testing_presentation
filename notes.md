# KSG Mutual Information with TDD, pFUnit, C++, and CUDA

---

## Tutorial Goal

Implement a **KSG Mutual Information estimator** in Fortran using
**pFUnit** and strict Test-Driven Development (TDD).

By the end:

* A correct KSG MI implementation
* A generalized correlation metric
* A C++ and CUDA-accelerated backend
* Without rewriting a single test

---

## Motivation — What is Mutual Information?

* Measures dependence between two random variables
* Detects linear and non-linear relationships
* Can be converted into a correlation metric:

$$
R = \sqrt{1 - \exp(-2,MI)}
$$

Example:

* $M_x = \sin(t)$
* $M_y = \cos(t)$

Pearson $\approx 0$
MI-based correlation $\approx 1$

---

## The KSG Estimator

$$
I^{(2)}(X,Y)
============

\psi(k)

* \frac{1}{k}
* \left\langle \psi(n_x) + \psi(n_y) \right\rangle

- \psi(N)
  $$

* Built entirely from nearest-neighbor geometry
* Uses the digamma function
* Operates point-by-point in joint space

---

## TDD Strategy

For every algorithm step:

1. Write the test
2. Stub the subroutine
3. Compile → test fails
4. Write minimal code → test passes
5. Refactor safely

We focus first on computing **ksg_counts** ($n_x, n_y$).

---

## Step 1 — Chebyshev Distance

**Test:** `test_max_norm_from_point`

```fortran
call max_norm_from_point(this%Mx, this%My, &
     this%Mx(1), this%My(1), dists)
```

**Implementation**

```fortran
dists(i) = max( abs(Mx(i) - xref_x), &
                abs(My(i) - xref_y) )
```

---

## Step 2 — k Smallest Distances (`k_argsort`)

**Test:** `test_k_argsort`

Initial implementation copies the array and repeatedly finds the minimum.
After the test passes, we refactor to maintain:

* `best_vals(k)`
* `best_idxs(k)`

Reducing memory from $O(n)$ to $O(k)$ safely under test protection.

---

## Step 3 — Marginal Radius from Neighbors

**Test:** `test_max_neighbor_distance`

Simple max loop over neighbor indices.

---

## Step 4 — Count Neighbors in Radius

**Test:** `test_count_neighbors_within_radius`

Simple loop counting values within a radius and subtracting self.

---

## Step 5 — Compose into Single Point KSG Count

**Test:** `test_complete_ksg_count`

`f90_ksg_count` calls only previously tested routines. No new math.

---

## Step 6 — All Points

**Test:** `test_ksg_counts_complete`

Loop over all points calling `f90_ksg_count`.

---

## Testing MI with Mathematical Truths

* Independent variables → MI ≈ 0
* Analytic Gaussian → MI ≈ 0.830366
* MI decreases as $k$ increases

These tests drive `calc_mutual_information`.

---

## Implementation of MI

```fortran
psi(i) = psi(i-1) + 1.0/(i-1)
call f90_ksg_counts(...)
mi = psi(k+1) + psi(N+1) - avg - 1.0/k
```

---

## Critical Numeric Lesson

Single precision passed geometric tests but failed Gaussian tests.
Double precision is required by the literature.

---

## Generalized Correlation

$$
R = \sqrt{1 - \exp(-2,MI)}
$$

Tests verify nonlinear, negative, and independent cases.

---

## Using pFUnit Fixtures

* Arrays declared once
* Populated per test in `setup`
* Tests focus on assertions

---

## Performance Problem

The Fortran implementation is fine for small analysis but becomes unusable for:

* thousands of variables
* thousands of samples

The hot path is:

```
f90_ksg_counts
```

---

## Porting the Hot Path to C++

Before going to CUDA, the algorithm is rewritten in **C++**:

* Ensures correct translation of indexing logic
* Verifies no bugs introduced in the Fortran–C interface
* Uses Boost UT to mirror the pFUnit tests
* Serial C++ version must pass all tests before CUDA work begins

This step isolates:

* Language boundary issues
* Pointer/index mistakes
* ABI and wrapper correctness

---

## Porting the Hot Path to CUDA

Only `ksg_counts` is parallelized.

Important realities of the CUDA version:

* Algorithm structure must change for parallel efficiency
* Threads compute counts for different sample points
* Must handle:

    * race conditions
    * shared memory issues
    * off-by-one indexing bugs
* Behavior must remain **identical** to Fortran and C++ versions

Without tests, this port would be extremely risky.

With tests:

* Some bugs caught by small geometric tests
* Others caught only by high-level MI property tests
* No bug survived the test suite

---

## Dependency Injection (Strategy Pattern)

We modify:

```fortran
procedure(ksg_counts_i), optional :: ksg_counts_strategy
```

in:

* `calc_mutual_information`
* `calc_generalized_correlation`

This allows us to inject:

* Fortran implementation
* C++ implementation
* CUDA implementation

No API changes. No test changes.

This is the **strategy pattern**, naturally forced by TDD.

---

## Parameterized Tests

Using pFUnit parameterized tests:

* Same test suite runs for:

    * Fortran backend
    * C++ backend
    * CUDA backend

```fortran
@testCase(testParameters={getKsgCountsStrategyParams()})
```

---

## Conditional CUDA Testing

* Preprocessor guards enable CUDA tests only when supported
* Same test suite works on CUDA and non-CUDA systems

---

## What the Tests Enabled

During CUDA development, tests caught:

* Race conditions
* Off-by-one indexing errors
* Subtle numeric differences

Some caught by low-level tests, others by MI property tests.

---

## Final Takeaway

We never wrote a giant KSG routine.

We wrote:

* small tests
* small subroutines
* mathematical property tests

From that we got:

* Correct Fortran implementation
* Safe refactors
* C++ port
* CUDA acceleration

Without rewriting a single test.
