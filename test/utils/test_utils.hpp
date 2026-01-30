#ifndef FORTRAN_TESTING_PRESENTATION_TEST_UTILS_HPP
#define FORTRAN_TESTING_PRESENTATION_TEST_UTILS_HPP

#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>

namespace test_utils {
    using real64 = double;

    // -----------------------------------------------------------------------------
    // Internal RNG seeded similarly to Fortran's system_clock seeding
    // -----------------------------------------------------------------------------
    std::mt19937_64 &rng();

    // -----------------------------------------------------------------------------
    // 1) Independent uniform arrays in [0,1)
    // -----------------------------------------------------------------------------
    void generate_independent_arrays(std::vector<real64> &a,
                                     std::vector<real64> &b, std::size_t n);

    // -----------------------------------------------------------------------------
    // 2) Correlated Gaussian arrays via Box–Muller transform
    // -----------------------------------------------------------------------------
    void correlated_gaussian_arrays(std::vector<real64> &x,
                                           std::vector<real64> &y,
                                           std::size_t n_point, real64 rho);

    // -----------------------------------------------------------------------------
    // 3) Pearson correlation
    // -----------------------------------------------------------------------------
     real64 pearson_corr(const std::vector<real64> &x,
                               const std::vector<real64> &y);

    // -----------------------------------------------------------------------------
    // 4) sin/cos paired arrays from uniform angle
    // -----------------------------------------------------------------------------
     void sin_cos_arrays(std::vector<real64> &x, std::vector<real64> &y,
                               std::size_t n);
} // namespace test_utils

extern "C" {
void c_generate_independent_arrays(double *a, double *b, int n);

void c_correlated_gaussian_arrays(double *x, double *y, int n_point,
                                  double rho);

double c_pearson_corr(const double *x, const double *y, int n);

void c_sin_cos_arrays(double *x, double *y, int n);


void c_generate_independent_arrays(double *a, double *b, int n) {
    std::vector<test_utils::real64> va, vb;
    test_utils::generate_independent_arrays(va, vb, n);

    std::copy(va.begin(), va.end(), a);
    std::copy(vb.begin(), vb.end(), b);
}

void c_correlated_gaussian_arrays(double *x, double *y, int n_point,
                                  double rho) {
    std::vector<test_utils::real64> vx, vy;
    test_utils::correlated_gaussian_arrays(vx, vy, n_point, rho);

    std::copy(vx.begin(), vx.end(), x);
    std::copy(vy.begin(), vy.end(), y);
}

double c_pearson_corr(const double *x, const double *y, int n) {
    std::vector<test_utils::real64> vx(x, x + n);
    std::vector<test_utils::real64> vy(y, y + n);

    return test_utils::pearson_corr(vx, vy);
}

void c_sin_cos_arrays(double *x, double *y, int n) {
    std::vector<test_utils::real64> vx, vy;
    test_utils::sin_cos_arrays(vx, vy, n);

    std::copy(vx.begin(), vx.end(), x);
    std::copy(vy.begin(), vy.end(), y);
}
}


#endif //FORTRAN_TESTING_PRESENTATION_TEST_UTILS_HPP
