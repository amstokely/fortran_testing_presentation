#ifndef FORTRAN_TESTING_PRESENTATION_TEST_UTILS_HPP
#define FORTRAN_TESTING_PRESENTATION_TEST_UTILS_HPP

#include <vector>
#include <random>
#include <span>
#include <stdexcept>

namespace test_utils {

// -----------------------------------------------------------------------------
// RNG
// -----------------------------------------------------------------------------
std::mt19937_64& rng();

// -----------------------------------------------------------------------------
// 1) Independent uniform arrays in [0,1)
// -----------------------------------------------------------------------------
void generate_independent_arrays(std::vector<double>& a,
                                 std::vector<double>& b,
                                 std::size_t n);

// -----------------------------------------------------------------------------
// 2) Correlated Gaussian arrays via Box–Muller transform
// -----------------------------------------------------------------------------
void correlated_gaussian_arrays(std::vector<double>& x,
                                std::vector<double>& y,
                                std::size_t n_point,
                                double rho);

// -----------------------------------------------------------------------------
// 3) Pearson correlation
// -----------------------------------------------------------------------------
double pearson_corr(const std::vector<double>& x,
                    const std::vector<double>& y);

// -----------------------------------------------------------------------------
// 4) sin/cos paired arrays from uniform angle
// -----------------------------------------------------------------------------
void sin_cos_arrays(std::vector<double>& x,
                    std::vector<double>& y,
                    std::size_t n);

// -----------------------------------------------------------------------------
// Throwing backend for testing
// -----------------------------------------------------------------------------
struct ThrowingKsgCounts {
    void operator()(std::span<const double>,
                    std::span<const double>,
                    int,
                    int,
                    std::span<int>,
                    std::span<int>) const
    {
        throw std::runtime_error("forced backend failure");
    }
};

} // namespace test_utils

// -----------------------------------------------------------------------------
// C ABI wrappers (declarations only)
// -----------------------------------------------------------------------------
extern "C" {

void c_generate_independent_arrays(double* a, double* b, int n);

void c_correlated_gaussian_arrays(double* x, double* y, int n_point, double rho);

double c_pearson_corr(const double* x, const double* y, int n);

void c_sin_cos_arrays(double* x, double* y, int n);

void c_cpp_throws_ksg_counts(const double* Mx, const double* My,
                             int n_points, int k,
                             int* mx_counts, int* my_counts, int* err);
}

#endif
