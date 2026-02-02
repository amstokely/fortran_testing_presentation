#ifndef MUTUAL_INFORMATION_HPP
#define MUTUAL_INFORMATION_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <span>

namespace ksg {

    // -------------------------------------------------------------------------
    // Max norm distance from reference point (Mx/My version)
    // -------------------------------------------------------------------------
    inline std::vector<double>
    max_norm_from_point(std::span<const double> Mx,
                        std::span<const double> My,
                        double xref_x,
                        double xref_y)
    {
        const std::size_t n = Mx.size();
        std::vector<double> dists(n);

        for (std::size_t i = 0; i < n; ++i) {
            dists[i] = std::max(
                std::abs(Mx[i] - xref_x),
                std::abs(My[i] - xref_y)
            );
        }

        return dists;
    }

    // -------------------------------------------------------------------------
    // k smallest indices (argsort)
    // -------------------------------------------------------------------------
    inline std::vector<std::size_t>
    k_argsort(std::span<const double> X, std::size_t k)
    {
        std::vector<std::size_t> idx(X.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::ranges::partial_sort(
            idx, idx.begin() + k,
            [&](const std::size_t a, const std::size_t b) {
                return X[a] < X[b];
            });

        idx.resize(k);
        return idx;
    }

    // -------------------------------------------------------------------------
    // Max neighbor distance in one dimension
    // -------------------------------------------------------------------------
    inline double
    max_neighbor_distance(std::span<const double> X,
                          double xref,
                          std::span<const std::size_t> idxs)
    {
        double max_dist = 0.0;
        for (auto i : idxs)
            max_dist = std::max(max_dist, std::abs(X[i] - xref));
        return max_dist;
    }

    // -------------------------------------------------------------------------
    // Count neighbors within radius (excluding self)
    // -------------------------------------------------------------------------
    inline std::size_t
    count_neighbors_within_radius(std::span<const double> X,
                                  double xref,
                                  double radius)
    {
        const auto count = std::ranges::count_if(X, [&](double v) {
            return std::abs(v - xref) <= radius;
        });
        return count - 1;
    }

    // -------------------------------------------------------------------------
    // Fast KSG neighbor counting (no allocations except idx)
    // -------------------------------------------------------------------------
    inline std::array<std::size_t, 2>
    cpp_ksg_count(std::span<const double> Mx,
              std::span<const double> My,
              std::size_t ref_idx,
              std::size_t k)
    {
        const std::size_t n = Mx.size();
        const double xref_x = Mx[ref_idx];
        const double xref_y = My[ref_idx];

        std::vector<std::size_t> idx(n);
        std::iota(idx.begin(), idx.end(), 0);

        auto max_norm = [&](std::size_t i) {
            return std::max(
                std::abs(Mx[i] - xref_x),
                std::abs(My[i] - xref_y)
            );
        };

        std::ranges::nth_element(idx, idx.begin() + (k + 1),
                                 [&](const std::size_t a, const std::size_t b) {
                                     return max_norm(a) < max_norm(b);
                                 });

        double max_dx = 0.0;
        double max_dy = 0.0;

        for (std::size_t t = 1; t <= k; ++t) {
            auto i = idx[t];
            max_dx = std::max(max_dx, std::abs(Mx[i] - xref_x));
            max_dy = std::max(max_dy, std::abs(My[i] - xref_y));
        }

        std::size_t count_x = 0;
        std::size_t count_y = 0;

        for (std::size_t i = 0; i < n; ++i) {
            if (std::abs(Mx[i] - xref_x) <= max_dx) ++count_x;
            if (std::abs(My[i] - xref_y) <= max_dy) ++count_y;
        }

        return {count_x - 1, count_y - 1};
    }

} // namespace ksg


// -----------------------------------------------------------------------------
// C interface for Fortran (Mx/My version)
// -----------------------------------------------------------------------------
extern "C" {
void c_cpp_ksg_count(const double* Mx,
                 const double* My,
                 int n_points,
                 int ref_idx,   // 0-based
                 int k,
                 int& nx,
                 int& ny);
void c_cpp_ksg_counts(const double* Mx,
                  const double* My,
                  int n_points,
                  int k,
                  int* nx,
                  int* ny);
}

#endif // MUTUAL_INFORMATION_HPP
