#ifndef MUTUAL_INFORMATION_HPP
#define MUTUAL_INFORMATION_HPP


#include <vector>
#include <span>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace ksg {
    // ------------------------------------------------------------
    // Row-major matrix accessor
    // ------------------------------------------------------------
    inline double Xat(const std::vector<double> &X, std::size_t n_dims,
                      std::size_t i, std::size_t j) {
        return X[i * n_dims + j];
    }

    // ------------------------------------------------------------
    // Max norm distance from reference point
    // ------------------------------------------------------------
    inline std::vector<double>
    max_norm_from_point(const std::vector<double> &X,
                        const std::span<const double> xref,
                        const std::size_t n_points, const std::size_t n_dims) {
        std::vector<double> dists(n_points);

        for (std::size_t i = 0; i < n_points; ++i) {
            double max_diff = 0.0;
            for (std::size_t j = 0; j < n_dims; ++j) {
                max_diff = std::max(max_diff,
                                    std::abs(Xat(X, n_dims, i, j) - xref[j]));
            }
            dists[i] = max_diff;
        }

        return dists;
    }

    // ------------------------------------------------------------
    // k smallest indices (argsort)
    // ------------------------------------------------------------
    inline std::vector<std::size_t>
    k_argsort(const std::vector<double> &X, const std::size_t k) {
        std::vector<std::size_t> idx(X.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::ranges::partial_sort(idx, idx.begin() + k,
                                  [&](const std::size_t a, const std::size_t b) {
                                      return X[a] < X[b];
                                  });
        idx.resize(k);
        return idx;
    }

    // ------------------------------------------------------------
    // Max neighbor distance in one dimension
    // ------------------------------------------------------------
    inline double
    max_neighbor_distance(std::span<const double> X, double xref,
                          std::span<const std::size_t> idxs) {
        double max_dist = 0.0;

        for (const auto i: idxs) max_dist = std::max(max_dist, std::abs(X[i] - xref));

        return max_dist;
    }

    // ------------------------------------------------------------
    // Count neighbors within radius (excluding self)
    // ------------------------------------------------------------
    inline std::size_t
    count_neighbors_within_radius(std::span<const double> X, const double xref,
                                  const double radius) {
        const std::size_t count = std::ranges::count_if(X, [&](const double v) {
            return std::abs(v - xref) <= radius;
        });

        return count - 1;
    }

    // ------------------------------------------------------------
    // Core KSG neighbor counting
    // ------------------------------------------------------------
    inline std::array<std::size_t, 2>
    ksg_count(const std::vector<double> &J, std::size_t n_points,
              std::size_t n_dims, std::size_t ref_idx, const std::size_t k) {
        // Extract reference row
        std::vector<double> xref(n_dims);
        for (std::size_t j = 0; j < n_dims; ++j)
            xref[j] = Xat(J, n_dims, ref_idx, j);

        // Compute max-norm distances
        const auto dists = max_norm_from_point(J, xref, n_points, n_dims);

        // k+1 nearest neighbors (includes self)
        auto neighbor_idxs = k_argsort(dists, k + 1);

        // Project columns without copying whole matrix
        std::vector<double> Xcol(n_points), Ycol(n_points);
        for (std::size_t i = 0; i < n_points; ++i) {
            Xcol[i] = Xat(J, n_dims, i, 0);
            Ycol[i] = Xat(J, n_dims, i, 1);
        }

        // Skip self
        const std::span<const std::size_t> neigh(neighbor_idxs.begin() + 1,
                                           neighbor_idxs.end());

        const double max_dist_x = max_neighbor_distance(Xcol, Xcol[ref_idx], neigh);

        const double max_dist_y = max_neighbor_distance(Ycol, Ycol[ref_idx], neigh);

        const std::size_t count_x = count_neighbors_within_radius(
            Xcol, Xcol[ref_idx], max_dist_x);

        const std::size_t count_y = count_neighbors_within_radius(
            Ycol, Ycol[ref_idx], max_dist_y);
        return {count_x, count_y};
    }
} // namespace ksg

extern "C" {
// J is row-major: n_points × n_dims
// counts must be size 2
void c_ksg_count(const double *J, int n_points, int n_dims, int ref_idx,
                 // 0-based!
                 int k, int *counts);
}

#endif // MUTUAL_INFORMATION_HPP
