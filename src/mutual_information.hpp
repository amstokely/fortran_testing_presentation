#ifndef MUTUAL_INFORMATION_HPP
#define MUTUAL_INFORMATION_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <span>

namespace ksg {
    inline std::vector<double>
    max_norm_from_point(const std::span<const double> Mx,
                        const std::span<const double> My, const double mx_ref,
                        const double my_ref) {
        const std::size_t n = Mx.size();
        std::vector<double> dists(n);

        for (std::size_t i = 0; i < n; ++i) {
            dists[i] = std::max(std::abs(Mx[i] - mx_ref),
                                std::abs(My[i] - my_ref));
        }

        return dists;
    }

    // -------------------------------------------------------------------------
    // k smallest indices (argsort)
    // -------------------------------------------------------------------------
    inline std::vector<std::size_t>
    k_argsort(const std::span<const double> X, const std::size_t k) {
        std::vector<std::size_t> idx(X.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::ranges::partial_sort(idx, idx.begin() + k,
                                  [&](const std::size_t a,
                                      const std::size_t b) {
                                      return X[a] < X[b];
                                  });

        idx.resize(k);
        return idx;
    }

    inline double
    max_neighbor_distance(const std::span<const double> X, const double xref,
                          std::span<const std::size_t> idxs) {
        double max_dist = 0.0;
        for (const auto i: idxs) max_dist = std::max(
                                     max_dist, std::abs(X[i] - xref));
        return max_dist;
    }

    inline std::size_t
    count_neighbors_within_radius(std::span<const double> X, const double xref,
                                  const double radius) {
        const auto count = std::ranges::count_if(X, [&](const double v) {
            return std::abs(v - xref) <= radius;
        });
        return count - 1;
    }

    inline std::array<std::size_t, 2>
    cpp_ksg_count(const std::span<const double> Mx,
                  const std::span<const double> My, const std::size_t ref_idx,
                  const std::size_t k) {
        const std::size_t n = Mx.size();
        const double mx_ref = Mx[ref_idx];
        const double my_ref = My[ref_idx];

        std::vector<std::size_t> idx(n);
        std::iota(idx.begin(), idx.end(), 0);

        auto max_norm = [&](std::size_t i) {
            return std::max(std::abs(Mx[i] - mx_ref), std::abs(My[i] - my_ref));
        };

        std::ranges::sort(idx, [&](const std::size_t a, const std::size_t b) {
            return max_norm(a) < max_norm(b);
        });
        double max_dmx = 0.0;
        double max_dmy = 0.0;

        for (std::size_t t = 1; t <= k; ++t) {
            const auto i = idx[t];
            max_dmx = std::max(max_dmx, std::abs(Mx[i] - mx_ref));
            max_dmy = std::max(max_dmy, std::abs(My[i] - my_ref));
        }

        std::size_t mx_count = 0;
        std::size_t my_count = 0;

        for (std::size_t i = 0; i < n; ++i) {
            if (std::abs(Mx[i] - mx_ref) <= max_dmx) ++mx_count;
            if (std::abs(My[i] - my_ref) <= max_dmy) ++my_count;
        }

        return {mx_count - 1, my_count - 1};
    }

    struct cpp_ksg_counts {
        void operator()(const std::span<const double> Mx,
                        const std::span<const double> My, const int n_points,
                        const int k, int *mx_counts, int *my_counts) const {
            for (int i = 0; i < n_points; ++i) {
                std::tie(mx_counts[i], my_counts[i]) = cpp_ksg_count(
                    Mx.subspan(0, n_points), My.subspan(0, n_points),
                    static_cast<std::size_t>(i), static_cast<std::size_t>(k));
            }
        }
    };

    template<typename KsgCountsStrategy = cpp_ksg_counts>
    void
    ksg_counts(const std::span<const double> Mx,
               const std::span<const double> My, const int n_points,
               const int k, int *mx_counts, int *my_counts) {
        KsgCountsStrategy{}(Mx, My, n_points, k, mx_counts, my_counts);
    }


#ifdef CUDA_SUPPORT
    struct cuda_ksg_counts {
        void operator()(std::span<const double> Mx, std::span<const double> My,
                        int n_points, int k, int *mx_counts,
                        int *my_counts) const;
    };
#endif
} // namespace ksg


// -----------------------------------------------------------------------------
// C interface for Fortran (Mx/My version)
// -----------------------------------------------------------------------------
extern "C" {
void c_cpp_ksg_counts(const double *Mx, const double *My, int n_points, int k,
                      int *mx_counts, int *my_counts);

#ifdef CUDA_SUPPORT
void c_cuda_ksg_counts(const double *Mx, const double *My, int n_points, int k,
                       int *mx_counts, int *my_counts);
#endif
}

#endif // MUTUAL_INFORMATION_HPP
