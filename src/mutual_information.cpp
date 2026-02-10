#include "mutual_information.hpp"
#include "cuda_mutual_information.cuh"
#include "config.hpp"
#include <stdexcept>

extern "C" {
void c_cpp_ksg_counts(
    const double *Mx, const double *My, const int n_points, const int k,
    int *mx_counts, int *my_counts
) {
    const std::span<const double> Mx_view{
        Mx, static_cast<size_t>(n_points)
    };
    const std::span<const double> My_view{
        My, static_cast<size_t>(n_points)
    };
    std::span<int> mx_view{mx_counts, static_cast<size_t>(n_points)};
    std::span<int> my_view{my_counts, static_cast<size_t>(n_points)};

    ksg::ksg_counts(Mx_view, My_view, n_points, k, mx_view, my_view);
}

void c_cuda_ksg_counts(
    const double *Mx, const double *My, const int n_points, const int k,
    int *mx_counts, int *my_counts
) {
    const std::span<const double> Mx_view{
        Mx, static_cast<std::size_t>(n_points)
    };
    const std::span<const double> My_view{
        My, static_cast<std::size_t>(n_points)
    };
    std::span<int> mx_view{mx_counts, static_cast<size_t>(n_points)};
    std::span<int> my_view{my_counts, static_cast<size_t>(n_points)};


    if constexpr (config::cuda_enabled) {
        ksg::ksg_counts<ksg::cuda_ksg_counts>(
            Mx_view, My_view, n_points, k, mx_view, my_view
        );
    } else {
        throw std::runtime_error("CUDA support not enabled");
    }
}
}
