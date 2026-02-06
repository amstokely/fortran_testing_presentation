#include "mutual_information.hpp"

extern "C" {
void c_cpp_ksg_counts(
    const double *Mx, const double *My, const int n_points, const int k, int *mx_counts, int *my_counts
) {
    const std::span Mx_view{Mx, static_cast<std::size_t>(n_points)};
    const std::span My_view{My, static_cast<std::size_t>(n_points)};
    ksg::ksg_counts(Mx_view, My_view, n_points, k, mx_counts, my_counts);
}

void c_cuda_ksg_counts(
    const double *Mx, const double *My, const int n_points, const int k, int *mx_counts, int *my_counts
) {
    const std::span Mx_view{Mx, static_cast<std::size_t>(n_points)};
    const std::span My_view{My, static_cast<std::size_t>(n_points)};
    ksg::ksg_counts<ksg::cuda_ksg_counts>(Mx_view, My_view, n_points, k, mx_counts, my_counts);
}
}
