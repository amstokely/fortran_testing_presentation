#include "mutual_information.hpp"

extern "C" {
void c_cpp_ksg_count(const double *Mx, const double *My, const int n_points,
                     const int ref_idx, const int k, int &nx, int &ny) {
    const std::span Mx_view{Mx, static_cast<std::size_t>(n_points)};
    const std::span My_view{My, static_cast<std::size_t>(n_points)};

    const auto result = ksg::cpp_ksg_count(Mx_view, My_view,
                                           static_cast<std::size_t>(ref_idx),
                                           static_cast<std::size_t>(k));
    nx = static_cast<int>(result[0]);
    ny = static_cast<int>(result[1]);
}

void c_cpp_ksg_counts(const double *Mx, const double *My, const int n_points,
                      const int k, int *nx, int *ny) {
    const std::span Mx_view{Mx, static_cast<std::size_t>(n_points)};
    const std::span My_view{My, static_cast<std::size_t>(n_points)};
    ksg::ksg_counts(Mx_view, My_view, n_points, k, nx, ny);
}

void c_cuda_ksg_counts(const double *Mx, const double *My, const int n_points,
                       const int k, int *nx, int *ny) {
    ksg::cuda_ksg_counts(Mx, My, n_points, k, nx, ny);
}
}
