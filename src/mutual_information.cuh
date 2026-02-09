#ifndef FORTRAN_TESTING_PRESENTATION_MUTUAL_INFORMATION_CUH
#define FORTRAN_TESTING_PRESENTATION_MUTUAL_INFORMATION_CUH

#include <span>

namespace ksg {
    struct cuda_ksg_counts {
        void operator()(
            std::span<const double> Mx, std::span<const double> My,
            int n_points, int k, std::span<int> mx_counts, std::span<int> my_counts
        ) const;
    };
} // namespace ksg

#endif
