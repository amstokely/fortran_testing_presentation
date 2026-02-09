#ifndef FORTRAN_TESTING_PRESENTATION_MUTUAL_INFORMATION_CUH
#define FORTRAN_TESTING_PRESENTATION_MUTUAL_INFORMATION_CUH
#include <span>

struct cuda_ksg_counts {
    void operator()(std::span<const double> Mx, std::span<const double> My,
                    int n_points, int k, int *mx_counts,
                    int *my_counts) const;
};

#endif //FORTRAN_TESTING_PRESENTATION_MUTUAL_INFORMATION_CUH