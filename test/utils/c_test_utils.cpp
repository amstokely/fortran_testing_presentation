#include "test_utils.hpp"
#include <algorithm>
#include <mutual_information.hpp>

extern "C" {

    void c_generate_independent_arrays(double* a, double* b, int n)
    {
        std::vector<double> va, vb;
        test_utils::generate_independent_arrays(va, vb, n);

        std::copy(va.begin(), va.end(), a);
        std::copy(vb.begin(), vb.end(), b);
    }

    void c_correlated_gaussian_arrays(double* x, double* y, int n_point, double rho)
    {
        std::vector<double> vx, vy;
        test_utils::correlated_gaussian_arrays(vx, vy, n_point, rho);

        std::copy(vx.begin(), vx.end(), x);
        std::copy(vy.begin(), vy.end(), y);
    }

    double c_pearson_corr(const double* x, const double* y, int n)
    {
        std::vector<double> vx(x, x + n);
        std::vector<double> vy(y, y + n);
        return test_utils::pearson_corr(vx, vy);
    }

    void c_sin_cos_arrays(double* x, double* y, int n)
    {
        std::vector<double> vx, vy;
        test_utils::sin_cos_arrays(vx, vy, n);

        std::copy(vx.begin(), vx.end(), x);
        std::copy(vy.begin(), vy.end(), y);
    }

    void c_cpp_throws_ksg_counts(const double* Mx, const double* My,
                                 int n_points, int k,
                                 int* mx_counts, int* my_counts, int* err)
    {
        *err = 0;

        try {

            const std::span<const double> Mx_view;
            const std::span<const double> My_view;
            std::span<int> mx_view;
            std::span<int> my_view;

            ksg::ksg_counts<test_utils::ThrowingKsgCounts>(Mx_view, My_view, n_points, k,
                                              mx_view, my_view);
        }
        catch (...) {
            *err = -1;
        }
    }

}
