#include <ut.hpp>
#include <vector>
#include <array>
#include <config.hpp>

#include "mutual_information.hpp"
#include "mutual_information.cuh"
#include "config.hpp"

using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "KSG helper algorithms"_test = [] {
        given("the same 8-point dataset as the Fortran fixture") = [] {
            const std::vector<double> Mx{1, 4, 6, 3, 5, 2, 7, 8};
            const std::vector<double> My{5, 6, 4, 8, 1, 7, 2, 3};

            when("computing max-norm distances from the reference point") = [&
            ] {
                const auto dists = ksg::max_norm_from_point(
                    Mx, My, Mx[0], My[0]);

                then("the distances match the Fortran reference values") = [&] {
                    constexpr std::array<double, 8> expected{
                        0, 3, 5, 3, 4, 2, 6, 7
                    };
                    for (std::size_t i = 0; i < 8; ++i)
                        expect(approx(dists[i], expected[i], 1e-12));
                };
            };

            when("finding the k smallest indices") = [] {
                const std::vector<double> X{0, 3, 5, 3, 4, 2, 6, 7};
                auto idxs = ksg::k_argsort(X, 4);

                then("the indices match the Fortran argsort") = [&] {
                    constexpr std::array<std::size_t, 4> expected{0, 5, 1, 3};
                    for (std::size_t i = 0; i < 4; ++i)
                        expect(idxs[i] == expected[i]);
                };
            };

            when("computing max neighbor distance in one dimension") = [&] {
                constexpr std::array<std::size_t, 3> neighbors{5, 1, 3};

                const auto dist =
                        ksg::max_neighbor_distance(My, 5.0, neighbors);

                then("the distance equals the Fortran result") = [&] {
                    expect(approx(dist, 3.0, 1e-12));
                };
            };

            when("counting neighbors within a radius") = [&] {
                const auto count = ksg::count_neighbors_within_radius(
                    My, 5.0, 3.0);

                then("the count matches the Fortran reference") = [&] {
                    expect(count == 6u);
                };
            };

            when("running the complete ksg_count algorithm") = [&] {
                auto mx_counts = std::make_unique<int[]>(8);
                auto my_counts = std::make_unique<int[]>(8);
                auto mx_view = std::span<int>{mx_counts.get(), 8};
                auto my_view = std::span<int>{my_counts.get(), 8};
                ksg::ksg_counts<ksg::cpp_ksg_counts>(
                    Mx, My, static_cast<int>(Mx.size()),
                    /*k*/ 3, mx_view, my_view);

                then("the first values in mx_counts and my_counts match the Fortran implementation") = [
                    &] {
                    expect(mx_counts[0] == 3u);
                    expect(my_counts[0] == 6u);
                };
            };
            if constexpr(config::cuda_enabled) {
                when("running the complete ksg_count cuda algorithm") = [&] {
                    auto mx_counts = std::make_unique<int[]>(8);
                    auto my_counts = std::make_unique<int[]>(8);
                    auto mx_view = std::span<int>{mx_counts.get(), 8};
                    auto my_view = std::span<int>{my_counts.get(), 8};
                    ksg::ksg_counts<ksg::cuda_ksg_counts>(
                        Mx, My, static_cast<int>(Mx.size()),
                        /*k*/ 3, mx_view, my_view);

                    then("the first values in mx_counts and my_counts match the Fortran implementation") = [
                        &] {
                        expect(mx_counts[0] == 3u);
                        expect(my_counts[0] == 6u);
                    };
                };
            }
        };
    };
}
