#include <ut.hpp>
#include <vector>
#include <array>

#include "mutual_information.hpp"

using namespace boost::ut;
using real64 = double;

int main() {

    // -------------------------------------------------------------------------
    // Shared test data (same as Fortran fixture setup)
    // -------------------------------------------------------------------------
    const std::vector<real64> Mx{
        1,4,6,3,5,2,7,8
    };
    const std::vector<real64> My{
        5,6,4,8,1,7,2,3
    };

    // Build row-major J (n_points x 2)
    std::vector<real64> J(8 * 2);
    for (std::size_t i = 0; i < 8; ++i) {
        J[i*2 + 0] = Mx[i];
        J[i*2 + 1] = My[i];
    }

    const std::array<real64,2> jref{1.0, 5.0};

    // -------------------------------------------------------------------------
    "max_norm_from_point"_test = [&] {
        constexpr std::array<real64,8> expected{
            0,3,5,3,4,2,6,7
        };

        const auto dists = ksg::max_norm_from_point(
            J, jref, 8, 2
        );

        for (std::size_t i = 0; i < 8; ++i)
            expect(approx(dists[i], expected[i], 1e-12));
    };

    // -------------------------------------------------------------------------
    "k_argsort"_test = [] {
        const std::vector<real64> X{
            0,3,5,3,4,2,6,7
        };

        constexpr std::array<std::size_t,4> expected{
            0,5,1,3   // zero-based version of [1,6,2,4]
        };

        auto idxs = ksg::k_argsort(X, 4);

        for (std::size_t i = 0; i < 4; ++i)
            expect(idxs[i] == expected[i]);
    };

    // -------------------------------------------------------------------------
    "max_neighbor_distance"_test = [&] {
        constexpr std::array<std::size_t,3> neighbors{
            5,1,3   // zero-based [6,2,4]
        };

        auto dist = ksg::max_neighbor_distance(
            My, 5.0, neighbors
        );

        expect(approx(dist, 3.0, 1e-12));
    };

    // -------------------------------------------------------------------------
    "count_neighbors_within_radius"_test = [&] {
        auto count = ksg::count_neighbors_within_radius(
            My, 5.0, 3.0
        );

        expect(count == 6u);
    };

    // -------------------------------------------------------------------------
    "complete_ksg_count"_test = [&] {
        auto count = ksg::ksg_count(
            J,
            /*n_points*/ 8,
            /*n_dims*/ 2,
            /*ref_idx*/ 0,  // Fortran ref 1 → C++ 0
            /*k*/ 3
        );

        expect(count[0] == 3u);
        expect(count[1] == 6u);
    };
};
