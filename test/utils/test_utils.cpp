#include "test_utils.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

namespace test_utils {

// -----------------------------------------------------------------------------
// RNG
// -----------------------------------------------------------------------------
std::mt19937_64& rng() {
    static std::mt19937_64 engine([] {
        const auto now =
            std::chrono::high_resolution_clock::now().time_since_epoch().count();

        std::seed_seq seq{
            static_cast<std::uint32_t>(now),
            static_cast<std::uint32_t>(now >> 32),
            0x9E3779B9u, 0x243F6A88u, 0xB7E15162u
        };
        return std::mt19937_64(seq);
    }());

    return engine;
}

// -----------------------------------------------------------------------------
// Independent uniform arrays
// -----------------------------------------------------------------------------
void generate_independent_arrays(std::vector<double>& a,
                                 std::vector<double>& b,
                                 std::size_t n)
{
    a.resize(n);
    b.resize(n);

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t i = 0; i < n; ++i) {
        a[i] = dist(rng());
        b[i] = dist(rng());
    }
}

// -----------------------------------------------------------------------------
// Correlated Gaussian arrays
// -----------------------------------------------------------------------------
void correlated_gaussian_arrays(std::vector<double>& x,
                                std::vector<double>& y,
                                std::size_t n_point,
                                double rho)
{
    constexpr double pi = 3.1415926535897932384626433832795;
    const double s = std::sqrt(1.0 - rho * rho);

    x.resize(n_point);
    y.resize(n_point);

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t i = 0; i < n_point; ++i) {
        double u1 = std::max(dist(rng()), 1e-12);
        const double u2 = dist(rng());

        const double r = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * pi * u2;

        const double gx = r * std::cos(theta);
        const double gy = r * std::sin(theta);

        x[i] = gx;
        y[i] = rho * gx + s * gy;
    }
}

// -----------------------------------------------------------------------------
// Pearson correlation
// -----------------------------------------------------------------------------
double pearson_corr(const std::vector<double>& x,
                    const std::vector<double>& y)
{
    const std::size_t n = x.size();

    const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
    const double my = std::accumulate(y.begin(), y.end(), 0.0) / n;

    double sx = 0.0, sy = 0.0, sxy = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sx += dx * dx;
        sy += dy * dy;
        sxy += dx * dy;
    }

    return sxy / std::sqrt(sx * sy);
}

// -----------------------------------------------------------------------------
// sin/cos arrays
// -----------------------------------------------------------------------------
void sin_cos_arrays(std::vector<double>& x,
                    std::vector<double>& y,
                    std::size_t n)
{
    constexpr double pi = 3.1415926535897932384626433832795;

    x.resize(n);
    y.resize(n);

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t i = 0; i < n; ++i) {
        const double theta = 2.0 * pi * dist(rng());
        x[i] = std::sin(theta);
        y[i] = std::cos(theta);
    }
}

} // namespace test_utils
