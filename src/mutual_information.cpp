#include "mutual_information.hpp"

// J is row-major: n_points × n_dims
  // counts must be size 2
  void c_ksg_count(const double* J, const int n_points, const int n_dims, const int ref_idx,   // 0-based!
                   const int k,
                   int* counts)
  {
      // Wrap raw memory in a vector view (no copy)
      std::vector<double> Jvec(J, J + static_cast<std::size_t>(n_points) * n_dims);

      const auto result = ksg::ksg_count(
          Jvec,
          static_cast<std::size_t>(n_points),
          static_cast<std::size_t>(n_dims),
          static_cast<std::size_t>(ref_idx),
          static_cast<std::size_t>(k)
      );

      counts[0] = static_cast<int>(result[0]);
      counts[1] = static_cast<int>(result[1]);
  }
