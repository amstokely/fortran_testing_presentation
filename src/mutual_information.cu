#include "mutual_information.hpp"


template<typename T> struct CudaArray {
    size_t bytes{};
    T *d_array = nullptr;

    explicit CudaArray(const size_t n_elements) : bytes(
        n_elements * sizeof(T)
    ), d_array(nullptr) {
        cudaMalloc(&d_array, bytes);
    }

    void to_device(const T *h_array) const {
        cudaMemcpy(d_array, h_array, bytes, cudaMemcpyHostToDevice);
    }

    void to_host(T *h_array) const {
        cudaMemcpy(h_array, d_array, bytes, cudaMemcpyDeviceToHost);
    }

    ~CudaArray() {
        cudaFree(d_array);
    }
};

struct CudaConfig {
    int block_dim_x{};
    int grid_dim_x{};
    int local_m_size{};

    explicit CudaConfig(const int n_points) {
        block_dim_x = 1024;
        grid_dim_x = (n_points + block_dim_x - 1) / block_dim_x;
        local_m_size = (n_points / block_dim_x) + 1;
    }
};


__device__ void warp_reduce_kernel(
    volatile int *s_a, const int l_tidx
) {
    s_a[l_tidx] += s_a[l_tidx + 32];
    s_a[l_tidx] += s_a[l_tidx + 16];
    s_a[l_tidx] += s_a[l_tidx + 8];
    s_a[l_tidx] += s_a[l_tidx + 4];
    s_a[l_tidx] += s_a[l_tidx + 2];
    s_a[l_tidx] += s_a[l_tidx + 1];
}

__device__ void min_warp_reduce_kernel(
    volatile double *s_min, volatile int *s_min_idx, const int l_tidx
) {
    if (s_min[l_tidx] > s_min[l_tidx + 32]) {
        s_min[l_tidx] = s_min[l_tidx + 32];
        s_min_idx[l_tidx] = s_min_idx[l_tidx + 32];
    }
    if (s_min[l_tidx] > s_min[l_tidx + 16]) {
        s_min[l_tidx] = s_min[l_tidx + 16];
        s_min_idx[l_tidx] = s_min_idx[l_tidx + 16];
    }
    if (s_min[l_tidx] > s_min[l_tidx + 8]) {
        s_min[l_tidx] = s_min[l_tidx + 8];
        s_min_idx[l_tidx] = s_min_idx[l_tidx + 8];
    }
    if (s_min[l_tidx] > s_min[l_tidx + 4]) {
        s_min[l_tidx] = s_min[l_tidx + 4];
        s_min_idx[l_tidx] = s_min_idx[l_tidx + 4];
    }
    if (s_min[l_tidx] > s_min[l_tidx + 2]) {
        s_min[l_tidx] = s_min[l_tidx + 2];
        s_min_idx[l_tidx] = s_min_idx[l_tidx + 2];
    }
    if (s_min[l_tidx] > s_min[l_tidx + 1]) {
        s_min[l_tidx] = s_min[l_tidx + 1];
        s_min_idx[l_tidx] = s_min_idx[l_tidx + 1];
    }
}

template<int LocalMSize> __global__ void ksg_counts_kernel(
    const double *mx, const double *my, const int k, const int n_points,
    int *mx_counts, int *my_counts
) {
    __shared__ int s_min_idx[1024];
    __shared__ volatile double s_min[1024];
    __shared__ double s_eps_mx[1];
    __shared__ double s_eps_my[1];
    __shared__ int s_mx_counts[1024];
    __shared__ int s_my_counts[1024];
    for (auto i = blockIdx.x; i < n_points; i += gridDim.x) {
        double l_mx[LocalMSize];
        double l_my[LocalMSize];
        double l_mx_i = mx[i];
        double l_my_i = my[i];
        const unsigned int l_tidx = threadIdx.x;
        unsigned int idx = l_tidx;
        s_mx_counts[l_tidx] = 0;
        s_my_counts[l_tidx] = 0;
        if (l_tidx == 0) {
            s_eps_mx[l_tidx] = 0.0;
            s_eps_my[l_tidx] = 0.0;
        }
        int l_m_idx = 0;
        while (idx < n_points) {
            l_mx[l_m_idx] = mx[idx];
            l_my[l_m_idx] = my[idx];
            idx += blockDim.x;
            l_m_idx++;
        }
        __syncthreads();
        for (int j = 0; j < k + 1; j++) {
            idx = l_tidx;
            auto l_min = static_cast<double>(RAND_MAX);
            int l_min_idx = 0;
            l_m_idx = 0;
            while (idx < n_points) {
                const double dX = fmax(
                    abs(l_mx_i - l_mx[l_m_idx]), abs(l_my_i - l_my[l_m_idx])
                );
                if (dX < l_min) {
                    l_min = dX;
                    l_min_idx = idx;
                }
                l_m_idx++;
                idx += blockDim.x;
            }
            s_min[l_tidx] = l_min;
            s_min_idx[l_tidx] = l_min_idx;
            __syncthreads();
            if (l_tidx < 512) {
                if (s_min[l_tidx] > s_min[l_tidx + 512]) {
                    s_min[l_tidx] = s_min[l_tidx + 512];
                    s_min_idx[l_tidx] = s_min_idx[l_tidx + 512];
                }
            }
            __syncthreads();
            if (l_tidx < 256) {
                if (s_min[l_tidx] > s_min[l_tidx + 256]) {
                    s_min[l_tidx] = s_min[l_tidx + 256];
                    s_min_idx[l_tidx] = s_min_idx[l_tidx + 256];
                }
            }
            __syncthreads();
            if (l_tidx < 128) {
                if (s_min[l_tidx] > s_min[l_tidx + 128]) {
                    s_min[l_tidx] = s_min[l_tidx + 128];
                    s_min_idx[l_tidx] = s_min_idx[l_tidx + 128];
                }
            }
            __syncthreads();
            if (l_tidx < 64) {
                if (s_min[l_tidx] > s_min[l_tidx + 64]) {
                    s_min[l_tidx] = s_min[l_tidx + 64];
                    s_min_idx[l_tidx] = s_min_idx[l_tidx + 64];
                }
            }
            __syncthreads();
            if (l_tidx < 32) {
                min_warp_reduce_kernel(s_min, s_min_idx, l_tidx);
            }
            __syncthreads();
            if ((s_min_idx[0] - (
                     (s_min_idx[0] / blockDim.x) * blockDim.x)) ==
                l_tidx) {
                if (j > 0) {
                    s_eps_mx[0] = fmax(
                        s_eps_mx[0], abs(
                            l_mx_i - l_mx[s_min_idx[0] / blockDim.x]
                        )
                    );
                    s_eps_my[0] = fmax(
                        s_eps_my[0], abs(
                            l_my_i - l_my[s_min_idx[0] / blockDim.x]
                        )
                    );
                }
                l_mx[s_min_idx[0] / blockDim.x] = static_cast<double>(
                    RAND_MAX);
                l_my[s_min_idx[0] / blockDim.x] = static_cast<double>(
                    RAND_MAX);
            }
            __syncthreads();
        }
        idx = l_tidx;
        l_m_idx = 0;
        while (idx < n_points) {
            s_mx_counts[l_tidx] +=
                    (abs(l_mx[l_m_idx] - l_mx_i) <= s_eps_mx[0]) ? 1 : 0;
            s_my_counts[l_tidx] +=
                    (abs(l_my[l_m_idx] - l_my_i) <= s_eps_my[0]) ? 1 : 0;
            l_m_idx++;
            idx += blockDim.x;
        }
        __syncthreads();
        if (l_tidx < 512) {
            s_mx_counts[l_tidx] += s_mx_counts[l_tidx + 512];
            s_my_counts[l_tidx] += s_my_counts[l_tidx + 512];
        }
        __syncthreads();
        if (l_tidx < 256) {
            s_mx_counts[l_tidx] += s_mx_counts[l_tidx + 256];
            s_my_counts[l_tidx] += s_my_counts[l_tidx + 256];
        }
        __syncthreads();
        if (l_tidx < 128) {
            s_mx_counts[l_tidx] += s_mx_counts[l_tidx + 128];
            s_my_counts[l_tidx] += s_my_counts[l_tidx + 128];
        }
        __syncthreads();
        if (l_tidx < 64) {
            s_mx_counts[l_tidx] += s_mx_counts[l_tidx + 64];
            s_my_counts[l_tidx] += s_my_counts[l_tidx + 64];
        }
        __syncthreads();

        if (l_tidx < 32) {
            warp_reduce_kernel(s_mx_counts, l_tidx);
            warp_reduce_kernel(s_my_counts, l_tidx);
        }
        if (l_tidx == 0) {
            mx_counts[i] = s_mx_counts[0] + k;
            my_counts[i] = s_my_counts[0] + k;
        }
        __syncthreads();
    }
}


template<int LocalMSize> void launch_kernel(
    const double *d_mx, const double *d_my, int k, const int n_points,
    int *d_mx_counts, int *d_my_counts, int grid_dim_x, int block_dim_x
) {
    ksg_counts_kernel<LocalMSize> <<<grid_dim_x, block_dim_x>>>(
        d_mx, d_my, k, n_points, d_mx_counts, d_my_counts
    );
}

template<int... I> void dispatch_kernel(
    int local_m_size, std::integer_sequence<int, I...>,
    const double *d_mx, const double *d_my, const int k,
    const int n_points, int *d_mx_counts, int *d_my_counts,
    const int grid_dim_x, const int block_dim_x
) {
    ((local_m_size == I + 1 ? (launch_kernel<I + 1>(
                                   d_mx, d_my, k, n_points, d_mx_counts,
                                   d_my_counts, grid_dim_x, block_dim_x
                               ), 0) : 0), ...);
}

namespace ksg {
    void cuda_ksg_counts(
        const double *mx, const double *my, const int n_points,
        const int k, int *mx_counts, int *my_counts
    ) {
        const auto cuda_mx = CudaArray<double>(n_points);
        const auto cuda_my = CudaArray<double>(n_points);
        const auto cuda_nx = CudaArray<int>(n_points);
        const auto cuda_ny = CudaArray<int>(n_points);

        cuda_mx.to_device(mx);
        cuda_my.to_device(my);


        const auto config = CudaConfig(n_points);

        dispatch_kernel(
            config.local_m_size, std::make_integer_sequence<int, 101>{},
            cuda_mx.d_array, cuda_my.d_array, k, n_points,
            cuda_nx.d_array, cuda_ny.d_array, config.grid_dim_x,
            config.block_dim_x
        );

        cuda_nx.to_host(mx_counts);
        cuda_ny.to_host(my_counts);
    }
}
