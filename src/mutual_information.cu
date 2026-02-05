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
    int blockSize{};
    int numBlocks{};
    int registerSize{};

    explicit CudaConfig(const int n_points) {
        blockSize = 1024;
        numBlocks = (n_points + blockSize - 1) / blockSize;
        registerSize = (n_points / blockSize) + 1;
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


template<int XRegisterSize> __global__ void ksg_counts_kernel(
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
        double r_mx[XRegisterSize];
        double r_my[XRegisterSize];
        double r_mx_i = mx[i];
        double r_my_i = my[i];
        const unsigned int l_tidx = threadIdx.x;
        unsigned int idx = l_tidx;
        s_mx_counts[l_tidx] = 0;
        s_my_counts[l_tidx] = 0;
        if (l_tidx == 0) {
            s_eps_mx[l_tidx] = 0.0;
            s_eps_my[l_tidx] = 0.0;
        }
        int r_idx = 0;
        while (idx < n_points) {
            r_mx[r_idx] = mx[idx];
            r_my[r_idx] = my[idx];
            idx += blockDim.x;
            r_idx++;
        }
        __syncthreads();
        for (int j = 0; j < k + 1; j++) {
            idx = l_tidx;
            auto l_min = static_cast<double>(RAND_MAX);
            int l_min_idx = 0;
            r_idx = 0;
            while (idx < n_points) {
                const double dX = fmax(
                    abs(r_mx_i - r_mx[r_idx]), abs(r_my_i - r_my[r_idx])
                );
                if (dX < l_min) {
                    l_min = dX;
                    l_min_idx = idx;
                }
                r_idx++;
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
                            r_mx_i - r_mx[s_min_idx[0] / blockDim.x]
                        )
                    );
                    s_eps_my[0] = fmax(
                        s_eps_my[0], abs(
                            r_my_i - r_my[s_min_idx[0] / blockDim.x]
                        )
                    );
                }
                r_mx[s_min_idx[0] / blockDim.x] = static_cast<double>(
                    RAND_MAX);
                r_my[s_min_idx[0] / blockDim.x] = static_cast<double>(
                    RAND_MAX);
            }
            __syncthreads();
        }
        idx = l_tidx;
        r_idx = 0;
        while (idx < n_points) {
            s_mx_counts[l_tidx] +=
                    (abs(r_mx[r_idx] - r_mx_i) <= s_eps_mx[0]) ? 1 : 0;
            s_my_counts[l_tidx] +=
                    (abs(r_my[r_idx] - r_my_i) <= s_eps_my[0]) ? 1 : 0;
            r_idx++;
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


template<int R> void launch_kernel(
    const double *d_Mx, const double *d_My, int k, const int n_points,
    int *d_nx, int *d_ny, int numBlocks, int blockSize
) {
    ksg_counts_kernel<R> <<<numBlocks, blockSize>>>(
        d_Mx, d_My, k, n_points, d_nx, d_ny
    );
}

template<int... I> void dispatch_kernel(
    int r, std::integer_sequence<int, I...>, const double *d_Mx,
    const double *d_My, int k, int n_points, int *d_nx, int *d_ny,
    const int numBlocks, const int blockSize
) {
    ((r == I + 1 ? (launch_kernel<I + 1>(
                        d_Mx, d_My, k, n_points, d_nx, d_ny, numBlocks,
                        blockSize
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
            config.registerSize, std::make_integer_sequence<int, 101>{},
            cuda_mx.d_array, cuda_my.d_array, k, n_points,
            cuda_nx.d_array, cuda_ny.d_array, config.numBlocks,
            config.blockSize
        );

        cuda_nx.to_host(mx_counts);
        cuda_ny.to_host(my_counts);
    }
}
