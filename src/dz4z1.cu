#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <cuda_runtime.h>
#include "util.h"

// __global__ void primes_kernel(unsigned int *results, int n_high)
// {
//     extern __shared__ int sdata[];
//     sdata[threadIdx.x] = 0;
//     for (int i = 3 + (blockIdx.x * blockDim.x + threadIdx.x) * 2; i <= n_high; i += blockDim.x * gridDim.x)
//     {
//         int prime = 1;
//         for (int j = 3; j < i; j += 2)
//         {
//             if (i % j == 0)
//             {
//                 prime = 0;
//                 break;
//             }
//         }
//         if (prime)
//         {
//             sdata[threadIdx.x]++;
//         }
//     }
//     __syncthreads();
//     if (threadIdx.x == 0)
//     {
//         unsigned primes = 0;
//         for (unsigned i = 0; i < blockDim.x; i++)
//         {
//             primes += sdata[i];
//         }
//         results[blockIdx.x] = primes;
//     }
// }

__global__ void primes_kernel(unsigned int *results, int n_high)
{
    extern __shared__ int sdata[];
    sdata[threadIdx.x] = 0;
    for (int i = 3 + threadIdx.x * 2; i <= n_high; i += blockDim.x * gridDim.x)
    {
        int prime = 1;
        for (int j = 3; j < i; j += 2)
        {
            if (i % j == 0)
            {
                prime = 0;
                break;
            }
        }
        sdata[threadIdx.x] += prime;
    }
    __syncthreads();
    if (threadIdx.x == 0)
    {
        unsigned primes = 0;
        for (unsigned i = 0; i < blockDim.x; i++)
        {
            primes += sdata[i];
        }
        results[blockIdx.x] = primes;
    }
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n = n_lo;

    while (n <= n_hi)
    {
        float wtime = 0;
        cudaEvent_t start, stop;

        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        int primes = 0;

        unsigned blocks = (n / 2) / NUM_OF_GPU_THREADS + ((n / 2) % NUM_OF_GPU_THREADS != 0);

        unsigned *gpuResults;
        cudaMalloc(&gpuResults, blocks * sizeof(int));
        unsigned *results = (unsigned *)malloc(blocks * sizeof(int));
        printf("blocks: %d --", blocks);

        primes_kernel<<<blocks, NUM_OF_GPU_THREADS, NUM_OF_GPU_THREADS * sizeof(int)>>>(gpuResults, n);

        cudaMemcpy(results, gpuResults, blocks * sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(gpuResults);

        if (n > 1)
            primes++;
        for (unsigned i = 0; i < blocks; i++)
        {
            printf("%d ", results[i]);
            primes += results[i];
        }

        free(results);

        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&wtime, start, stop);

        printf("  %8d  %8d  %14f\n", n, primes, wtime / 1000.);
        n = n * n_factor;
    }
}

int main(int argc, char *argv[])
{
    int lo;
    int hi;
    int factor;

    if (argc != 4)
    {
        lo = 1;
        hi = 131072;
        factor = 2;
    }
    else
    {
        lo = atoi(argv[1]);
        hi = atoi(argv[2]);
        factor = atoi(argv[3]);
    }

    printf("TEST: lo=%d, hi=%d, factor=%d\n", lo, hi, factor);
    test(lo, hi, factor);

    return 0;
}