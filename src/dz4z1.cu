#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <cuda_runtime.h>
#include "util.h"

__global__ void primes_kernel(unsigned int *results, int n_high)
{
    __shared__ int sdata[NUM_OF_GPU_THREADS];
    sdata[threadIdx.x] = 0;
    int i = 3 + (blockIdx.x * blockDim.x + threadIdx.x) * 2;
    if (i <= n_high)
    {
        int prime = 1;
        for (int j = 3; j < i; j += 2)
        {
            if (i % j == 0)
            {
                prime = 0;
            }
        }
        sdata[threadIdx.x] = prime;
    }
    __syncthreads();
    for (unsigned i = blockDim.x >> 1; i > 0; i >>= 1)
    {
        if (threadIdx.x < i)
        {
            sdata[threadIdx.x] += sdata[threadIdx.x + i];
        }
        __syncthreads();
    }
    results[blockIdx.x] = sdata[0];
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

        unsigned blocks = (n / 2) / NUM_OF_GPU_THREADS + (((n / 2) % NUM_OF_GPU_THREADS) != 0);

        unsigned *gpuResults;
        cudaMalloc(&gpuResults, blocks * sizeof(int));
        unsigned *results = (unsigned *)malloc(blocks * sizeof(int));

        primes_kernel<<<blocks, NUM_OF_GPU_THREADS>>>(gpuResults, n);

        cudaMemcpy(results, gpuResults, blocks * sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(gpuResults);

        if (n > 1)
            primes++;
        for (unsigned i = 0; i < blocks; i++)
        {
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
