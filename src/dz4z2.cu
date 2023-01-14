#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "util.h"

const float a = 3.0;
const float b = 2.0;
const float c = 1.0;
const float h = 0.001;
const int ni = 16;
const int nj = 11;
const int nk = 6;
float paramData[ni * nj * nk * 3];
float wExact[ni * nj * nk];
float results[ni * nj * nk];

__device__ float potential(float x, float y, float z)
{
    return 2.0 * ((x * x) / 81 + (y * y) / 16 + (z * z)) + 1.3611111111111112;
}

__global__ void feynman(const int N, float *params, float *results)
{
    __shared__ int sdata[NUM_OF_GPU_THREADS];
    curandState rand;
    curand_init(123456789, threadIdx.x, 0, &rand);
    float x = params[blockIdx.x * 3 + 0];
    float y = params[blockIdx.x * 3 + 1];
    float z = params[blockIdx.x * 3 + 2];
    float wt = 0;
    if (threadIdx.x == 0)
    {
        results[blockIdx.x] = 0;
    }
    for (int trial = threadIdx.x; trial < N; trial += blockDim.x)
    {
        float x1 = x;
        float x2 = y;
        float x3 = z;
        float w = 1.0;
        float chk = 0.0;
        float oldPotential = potential(x1, x2, x3);
        while (chk < 1.0)
        {
            float ut = curand_uniform(&rand);
            float us;
            float dx;
            if (ut < 0.3333333333333333)
            {
                us = curand_uniform(&rand) - 0.5;
                if (us < 0.0)
                    dx = -0.05477225575051661;
                else
                    dx = 0.05477225575051661;
            }
            else
                dx = 0.0;

            float dy;
            ut = curand_uniform(&rand);
            if (ut < 0.3333333333333333)
            {
                us = curand_uniform(&rand) - 0.5;
                if (us < 0.0)
                    dy = -0.05477225575051661;
                else
                    dy = 0.05477225575051661;
            }
            else
                dy = 0.0;

            float dz;
            ut = curand_uniform(&rand);
            if (ut < 0.3333333333333333)
            {
                us = curand_uniform(&rand) - 0.5;
                if (us < 0.0)
                    dz = -0.05477225575051661;
                else
                    dz = 0.05477225575051661;
            }
            else
                dz = 0.0;

            x1 += dx;
            x2 += dy;
            x3 += dz;

            float newPotential = potential(x1, x2, x3);
            float we = (1.0 - h * oldPotential) * w;
            w = w - 0.5 * h * (newPotential * we + oldPotential * w);
            chk = (x1 * x1) / 9.0 + (x2 * x2) / 4.0 + (x3 * x3);
            oldPotential = newPotential;
        }
        wt += w;
    }
    sdata[threadIdx.x] = wt;
    __syncthreads();
    float sdataWt = 0.0;
    if (threadIdx.x == 0)
    {
        for (int i = 0; i < NUM_OF_GPU_THREADS; ++i)
        {
            sdataWt += sdata[i];
        }
    }
    atomicAdd(&(results[blockIdx.x]), wt);
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Invalid number of arguments passed.\n");
        return 1;
    }
    const int N = atoi(argv[1]);

    printf("TEST: N=%d, num_threads=%d\n", N, NUM_OF_GPU_THREADS);

    // Time recording
    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    float wtime = 0;

    // Pre-calculate certain parameters
    int n_inside = 0;
    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            for (int k = 0; k < nk; k++)
            {
                float x = ((float)(ni - i + 1) * (-a) + (float)(i)*a) / (float)(ni - 1);
                float y = ((float)(nj - j + 1) * (-b) + (float)(j)*b) / (float)(nj - 1);
                float z = ((float)(nk - k + 1) * (-c) + (float)(k)*c) / (float)(nk - 1);
                float chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);

                if (1.0 < chk)
                {
                    continue;
                }

                float w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);
                paramData[n_inside * 3 + 0] = x;
                paramData[n_inside * 3 + 1] = y;
                paramData[n_inside * 3 + 2] = z;
                wExact[n_inside] = w_exact;
                ++n_inside;
            }
        }
    }

    // Allocate required GPU memory and copy input data
    float *gpuParams;
    float *gpuResults;
    cudaMalloc(&gpuParams, n_inside * 3 * sizeof(float));
    cudaMalloc(&gpuResults, n_inside * sizeof(float));
    cudaMemcpy(gpuParams, paramData, n_inside * 3 * sizeof(float), cudaMemcpyHostToDevice);

    // Run the kernel
    feynman<<<n_inside, NUM_OF_GPU_THREADS>>>(N, gpuParams, gpuResults);

    // Copy kernel results and free GPU data
    cudaMemcpy(results, gpuResults, n_inside * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(gpuParams);
    cudaFree(gpuResults);

    // Do block reduction
    float err = 0.0;
    for (int i = 0; i < n_inside; ++i)
    {
        err += pow(wExact[i] - (results[i] / (float)(N)), 2);
    }
    err = sqrt(err / (float)(n_inside));

    // Stop timers
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&wtime, start, stop);

    // Output results
    printf("%d    %lf    %lf\n", N, err, wtime / 1000.);
    printf("TEST END\n");

    return 0;
}
