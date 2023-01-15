#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include <cuda_runtime.h>

#define MM 15
#define NPART 13500
#define NPART3 40500
#define DEN 0.83134
#define SIDE 25.323178916891326
#define TREF 0.722
#define RCOFF 3.75
#define H 0.064
#define IREP 10
#define ISTOP 20
#define MOVEMX 20
#define A 1.688211927792755
#define HSQ 0.004096
#define HSQ2 0.002048
#define TSCALE 0.0011852729831839394
#define VAVER 0.19599338849393194
#define SIDEH 12.661589458445663
#define RCOFFS 14.0625

double epot;
double vir;
double count;

void fcc(double x[])
{
    int ijk = 0;

    for (int lg = 0; lg < 2; lg++)
        for (int i = 0; i < MM; i++)
            for (int j = 0; j < MM; j++)
                for (int k = 0; k < MM; k++)
                {
                    x[ijk] = i * A + lg * A * 0.5;
                    x[ijk + 1] = j * A + lg * A * 0.5;
                    x[ijk + 2] = k * A;
                    ijk += 3;
                }

    for (int lg = 1; lg < 3; lg++)
        for (int i = 0; i < MM; i++)
            for (int j = 0; j < MM; j++)
                for (int k = 0; k < MM; k++)
                {
                    x[ijk] = i * A + (2 - lg) * A * 0.5;
                    x[ijk + 1] = j * A + (lg - 1) * A * 0.5;
                    x[ijk + 2] = k * A + A * 0.5;
                    ijk += 3;
                }
}

void mxwell(double vh[])
{
    double ekin = 0.0;
    double sp = 0.0;

    srand48(4711);

    for (int i = 0; i < NPART3; i += 2)
    {
        double v1;
        double v2;
        double s = 2.0;
        while (s >= 1.0)
        {
            v1 = 2.0 * drand48() - 1.0;
            v2 = 2.0 * drand48() - 1.0;
            s = v1 * v1 + v2 * v2;
        }
        double r = sqrt(-2.0 * log(s) / s);
        vh[i] = v1 * r;
        vh[i + 1] = v2 * r;
    }

    for (int i = 0; i < NPART3; i += 3)
        sp += vh[i];
    sp /= (double)NPART;
    for (int i = 0; i < NPART3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    sp = 0.0;
    for (int i = 1; i < NPART3; i += 3)
        sp += vh[i];
    sp /= (double)NPART;
    for (int i = 1; i < NPART3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    sp = 0.0;
    for (int i = 2; i < NPART3; i += 3)
        sp += vh[i];
    sp /= (double)NPART;
    for (int i = 2; i < NPART3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    double sc = H * sqrt(TREF / (TSCALE * ekin));
    for (int i = 0; i < NPART3; i++)
        vh[i] *= sc;
}

void domove(double x[], double vh[], double f[])
{
    for (int i = 0; i < NPART3; i++)
    {
        x[i] += vh[i] + f[i];
        if (x[i] < 0.0)
            x[i] += SIDE;
        if (x[i] > SIDE)
            x[i] -= SIDE;
        vh[i] += f[i];
        f[i] = 0.0;
    }
}

__global__ void forces(double *x, double *f, double *vir, double *epot)
{
    __shared__ double virs[NUM_OF_GPU_THREADS];
    __shared__ double epots[NUM_OF_GPU_THREADS];
    double myVir = 0.0;
    double myEpot = 0.0;
    int i = blockIdx.x * 3;
    for (int j = i + (threadIdx.x + 1) * 3; j < NPART3; j += 3 * blockDim.x)
    {
        double xx = x[i] - x[j];
        double yy = x[i + 1] - x[j + 1];
        double zz = x[i + 2] - x[j + 2];
        if (xx < -SIDEH)
            xx += SIDE;
        if (xx > SIDEH)
            xx -= SIDE;
        if (yy < -SIDEH)
            yy += SIDE;
        if (yy > SIDEH)
            yy -= SIDE;
        if (zz < -SIDEH)
            zz += SIDE;
        if (zz > SIDEH)
            zz -= SIDE;
        double rd = xx * xx + yy * yy + zz * zz;

        if (rd <= RCOFFS)
        {
            double rrd = 1.0 / rd;
            double rrd2 = rrd * rrd;
            double rrd3 = rrd2 * rrd;
            double rrd4 = rrd2 * rrd2;
            double rrd6 = rrd2 * rrd4;
            double rrd7 = rrd6 * rrd;
            myEpot += rrd6 - rrd3;
            double r148 = rrd7 - 0.5 * rrd4;
            myVir -= rd * r148;
            double forcex = xx * r148;
            double forcey = yy * r148;
            double forcez = zz * r148;
            atomicAdd(&f[i], forcex);
            atomicAdd(&f[j], -forcex);
            atomicAdd(&f[i + 1], forcey);
            atomicAdd(&f[j + 1], -forcey);
            atomicAdd(&f[i + 2], forcez);
            atomicAdd(&f[j + 2], -forcez);
        }
    }
    virs[threadIdx.x] = myVir;
    epots[threadIdx.x] = myEpot;
    __syncthreads();
    if (threadIdx.x == 0)
    {
        double allVirs = 0.0;
        double allEpots = 0.0;
        for (int i = 0; i < NUM_OF_GPU_THREADS; ++i)
        {
            allVirs += virs[i];
            allEpots += epots[i];
        }
        atomicAdd(vir, allVirs);
        atomicAdd(epot, allEpots);
    }
}

double mkekin(double f[], double vh[])
{
    double sum = 0.0;

    for (int i = 0; i < NPART3; i++)
    {
        f[i] *= HSQ2;
        vh[i] += f[i];
        sum += vh[i] * vh[i];
    }

    return sum / HSQ;
}

void prnout(int move, double ekin, double epot, double vir, double vel, double count)
{
    double ek = 24.0 * ekin;
    epot *= 4.0;
    double etot = ek + epot;
    double temp = TSCALE * ekin;
    double pres = DEN * 16.0 * (ekin - vir) / (double)NPART;
    vel /= (double)NPART;
    double rp = (count / (double)NPART) * 100.0;
    printf(" %6d%12.4f%12.4f%12.4f%10.4f%10.4f%10.4f%6.1f\n", move, ek, epot, etot, temp, pres, vel, rp);
}

double velavg(double vh[])
{
    double vaverh = VAVER * H;
    double vel = 0.0;

    count = 0.0;
    for (int i = 0; i < NPART3; i += 3)
    {
        double sq = sqrt(vh[i] * vh[i] + vh[i + 1] * vh[i + 1] + vh[i + 2] * vh[i + 2]);
        if (sq > vaverh)
            count++;
        vel += sq;
    }
    vel /= H;

    return vel;
}

int main(void)
{
    double x[NPART3];
    double vh[NPART3];
    double f[NPART3];

    printf("TEST: num_threads=%d\n", NUM_OF_GPU_THREADS);

    fcc(x);
    mxwell(vh);
    for (int i = 0; i < NPART3; ++i)
        f[i] = 0.0;

    double time = 0.0;

    // Allocate required GPU memory
    double *gpuX;
    double *gpuF;
    double *gpuVir;
    double *gpuEpot;
    cudaMalloc(&gpuX, NPART3 * sizeof(double));
    cudaMalloc(&gpuF, NPART3 * sizeof(double));
    cudaMalloc(&gpuVir, sizeof(double));
    cudaMalloc(&gpuEpot, sizeof(double));

    for (int move = 1; move <= MOVEMX; ++move)
    {
        // Time recording
        cudaEvent_t start;
        cudaEvent_t stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        float wtime = 0;

        // Initial setup
        domove(x, vh, f);

        // Copy the data to GPU
        cudaMemcpy(gpuX, x, NPART3 * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemset(gpuF, 0, NPART3 * sizeof(double));
        cudaMemset(gpuVir, 0, sizeof(double));
        cudaMemset(gpuEpot, 0, sizeof(double));

        // Call the kernel
        forces<<<NPART, NUM_OF_GPU_THREADS>>>(gpuX, gpuF, gpuVir, gpuEpot);

        // Copy the data back to CPU
        cudaMemcpy(f, gpuF, NPART3 * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&vir, gpuVir, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&epot, gpuEpot, sizeof(double), cudaMemcpyDeviceToHost);

        // Continue with the iteration
        double ekin = mkekin(f, vh);
        double vel = velavg(vh);
        if (move < ISTOP && fmod(move, IREP) == 0)
        {
            double sc = sqrt(TREF / (TSCALE * ekin));
            for (int i = 0; i < NPART3; i++)
                vh[i] *= sc;
            ekin = TREF / TSCALE;
        }

        // Stop timers
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&wtime, start, stop);
        time += wtime / 1000.;
        prnout(move, ekin, epot, vir, vel, count);
    }

    // Free required GPU memory
    cudaFree(gpuX);
    cudaFree(gpuF);

    // Print results
    printf("%d    %lf\n", MOVEMX, time);
    return 0;
}
