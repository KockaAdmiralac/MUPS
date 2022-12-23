#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "util.h"

// Number of particles
#define MM 15
#define NPART 4 * MM *MM *MM

// Block partitioning
#define BLOCK_SIZE 450

// Tags
enum tags
{
    END_TAG = 1000,
    BLOCK_START_TAG,
    VIR_TAG,
    EPOT_TAG,
    EKIN_SUM_TAG,
    VEL_TAG,
    COUNT_TAG,
    SCALE_TAG,
};

// Constants
const double c_den = 0.83134;
const double c_tref = 0.722;
double c_rcoff = (double)MM / 4.0;
const double c_h = 0.064;
const int c_irep = 10;
const int c_istop = 20;
const int c_movemx = 20;
const double c_hsq = c_h * c_h;
const double c_hsq2 = c_hsq * 0.5;
const double c_tscale = 16.0 / ((double)NPART - 1.0);

// FCC lattice
void fcc(double x[], const int mm, const double a)
{
    int ijk = 0;

    for (int lg = 0; lg < 2; lg++)
        for (int i = 0; i < mm; i++)
            for (int j = 0; j < mm; j++)
                for (int k = 0; k < mm; k++)
                {
                    x[ijk] = i * a + lg * a * 0.5;
                    x[ijk + 1] = j * a + lg * a * 0.5;
                    x[ijk + 2] = k * a;
                    ijk += 3;
                }

    for (int lg = 1; lg < 3; lg++)
        for (int i = 0; i < mm; i++)
            for (int j = 0; j < mm; j++)
                for (int k = 0; k < mm; k++)
                {
                    x[ijk] = i * a + (2 - lg) * a * 0.5;
                    x[ijk + 1] = j * a + (lg - 1) * a * 0.5;
                    x[ijk + 2] = k * a + a * 0.5;
                    ijk += 3;
                }
}

void dfill(double a[], const int n, const double val, const int ia)
{
    for (int i = 0; i < (n - 1) * ia + 1; i += ia)
    {
        a[i] = val;
    }
}

// Change in position
// Changes: x, vh, f
// Needs: x, vh, f
void domove(double x[], double vh[], double f[], const int n3, const double side)
{
    for (int i = 0; i < n3; i++)
    {
        x[i] += vh[i] + f[i];
        // Periodic boundary conditions
        if (x[i] < 0.0)
            x[i] += side;
        if (x[i] > side)
            x[i] -= side;
        // Partial velocity updates
        vh[i] += f[i];
        // Initialise forces for the next iteration
        f[i] = 0.0;
    }
}

// Sample Maxwell distribution at temperature tref
// Used only in initialization of vh
void mxwell(double vh[], const int n3)
{
    int npart = n3 / 3;
    double ekin = 0.0;
    double sp = 0.0;

    srand48(4711);
    double tscale = 16.0 / ((double)npart - 1.0);

    for (int i = 0; i < n3; i += 2)
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

    for (int i = 0; i < n3; i += 3)
        sp += vh[i];
    sp /= (double)npart;
    for (int i = 0; i < n3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    sp = 0.0;
    for (int i = 1; i < n3; i += 3)
        sp += vh[i];
    sp /= (double)npart;
    for (int i = 1; i < n3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    sp = 0.0;
    for (int i = 2; i < n3; i += 3)
        sp += vh[i];
    sp /= (double)npart;
    for (int i = 2; i < n3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    double sc = c_h * sqrt(c_tref / (tscale * ekin));
    for (int i = 0; i < n3; i++)
        vh[i] *= sc;
}

void forces(int start, int npart, int fsize, double x[], double f[], double side, double *vir_o, double *epot_o)
{
    double vir = 0.0;
    double epot = 0.0;
    double sideh = 0.5 * side;
    double rcoffs = c_rcoff * c_rcoff;

    for (int i = 0; i < npart * 3; i += 3)
    {
        for (int j = i + 3; j < npart * 3; j += 3)
        {
            double xx = x[i] - x[j];
            double yy = x[i + 1] - x[j + 1];
            double zz = x[i + 2] - x[j + 2];
            if (xx < -sideh)
                xx += side;
            if (xx > sideh)
                xx -= side;
            if (yy < -sideh)
                yy += side;
            if (yy > sideh)
                yy -= side;
            if (zz < -sideh)
                zz += side;
            if (zz > sideh)
                zz -= side;
            double rd = xx * xx + yy * yy + zz * zz;

            if (rd <= rcoffs)
            {
                double rrd = 1.0 / rd;
                double rrd2 = rrd * rrd;
                double rrd3 = rrd2 * rrd;
                double rrd4 = rrd2 * rrd2;
                double rrd6 = rrd2 * rrd4;
                double rrd7 = rrd6 * rrd;
                epot += rrd6 - rrd3;
                double r148 = rrd7 - 0.5 * rrd4;
                vir -= rd * r148;
                double forcex = xx * r148;
                double forcey = yy * r148;
                double forcez = zz * r148;
                if (i >= start && i < start + fsize)
                {
                    f[i - start] += forcex;
                    f[i + 1 - start] += forcey;
                    f[i + 2 - start] += forcez;
                }
                if (j >= start && j < start + fsize)
                {
                    f[j - start] -= forcex;
                    f[j + 1 - start] -= forcey;
                    f[j + 2 - start] -= forcez;
                }
            }
        }
    }
    if (epot_o != NULL)
        *epot_o = epot;
    if (vir_o != NULL)
        *vir_o = vir;
}

double mkekin(const int npart, double f[], double vh[])
{
    double sum = 0.0;

    for (int i = 0; i < 3 * npart; i++)
    {
        f[i] *= c_hsq2;
        vh[i] += f[i];
        sum += vh[i] * vh[i];
    }

    return sum;
}

double velavg(const int npart, double vh[], const double vaver, double *count_o)
{
    double vaverh = vaver * c_h;
    double vel = 0.0;

    double count = 0.0;
    for (int i = 0; i < npart * 3; i += 3)
    {
        double sq = sqrt(vh[i] * vh[i] + vh[i + 1] * vh[i + 1] + vh[i + 2] * vh[i + 2]);
        if (sq > vaverh)
            count++;
        vel += sq;
    }

    // vel /= c_h;
    if (count_o != NULL)
        *count_o = count;

    return vel;
}

int main(int argc, char *argv)
{
    const double c_side = pow((double)NPART / c_den, 0.3333333);
    const double c_a = c_side / (double)MM;
    const double c_vaver = 1.13 * sqrt(c_tref / 24.0);

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER)
    {
        // MPI size
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        double x[NPART * 3];
        // Temporary vh and f for initializing x
        double *vh_t = malloc(NPART * 3 * sizeof(double));
        double *f_t = malloc(NPART * 3 * sizeof(double));

        // Initialize x
        fcc(x, MM, c_a);
        mxwell(vh_t, 3 * NPART);
        dfill(f_t, 3 * NPART, 0.0, 1);
        domove(x, vh_t, f_t, 3 * NPART, c_side);

        free(vh_t);
        free(f_t);

        // Send initial x
        MPI_Bcast(x, NPART * 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

        // Time measuring
        double time = 0.0;

        unsigned excess = 0;
        for (int move = 1; move <= c_movemx; ++move)
        {
            // Time measuring
            double start = MPI_Wtime();

            // Aggregated results
            double vir_sum = 0;
            double epot_sum = 0;
            double ekin = 0;
            double vel_sum = 0;
            double count_sum = 0;

            // Block distribution
            unsigned total_blocks = NPART / BLOCK_SIZE + (NPART % BLOCK_SIZE != 0);
            unsigned sent = 0;
            unsigned processed = 0;

            int scale = (move < c_istop && fmod(move, c_irep) == 0);

            // Send everyone a block to process
            while (sent < size - 1)
            {
                int block = sent * BLOCK_SIZE;
                if (sent < total_blocks)
                {
                    MPI_Send(&block, 1, MPI_INT, sent + 1, BLOCK_START_TAG, MPI_COMM_WORLD);
                    sent++;
                }
                else
                {
                    MPI_Send(&block, 1, MPI_INT, sent + 1, END_TAG, MPI_COMM_WORLD);
                    excess++;
                }
            }
            while (processed < total_blocks)
            {
                int block;
                MPI_Status status;
                MPI_Recv(&block, 1, MPI_LONG_INT, MPI_ANY_SOURCE, BLOCK_START_TAG, MPI_COMM_WORLD, &status);
                int worker = status.MPI_SOURCE;

                // Receive vir and epot
                double vir;
                double epot;
                MPI_Recv(&vir, 1, MPI_DOUBLE, worker, VIR_TAG, MPI_COMM_WORLD, &status);
                vir_sum += vir;
                MPI_Recv(&epot, 1, MPI_DOUBLE, worker, EPOT_TAG, MPI_COMM_WORLD, &status);
                epot_sum += epot;

                // Receive ekin sum
                double ekin_sum;
                MPI_Recv(&ekin_sum, 1, MPI_DOUBLE, worker, EPOT_TAG, MPI_COMM_WORLD, &status);
                ekin += ekin_sum;

                // Receive vel and count (velavg)
                double vel;
                double count;
                MPI_Recv(&vel, 1, MPI_DOUBLE, worker, VEL_TAG, MPI_COMM_WORLD, &status);
                vel_sum += vel;
                MPI_Recv(&count, 1, MPI_DOUBLE, worker, COUNT_TAG, MPI_COMM_WORLD, &status);
                count_sum += count;

                processed++;
            }

            // TODO Synchronization

            ekin /= c_hsq;
            vel_sum /= c_h;

            // Send scaling
            // Send entire x

            time += MPI_Wtime() - start;
            prnout(move, ekin, epot_sum, c_tscale, vir_sum, vel_sum, count_sum, NPART, c_den);
        }

        printf("%d    %lf\n", c_movemx, time);
    }
    else
    {
        // Entire x, needed for forces
        double x[3 * NPART];
        // Blocks of f and vh
        double f[3 * BLOCK_SIZE];
        double vh[3 * BLOCK_SIZE];

        // Initialize vh
        mxwell(vh, 3 * BLOCK_SIZE);

        // Initialize f
        dfill(f, 3 * BLOCK_SIZE, 0.0, 1);

        // Receive initial x
        MPI_Bcast(x, NPART * 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

        int needed = 1;
        while (needed)
        {
            int block;
            MPI_Status status;
            MPI_Recv(&block, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == BLOCK_START_TAG)
            {
                int start = block;

                double epot;
                double vir;
                forces(start, 3 * NPART, 3 * BLOCK_SIZE, x, f, c_side, &vir, &epot);
                // f <- x, -> vir, epot

                double sum;
                sum = mkekin(BLOCK_SIZE, f, vh);
                // vh <- f, -> sum

                double count;
                double vel;
                vel = velavg(BLOCK_SIZE, vh, c_vaver, &count);
                // vh, -> vel, count

                
                // Send block
                MPI_Send(&start, 1, MPI_INT, MASTER, BLOCK_START_TAG, MPI_COMM_WORLD);
                // Send vir (forces)
                MPI_Send(&vir, 1, MPI_DOUBLE, MASTER, VIR_TAG, MPI_COMM_WORLD);
                // Send epot (forces)
                MPI_Send(&epot, 1, MPI_DOUBLE, MASTER, EPOT_TAG, MPI_COMM_WORLD);
                // Send sum (mkekin)
                MPI_Send(&sum, 1, MPI_DOUBLE, MASTER, EKIN_SUM_TAG, MPI_COMM_WORLD);
                // Send vel (velavg)
                MPI_Send(&vel, 1, MPI_DOUBLE, MASTER, VEL_TAG, MPI_COMM_WORLD);
                // Send count (velavg)
                MPI_Send(&count, 1, MPI_DOUBLE, MASTER, COUNT_TAG, MPI_COMM_WORLD);
            }
            else if (status.MPI_TAG == END_TAG)
            {
                needed = 0;
            }
        }
    }

    MPI_Finalize();
    return 0;
}