#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include "util.h"

double potential(double a, double b, double c, double x, double y, double z)
{
    return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

double r8_uniform_01(int *seed)
{
    int k = *seed / 127773;

    *seed = 16807 * (*seed - k * 127773) - k * 2836;

    if (*seed < 0)
    {
        *seed = *seed + 2147483647;
    }
    return (double)(*seed) * 4.656612875E-10;
}

double feynman(const int N)
{
    const double a = 3.0;
    const double b = 2.0;
    const double c = 1.0;
    const double h = 0.001;
    const double stepsz = sqrt(3 * h);
    const int ni = 16;
    const int nj = 11;
    const int nk = 6;
    int n_inside = 0;
    double err = 0.0;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int seed = 123456789 + rank;
    for (int ijk = rank; ijk < ni * nj * nk; ijk += size)
    {
        int k = ijk % nk;
        int j = (ijk / nk) % nj;
        int i = (ijk / nk / nj) % ni;
        double x = ((double)(ni - i + 1) * (-a) + (double)(i)*a) / (double)(ni - 1);
        double y = ((double)(nj - j + 1) * (-b) + (double)(j)*b) / (double)(nj - 1);
        double z = ((double)(nk - k + 1) * (-c) + (double)(k)*c) / (double)(nk - 1);
        double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);
        double w_exact = 0.0;
        double wt = 0.0;

        if (1.0 < chk)
        {
            continue;
        }

        ++n_inside;

        w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);

        for (int trial = 0; trial < N; trial++)
        {
            double x1 = x;
            double x2 = y;
            double x3 = z;
            double w = 1.0;
            chk = 0.0;
            while (chk < 1.0)
            {
                double ut = r8_uniform_01(&seed);
                double us;

                double dx;
                if (ut < 1.0 / 3.0)
                {
                    us = r8_uniform_01(&seed) - 0.5;
                    if (us < 0.0)
                        dx = -stepsz;
                    else
                        dx = stepsz;
                }
                else
                    dx = 0.0;

                double dy;
                ut = r8_uniform_01(&seed);
                if (ut < 1.0 / 3.0)
                {
                    us = r8_uniform_01(&seed) - 0.5;
                    if (us < 0.0)
                        dy = -stepsz;
                    else
                        dy = stepsz;
                }
                else
                    dy = 0.0;

                double dz;
                ut = r8_uniform_01(&seed);
                if (ut < 1.0 / 3.0)
                {
                    us = r8_uniform_01(&seed) - 0.5;
                    if (us < 0.0)
                        dz = -stepsz;
                    else
                        dz = stepsz;
                }
                else
                    dz = 0.0;

                double vs = potential(a, b, c, x1, x2, x3);
                x1 = x1 + dx;
                x2 = x2 + dy;
                x3 = x3 + dz;

                double vh = potential(a, b, c, x1, x2, x3);

                double we = (1.0 - h * vs) * w;
                w = w - 0.5 * h * (vh * we + vs * w);

                chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
            }
            wt += w;
        }
        err += pow(w_exact - (wt / (double)(N)), 2);
    }
    double mpi_err = 0.0;
    int mpi_n_inside = 0;
    MPI_Reduce(&err, &mpi_err, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&n_inside, &mpi_n_inside, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    if (rank == MASTER) {
        return sqrt(mpi_err / (double)(mpi_n_inside));
    }
    return 0.0;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Invalid number of arguments passed.\n");
        return 1;
    }
    const int N = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == MASTER)
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        printf("TEST: N=%d, num_procs=%d\n", N, size);
    }
    double wtime = MPI_Wtime();
    double err = feynman(N);
    wtime = MPI_Wtime() - wtime;
    if (rank == MASTER)
    {
        printf("%d    %lf    %lf\n", N, err, wtime);
    }
    MPI_Finalize();

    return 0;
}
