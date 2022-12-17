#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define MASTER 0

int prime_number(int n)
{
    int total = 0;
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (int i = 2 + rank; i <= n; i += size)
    {
        int prime = 1;
        for (int j = 2; j < i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total = total + prime;
    }
    int mpi_total = 0;
    MPI_Reduce(&total, &mpi_total, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    return mpi_total;
}

int main(int argc, char *argv[])
{
    int factor;
    int hi;
    int lo;
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
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == MASTER)
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        printf("TEST: lo=%d, hi=%d, factor=%d, num_procs=%d\n", lo, hi, factor, size);
    }
    int n = lo;
    while (n <= hi)
    {
        double ctime = MPI_Wtime();
        int primes = prime_number(n);
        ctime = MPI_Wtime() - ctime;
        if (rank == MASTER)
        {
            printf("  %8d  %8d  %14f\n", n, primes, ctime);
        }
        n *= factor;
    }
    MPI_Finalize();
    return 0;
}
