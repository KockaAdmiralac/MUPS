#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include "util.h"

int prime_number_z1(int n)
{
    int total = 0;

#pragma omp parallel default(none) shared(n) reduction(+ \
                                                       : total)
    {
        int num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();
        for (int i = 2 + thread_id; i <= n; i += num_threads)
        {
            int found = 1;
            for (int j = 2; j < i; j++)
            {
                if ((i % j) == 0)
                {
                    found = 0;
                    break;
                }
            }
            total += found;
        }
    }
    return total;
}

int prime_number_z2(int n)
{
    int total = 0;

#pragma omp parallel default(none) shared(n) reduction(+ \
                                                       : total)
    {
#pragma omp for schedule(static, 1)
        for (int i = 2; i <= n; ++i)
        {
            int found = 1;
            for (int j = 2; j < i; j++)
            {
                if ((i % j) == 0)
                {
                    found = 0;
                    break;
                }
            }
            total += found;
        }
    }
    return total;
}

void test(int (*func)(int), int n_lo, int n_hi, int n_factor)
{
    int n = n_lo;

    while (n <= n_hi)
    {
        double wtime = omp_get_wtime();
        int primes = func(n);
        wtime = omp_get_wtime() - wtime;
        printf("  %8d  %8d  %14f\n", n, primes, wtime);
        n = n * n_factor;
    }
}

int (*FUNCS[])(int) = {prime_number_z1, prime_number_z2};

int main(int argc, char *argv[])
{
    int func = 0;
    int lo;
    int hi;
    int factor;

    if (argc != 5)
    {
        func = 0;
        lo = 1;
        hi = 131072;
        factor = 2;
    }
    else
    {
        func = atoi(argv[1]);
        lo = atoi(argv[2]);
        hi = atoi(argv[3]);
        factor = atoi(argv[4]);
    }

    printf("TEST: func=%d, lo=%d, hi=%d, factor=%d, num_threads=%ld\n", func, lo, hi, factor, get_num_threads());
    test(FUNCS[func], lo, hi, factor);

    return 0;
}
