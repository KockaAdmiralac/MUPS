#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/param.h>
#include <omp.h>
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

// print na stdout upotrebiti u validaciji paralelnog resenja
int main(int argc, char **argv)
{
    const double a = 3.0;
    const double b = 2.0;
    const double c = 1.0;
    const double h = 0.001;
    const double stepsz = sqrt(3 * h);

    int seed = 123456789;

    if (argc < 2)
    {
        printf("Invalid number of arguments passed.\n");
        return 1;
    }
    const int N = atoi(argv[1]);

    int ni;
    int nj;
    int nk;
    if (a == MIN(MIN(a, b), c))
    {
        ni = 6;
        nj = 1 + ceil(b / a) * (ni - 1);
        nk = 1 + ceil(c / a) * (ni - 1);
    }
    else if (b == MIN(MIN(a, b), c))
    {
        nj = 6;
        ni = 1 + ceil(a / b) * (nj - 1);
        nk = 1 + ceil(c / b) * (nj - 1);
    }
    else
    {
        nk = 6;
        ni = 1 + ceil(a / c) * (nk - 1);
        nj = 1 + ceil(b / c) * (nk - 1);
    }

    double err = 0.0;
    int n_inside = 0;

    printf("TEST: N=%d, num_threads=%ld\n", N, get_num_threads());
    double wtime = omp_get_wtime();

#pragma omp parallel default(none) shared(a, b, c, h, stepsz, ni, nj, nk, N) firstprivate(seed) reduction(+                  \
                                                                                                          : err) reduction(+ \
                                                                                                                           : n_inside)
    {
        seed += omp_get_thread_num();
#pragma omp for collapse(3)
        for (int i = 1; i <= ni; i++)
        {
            for (int j = 1; j <= nj; j++)
            {
                for (int k = 1; k <= nk; k++)
                {
                    double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
                    double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
                    double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);
                    double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);

                    if (1.0 < chk)
                    {
#ifdef DEBUG
                        printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
                               x, y, z, 1.0, 1.0, 0.0, 0);
#endif
                        continue;
                    }

                    ++n_inside;

                    double w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);
                    double wt = 0.0;
#ifdef DEBUG
                    int steps = 0;
#endif
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

#ifdef DEBUG
                            ++steps;
#endif

                            double vh = potential(a, b, c, x1, x2, x3);

                            double we = (1.0 - h * vs) * w;
                            w = w - 0.5 * h * (vh * we + vs * w);

                            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
                        }
                        wt = wt + w;
                    }
                    wt = wt / (double)(N);

                    err += pow(w_exact - wt, 2);

#ifdef DEBUG
                    printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
                           x, y, z, wt, w_exact, fabs(w_exact - wt), steps / N);
#endif
                }
            }
        }
    }
    err = sqrt(err / (double)(n_inside));
    wtime = omp_get_wtime() - wtime;

    printf("RMS absolute error in solution = %e\n", err);
    printf("Time: %lf\n", wtime);
    printf("TEST END\n");

    return 0;
}
