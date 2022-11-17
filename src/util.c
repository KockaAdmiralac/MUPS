#include "util.h"
#include <unistd.h>
#include <stdlib.h>

long get_num_threads(void)
{
    long num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    char* num_threads_str = getenv("OMP_NUM_THREADS");
    if (num_threads_str != NULL)
    {
        num_threads = atol(num_threads_str);
    }
    return num_threads;
}
