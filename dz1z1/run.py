#!/usr/bin/env python
from subprocess import Popen, PIPE
from typing import List
from os import environ as env
from sys import exit, stderr
from matplotlib import pyplot as plt
import numpy as np

Result = List[List[str]]

ARGS = [
    [1, 131072, 2],
    [5, 500000, 10],
    [1, 65536, 4]
]

NUM_FUNCS = 2

THREADS = [1, 2, 4, 8]

WIDTH = 0.9

def test(func_num: int, args: List[int], num_threads: int) -> Result:
    process_env = env.copy()
    process_env['OMP_NUM_THREADS'] = str(num_threads)
    process_args = ['./prime', str(func_num)] + [str(arg) for arg in args]
    process = Popen(process_args, env=process_env, stdout=PIPE)
    if process.wait() != 0 or not process.stdout:
        return []
    results = []
    log_filename = f'{" ".join(process_args + [str(num_threads)])}.log'
    with open(log_filename, 'w', encoding='utf-8') as log_file:
        for line in process.stdout:
            line = line.decode('utf-8')
            log_file.write(line)
            if not line.startswith('TEST'):
                results.append(line.split())
    return results

def check_same(result1: Result, result2: Result) -> bool:
    primes1 = [int(row[1]) for row in result1]
    primes2 = [int(row[1]) for row in result2]
    return primes1 == primes2

def get_y_axis(result: Result, seq_result: Result) -> List[float]:
    return [float(seq_result[idx][2]) / float(row[2]) for idx, row in enumerate(result)]

def get_x_axis(result: Result) -> List[int]:
    return [int(row[0]) for row in result]

def main():
    for func_num in range(NUM_FUNCS):
        for arg_idx, args in enumerate(ARGS):
            seq_results = []
            x_axis = np.array([])
            x_labels = []
            plt.figure(figsize=(15, 6))
            for num_threads in THREADS:
                if num_threads == 1:
                    seq_results = test(func_num, args, num_threads)
                    x_labels = get_x_axis(seq_results)
                    x_axis = np.arange(len(x_labels)) * WIDTH
                else:
                    results = test(func_num, args, num_threads)
                    if len(results) == 0:
                        print('An error occurred while getting results for ', func_num, args, num_threads, file=stderr)
                        exit(1)
                    if not check_same(seq_results, results):
                        print('Results mismatch for ', func_num, args, num_threads, seq_results, results, file=stderr)
                        exit(2)
                    speedups = get_y_axis(results, seq_results)
                    bar_width = WIDTH / (len(THREADS) - 1)
                    x_my = x_axis - (WIDTH / 2) + (THREADS.index(num_threads) - 1) * bar_width + (bar_width / 2)
                    bar = plt.bar(x_my, speedups, label=f'threads={num_threads}', width=bar_width)
                    plt.bar_label(bar, [round(speedup, 1) for speedup in speedups])
            plt.title(f'Results for function index {func_num} and arguments {args}')
            plt.xlabel('$N$')
            plt.ylabel('Speedup')
            plt.xticks(x_axis, x_labels)
            plt.legend()
            plt.savefig(f'results-{func_num}-{arg_idx}.svg')


if __name__ == '__main__':
    main()
