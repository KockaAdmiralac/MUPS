#!/usr/bin/env python
from subprocess import Popen, PIPE
from typing import Any, Dict, List
from os import environ as env
from os.path import dirname, realpath, join
from sys import argv, exit, stderr
from matplotlib import pyplot as plt
import numpy as np

Result = List[List[str]]

SCRIPT_DIR = dirname(realpath(__file__))
BUILD_DIR = join(SCRIPT_DIR, 'gen')

TESTS = {
    'prime': {
        'type': 'omp',
        'args': [
            [1, 131072, 2],
            [5, 500000, 10],
            [1, 65536, 4]
        ],
        'funcs': 2,
        'x': lambda result: [int(row[0]) for row in result],
        'y': lambda result, seq_result: [max(float(seq_result[idx][2]), 0.0000001) / max(float(row[2]), 0.0000001) for idx, row in enumerate(result)],
        'same': lambda result1, result2: [int(row[1]) for row in result1] == [int(row[1]) for row in result2],
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman': {
        'type': 'omp',
        'args': [[1000], [5000], [10000], [20000]],
        'funcs': 2,
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= 0.01),
        'threads': [1, 2, 4, 8, 16]
    },
    'moldyn': {
        'type': 'omp',
        'x': lambda result: [int(result[-1][0])],
        'y': lambda result, seq_result: [max(float(seq_result[-1][1]), 0.0000001) / max(float(result[-1][1]), 0.0000001)],
        'same': lambda result1, result2: [float(row[1]) for row in result1[:-1]] == [float(row[1]) for row in result2[:-1]],
        'threads': [1, 2, 4, 8, 16]
    },
    'prime-mpi': {
        'type': 'mpi',
        'args': [
            [1, 131072, 2],
            [5, 500000, 10],
            [1, 65536, 4]
        ],
        'x': lambda result: [int(row[0]) for row in result],
        'y': lambda result, seq_result: [max(float(seq_result[idx][2]), 0.0000001) / max(float(row[2]), 0.0000001) for idx, row in enumerate(result)],
        'same': lambda result1, result2: [int(row[1]) for row in result1] == [int(row[1]) for row in result2],
        'threads': [1, 4]
    },
    'feynman-mpi': {
        'type': 'mpi',
        'args': [[1000], [5000], [10000], [20000]],
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= 0.01),
        'threads': [1, 4]
    },
    'moldyn-mpi': {
        'type': 'mpi',
        'x': lambda result: [int(result[-1][0])],
        'y': lambda result, seq_result: [max(float(seq_result[-1][1]), 0.0000001) / max(float(result[-1][1]), 0.0000001)],
        'same': lambda result1, result2: [float(row[1]) for row in result1[:-1]] == [float(row[1]) for row in result2[:-1]],
        'threads': [1, 4]
    },
    'moldyn-mw-mpi': {
        'type': 'mpi',
        'x': lambda result: [int(result[-1][0])],
        'y': lambda result, seq_result: [max(float(seq_result[-1][1]), 0.0000001) / max(float(result[-1][1]), 0.0000001)],
        'same': lambda result1, result2: [float(row[1]) for row in result1[:-1]] == [float(row[1]) for row in result2[:-1]],
        'threads': [1, 4]
    },
    'prime-cuda': {
        'type': 'cuda',
        'args': [
            [1, 131072, 2],
            [5, 500000, 10],
            [1, 65536, 4]
        ],
        'x': lambda result: [int(row[0]) for row in result],
        'y': lambda result, seq_result: [max(float(seq_result[idx][2]), 0.0000001) / max(float(row[2]), 0.0000001) for idx, row in enumerate(result)],
        'same': lambda result1, result2: [int(row[1]) for row in result1] == [int(row[1]) for row in result2],
        'threads': [1, 2]
    }
}

WIDTH = 1.0

def run_test(func_num: int, test_type: str, exe_name: str, args: List[int], num_threads: int) -> Result:
    process_env = env.copy()
    stringified_args = [str(arg) for arg in args]
    if test_type == 'omp':
        process_env['OMP_NUM_THREADS'] = str(num_threads)
        process_args = [f'{BUILD_DIR}/{exe_name}', str(func_num)]
        log_filename = ' '.join(process_args + stringified_args + [str(num_threads)])
    elif test_type == 'mpi':
        process_args = ['mpiexec', '-np', str(num_threads), f'{BUILD_DIR}/{exe_name}']
        log_filename = f'{exe_name} {" ".join(stringified_args)} {num_threads}'
    elif test_type == 'cuda':
        exe_path = f'{BUILD_DIR}/{exe_name}'
        log_filename = ' '.join([exe_path] + stringified_args + [str(num_threads)])
        if num_threads == 1:
            process_env['OMP_NUM_THREADS'] = '1'
            exe_path = exe_path.replace('-cuda', '')
            stringified_args = ['0'] + stringified_args
        process_args = [exe_path]
    else:
        raise BaseException('Unknown test type.')
    process_args += stringified_args
    process = Popen(process_args, env=process_env, stdout=PIPE)
    if process.wait() != 0 or not process.stdout:
        return []
    results = []
    log_filename = join(BUILD_DIR, f'{log_filename}.log')
    with open(log_filename, 'w', encoding='utf-8') as log_file:
        for line in process.stdout:
            line = line.decode('utf-8')
            log_file.write(line)
            if not line.startswith('TEST'):
                results.append(line.split())
    return results

def run_tests(test_name: str, test_data: Dict[str, Any]):
    print('Running', test_name, 'tests')
    test_data = TESTS[test_name]
    num_funcs = test_data['funcs'] if 'funcs' in test_data else 1
    test_args = test_data['args'] if 'args' in test_data else [[]]
    get_x_axis = test_data['x']
    get_y_axis = test_data['y']
    check_same = test_data['same']
    test_type = test_data['type']
    threads = test_data['threads']
    for func_num in range(num_funcs):
        for arg_idx, args in enumerate(test_args):
            seq_results = []
            x_axis = np.array([])
            x_labels = []
            plt.figure(figsize=(15, 6))
            for num_threads in threads:
                print('Running test with function', func_num, 'arguments', args, 'and', num_threads, 'threads')
                if num_threads == threads[0]:
                    seq_results = run_test(func_num, test_type, test_name, args, num_threads)
                    x_labels = get_x_axis(seq_results)
                    x_axis = np.arange(len(x_labels)) * WIDTH
                else:
                    results = run_test(func_num, test_type, test_name, args, num_threads)
                    if len(results) == 0:
                        print('An error occurred while getting results for ', func_num, args, num_threads, file=stderr)
                        exit(1)
                    if not check_same(seq_results, results):
                        print('Results mismatch for ', func_num, args, num_threads, seq_results, results, file=stderr)
                        print('Test FAILED')
                        exit(2)
                    speedups = get_y_axis(results, seq_results)
                    bar_width = WIDTH / (len(threads) - 1)
                    x_my = x_axis - (WIDTH / 2) + (threads.index(num_threads) - 1) * bar_width + (bar_width / 2)
                    bar = plt.bar(x_my, speedups, label=f'threads={num_threads}', width=bar_width)
                    plt.bar_label(bar, [round(speedup, 1) for speedup in speedups])
            plt.title(f'Results for function index {func_num} and arguments {args}')
            plt.xlabel('$N$')
            plt.ylabel('Speedup')
            plt.xticks(x_axis, x_labels)
            plt.legend()
            plt.savefig(join(BUILD_DIR, f'results-{test_name}-{func_num}-{arg_idx}.svg'))
    print('Test PASSED')

def main():
    if len(argv) > 1:
        test_name = argv[1]
        if test_name not in TESTS:
            print('Invalid test name.')
            exit(3)
        run_tests(test_name, TESTS[test_name])
    else:
        for test_name, test_data in TESTS.items():
            run_tests(test_name, test_data)


if __name__ == '__main__':
    main()
