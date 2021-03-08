#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import argparse
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from tests_runner import BlockLengthTestRes, BlockLengthTestsRes
from collections import OrderedDict as odict
from statistics import mean, stdev

def main():
    parser = argparse.ArgumentParser(description="ECRAD block tests results plotter")
    parser.add_argument("--block-len-test-file-path", type=str, default='block-len-test-results.bin', required=True)
    args = parser.parse_args()
    f_path = args.block_len_test_file_path
    with open(f_path, 'rb') as f:
        all_res = pickle.load(f)
    sizes = [sz for sz in range(all_res.min_block_size, all_res.max_block_size + 1)]
    n_sizes = len(sizes)

    def size_idx(sz: int, min_size: int = all_res.min_block_size):
        return sz - min_size

    timings_by_type = odict()
    timings_by_type_omp = odict()
    timings_stddev_by_type_omp = odict()
    for solver_type in all_res.solver_types:
        timings_by_type[solver_type] = [None] * n_sizes
        timings_by_type_omp[solver_type] = [None] * n_sizes
        timings_stddev_by_type_omp[solver_type] = [None] * n_sizes
    for res in all_res.results:
        ave_time = mean(res.timings)
        std_dev = stdev(res.timings)
        if not res.omp_enabled:
            timings_by_type[res.solver_name][size_idx(res.block_size)] = ave_time / 1000.
        else:
            timings_by_type_omp[res.solver_name][size_idx(res.block_size)] = ave_time / 1000.
            timings_stddev_by_type_omp[res.solver_name][size_idx(res.block_size)] = std_dev / 1000.
    #Plot
    fig, ax = plt.subplots()  # Create a figure and an axes.
    for solver_type, timings in timings_by_type.items():
        ax.plot(sizes, timings, label=solver_type)
    ax.set_xlabel('# columns per block')
    ax.set_ylabel('runtime, ms')
    ax.set_title(f'{all_res.cpu_model_name}, single thread')
    ax.legend()
    ax.set_ylim(ymin=0)
    fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid()

    fig_omp, ax_omp = plt.subplots()  # Create a figure and an axes.
    for solver_type, timings in timings_by_type_omp.items():
        ax_omp.plot(sizes, timings, label=solver_type)
        std_dev = timings_stddev_by_type_omp[solver_type]
        timings_add_stddev = [timings[i] + 3 * std_dev[i] for i in range(len(timings))]
        timings_sub_stddev = [timings[i] - 3 * std_dev[i] for i in range(len(timings))]
        last_color = fig_omp.gca().lines[-1].get_color()
        ax_omp.plot(sizes, timings_add_stddev, '--', color=last_color)
        ax_omp.plot(sizes, timings_sub_stddev, '--', color=last_color)
    ax_omp.set_xlabel('# columns per block')
    ax_omp.set_ylabel('runtime, ms')
    ax_omp.set_title(f'{all_res.cpu_model_name}, openmp')
    ax_omp.legend()
    ax_omp.set_ylim(ymin=0)
    fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid()

    plt.show()
    return True


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
