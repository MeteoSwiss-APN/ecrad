#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import argparse
import pickle
import sys
from tests_runner import BlockLengthTestRes, BlockLengthTestsRes


def main():
    parser = argparse.ArgumentParser(description="Collect ECRAD block tests results into 1 file")
    parser.add_argument("--input-dir", type=str, required=True)
    parser.add_argument("--min-size", type=int, default=1)
    parser.add_argument("--max-size", type=int, default=32)
    parser.add_argument("--out-path", type=str, required=True)
    args = parser.parse_args()
    input_dir_path = args.input_dir
    min_size = args.min_size
    max_size = args.max_size
    out_path = args.out_path

    cpu_model_name = None
    num_runs = None
    test_cases = None
    solver_types = None
    results = []

    for i in range(min_size, max_size + 1):
        in_res_file_path = os.path.join(input_dir_path, f'result-{i}.bin')
        assert os.path.exists(in_res_file_path)
        all_res = None
        with open(in_res_file_path, 'rb') as f:
            all_res = pickle.load(f)
        if i == 1:
            cpu_model_name = all_res.cpu_model_name
            num_runs = all_res.num_runs
            solver_types = all_res.solver_types
            test_cases = all_res.test_cases
        else:
            assert cpu_model_name == all_res.cpu_model_name
            assert num_runs == all_res.num_runs
            assert solver_types == all_res.solver_types
            assert test_cases == all_res.test_cases
        assert all_res.min_block_size == i
        assert all_res.max_block_size == i
        results.extend(all_res.results)
    all_res = BlockLengthTestsRes(cpu_model_name,
                                  min_size,
                                  max_size,
                                  num_runs,
                                  test_cases,
                                  solver_types,
                                  results)
    with open(out_path, 'wb') as f:
        s = pickle.dumps(all_res)
        f.write(s)
    return True


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
