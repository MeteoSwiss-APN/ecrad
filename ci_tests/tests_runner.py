#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import sys
assert sys.version_info[0] >= 3 and sys.version_info[1] >= 6, 'Python >= 3.6 is required'
import config
from typing import Optional, NamedTuple, List
import subprocess, os, tempfile, argparse
from contextlib import contextmanager
import xarray as xr
import numpy as np
sys.path.append(config.DRIVER_REPORT_DIR)
from driver_report_serialization import driver_report as DriverReport, TestStep as DriverReportStep, \
    parse as parse_driver_report
sys.path.append(config.CI_TESTS_DIR)
from cpu_reg_test_report_serialization import cpu_regression_test_report as CPURegTestReport, \
    TestStep as CPURegTestReportStep, parse as parse_cpu_reg_test_report
from cpuinfo import get_cpu_info
from statistics import mean, stdev

THIS_DIR_PATH = os.path.dirname(os.path.realpath(__file__))
CPU_MODEL_NAME = get_cpu_info()['brand_raw']
TMPFS_DIR = '/dev/shm'


class TerminalColor:
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    RESET = '\033[m'


def green_text(s: str):
    return TerminalColor.GREEN + s + TerminalColor.RESET


def red_text(s: str):
    return TerminalColor.RED + s + TerminalColor.RESET


def file_exists(path: str):
    return os.path.exists(path) and os.path.isfile(path)


@contextmanager
def prepare_dir(dir_path=None, parent_dir=None):
    if dir_path is not None:
        os.makedirs(dir_path, exist_ok=True)
        yield dir_path
    else:
        d = tempfile.TemporaryDirectory(dir=parent_dir)
        yield d.name
        d.cleanup()


class ToleranceTestcaseCFG(NamedTuple):
    openmp_enabled: bool
    block_size: int


DEFAULT_TOL_TEST_CASE_CFG = ToleranceTestcaseCFG(True, 32)


class CPURegressionTestcaseCFG(NamedTuple):
    openmp_enabled: bool
    block_size: int
    num_runs: int
    max_time_rel_tol: float
    avg_time_rel_tolerance: float
    std_dev_time_rel_tolerance: float
    min_relevant_time_us: int


DEFAULT_CPU_REG_TEST_CASE_CFG = CPURegressionTestcaseCFG(False, 32, 10, 0.2, 0.1, 0.1, 50000)


def run_testcase(name: str, solver_type: str, cfg: ToleranceTestcaseCFG, data_multiplier: int = 1,
                 output_file_path: Optional[str] = None, output_report_path: Optional[str] = None,
                 reference_file_path: Optional[str] = None, cmp_abs_tolerance: float = 0.0001,
                 output_cfg_path: Optional[str] = None, working_dir: Optional[str] = None) -> bool:
    in_working_dir = working_dir
    in_file_path = os.path.join(config.TEST_DATA_INPUT_DIR, name + '.nc')
    with prepare_dir(in_working_dir, TMPFS_DIR) as working_dir:
        basename = f'{name}_{solver_type}_{data_multiplier}'
        if output_file_path is None:
            output_file_path = os.path.join(working_dir, f'{basename}.nc')
        if output_report_path is None:
            output_report_path = os.path.join(working_dir, f'{basename}_report.xml')
        if output_cfg_path is None:
            output_cfg_path = os.path.join(working_dir, f'{basename}_config.nam')
        assert file_exists(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH), \
            f'Input cfg template "{config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH}" not found'
        cfg_file_txt = None
        with open(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH, 'r') as f:
            cfg_file_txt = f.read()
        cfg_file_txt = cfg_file_txt.format(openmp_enabled=cfg.openmp_enabled,
                                           block_size=cfg.block_size,
                                           output_report_path=output_report_path,
                                           output_cfg_path=output_cfg_path,
                                           cfg_dir_path=config.RADIATION_CFG_DIR,
                                           sw_solver_type=solver_type,
                                           lw_solver_type=solver_type)
        with open(output_cfg_path, 'w') as f:
            f.write(cfg_file_txt)
        args = (config.DRIVER_BIN, output_cfg_path, in_file_path, output_file_path)
        p = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=working_dir)

        def on_fail(reason: str, p=p, args=args):
            print(f'ECRAD driver run failed! Reason: {reason}')
            print(f'args: {args}')
            stderr = p.stderr.decode('utf-8').lower()
            print(f'stderr: {stderr}')

        if p.returncode != 0:
            on_fail('return code not zero')
            return False
        if not file_exists(output_file_path):
            on_fail('output file missing')
            return False
        if not file_exists(output_report_path):
            on_fail('report missing')
            return False
        if reference_file_path is not None:
            if not file_exists(reference_file_path):
                print(f'Reference file {reference_file_path} not found')
                return False
            if not nc_equal(reference_file_path, output_file_path, cmp_abs_tolerance):
                print(f'Output file {output_file_path} is not equal to {reference_file_path} within absolute tolerance'+
                      f' {cmp_abs_tolerance}')
                return False
        return True


def get_step_duration(step: DriverReportStep) -> int:
    return step.get_end_timestamp_us() - step.get_start_timestamp_us()


def run_cpu_regression_testcase(name: str, solver_type: str, out_report_path: str, cfg: CPURegressionTestcaseCFG,
                                data_multiplier: int = 1, reference_file_path: Optional[str] = None,
                                cmp_abs_tolerance: float = 0.0001, working_dir: Optional[str] = None) -> bool:
    in_working_dir = working_dir
    in_file_path = os.path.join(config.TEST_DATA_INPUT_DIR, name + '.nc')
    with prepare_dir(in_working_dir, TMPFS_DIR) as working_dir:
        basename = f'{name}_{solver_type}_{data_multiplier}'
        output_file_path = os.path.join(working_dir, f'{basename}.nc')
        driver_report_path = os.path.join(working_dir, f'{basename}_drv_report.xml')
        driver_cfg_path = os.path.join(working_dir, f'{basename}_config.nam')
        assert file_exists(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH), \
            f'Input cfg template "{config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH}" not found'
        cfg_file_txt = None
        with open(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH, 'r') as f:
            cfg_file_txt = f.read()
        cfg_file_txt = cfg_file_txt.format(openmp_enabled=cfg.openmp_enabled,
                                           block_size=cfg.block_size,
                                           output_report_path=driver_report_path,
                                           output_cfg_path=driver_cfg_path,
                                           cfg_dir_path=config.RADIATION_CFG_DIR,
                                           sw_solver_type=solver_type,
                                           lw_solver_type=solver_type)
        with open(driver_cfg_path, 'w') as f:
            f.write(cfg_file_txt)
        args = (config.DRIVER_BIN, driver_cfg_path, in_file_path, output_file_path)
        step_name_by_idx = []
        step_durations_by_idx = []
        for i in range(0, cfg.num_runs):
            p = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=working_dir)

            def on_fail(reason: str, p=p, args=args):
                print(f'ECRAD driver run failed! Reason: {reason}')
                print(f'args: {args}')
                if p.stderr is not None:
                    stderr = p.stderr.decode('utf-8').lower()
                    print(f'stderr: {stderr}')

            if p.returncode != 0:
                on_fail('return code not zero')
                return False
            if not file_exists(output_file_path):
                on_fail('output file missing')
                return False
            if not file_exists(driver_report_path):
                on_fail('report missing')
                return False
            if reference_file_path is not None:
                if not file_exists(reference_file_path):
                    print(f'Reference file {reference_file_path} not found')
                    return False
                if not nc_equal(reference_file_path, output_file_path, cmp_abs_tolerance):
                    print(f'Output file {output_file_path} is not equal to {reference_file_path} within absolute tolerance'+
                          f' {cmp_abs_tolerance}')
                    return False

            drv_report = parse_driver_report(driver_report_path)
            if i != 0:
                if len(drv_report.get_step()) != len(step_name_by_idx):
                    on_fail('Run {i}: different number of steps')
                    return False

            for si in range(0, len(drv_report.get_step())):
                step = drv_report.get_step()[si]
                step_name = step.get_name()
                step_result = step.get_result()
                if not step_result:
                    on_fail(f'Run {i}: Step "{step_name}" failed')
                    return False
                step_duration = get_step_duration(step)
                if i == 0:
                    step_name_by_idx.append(step_name)
                    step_durations_by_idx.append([])
                if step_name != step_name_by_idx[si]:
                    on_fail(f'Run {i}: Step "{step_name}" unexpected')
                    return False
                step_durations_by_idx[si].append(step_duration)
        num_steps = len(step_durations_by_idx)
        step_avg_duration_by_idx = [[] for i in range(num_steps)]
        step_max_duration_by_idx = [[] for i in range(num_steps)]
        step_stddev_duration_by_idx = [[] for i in range(num_steps)]
        for i in range(0, num_steps):
            step_durations = step_durations_by_idx[i]
            step_avg_duration_by_idx[i] = mean(step_durations)
            step_max_duration_by_idx[i] = max(step_durations)
            step_stddev_duration_by_idx[i] = stdev(step_durations)
        reg_report = CPURegTestReport(CPU_MODEL_NAME, cfg.openmp_enabled)
        for i in range(0, num_steps):
            reg_report.add_step(CPURegTestReportStep(step_name_by_idx[i], step_avg_duration_by_idx[i],
                                                     step_stddev_duration_by_idx[i], step_max_duration_by_idx[i]))
        with open(out_report_path, 'w') as f:
            reg_report.export(f, 0)
    return True


def tol_test_case_ref_file_path(test_case_name: str, solver_type: str):
    return os.path.join(config.TOLERANCE_TEST_DATA_REFERENCE_DIR, f'{test_case_name}_{solver_type}_result.nc')


def cpu_reg_test_case_ref_report_path(test_case_name: str, solver_type: str):
    return os.path.join(config.CPU_REGRESSION_TEST_DATA_REFERENCE_DIR, f'{test_case_name}_{solver_type}_report.xml')


def reset_ci_tests():
    reset_ci_tol_tests()
    reset_ci_cpu_regression_tests()


def reset_ci_tol_tests():
    print("Resetting CI tolerance tests")
    for test_case_name in config.TEST_CASES:
        for solver_type in config.SOLVER_TYPES:
            out_file_path = tol_test_case_ref_file_path(test_case_name, solver_type)
            print(f'\tGenerating tolerance reference file for test case: {test_case_name} solver: {solver_type}')
            run_testcase(test_case_name, solver_type, DEFAULT_TOL_TEST_CASE_CFG, output_file_path=out_file_path)


def reset_ci_cpu_regression_tests():
    print("Resetting CI CPU regression tests")
    for test_case_name in config.TEST_CASES:
        for solver_type in config.SOLVER_TYPES:
            output_report_path = cpu_reg_test_case_ref_report_path(test_case_name, solver_type)
            print(f'\tGenerating reference report file for test case: {test_case_name} solver: {solver_type}')
            run_cpu_regression_testcase(test_case_name, solver_type, output_report_path, DEFAULT_CPU_REG_TEST_CASE_CFG)


def run_ci_tests() -> bool:
    all_res = True
    tol_res = run_ci_tol_tests()
    if not tol_res:
        all_res = False
        print(red_text('CI tolerance tests failed'))
    cpu_reg_res = run_ci_cpu_regression_tests()
    if not cpu_reg_res:
        all_res = False
        print(red_text('CI tolerance tests failed'))
    return all_res


def run_ci_cpu_regression_tests() -> bool:
    print('CPU regression tests:')
    all_tests_res = True
    with prepare_dir(parent_dir=TMPFS_DIR) as working_dir:
        for test_case_name in config.TEST_CASES:
            for solver_type in config.SOLVER_TYPES:
                ref_file_path = tol_test_case_ref_file_path(test_case_name, solver_type)
                ref_report_path = cpu_reg_test_case_ref_report_path(test_case_name, solver_type)
                report_path = os.path.join(working_dir, f'{test_case_name}_{solver_type}.xml')
                print(f'\t{test_case_name} solver: {solver_type} {green_text("STARTED")}')
                res = run_cpu_regression_testcase(test_case_name, solver_type, report_path, DEFAULT_CPU_REG_TEST_CASE_CFG,
                                            reference_file_path=ref_file_path)
                if not res:
                    all_tests_res = False
                else:
                    res = cpu_reg_report_equal(ref_report_path, report_path, DEFAULT_CPU_REG_TEST_CASE_CFG)
                    if not res:
                        all_tests_res = False
                print(f'\t{test_case_name} solver: {solver_type} {green_text("SUCCESS") if res else red_text("FAILED")}')
    return all_tests_res


def run_ci_tol_tests() -> bool:
    print('Tolerance tests:')
    all_tests_res = True
    for test_case_name in config.TEST_CASES:
        for solver_type in config.SOLVER_TYPES:
            ref_file_path = tol_test_case_ref_file_path(test_case_name, solver_type)
            print(f'\t{test_case_name} solver: {solver_type} {green_text("STARTED")}')
            res = run_testcase(test_case_name, solver_type, DEFAULT_TOL_TEST_CASE_CFG, reference_file_path=ref_file_path)
            if not res:
                all_tests_res = False
            print(f'\t{test_case_name} solver: {solver_type} {green_text("SUCCESS") if res else red_text("FAILED")}')
    return all_tests_res


def nc_equal(ref_filename: str, in_filename: str, abs_threshold: float=0.0001) -> bool:
    ref = xr.open_dataset(ref_filename)
    input = xr.open_dataset(in_filename)
    # Compare dimensions
    ref_dim_names = sorted(ref.dims.keys())
    in_dim_names = sorted(input.dims.keys())
    if ref_dim_names != in_dim_names:
        return False
    for dim_name in ref_dim_names:
        if ref.dims[dim_name] != input.dims[dim_name]:
            return False
    # Compare data
    ref_var_names = sorted(ref.data_vars.keys())
    in_var_names = sorted(input.data_vars.keys())
    if ref_var_names != in_var_names:
        return False
    for var_name in ref_var_names:
        ref_var_data = ref.data_vars[var_name]
        in_var_data = input.data_vars[var_name]
        if ref_var_data.name != in_var_data.name:
            return False
        if ref_var_data.dims != in_var_data.dims:
            return False
        if ref_var_data.sizes != in_var_data.sizes:
            return False
        if not np.allclose(ref_var_data, in_var_data, atol=abs_threshold, equal_nan=True):
            return False
    return True


def rel_diff(val: int, ref_val: int) -> float:
    rel_change = float(ref_val - val) / float(ref_val)
    return rel_change


def cpu_reg_report_equal(reference_report_path: str, report_path: str, report_cfg: CPURegressionTestcaseCFG) -> bool:
    report = parse_cpu_reg_test_report(report_path)
    ref_report = parse_cpu_reg_test_report(reference_report_path)
    cpu = report.get_cpu_model()
    ref_cpu = ref_report.get_cpu_model()
    if cpu != ref_cpu:
        print('Warning! Testcase was run on a different system, comparison may not be valid:\n' +
              f'\tExpected: {ref_cpu}, Got: {cpu}')
    if len(report.get_step()) != len(ref_report.get_step()):
        return False
    num_steps = len(report.get_step())
    for i in range(0, num_steps):
        step = report.get_step()[i]
        ref_step = ref_report.get_step()[i]
        if step.get_name() != ref_step.get_name():
            return False
        step_name = ref_step.get_name()
        step_mean = step.get_avg_duration_us()
        ref_step_mean = ref_step.get_avg_duration_us()
        if step_mean < report_cfg.min_relevant_time_us and ref_step_mean < report_cfg.min_relevant_time_us:
            continue  # There is no point in comparison if durations are too short
        step_mean_rel_diff = rel_diff(step_mean, ref_step_mean)
        abs_step_mean_rel_diff = abs(step_mean_rel_diff)
        if abs_step_mean_rel_diff >= report_cfg.avg_time_rel_tolerance:
            if step_mean_rel_diff < 0:
                print(red_text(f'Step "{step_name}" is now {int(abs_step_mean_rel_diff * 100.)}% slower on average'))
                return False
            else:
                print(green_text(f'Step "{step_name}" is now {int(abs_step_mean_rel_diff * 100.)}% faster on average'))
        step_max = step.get_max_duration_us()
        ref_step_max = ref_step.get_max_duration_us()
        step_max_rel_diff = rel_diff(step_max, ref_step_max)
        abs_step_max_rel_diff = abs(step_max_rel_diff)
        if abs_step_max_rel_diff >= report_cfg.max_time_rel_tol:
            if step_max_rel_diff < 0:
                print(red_text(f'Step "{step_name}" is now {int(abs_step_max_rel_diff * 100.)}% slower on maximum'))
                return False
            else:
                print(green_text(f'Step "{step_name}" is now {int(abs_step_max_rel_diff * 100.)}% faster on maximum'))
        step_stddev = step.get_std_dev_duration_us()
        step_rel_stddev = float(step_stddev) / float(ref_step_mean)
        if step_rel_stddev >= report_cfg.std_dev_time_rel_tolerance:
            print(red_text(f'Step "{step_name}" duration now varies too much between runs. '
                           f'Stddev exceeds allowed {report_cfg.std_dev_time_rel_tolerance} from mean '))
            return False
    return True


def report_equal(reference_report_path: str, report_path: str, time_rel_tol: float = 0.1,
                 time_min_duration_us=5000) -> bool:
    report = parse_driver_report(report_path)
    ref_report = parse_driver_report(reference_report_path)
    if len(report.get_step()) != len(ref_report.get_step()):
        return False
    num_steps = len(report.get_step())
    for i in range(0, num_steps):
        step = report.get_step()[i]
        ref_step = ref_report.get_step()[i]
        if step.get_name() != ref_step.get_name():
            return False
        step_name = ref_step.get_name()
        step_duration = get_step_duration(step)
        ref_step_duration = get_step_duration(ref_step)
        if step_duration < time_min_duration_us and ref_step_duration < time_min_duration_us:
            continue  # There is no point in comparison if durations are too short
        rel_change = float(ref_step_duration - step_duration) / float(ref_step_duration)
        rel_change_abs = abs(rel_change)
        if rel_change_abs >= time_rel_tol:
            if rel_change < 0:
                print(red_text(f'Step "{step_name}" is now {int(rel_change_abs * 100.)}% slower'))
                return False
            else:
                print(green_text(f'Step "{step_name}" is now {int(rel_change_abs * 100.)}% faster'))
    return True


def multiply_dataset(in_filename: str, out_filename: str, n_times: int):
    """Enlarges dataset by copying it n_times"""
    assert False, 'Implementation not finished'
    assert n_times > 0
    data = xr.open_dataset(in_filename)
    out = xr.Dataset()

    new_dims = {key: val for key, val in data.dims.items() if key != "column"}

    out = out.expand_dims(new_dims)
    out = out.expand_dims({"column": n})

    lon = data["longitude"][:]
    lat = data["latitude"][:]

    out["longitude"] = (("column",), np.linspace(np.min(lon), np.max(lon), n))
    out["latitude"] = (("column",), np.linspace(np.min(lat), np.max(lat), n))

    n_repeat = n // data.dims["column"] + 1

    for vn in data.data_vars:
        var = data[vn]
        if vn not in ["longitude", "latitude"] and "column" in var.dims:
            colax = var.dims.index("column")
            out[vn] = (var.dims, np.repeat(var[:], n_repeat, axis=colax)[:n])

    out.to_netcdf(out_filename)


def run_timing_testcase(name: str, solver_type: str, block_size: int, omp_enabled: bool, num_runs: int,
                        out_timings: List[int], timing_step_name: str = 'computation',
                        reference_file_path: Optional[str] = None,  cmp_abs_tolerance: float = 0.0001,
                        working_dir: Optional[str] = None) -> bool:
    in_working_dir = working_dir
    in_file_path = os.path.join(config.TEST_DATA_INPUT_DIR, name + '.nc')
    with prepare_dir(in_working_dir, TMPFS_DIR) as working_dir:
        basename = f'{name}_{solver_type}'
        output_file_path = os.path.join(working_dir, f'{basename}.nc')
        driver_report_path = os.path.join(working_dir, f'{basename}_drv_report.xml')
        driver_cfg_path = os.path.join(working_dir, f'{basename}_config.nam')
        assert file_exists(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH), \
            f'Input cfg template "{config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH}" not found'
        cfg_file_txt = None
        with open(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH, 'r') as f:
            cfg_file_txt = f.read()
        cfg_file_txt = cfg_file_txt.format(openmp_enabled=omp_enabled,
                                           block_size=block_size,
                                           output_report_path=driver_report_path,
                                           output_cfg_path=driver_cfg_path,
                                           cfg_dir_path=config.RADIATION_CFG_DIR,
                                           sw_solver_type=solver_type,
                                           lw_solver_type=solver_type)
        with open(driver_cfg_path, 'w') as f:
            f.write(cfg_file_txt)
        args = (config.DRIVER_BIN, driver_cfg_path, in_file_path, output_file_path)
        step_name_by_idx = []
        step_durations_by_idx = []
        for i in range(0, num_runs):
            p = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=working_dir)

            def on_fail(reason: str, p=p, args=args):
                print(f'ECRAD driver run failed! Reason: {reason}')
                print(f'args: {args}')
                if p.stderr is not None:
                    stderr = p.stderr.decode('utf-8').lower()
                    print(f'stderr: {stderr}')

            if p.returncode != 0:
                on_fail('return code not zero')
                return False
            if not file_exists(output_file_path):
                on_fail('output file missing')
                return False
            if not file_exists(driver_report_path):
                on_fail('report missing')
                return False
            if reference_file_path is not None:
                if not file_exists(reference_file_path):
                    print(f'Reference file {reference_file_path} not found')
                    return False
                if not nc_equal(reference_file_path, output_file_path, cmp_abs_tolerance):
                    print(f'Output file {output_file_path} is not equal to {reference_file_path} within absolute tolerance'+
                          f' {cmp_abs_tolerance}')
                    return False

            drv_report = parse_driver_report(driver_report_path)
            if i != 0:
                if len(drv_report.get_step()) != len(step_name_by_idx):
                    on_fail('Run {i}: different number of steps')
                    return False

            for si in range(0, len(drv_report.get_step())):
                step = drv_report.get_step()[si]
                step_name = step.get_name()
                step_result = step.get_result()
                if not step_result:
                    on_fail(f'Run {i}: Step "{step_name}" failed')
                    return False
                step_duration = get_step_duration(step)
                if i == 0:
                    step_name_by_idx.append(step_name)
                    step_durations_by_idx.append([])
                if step_name != step_name_by_idx[si]:
                    on_fail(f'Run {i}: Step "{step_name}" unexpected')
                    return False
                step_durations_by_idx[si].append(step_duration)
        out_timings.clear()
        step_name_idx = step_name_by_idx.index(timing_step_name)
        out_timings.extend(step_durations_by_idx[step_name_idx])
    return True


class BlockLengthTestRes(NamedTuple):
    block_size: int
    test_case_name: str
    solver_name: str
    omp_enabled: bool
    timings: List[int]


class BlockLengthTestsRes(NamedTuple):
    cpu_model_name: str
    min_block_size: int
    max_block_size: int
    num_runs: int
    test_cases: List[str]
    solver_types: List[str]
    results: List[BlockLengthTestRes]

import pickle

def run_block_length_tests(out_file_path: str, min_block_size: int = 1, max_block_size: int = 2, n_runs: int = 1):
    print('Block length tests:')
    all_tests_res = True
    res_data = []
    for block_size in range(min_block_size, max_block_size + 1):
        for test_case_name in config.TEST_CASES:
            for solver_type in config.SOLVER_TYPES:
                print(f'\t{test_case_name} solver: {solver_type} {green_text("STARTED")}')
                for omp_enabled in (False, True):
                    for run_idx in range(1, n_runs + 1):
                        timings = []
                        res = run_timing_testcase(test_case_name, solver_type, block_size, omp_enabled, n_runs, timings)
                        if res:
                            res_data.append(BlockLengthTestRes(block_size, test_case_name, solver_type, omp_enabled, timings))
                        else:
                            all_tests_res = False
                print(f'\t{test_case_name} solver: {solver_type} {green_text("SUCCESS") if res else red_text("FAILED")}')
    all_res = BlockLengthTestsRes(cpu_model_name=CPU_MODEL_NAME,
                                  min_block_size=min_block_size,
                                  max_block_size=max_block_size,
                                  num_runs=n_runs,
                                  test_cases=config.TEST_CASES,
                                  solver_types=config.SOLVER_TYPES,
                                  results=res_data)
    with open(out_file_path, 'wb') as f:
        s = pickle.dumps(all_res)
        f.write(s)


def main() -> bool:
    parser = argparse.ArgumentParser(description="ECRAD automated tests runner")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--reset-ci-tests", help="Generate reference results for ci testcases", action='store_true')
    group.add_argument("--reset-ci-tol-tests", help="Generate reference results for ci tolerance testcases", action='store_true')
    group.add_argument("--reset-cpu-reg-tests", help="Generate reference results for CPU regression tests", action='store_true')
    group.add_argument("--run-ci-tests", help="Run continuous integration tests", action='store_true')
    group.add_argument("--run-ci-tol-tests", help="Run continuous integration tolerance tests", action='store_true')
    group.add_argument("--run-cpu-reg-tests", help="Run CPU regression tests", action='store_true')
    group.add_argument("--run-block-len-test", help="Find optimal block length", action='store_true')
    parser.add_argument("--block-len-test-min-size", type=int, default=1)
    parser.add_argument("--block-len-test-max-size", type=int, default=64)
    parser.add_argument("--block-len-test-num-runs", type=int, default=10)
    parser.add_argument("--block-len-test-out-file-name", type=str, default='block-len-test-results.bin')
    args = parser.parse_args()
    if args.reset_ci_tests:
        return reset_ci_tests()
    elif args.reset_ci_tol_tests:
        return reset_ci_tol_tests()
    elif args.reset_cpu_reg_tests:
        return reset_ci_cpu_regression_tests()
    elif args.run_ci_tests:
        return run_ci_tests()
    elif args.run_ci_tol_tests:
        return run_ci_tol_tests()
    elif args.run_cpu_reg_tests:
        return run_ci_cpu_regression_tests()
    elif args.run_block_len_test:
        min_size = args.block_len_test_min_size
        max_size = args.block_len_test_max_size
        num_runs = args.block_len_test_num_runs
        out_file_path = args.block_len_test_out_file_name
        if not os.path.isabs(out_file_path):
            out_file_path = os.path.join(THIS_DIR_PATH, out_file_path)
        run_block_length_tests(out_file_path, min_size, max_size, num_runs)
    return True


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
