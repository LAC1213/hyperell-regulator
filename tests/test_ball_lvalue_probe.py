"""Optional integration checks for the compiled, values-only L-value probe.

Set HYPERELL_BALL_LVALUE_PROBE and HYPERELL_BALL_LVALUE_FACTORS to enable.
The factors must cover primes <= 53443 for the standard motive of 249.a.249.1.
HYPERELL_BALL_LVALUE_FULL_CONTROL optionally enables the upstream control run.
"""

from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess

import pytest


@pytest.fixture
def probe_paths():
    binary = os.environ.get('HYPERELL_BALL_LVALUE_PROBE')
    factors = os.environ.get('HYPERELL_BALL_LVALUE_FACTORS')
    if not binary or not factors:
        pytest.skip('optional compiled L-value probe and prime-factor input not configured')
    return str(Path(binary).resolve()), str(Path(factors).resolve())


def run_probe(binary, factors, directory):
    return subprocess.run([binary, '80', '62001', factors], cwd=directory,
                          text=True, capture_output=True, timeout=60)


def dyadic_interval(output, label):
    line = next(row for row in output.splitlines() if row.startswith(label + ' '))
    _, lo, hi, exponent = line.split()
    scale = Fraction(2) ** int(exponent)
    return int(lo) * scale, int(hi) * scale


def test_low_target_rejects_insufficient_final_accuracy(probe_paths, tmp_path):
    binary, factors = probe_paths
    result = run_probe(binary, factors, tmp_path)
    assert result.returncode == 9, result.stderr
    assert 'computation values_only' in result.stdout
    assert 'L_second_derivative ' not in result.stdout
    recorded = json.loads((Path(__file__).parents[1] / 'reports/ball_lvalue_249_80.json').read_text())
    assert recorded['status'] == 'withdrawn_error_bounds'


def test_error_bound_regression_binary(tmp_path):
    binary = os.environ.get('HYPERELL_BALL_ERROR_TEST')
    if not binary:
        pytest.skip('optional compiled error-bound regression not configured')
    result = subprocess.run([str(Path(binary).resolve()), '400'], cwd=tmp_path,
                            text=True, capture_output=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'independent FLINT enclosure' in result.stdout
    assert 'all checked Fourier entries' in result.stdout


def test_incomplete_prime_input_is_rejected(probe_paths, tmp_path):
    binary, _ = probe_paths
    incomplete = tmp_path / 'incomplete.txt'
    incomplete.write_text('2 1\n')
    result = run_probe(binary, str(incomplete), tmp_path)
    assert result.returncode == 5
    assert 'L_second_derivative ' not in result.stdout


def test_full_pipeline_still_reports_zero_count_warning(probe_paths, tmp_path):
    binary = os.environ.get('HYPERELL_BALL_LVALUE_FULL_CONTROL')
    if not binary:
        pytest.skip('optional full-pipeline regression control not configured')
    _, factors = probe_paths
    result = run_probe(str(Path(binary).resolve()), factors, tmp_path)
    assert result.returncode == 6
    assert 'Found too many zeros' in result.stderr
    assert 'L_second_derivative ' not in result.stdout
