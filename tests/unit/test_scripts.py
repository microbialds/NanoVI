"""Smoke tests for the executables in bin/ that the Nextflow modules call by name (F38).

Nextflow puts bin/ on PATH in every task, so each script must be executable,
start with a python3 shebang and import its sibling modules when run from
another directory.
"""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from .synthetic import synthetic_reads

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))

SCRIPTS = {
    'read_likelihood.py': '--alignments',
    'estimate_abundance.py': '--scores',
    'write_abundance_table.py': '--read-counts',
    'generate_phyloseq_tables.py': '--output-dir',
    'main.py': 'build-database',
}


@pytest.mark.parametrize('script', sorted(SCRIPTS) + ['python_versions.py'])
def test_script_is_executable_with_python3_shebang(script):
    path = os.path.join(BIN, script)
    assert os.access(path, os.X_OK)
    with open(path) as f:
        assert f.readline().strip() == '#!/usr/bin/env python3'


@pytest.mark.parametrize('script,option', sorted(SCRIPTS.items()))
def test_help_runs_from_another_directory(script, option, tmp_path):
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, script), '--help'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert option in result.stdout


def test_python_versions_prints_a_versions_block(tmp_path):
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'python_versions.py'), 'NANOVI:STEP'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0] == '"NANOVI:STEP":'
    assert lines[1].startswith('    python: "3.')


# ─── estimate_abundance.py (F20) ─────────────────────────────────────────────

def run_estimate(tmp_path, *options):
    lp = synthetic_reads({1: 700, 2: 250, 3: 45, 4: 5}, seed=3)
    sys.path.insert(0, BIN)
    from likelihood import read_scores_to_arrays
    offsets, tax_ids, scores = read_scores_to_arrays(lp)
    npz = tmp_path / 'read_scores.npz'
    np.savez_compressed(npz, offsets=offsets, tax_ids=tax_ids, scores=scores)
    out = tmp_path / 'abundance.json'
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', str(npz),
         '--output', str(out), *options],
        cwd=tmp_path, capture_output=True, text=True,
    )
    return result, (json.loads(out.read_text()) if out.exists() else None), lp


def test_estimate_abundance_defaults_match_the_bayesian_wrapper(tmp_path):
    sys.path.insert(0, BIN)
    from variational_inference_bayesian import variational_inference_iterations
    result, data, lp = run_estimate(tmp_path)
    assert result.returncode == 0, result.stderr
    assert set(data) == {'abundance', 'below_detection_count'}
    lp = json.loads(json.dumps(lp))  # integer IDs as the JSON round trip gives them
    db_ids = list({tid for read in lp.values() for tid in read[0]})
    freq_full, freq_thresh, below = variational_inference_iterations(lp, db_ids, 0.0001, 0.001)
    assert data['abundance'] == {str(k): v for k, v in (freq_thresh or freq_full).items()}
    assert data['below_detection_count'] == below


def test_estimate_abundance_em_and_floor(tmp_path):
    result, data, _ = run_estimate(tmp_path, '--estimator', 'em', '--min-abundance', '0.02')
    assert result.returncode == 0, result.stderr
    # species 4 (5 of 1,000 reads) is below the 2 % floor; its reads go below detection
    assert set(data['abundance']) == {'1', '2', '3'}
    assert sum(data['abundance'].values()) == pytest.approx(1.0)
    assert data['below_detection_count'] > 0


def test_estimate_abundance_credible_intervals(tmp_path):
    result, data, _ = run_estimate(tmp_path, '--credible-intervals')
    assert result.returncode == 0, result.stderr
    assert set(data['posterior']) == set(data['abundance'])
    assert set(data['posterior']['1']) == {'alpha', 'expected_reads', 'shared_reads'}


def test_estimate_abundance_rejects_intervals_with_em(tmp_path):
    result, data, _ = run_estimate(tmp_path, '--estimator', 'em', '--credible-intervals')
    assert result.returncode == 2
    assert '--credible-intervals needs --estimator vb' in result.stderr


# ─── Estimator options removed after the Phase 7a benchmark ─────────────────

def write_npz(path, per_read):
    sys.path.insert(0, BIN)
    from likelihood import read_scores_to_arrays
    offsets, tax_ids, scores = read_scores_to_arrays(per_read)
    np.savez_compressed(path, offsets=offsets, tax_ids=tax_ids, scores=scores)
    return str(path)


def estimate(tmp_path, npz, *options):
    out = tmp_path / 'abundance.json'
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', npz, '--output', str(out), *options],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(out.read_text())


@pytest.mark.parametrize('estimator', ['vb', 'em'])
def test_species_below_the_pruning_threshold_are_pruned(tmp_path, estimator):
    # 2,005 reads: the pruning threshold is 10/2,005; species 2 has 5 reads of its own.
    # --no-pruning, which kept it, lost in the Phase 7a benchmark (D31)
    per_read = {f'a{i}': ([1], [-1.0]) for i in range(2000)}
    per_read.update({f'b{i}': ([2], [-1.0]) for i in range(5)})
    npz = write_npz(tmp_path / 's.npz', per_read)
    data = estimate(tmp_path, npz, '--estimator', estimator, '--min-abundance', '0')
    assert set(data['abundance']) == {'1'} and data['below_detection_count'] == 5
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', npz,
         '--output', str(tmp_path / 'o.json'), '--no-pruning'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert 'unrecognized arguments: --no-pruning' in result.stderr


def test_estimate_abundance_has_no_unknown_component(tmp_path):
    # The unknown-taxon prototype was dropped after the Phase 7a benchmark (D22)
    npz = write_npz(tmp_path / 's.npz', {'a': ([1], [-1.0])})
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', npz,
         '--output', str(tmp_path / 'o.json'), '--unknown-component'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert 'unrecognized arguments: --unknown-component' in result.stderr
