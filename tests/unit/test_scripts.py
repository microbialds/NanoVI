"""Smoke tests for the executables in bin/ that the Nextflow modules call by name (F38).

Nextflow puts bin/ on PATH in every task, so each script must be executable,
start with a python3 shebang and import its sibling modules when run from
another directory.
"""
import os
import subprocess
import sys

import pytest

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))

SCRIPTS = {
    'cigar_profile.py': '--sam',
    'read_scores.py': '--cigar-json',
    'estimate_abundance.py': '--logp',
    'write_abundance_table.py': '--output-prefix',
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
