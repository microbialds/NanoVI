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
    'identical_16s.py': '--cross-genus',
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
    assert set(data) == {'abundance', 'below_detection_count', 'read_summary'}
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


@pytest.mark.parametrize('estimator,keys', [('vb', {'alpha', 'expected_reads', 'shared_reads'}),
                                            ('em', {'expected_reads', 'shared_reads'})])
def test_estimate_abundance_writes_the_read_summary(tmp_path, estimator, keys):
    # Always written, for shared_read_fraction (D27); alpha only from VB, for the intervals
    result, data, _ = run_estimate(tmp_path, '--estimator', estimator)
    assert result.returncode == 0, result.stderr
    assert set(data['read_summary']) == set(data['abundance'])
    assert set(data['read_summary']['1']) == keys


def test_estimate_abundance_no_longer_takes_credible_intervals(tmp_path):
    # The interval columns are chosen in write_abundance_table.py (WRITE_OUTPUT)
    result, data, _ = run_estimate(tmp_path, '--credible-intervals')
    assert result.returncode == 2
    assert 'unrecognized arguments: --credible-intervals' in result.stderr


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


# ─── write_abundance_table.py --identical-16s (D24) ─────────────────────────

def test_write_abundance_table_identical_column(tmp_path):
    (tmp_path / 'abundance.json').write_text(json.dumps({'abundance': {'1': 0.7, '2': 0.3}, 'below_detection_count': 0}))
    (tmp_path / 'counts.json').write_text(json.dumps({'assigned_count': 10, 'unassigned_count': 2}))
    (tmp_path / 'taxonomy.tsv').write_text(
        'tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n'
        '1\tStaphylococcus aureus\tStaphylococcus\tF\tO\tC\tP\tBacteria\n'
        '2\tStaphylococcus argenteus\tStaphylococcus\tF\tO\tC\tP\tBacteria\n')
    (tmp_path / 'identical_16s.tsv').write_text(
        'tax_id\tspecies\tgenus\tgroup\tgroup_size\tidentical_16s_species\n'
        '2\tStaphylococcus argenteus\tStaphylococcus\tG1\t2\tStaphylococcus aureus\n'
        '1\tStaphylococcus aureus\tStaphylococcus\tG1\t2\tStaphylococcus argenteus\n')
    argv = [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'abundance.json',
            '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv']
    (tmp_path / 'empty.tsv').write_text('')
    for table, expected in (('identical_16s.tsv', ['Staphylococcus aureus', 'Staphylococcus argenteus']), ('', ['', '']),
                            ('empty.tsv', ['', ''])):
        result = subprocess.run(argv + ['--output-prefix', 'out', '--identical-16s', table],
                                cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        lines = [line.split('\t') for line in (tmp_path / 'out.tsv').read_text().splitlines()]
        assert lines[0][-1] == 'identical_16s_species'
        # rows in name order: S. argenteus, then S. aureus
        assert [line[-1] for line in lines[1:3]] == expected


# ─── shared_read_fraction by default (D27) ──────────────────────────────────

def test_vb_and_em_give_the_same_shared_read_fraction(tmp_path):
    # Species 1 and 2 have 400 reads of their own each and tie on 200 more: each
    # gets about 100 of the tied reads, a shared fraction of 100 / 500
    per_read = {f'a{i}': ([1], [-1.0]) for i in range(400)}
    per_read.update({f'b{i}': ([2], [-1.0]) for i in range(400)})
    per_read.update({f't{i}': ([1, 2], [-1.0, -1.0]) for i in range(200)})
    npz = write_npz(tmp_path / 's.npz', per_read)
    fractions = {}
    for estimator in ('vb', 'em'):
        data = estimate(tmp_path, npz, '--estimator', estimator)
        fractions[estimator] = {t: v['shared_reads'] / v['expected_reads'] for t, v in data['read_summary'].items()}
    assert fractions['vb'] == pytest.approx({'1': 0.2, '2': 0.2}, abs=1e-3)
    assert fractions['em'] == pytest.approx(fractions['vb'], abs=1e-3)


def write_table(tmp_path, data, *options):
    (tmp_path / 'abundance.json').write_text(json.dumps(data))
    (tmp_path / 'counts.json').write_text(json.dumps({'assigned_count': 1000, 'unassigned_count': 0}))
    (tmp_path / 'taxonomy.tsv').write_text(
        'tax_id\tspecies\tgenus\n1\tAlpha one\tAlpha\n2\tAlpha two\tAlpha\n')
    return subprocess.run(
        [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'abundance.json',
         '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv', '--output-prefix', 'out', *options],
        cwd=tmp_path, capture_output=True, text=True,
    )


@pytest.mark.parametrize('estimator', ['vb', 'em'])
def test_write_abundance_table_interval_columns_only_on_request(tmp_path, estimator):
    per_read = {f'a{i}': ([1], [-1.0]) for i in range(600)}
    per_read.update({f't{i}': ([1, 2], [-1.0, -1.0]) for i in range(100)})
    per_read.update({f'b{i}': ([2], [-1.0]) for i in range(300)})
    data = estimate(tmp_path, write_npz(tmp_path / 's.npz', per_read), '--estimator', estimator)
    result = write_table(tmp_path, data)
    assert result.returncode == 0, result.stderr
    header = (tmp_path / 'out.tsv').read_text().splitlines()[0].split('\t')
    assert header[-1] == 'shared_read_fraction' and 'ci_lower' not in header
    result = write_table(tmp_path, data, '--credible-intervals')
    if estimator == 'vb':
        assert result.returncode == 0, result.stderr
        header = (tmp_path / 'out.tsv').read_text().splitlines()[0].split('\t')
        assert header[-3:] == ['shared_read_fraction', 'ci_lower', 'ci_upper']
    else:
        assert result.returncode == 2
        assert '--credible-intervals needs abundances from estimate_abundance.py --estimator vb' in result.stderr


# ─── A sample with no scored read (F45) ─────────────────────────────────────

@pytest.mark.parametrize('estimator', ['vb', 'em'])
def test_estimate_abundance_with_no_scored_read(tmp_path, estimator):
    npz = write_npz(tmp_path / 'S1_read_scores.npz', {})
    out = tmp_path / 'abundance.json'
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', npz, '--output', str(out),
         '--estimator', estimator],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert 'S1_read_scores.npz' in result.stderr
    assert json.loads(out.read_text()) == {'abundance': {}, 'below_detection_count': 0, 'read_summary': {}}


@pytest.mark.parametrize('options', [[], ['--counts'], ['--counts', '--credible-intervals', '--binomial-names']])
def test_write_abundance_table_with_no_scored_read(tmp_path, options):
    # Only the Below detection (0) and Unassigned rows; every read is unassigned
    (tmp_path / 'abundance.json').write_text(json.dumps({'abundance': {}, 'below_detection_count': 0, 'read_summary': {}}))
    (tmp_path / 'counts.json').write_text(json.dumps({'assigned_count': 0, 'unassigned_count': 30}))
    (tmp_path / 'taxonomy.tsv').write_text('tax_id\tspecies\tgenus\n1\tAlpha one\tAlpha\n')
    (tmp_path / 'identical_16s.tsv').write_text('tax_id\tspecies\tgenus\tgroup\tgroup_size\tidentical_16s_species\n')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'abundance.json',
         '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv', '--output-prefix', 'S1_rel-abundance',
         '--identical-16s', 'identical_16s.tsv', *options],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    import pandas as pd
    table = pd.read_csv(tmp_path / 'S1_rel-abundance.tsv', sep='\t', keep_default_na=False)
    assert table['species'].tolist() == ['Below detection', 'Unassigned']
    assert table['abundance'].tolist() == [0, 0]
    if '--counts' in options:
        assert table['estimated_counts'].tolist() == [0, 30]
    if '--credible-intervals' in options:
        assert table['ci_lower'].tolist() == ['', '']

    # The phyloseq tables keep the sample, with zero counts beside a sample with reads
    other = tmp_path / 'S2_rel-abundance.tsv'
    other.write_text('species\tabundance\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\testimated_counts\n'
                     'Alpha one\t1.0\tAlpha\t\t\t\t\t\t10\n'
                     'Below detection\t0\t\t\t\t\t\t\t0\nUnassigned\t0\t\t\t\t\t\t\t2\n')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'generate_phyloseq_tables.py'), 'S1_rel-abundance.tsv', 'S2_rel-abundance.tsv',
         '--output-dir', 'phyloseq'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    abundance = pd.read_csv(tmp_path / 'phyloseq' / 'otu_table_abundance.tsv', sep='\t', index_col=0)
    assert abundance['S1'].tolist() == [0, 0, 0]
    metadata = pd.read_csv(tmp_path / 'phyloseq' / 'sample_metadata.tsv', sep='\t', index_col=0)
    assert metadata.loc['S1', 'n_taxa_detected'] == 0
    if '--counts' in options:
        counts = pd.read_csv(tmp_path / 'phyloseq' / 'otu_table_counts.tsv', sep='\t', index_col=0)
        assert counts.loc['Unassigned', 'S1'] == 30 and counts.loc['Alpha one', 'S1'] == 0
        assert metadata.loc['S1', 'assigned_counts'] == 0 and metadata.loc['S1', 'unassigned_counts'] == 30


# ─── shared_read_fraction per reported row (F46) ────────────────────────────

@pytest.mark.parametrize('estimator', ['vb', 'em'])
@pytest.mark.parametrize('binomial', [False, True])
def test_reads_split_within_one_row_are_not_shared(tmp_path, estimator, binomial):
    # Tax IDs 1 and 2 are strains of one species that split 400 reads 50/50; with
    # --binomial-names (or, without it, because they have the same name) they are one
    # row. Species 3 is a sibling: 200 reads tie between it and both strains and stay shared
    names = ['Alpha one', 'Alpha one str. B' if binomial else 'Alpha one', 'Alpha three']
    (tmp_path / 'taxonomy.tsv').write_text(
        'tax_id\tspecies\tgenus\n' + ''.join(f'{t}\t{n}\tAlpha\n' for t, n in zip([1, 2, 3], names)))
    per_read = {f's{i}': ([1, 2], [-1.0, -1.0]) for i in range(400)}
    per_read.update({f'c{i}': ([3], [-1.0]) for i in range(400)})
    per_read.update({f't{i}': ([1, 2, 3], [-1.0, -1.0, -1.0]) for i in range(200)})
    options = ['--estimator', estimator, '--taxonomy', 'taxonomy.tsv'] + (['--binomial-names'] if binomial else [])
    data = estimate(tmp_path, write_npz(tmp_path / 's.npz', per_read), *options)
    result = write_table_with(tmp_path, data, '--binomial-names' if binomial else None)
    assert result.returncode == 0, result.stderr
    import pandas as pd
    table = pd.read_csv(tmp_path / 'out.tsv', sep='\t').set_index('species')
    # Only the tied reads are shared: about 100 of the merged row's 500 and of species 3's 500
    assert table.loc['Alpha one', 'shared_read_fraction'] == pytest.approx(0.2, abs=0.02)
    assert table.loc['Alpha three', 'shared_read_fraction'] == pytest.approx(0.2, abs=0.02)


def write_table_with(tmp_path, data, option):
    (tmp_path / 'abundance.json').write_text(json.dumps(data))
    (tmp_path / 'counts.json').write_text(json.dumps({'assigned_count': 1000, 'unassigned_count': 0}))
    return subprocess.run(
        [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'abundance.json',
         '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv', '--output-prefix', 'out']
        + ([option] if option else []),
        cwd=tmp_path, capture_output=True, text=True,
    )
