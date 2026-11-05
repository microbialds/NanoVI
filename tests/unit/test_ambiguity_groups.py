"""Ambiguity groups: rows of the abundance table the sample cannot separate (F57, D39)."""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from .test_scripts import BIN, estimate, write_npz

sys.path.insert(0, BIN)
from ambiguity_groups import ambiguity_groups, link_groups, pairwise_shared, tax_families  # noqa: E402

# Four reads over three reported columns (-1: padding)
IDX = np.array([[0, 1], [0, -1], [0, 2], [1, 2]])
PHI = np.array([[0.5, 0.5], [1.0, 0.0], [0.995, 0.005], [0.9, 0.1]])
ONE_FAMILY = {11: 'Fam1', 12: 'Fam1', 13: 'Fam1'}


def test_pairwise_shared_counts_the_reads_split_with_each_partner():
    expected, shared = pairwise_shared(IDX, PHI, [0, 1, 2], 3)
    assert expected == pytest.approx([2.495, 1.4, 0.105])
    # shared[a, b]: a's responsibility on the reads that b takes at least 0.01 of
    assert shared == pytest.approx(np.array([[0.0, 0.5, 0.0],
                                             [0.5, 0.0, 0.9],
                                             [0.005, 0.1, 0.0]]))


def test_pairwise_shared_judges_the_split_among_reported_rows():
    # A read whose responsibility is mostly on species not reported: its 0.02 on column 1
    # is 0.2 of what the reported rows take, so it is a split
    expected, shared = pairwise_shared(np.array([[0, 1]]), np.array([[0.08, 0.02]]), [0, 1], 2)
    assert expected == pytest.approx([0.08, 0.02])
    assert shared[0, 1] == pytest.approx(0.08) and shared[1, 0] == pytest.approx(0.02)


def test_link_groups_uses_the_larger_share_and_chains_links():
    expected = np.array([100.0, 10.0, 10.0, 50.0])
    shared = np.zeros((4, 4))
    shared[0, 1], shared[1, 0] = 6.0, 6.0   # f_0(1) = 0.06, f_1(0) = 0.6
    shared[1, 2], shared[2, 1] = 3.0, 3.0   # 0.3 both ways
    shared[0, 3], shared[3, 0] = 2.0, 2.0   # 0.02 and 0.04
    groups, max_shared = link_groups(expected, shared, ['Fam1'] * 4, 0.25)
    assert groups == [0, 0, 0, 1]
    assert max_shared == [pytest.approx(0.6), None]
    groups, max_shared = link_groups(expected, shared, ['Fam1'] * 4, 0.7)
    assert groups == [0, 1, 2, 3] and max_shared == [None] * 4


def test_link_groups_links_only_rows_of_one_family():
    expected = np.array([10.0, 10.0, 10.0, 10.0])
    shared = np.full((4, 4), 5.0)   # f = 0.5 between every pair
    np.fill_diagonal(shared, 0.0)
    groups, max_shared = link_groups(expected, shared, ['Fam1', 'Fam2', 'Fam1', 'Fam2'])
    assert groups == [0, 1, 0, 1] and max_shared == [pytest.approx(0.5)] * 2
    # A row with no family links to nothing, not even to another row with none
    groups, max_shared = link_groups(expected, shared, ['Fam1', None, 'Fam1', None])
    assert groups == [0, 1, 0, 2] and max_shared == [pytest.approx(0.5), None, None]


def test_tax_families_are_none_where_the_taxonomy_has_no_family():
    taxonomy = pd.DataFrame({'species': ['A a', 'B b', 'C c', 'D d', 'D e'], 'family': ['Fam1 ', '', ' ', None, 'Fam9']},
                            index=pd.Index(['1', '2', '3', '4', '1'], name='tax_id'))
    assert tax_families([1, 2, 3, 4, 5], taxonomy) == {1: ('Fam1 ',), 2: None, 3: None, 4: None, 5: None}
    assert tax_families([1], taxonomy.drop(columns='family')) == {1: None}


def test_tax_families_tell_apart_one_family_name_in_two_orders():
    taxonomy = pd.DataFrame({'family': ['Family_XI', 'Family_XI', 'Family_XI'], 'order': ['Ord1', 'Ord2', 'Ord1']},
                            index=pd.Index(['1', '2', '3'], name='tax_id'))
    families = tax_families([1, 2, 3], taxonomy)
    assert families == {1: ('Family_XI', 'Ord1'), 2: ('Family_XI', 'Ord2'), 3: ('Family_XI', 'Ord1')}
    expected, shared = np.full(3, 10.0), np.full((3, 3), 5.0)
    np.fill_diagonal(shared, 0.0)
    assert link_groups(expected, shared, [families[t] for t in (1, 2, 3)])[0] == [0, 1, 0]


def test_group_shared_reads_count_only_splits_with_other_groups():
    groups = ambiguity_groups(IDX, PHI, [11, 12, 13], ONE_FAMILY, min_shared=0.3)
    # f_1(0) = 0.5 / 1.4 links columns 0 and 1; f_2(1) = 0.1 / 0.105 links column 2 too
    assert [g['tax_ids'] for g in groups] == [[11, 12, 13]]
    assert groups[0]['shared_reads'] == 0.0
    groups = ambiguity_groups(IDX, PHI, [11, 12, 13], ONE_FAMILY, min_shared=0.99)
    assert [g['tax_ids'] for g in groups] == [[11], [12], [13]]
    # Read 3 replaced by a read of column 2's own, and read 2 split 0.98 / 0.02: f_1(0) = 1
    # links columns 0 and 1; f_0(2) = 0.98 / 2.48 and f_2(0) = 0.02 / 1.02 stay below 0.45
    idx = np.array([[0, 1], [0, -1], [0, 2], [2, -1]])
    phi = np.array([[0.5, 0.5], [1.0, 0.0], [0.98, 0.02], [1.0, 0.0]])
    groups = ambiguity_groups(idx, phi, [11, 12, 13], ONE_FAMILY, min_shared=0.45)
    assert [g['tax_ids'] for g in groups] == [[11, 12], [13]]
    # Read 2 is split between the groups; read 0, split within the first, is its own
    assert groups[0]['shared_reads'] == pytest.approx(0.98)
    assert groups[1]['shared_reads'] == pytest.approx(0.02)
    assert groups[0]['max_shared'] == pytest.approx(1.0) and groups[1]['max_shared'] is None


def test_tax_ids_of_one_table_row_are_one_row_before_linking():
    # Two strains reported in one row (F46) split every read: the row is a group of one
    idx, phi = np.array([[0, 1]] * 4), np.array([[0.5, 0.5]] * 4)
    families = {1: 'Fam1', 2: 'Fam1'}
    groups = ambiguity_groups(idx, phi, [1, 2], families, row_keys={1: ('Alpha one',), 2: ('Alpha one',)})
    assert groups == [{'tax_ids': [1, 2], 'shared_reads': 0.0, 'max_shared': None}]
    groups = ambiguity_groups(idx, phi, [1, 2], families)
    assert groups == [{'tax_ids': [1, 2], 'shared_reads': 0.0, 'max_shared': pytest.approx(1.0)}]


# ─── estimate_abundance.py --groups and write_abundance_table.py --groups-output ───

TAXONOMY = {1: ('Alpha one', 'Alpha', 'Fam1'), 2: ('Alpha two', 'Alpha', 'Fam1'), 3: ('Beta three', 'Beta', 'Fam1'),
            4: ('Gamma four', 'Gamma', 'Fam2'), 5: ('Alpha five', 'Alpha', 'Fam1')}


def write_inputs(tmp_path):
    (tmp_path / 'taxonomy.tsv').write_text('tax_id\tspecies\tgenus\tfamily\n' + ''.join(
        f'{t}\t{s}\t{g}\t{f}\n' for t, (s, g, f) in TAXONOMY.items()))
    # Species 2 shares an identical 16S with 1, and 5 (never reported) with 1
    (tmp_path / 'identical_16s.tsv').write_text(
        'tax_id\tspecies\tgenus\tgroup\tgroup_size\tidentical_16s_species\n'
        '5\tAlpha five\tAlpha\tG1\t3\tAlpha one\n1\tAlpha one\tAlpha\tG1\t3\tAlpha five; Alpha two\n'
        '2\tAlpha two\tAlpha\tG1\t3\tAlpha one\n')
    # Species 1 and 2 tie on 300 reads and species 2 has 20 of its own: 0.41 of its 34
    # expected reads (and of species 1's 690) are on the tied reads; 3 and 4 are on their
    # own, but for 5 reads 1 / 3
    per_read = {f'a{i}': ([1], [-1.0]) for i in range(400)}
    per_read.update({f't{i}': ([1, 2], [-1.0, -1.0]) for i in range(300)})
    per_read.update({f'b{i}': ([2], [-1.0]) for i in range(20)})
    per_read.update({f'c{i}': ([3], [-1.0]) for i in range(300)})
    per_read.update({f'x{i}': ([1, 3], [-1.0, -1.5]) for i in range(5)})
    per_read.update({f'd{i}': ([4], [-1.0]) for i in range(175)})
    (tmp_path / 'counts.json').write_text(json.dumps({'assigned_count': 1200, 'unassigned_count': 7}))
    return write_npz(tmp_path / 's.npz', per_read)


def write_groups(tmp_path, *options):
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'abundance.json',
         '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv', '--output-prefix', 'S1_rel-abundance',
         '--groups-output', 'S1_groups', *options],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    read = dict(sep='\t', keep_default_na=False)
    return pd.read_csv(tmp_path / 'S1_rel-abundance.tsv', **read), pd.read_csv(tmp_path / 'S1_groups.tsv', **read)


@pytest.mark.parametrize('estimator', ['vb', 'em'])
def test_groups_leave_the_abundances_unchanged_and_partition_the_rows(tmp_path, estimator):
    npz = write_inputs(tmp_path)
    without = estimate(tmp_path, npz, '--estimator', estimator, '--taxonomy', 'taxonomy.tsv')
    data = estimate(tmp_path, npz, '--estimator', estimator, '--taxonomy', 'taxonomy.tsv', '--groups')
    assert {k: v for k, v in data.items() if k != 'groups'} == without
    assert sorted(t for g in data['groups'] for t in g['tax_ids']) == sorted(int(t) for t in data['abundance'])
    assert [sorted(g['tax_ids']) for g in data['groups'] if len(g['tax_ids']) > 1] == [[1, 2]]


@pytest.mark.parametrize('estimator', ['vb', 'em'])
def test_groups_table_sums_the_species_rows(tmp_path, estimator):
    npz = write_inputs(tmp_path)
    estimate(tmp_path, npz, '--estimator', estimator, '--taxonomy', 'taxonomy.tsv', '--groups')
    species, groups = write_groups(tmp_path, '--counts', '--identical-16s', 'identical_16s.tsv')
    assert list(groups.columns) == ['group', 'members', 'n_members', 'abundance', 'lowest_common_rank',
                                    'lowest_common_taxon', 'estimated_counts', 'identical_16s_relatives',
                                    'max_shared', 'shared_read_fraction']
    assert groups['group'].tolist() == ['G1', 'G2', 'G3', 'Below detection', 'Unassigned']
    pair = groups.iloc[0]
    assert pair['members'] == 'Alpha one; Alpha two' and pair['n_members'] == '2'
    by_name = species.set_index('species')
    assert pair['abundance'] == pytest.approx(by_name.loc[['Alpha one', 'Alpha two'], 'abundance'].sum(), abs=0)
    assert pair['estimated_counts'] == by_name.loc[['Alpha one', 'Alpha two'], 'estimated_counts'].sum()
    assert (pair['lowest_common_rank'], pair['lowest_common_taxon']) == ('genus', 'Alpha')
    # The reported partner is a member, not a relative
    assert pair['identical_16s_relatives'] == 'Alpha five'
    assert float(pair['max_shared']) > 0.25
    # The tied reads are the group's own now; only the reads split with Beta three are shared
    assert float(pair['shared_read_fraction']) < 0.02
    assert float(pair['shared_read_fraction']) < float(by_name.loc['Alpha two', 'shared_read_fraction'])
    singles = groups.iloc[1:3].set_index('members')
    assert list(singles.index) == ['Beta three', 'Gamma four']
    assert singles['lowest_common_rank'].tolist() == ['species', 'species']
    assert singles['max_shared'].tolist() == ['', '']
    assert singles.loc['Beta three', 'shared_read_fraction'] == by_name.loc['Beta three', 'shared_read_fraction']
    # Like the species table, the group table adds up to 1 and to the filtered reads
    assert groups['abundance'].sum() == pytest.approx(1.0)
    assert groups['estimated_counts'].sum() == species['estimated_counts'].sum() == 1207
    assert groups.iloc[3:][['members', 'abundance', 'estimated_counts']].values.tolist() == \
        species.iloc[-2:][['species', 'abundance', 'estimated_counts']].values.tolist()


def test_groups_table_intervals_come_from_the_summed_alpha(tmp_path):
    npz = write_inputs(tmp_path)
    estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups')
    species, groups = write_groups(tmp_path, '--credible-intervals')
    by_name = species.set_index('species')
    singles = groups.set_index('members').loc[['Beta three', 'Gamma four']]
    for col in ('ci_lower', 'ci_upper'):
        assert singles[col].tolist() == by_name.loc[['Beta three', 'Gamma four'], col].tolist()
    pair = groups.iloc[0]
    assert float(pair['ci_lower']) < pair['abundance'] < float(pair['ci_upper'])
    # Without --identical-16s there is no relatives column; without --counts no counts
    assert 'identical_16s_relatives' not in groups.columns and 'estimated_counts' not in groups.columns


def test_groups_with_binomial_names_group_table_rows(tmp_path):
    # 'Alpha two str. B' is cut to 'Alpha two'; the strains of one row never form a pair
    npz = write_inputs(tmp_path)
    taxonomy = (tmp_path / 'taxonomy.tsv').read_text().replace('Gamma four', 'Alpha one str. B') \
        .replace('Gamma\tFam2', 'Alpha\tFam1')
    (tmp_path / 'taxonomy.tsv').write_text(taxonomy)
    data = estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups', '--binomial-names')
    _, groups = write_groups(tmp_path, '--binomial-names')
    assert sorted(sorted(g['tax_ids']) for g in data['groups']) == [[1, 2, 4], [3]]
    assert groups['members'].tolist()[:2] == ['Alpha one; Alpha two', 'Beta three']


@pytest.mark.parametrize('family', ['Fam3', ''])
def test_groups_link_only_species_of_one_family(tmp_path, family):
    # Alpha two in another family, or with a blank one, no longer joins Alpha one
    npz = write_inputs(tmp_path)
    taxonomy = (tmp_path / 'taxonomy.tsv').read_text().replace('Alpha two\tAlpha\tFam1', f'Alpha two\tAlpha\t{family}')
    (tmp_path / 'taxonomy.tsv').write_text(taxonomy)
    data = estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups')
    assert all(len(g['tax_ids']) == 1 for g in data['groups'])


def test_groups_without_a_family_column_link_nothing_and_say_so(tmp_path):
    npz = write_inputs(tmp_path)
    (tmp_path / 'taxonomy.tsv').write_text('tax_id\tspecies\tgenus\n' + ''.join(
        f'{t}\t{s}\t{g}\n' for t, (s, g, _) in TAXONOMY.items()))
    result = subprocess.run([sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', npz,
                             '--output', 'abundance.json', '--taxonomy', 'taxonomy.tsv', '--groups'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'no family column' in result.stderr
    groups = json.loads((tmp_path / 'abundance.json').read_text())['groups']
    assert all(len(g['tax_ids']) == 1 for g in groups)


def test_groups_lowest_common_rank_skips_ranks_the_taxonomy_lacks(tmp_path):
    # No genus column: the table fills the genus with a placeholder, which is not a common taxon
    npz = write_inputs(tmp_path)
    (tmp_path / 'taxonomy.tsv').write_text('tax_id\tspecies\tfamily\n' + ''.join(
        f'{t}\t{s}\t{f}\n' for t, (s, _, f) in TAXONOMY.items()))
    estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups')
    _, groups = write_groups(tmp_path)
    pair = groups.iloc[0]
    assert pair['members'] == 'Alpha one; Alpha two'
    assert (pair['lowest_common_rank'], pair['lowest_common_taxon']) == ('family', 'Fam1')


def test_groups_with_a_species_missing_from_the_taxonomy(tmp_path):
    # Alpha two is an 'Unknown taxid' row with no family: it links to nothing
    npz = write_inputs(tmp_path)
    taxonomy = ''.join(line for line in (tmp_path / 'taxonomy.tsv').read_text().splitlines(keepends=True)
                       if not line.startswith('2\t'))
    (tmp_path / 'taxonomy.tsv').write_text(taxonomy)
    data = estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups')
    assert sorted(t for g in data['groups'] for t in g['tax_ids']) == sorted(int(t) for t in data['abundance'])
    assert all(len(g['tax_ids']) == 1 for g in data['groups'])


def test_groups_with_no_scored_read(tmp_path):
    write_inputs(tmp_path)
    data = estimate(tmp_path, write_npz(tmp_path / 'empty.npz', {}), '--taxonomy', 'taxonomy.tsv', '--groups')
    assert data['groups'] == []
    _, groups = write_groups(tmp_path, '--counts', '--identical-16s', '')
    assert groups['group'].tolist() == ['Below detection', 'Unassigned']
    assert groups['estimated_counts'].tolist() == [0, 7]


def test_groups_options_are_checked(tmp_path):
    npz = write_inputs(tmp_path)
    argv = [sys.executable, os.path.join(BIN, 'estimate_abundance.py'), '--scores', npz, '--output', 'a.json']
    result = subprocess.run(argv + ['--groups'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 2 and '--groups needs --taxonomy' in result.stderr
    result = subprocess.run(argv + ['--groups', '--taxonomy', 'taxonomy.tsv', '--group-min-shared', '0'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 2 and 'above 0 and at most 1' in result.stderr
    # A JSON written without --groups has no groups to write
    estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'abundance.json',
         '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv', '--output-prefix', 'out',
         '--groups-output', 'out_groups'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 2 and 'needs abundances from estimate_abundance.py --groups' in result.stderr


def test_group_shared_reads_equal_the_read_summary_with_groups_as_rows():
    # The sparse computation gives what posterior_read_summary gives with the groups as table rows
    from variational_inference_bayesian import posterior_read_summary
    rng = np.random.default_rng(5)
    idx = rng.integers(-1, 6, size=(300, 4))
    phi = rng.dirichlet(np.ones(4), size=300)
    ids = [11, 12, 13, 14, 15, 16]
    for min_shared in (0.05, 0.25, 0.6):
        groups = ambiguity_groups(idx, phi, ids, dict.fromkeys(ids, 'Fam1'), min_shared=min_shared)
        group_of = {sid: g for g, entry in enumerate(groups) for sid in entry['tax_ids']}
        summary = posterior_read_summary(None, idx, phi, ids, groups=[group_of[s] for s in ids])
        for g, entry in enumerate(groups):
            assert entry['shared_reads'] == pytest.approx(sum(summary[s]['shared_reads'] for s in entry['tax_ids']),
                                                          rel=1e-12, abs=1e-12)


def test_groups_of_one_row_of_several_tax_ids_match_the_species_row(tmp_path):
    # With --binomial-names three strains are one row; its group of one has the row's exact values
    npz = write_inputs(tmp_path)
    taxonomy = (tmp_path / 'taxonomy.tsv').read_text().replace('Gamma four\tGamma\tFam2', 'Beta three str. X\tBeta\tFam1')
    (tmp_path / 'taxonomy.tsv').write_text(taxonomy)
    estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups', '--binomial-names')
    species, groups = write_groups(tmp_path, '--binomial-names', '--credible-intervals')
    row = species.set_index('species').loc['Beta three']
    group = groups.set_index('members').loc['Beta three']
    for col in ('abundance', 'shared_read_fraction', 'ci_lower', 'ci_upper'):
        assert group[col] == row[col]


def test_a_tax_id_missing_from_the_read_summary_counts_zero(tmp_path):
    # As in the species table, which writes 0 for it
    npz = write_inputs(tmp_path)
    data = estimate(tmp_path, npz, '--taxonomy', 'taxonomy.tsv', '--groups')
    del data['read_summary']['4']
    (tmp_path / 'abundance.json').write_text(json.dumps(data))
    species, groups = write_groups(tmp_path)
    assert groups.set_index('members').loc['Gamma four', 'shared_read_fraction'] == \
        species.set_index('species').loc['Gamma four', 'shared_read_fraction']
