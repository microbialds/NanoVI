"""16S copy numbers and the cell-level abundance table (F58, D40)."""
import gzip
import json
import os
import subprocess
import sys

import pandas as pd
import pytest

from .test_scripts import BIN

sys.path.insert(0, BIN)
import build_gtdb_db  # noqa: E402
import copy_numbers  # noqa: E402
from copy_numbers import fill_copy_numbers, gtdb_metadata_values, read_copy_numbers  # noqa: E402

RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


def tax_row(tax_id, species, genus, family, order='O', cls='C', phylum='P'):
    return {'tax_id': str(tax_id), 'species': species, 'genus': genus, 'family': family, 'order': order,
            'class': cls, 'phylum': phylum, 'superkingdom': 'Bacteria'}


TAXONOMY = [
    tax_row(1, 'Alpha one', 'Alpha', 'F1'),
    tax_row(2, 'Alpha two', 'Alpha', 'F1'),
    tax_row(3, 'Alpha three', 'Alpha', 'F1'),      # no value: genus median of 4 and 5
    tax_row(4, 'Beta one', 'Beta', 'F1'),          # no value, genus has none: family median of 4, 5
    tax_row(5, 'Gamma one', 'Gamma', 'F2'),        # no value anywhere in F2
    tax_row(6, 'Alpha other', 'Alpha', 'F9'),      # an Alpha of another family is another genus
    tax_row(7, 'Unplaced one', '', ''),            # blank genus and family: none
    tax_row(8, 'Delta one', '', 'F1'),             # blank genus: family median
]


def by_tax(rows):
    return {r['tax_id']: (r['copy_number'], r['level']) for r in rows}


def test_fill_uses_the_genus_then_the_family_median():
    rows = fill_copy_numbers(TAXONOMY, {'Alpha one': (4.0, 3), 'Alpha two': (5.0, 1)}, 'test')
    assert by_tax(rows) == {
        '1': (4.0, 'species'), '2': (5.0, 'species'), '3': (4.5, 'genus'), '4': (4.5, 'family'),
        '5': (None, 'none'), '6': (None, 'none'), '7': (None, 'none'), '8': (4.5, 'family'),
    }
    assert [r['n_genomes'] for r in rows[:3]] == [3, 1, '']
    assert {r['source'] for r in rows} == {'test'}


def test_fill_takes_each_species_value_once_and_never_a_filled_one():
    taxonomy = [tax_row(1, 'A a', 'A', 'F'), tax_row(2, 'A b', 'A', 'F'), tax_row(3, 'A c', 'A', 'F'),
                tax_row(4, 'B a', 'B', 'F'), tax_row(5, 'B b', 'B', 'F')]
    rows = fill_copy_numbers(taxonomy, {'A a': (2.0, 1), 'A b': (2.0, 9), 'B a': (10.0, 1)}, 'test')
    # B b: genus B's one value; A c: genus A's; the family median is not needed
    assert by_tax(rows)['5'] == (10.0, 'genus') and by_tax(rows)['3'] == (2.0, 'genus')


def test_fill_by_tax_id_for_rrndb():
    rows = fill_copy_numbers(TAXONOMY[:3], {'1': (7.0, '12')}, 'rrndb', by='tax_id')
    assert by_tax(rows) == {'1': (7.0, 'species'), '2': (7.0, 'genus'), '3': (7.0, 'genus')}


def write_metadata(path, rows):
    header = 'accession\tgtdb_taxonomy\tncbi_assembly_level\tother\tssu_count\n'
    text = header + ''.join('\t'.join(r) + '\n' for r in rows)
    with gzip.open(path, 'wt') as fh:
        fh.write(text)


LIN = 'd__Bacteria;p__P;c__C;o__O;f__F;g__{g};s__{s}'
METADATA = [
    ('RS_GCF_1', LIN.format(g='Alpha_A', s='Alpha_A one'), 'Complete Genome', 'x', '6'),
    ('RS_GCF_2', LIN.format(g='Alpha_A', s='Alpha_A one'), 'Complete Genome', 'x', '8'),
    ('RS_GCF_3', LIN.format(g='Alpha_A', s='Alpha_A one'), 'Contig', 'x', '1'),       # a draft: not used
    ('RS_GCF_4', LIN.format(g='Alpha', s='Alpha one'), 'Complete Genome', 'x', '3'),
    ('RS_GCF_5', LIN.format(g='Alpha', s='Alpha one'), 'Complete Genome', 'x', '0'),  # no 16S found: skipped
    ('RS_GCF_6', LIN.format(g='Beta', s='Beta one'), 'Chromosome', 'x', '9'),
    ('RS_GCF_7', LIN.format(g='Beta', s='Beta one'), 'Complete Genome', 'x', 'none'),
]


def test_gtdb_metadata_uses_complete_genomes_only(tmp_path):
    path = tmp_path / 'bac120_metadata_r0.tsv.gz'
    write_metadata(path, METADATA)
    kept = gtdb_metadata_values([path], lambda lin: build_gtdb_db.parse_gtdb_lineage(lin)['species'])
    assert kept == {'Alpha_A one': (7.0, 2), 'Alpha one': (3.0, 1)}
    merged = gtdb_metadata_values([path], lambda lin: build_gtdb_db.parse_gtdb_lineage(lin, True)['species'])
    assert merged == {'Alpha one': (6.0, 3)}


def test_gtdb_metadata_must_have_the_columns(tmp_path):
    path = tmp_path / 'm.tsv'
    path.write_text('accession\tssu_count\nRS_1\t3\n')
    with pytest.raises(ValueError, match='not a GTDB metadata file'):
        gtdb_metadata_values([path], lambda lin: lin)


@pytest.mark.parametrize('text,message', [
    ('species\tcopy_number\nA a\t3\nA a\t4\n', 'more than once'),
    ('species\tcopy_number\nA a\t0\n', 'must be a positive finite number'),
    ('species\tcopy_number\nA a\tinf\n', 'must be a positive finite number'),
    ('species\tcopy_number\nA a\tmany\n', 'not a number'),
    ('name\tcopies\nA a\t3\n', "needs 'species' and 'copy_number'"),
])
def test_table_values_are_checked(tmp_path, text, message):
    path = tmp_path / 't.tsv'
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        copy_numbers.table_values(path)


def test_table_values_skip_blank_values_and_keep_extra_columns(tmp_path):
    path = tmp_path / 't.tsv'
    path.write_text('species\tcopy_number\tnote\nA a\t3\tfrom X\nA b\t\tunknown\n')
    assert copy_numbers.table_values(path) == {'A a': (3.0, '')}


def write_taxonomy(path, rows):
    pd.DataFrame(rows)[['tax_id'] + RANKS].to_csv(path, sep='\t', index=False)


def test_build_copy_numbers_writes_and_removes_the_file(tmp_path):
    write_taxonomy(tmp_path / 'taxonomy.tsv', TAXONOMY[:3])
    table = tmp_path / 'source.tsv'
    table.write_text('species\tcopy_number\nAlpha one\t4\nAlpha two\t7\n')
    copy_numbers.build_copy_numbers(str(tmp_path), tmp_path / 'taxonomy.tsv', table=str(table))
    out = pd.read_csv(tmp_path / 'copy_numbers.tsv', sep='\t', dtype=str, keep_default_na=False)
    assert list(out.columns) == copy_numbers.COPY_NUMBER_COLUMNS
    assert out.values.tolist() == [['1', 'Alpha one', '4', 'species', '', 'table:source.tsv'],
                                   ['2', 'Alpha two', '7', 'species', '', 'table:source.tsv'],
                                   ['3', 'Alpha three', '5.5', 'genus', '', 'table:source.tsv']]
    assert read_copy_numbers(tmp_path / 'copy_numbers.tsv')['3'] == ('Alpha three', 5.5, 'genus')
    # A build without a source removes the earlier build's file
    copy_numbers.build_copy_numbers(str(tmp_path), tmp_path / 'taxonomy.tsv')
    assert not (tmp_path / 'copy_numbers.tsv').exists()
    with pytest.raises(ValueError, match='one copy-number source'):
        copy_numbers.build_copy_numbers(str(tmp_path), tmp_path / 'taxonomy.tsv', table=str(table),
                                        rrndb_pantaxa=str(table))


def ssu_fasta(path):
    seq_a, seq_b = 'ACGT' * 300, 'TTGCA' * 240
    records = [(f'>RS_GCF_1~c1 {LIN.format(g="Alpha_A", s="Alpha_A one")} [ssu_len=1200]', seq_a),
               (f'>RS_GCF_4~c4 {LIN.format(g="Alpha", s="Alpha one")} [ssu_len=1200]', seq_b),
               (f'>RS_GCF_6~c6 {LIN.format(g="Beta", s="Beta one")} [ssu_len=1200]', seq_b[::-1])]
    path.write_text(''.join(f'{h}\n{s}\n' for h, s in records))


def test_build_gtdb_db_writes_copy_numbers_from_the_metadata(tmp_path):
    ssu, meta, out = tmp_path / 'ssu.fna', tmp_path / 'bac120_metadata_r0.tsv.gz', tmp_path / 'db'
    ssu_fasta(ssu)
    write_metadata(meta, METADATA)
    result = subprocess.run([sys.executable, os.path.join(BIN, 'build_gtdb_db.py'), '--ssu', str(ssu),
                             '--output-dir', str(out), '--gtdb-metadata', str(meta)],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    cn = pd.read_csv(out / 'copy_numbers.tsv', sep='\t', dtype=str, keep_default_na=False)
    assert cn[['species', 'copy_number', 'level', 'n_genomes']].values.tolist() == [
        ['Alpha one', '3', 'species', '1'], ['Alpha_A one', '7', 'species', '2'],
        # Beta one has only a Chromosome assembly and a blank count: the family median of 3 and 7
        ['Beta one', '5', 'family', '']]
    assert set(cn['source']) == {'gtdb_complete_genomes:bac120_metadata_r0.tsv.gz'}
    # Rebuilt without the metadata: no copy numbers
    build_gtdb_db.build_database(str(ssu), 'db_gtdb', str(out))
    assert not (out / 'copy_numbers.tsv').exists()


def test_build_gtdb_db_takes_one_source(tmp_path):
    result = subprocess.run([sys.executable, os.path.join(BIN, 'build_gtdb_db.py'), '--ssu', 'x.fna',
                             '--gtdb-metadata', 'm.tsv', '--copy-numbers', 't.tsv'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0 and 'not allowed with' in result.stderr


def build_database_inputs(tmp_path):
    seq_a, seq_b = 'ACGT' * 300, 'TTGCA' * 240
    (tmp_path / 'seqs.fa').write_text(f'>a1\n{seq_a}\n>b1\n{seq_b}\n')
    (tmp_path / 'seq2tax.tsv').write_text('a1\t287\nb1\t1001\n')


def test_build_database_takes_a_copy_number_table(tmp_path):
    build_database_inputs(tmp_path)
    (tmp_path / 'seq2tax.tsv').write_text('a1\t10\nb1\t20\n')
    (tmp_path / 'lineages.tsv').write_text(
        'tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n'
        '10\tAlpha one\tAlpha\tF1\tO\tC\tP\tBacteria\n20\tAlpha two\tAlpha\tF1\tO\tC\tP\tBacteria\n')
    (tmp_path / 'cn.tsv').write_text('species\tcopy_number\nAlpha one\t6\n')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'main.py'), 'build-database', 'mydb', '--sequences', 'seqs.fa',
         '--seq2tax', 'seq2tax.tsv', '--taxonomy-list', 'lineages.tsv', '--copy-numbers', 'cn.tsv'],
        cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert read_copy_numbers(tmp_path / 'mydb' / 'copy_numbers.tsv') == {
        '10': ('Alpha one', 6.0, 'species'), '20': ('Alpha two', 6.0, 'genus')}


def test_build_database_looks_rrndb_up_by_ncbi_taxid(tmp_path, mock_nodes_dict, mock_names_dict):
    build_database_inputs(tmp_path)
    dump = tmp_path / 'taxdump'
    dump.mkdir()
    nodes = dict(mock_nodes_dict, **{'1001': ('286', 'species')})
    names = dict(mock_names_dict, **{'1001': 'Pseudomonas putida'})
    (dump / 'nodes.dmp').write_text(''.join(f'{t}\t|\t{p}\t|\t{r}\t|\n' for t, (p, r) in nodes.items()))
    (dump / 'names.dmp').write_text(''.join(f'{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n' for t, n in names.items()))
    (tmp_path / 'pantaxa.tsv').write_text(
        'taxid\trank\tname\tchildcount\tmin\tmax\tmode\tmedian\tmean\tstddev\tsum16slist\n'
        '287\tspecies\tPseudomonas aeruginosa\t3\t4\t4\t4\t4\t4.00\t0.00\t[4, 4, 4]\n'
        '286\tgenus\tPseudomonas\t9\t1\t7\t4\t5\t4.5\t1\t[]\n')
    run = lambda *extra: subprocess.run(
        [sys.executable, os.path.join(BIN, 'main.py'), 'build-database', 'mydb', '--sequences', 'seqs.fa',
         '--seq2tax', 'seq2tax.tsv', *extra, '--rrndb-pantaxa', 'pantaxa.tsv'],
        cwd=tmp_path, capture_output=True, text=True)
    result = run('--ncbi-taxonomy', str(dump))
    assert result.returncode == 0, result.stderr
    cn = pd.read_csv(tmp_path / 'mydb' / 'copy_numbers.tsv', sep='\t', dtype=str, keep_default_na=False)
    # P. putida has no species row: the genus median of the database's species values (not rrnDB's genus row)
    assert cn[['tax_id', 'copy_number', 'level', 'n_genomes']].values.tolist() == [
        ['287', '4', 'species', '3'], ['1001', '4', 'genus', '']]
    (tmp_path / 'lineages.tsv').write_text('tax_id\tspecies\n287\tPseudomonas aeruginosa\n')
    result = run('--taxonomy-list', 'lineages.tsv')
    assert result.returncode != 0 and '--rrndb-pantaxa needs --ncbi-taxonomy' in result.stderr


# ─── write_abundance_table.py --copy-numbers (the cell-level table) ──────────

def write_inputs(tmp_path, freq, taxonomy, cn_rows, groups=None, below=5, unassigned=7):
    data = {'abundance': freq, 'below_detection_count': below,
            'read_summary': {t: {'expected_reads': 100 * v, 'shared_reads': 0.0} for t, v in freq.items()}}
    if groups is not None:
        data['groups'] = groups
    (tmp_path / 'ab.json').write_text(json.dumps(data))
    (tmp_path / 'counts.json').write_text(json.dumps({'assigned_count': 100, 'unassigned_count': unassigned}))
    write_taxonomy(tmp_path / 'taxonomy.tsv', taxonomy)
    copy_numbers.write_copy_numbers(cn_rows, tmp_path / 'copy_numbers.tsv')


def run_table(tmp_path, *extra, cell=True):
    args = [sys.executable, os.path.join(BIN, 'write_abundance_table.py'), '--abundance', 'ab.json',
            '--read-counts', 'counts.json', '--taxonomy', 'taxonomy.tsv', '--output-prefix', 's_rel', '--counts',
            *extra]
    if cell:
        args += ['--copy-numbers', 'copy_numbers.tsv', '--cell-output', 's_cell']
    return subprocess.run(args, cwd=tmp_path, capture_output=True, text=True)


def cn_row(tax_id, species, cn, level='species'):
    return {'tax_id': str(tax_id), 'species': species, 'copy_number': cn, 'level': level, 'n_genomes': '',
            'source': 'test'}


def read_cell(path):
    lines = path.read_text().splitlines()
    return lines[:2], pd.read_csv(path, sep='\t', comment='#')


def test_cell_table_divides_by_copy_number_and_leaves_out_what_has_none(tmp_path):
    taxonomy = [tax_row(1, 'Alpha one', 'Alpha', 'F1'), tax_row(2, 'Beta one', 'Beta', 'F1'),
                tax_row(3, 'Gamma one', 'Gamma', 'F2')]
    write_inputs(tmp_path, {'1': 0.5, '2': 0.3, '3': 0.2}, taxonomy,
                 [cn_row(1, 'Alpha one', 10.0), cn_row(2, 'Beta one', 2.5, 'genus'), cn_row(3, 'Gamma one', None, 'none')])
    result = run_table(tmp_path)
    assert result.returncode == 0, result.stderr
    comments, cell = read_cell(tmp_path / 's_cell.tsv')
    assert comments[0].startswith('# Cell-level abundance')
    assert comments[1] == ('# Not included: Below detection (5 reads) and Unassigned (7 reads), '
                           'and 1 species without a copy number (read-level abundance 0.2).')
    assert list(cell.columns) == copy_numbers.CELL_COLUMNS + RANKS[1:]
    assert list(cell['species']) == ['Alpha one', 'Beta one']
    # 0.5 / 10 = 0.05 and 0.3 / 2.5 = 0.12, renormalised
    assert cell['abundance'].tolist() == pytest.approx([0.05 / 0.17, 0.12 / 0.17])
    assert cell['read_abundance'].tolist() == [0.5, 0.3]
    assert cell['copy_number'].tolist() == [10.0, 2.5]
    assert cell['copy_number_level'].tolist() == ['species', 'genus']


def test_species_table_does_not_change(tmp_path):
    taxonomy = [tax_row(1, 'Alpha one', 'Alpha', 'F1'), tax_row(2, 'Beta one', 'Beta', 'F1')]
    write_inputs(tmp_path, {'1': 0.6, '2': 0.4}, taxonomy, [cn_row(1, 'Alpha one', 3.0), cn_row(2, 'Beta one', 5.0)])
    assert run_table(tmp_path, cell=False).returncode == 0
    before = (tmp_path / 's_rel.tsv').read_bytes()
    assert run_table(tmp_path).returncode == 0
    assert (tmp_path / 's_rel.tsv').read_bytes() == before


def test_a_merged_row_takes_the_cell_weighted_mean(tmp_path):
    # --binomial-names merges the two strains into one row; the third tax ID has no copy number
    taxonomy = [tax_row(1, 'Alpha one A', 'Alpha', 'F1'), tax_row(2, 'Alpha one B', 'Alpha', 'F1'),
                tax_row(3, 'Alpha one C', 'Alpha', 'F1'), tax_row(4, 'Beta one', 'Beta', 'F1')]
    write_inputs(tmp_path, {'1': 0.4, '2': 0.2, '3': 0.1, '4': 0.3}, taxonomy,
                 [cn_row(1, 'Alpha one A', 8.0), cn_row(2, 'Alpha one B', 2.0, 'genus'),
                  cn_row(3, 'Alpha one C', None, 'none'), cn_row(4, 'Beta one', 3.0)])
    result = run_table(tmp_path, '--binomial-names')
    assert result.returncode == 0, result.stderr
    comments, cell = read_cell(tmp_path / 's_cell.tsv')
    alpha_cells = 0.4 / 8 + 0.2 / 2
    assert cell['species'].tolist() == ['Alpha one', 'Beta one']
    assert cell['copy_number'].tolist() == pytest.approx([0.6 / alpha_cells, 3.0])
    assert cell['copy_number_level'].tolist() == ['genus', 'species']
    assert cell['abundance'].tolist() == pytest.approx([alpha_cells / (alpha_cells + 0.1), 0.1 / (alpha_cells + 0.1)])
    assert cell['read_abundance'].tolist() == pytest.approx([0.7, 0.3])
    assert comments[1].endswith('and 0 species without a copy number (read-level abundance 0); tax IDs without one '
                                'in rows listed below: read-level abundance 0.1, left out of their row\'s cell share '
                                'and copy number.')


def test_groups_get_the_members_cells_and_copy_number(tmp_path):
    taxonomy = [tax_row(1, 'Alpha one', 'Alpha', 'F1'), tax_row(2, 'Alpha two', 'Alpha', 'F1'),
                tax_row(3, 'Beta one', 'Beta', 'F2')]
    groups = [{'tax_ids': ['1', '2'], 'shared_reads': 10.0, 'max_shared': 0.5},
              {'tax_ids': ['3'], 'shared_reads': 0.0, 'max_shared': 0.0}]
    write_inputs(tmp_path, {'1': 0.5, '2': 0.25, '3': 0.25}, taxonomy,
                 [cn_row(1, 'Alpha one', 5.0), cn_row(2, 'Alpha two', 1.0), cn_row(3, 'Beta one', 5.0, 'family')],
                 groups=groups)
    result = run_table(tmp_path, '--groups-output', 's_groups')
    assert result.returncode == 0, result.stderr
    g = pd.read_csv(tmp_path / 's_groups.tsv', sep='\t')
    _, cell = read_cell(tmp_path / 's_cell.tsv')
    assert g.columns[-3:].tolist() == ['cell_abundance', 'copy_number', 'copy_number_level']
    assert g['cell_abundance'][:2].tolist() == pytest.approx(
        [cell['abundance'][:2].sum(), cell['abundance'][2]])
    # (0.5 + 0.25) / (0.5 / 5 + 0.25 / 1): the cell-weighted mean of 5 and 1
    assert g['copy_number'][:2].tolist() == pytest.approx([0.75 / 0.35, 5.0])
    assert g['copy_number_level'][:2].tolist() == ['species', 'family']
    assert g['cell_abundance'][2:].isna().all()
    # Without copy numbers the group table has no cell columns
    assert run_table(tmp_path, '--groups-output', 's_groups', cell=False).returncode == 0
    assert 'cell_abundance' not in pd.read_csv(tmp_path / 's_groups.tsv', sep='\t').columns


def test_copy_numbers_of_another_database_are_refused(tmp_path):
    taxonomy = [tax_row(1, 'Alpha one', 'Alpha', 'F1')]
    write_inputs(tmp_path, {'1': 1.0}, taxonomy, [cn_row(1, 'Zeta one', 3.0)])
    result = run_table(tmp_path)
    assert result.returncode != 0
    assert "copy_numbers.tsv does not match taxonomy.tsv: tax_id 1 is 'Zeta one'" in result.stderr


def test_a_sample_with_no_species_gets_an_empty_cell_table(tmp_path):
    write_inputs(tmp_path, {}, [tax_row(1, 'Alpha one', 'Alpha', 'F1')], [cn_row(1, 'Alpha one', 3.0)],
                 below=0, unassigned=12)
    result = run_table(tmp_path)
    assert result.returncode == 0, result.stderr
    comments, cell = read_cell(tmp_path / 's_cell.tsv')
    assert len(cell) == 0 and list(cell.columns)[:5] == copy_numbers.CELL_COLUMNS
    assert 'Unassigned (12 reads)' in comments[1]


def test_copy_numbers_and_cell_output_go_together(tmp_path):
    write_inputs(tmp_path, {'1': 1.0}, [tax_row(1, 'Alpha one', 'Alpha', 'F1')], [cn_row(1, 'Alpha one', 3.0)])
    result = run_table(tmp_path, '--copy-numbers', 'copy_numbers.tsv', cell=False)
    assert result.returncode != 0 and 'go together' in result.stderr


# ─── Edge cases (review of the branch) ───────────────────────────────────────

def test_blank_na_and_padded_species_names_pass_the_database_check(tmp_path):
    # copy_numbers.tsv is written from the taxonomy by build_copy_numbers, as the builders do
    write_inputs(tmp_path, {'1': 0.4, '2': 0.3, '3': 0.2, '4': 0.1}, [tax_row(1, 'x', 'X', 'F')], [])
    (tmp_path / 'taxonomy.tsv').write_text(
        'tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n'
        '1\tAlpha one\tAlpha\tF1\tO\tC\tP\tBacteria\n'
        '2\t\tAlpha\tF1\tO\tC\tP\tBacteria\n'
        '3\tNA\tAlpha\tF1\tO\tC\tP\tBacteria\n'
        '4\tBeta one \tAlpha\tF1\tO\tC\tP\tBacteria\n')
    (tmp_path / 'source.tsv').write_text('species\tcopy_number\nAlpha one\t4\nBeta one\t2\n')
    with open(os.devnull, 'w') as quiet:
        stdout, sys.stdout = sys.stdout, quiet
        try:
            copy_numbers.build_copy_numbers(str(tmp_path), tmp_path / 'taxonomy.tsv', table=str(tmp_path / 'source.tsv'))
        finally:
            sys.stdout = stdout
    result = run_table(tmp_path)
    assert result.returncode == 0, result.stderr
    _, cell = read_cell(tmp_path / 's_cell.tsv')
    assert len(cell) == 4 and cell['copy_number'].sum() == pytest.approx(4 + 3 + 3 + 2)


def test_a_tax_id_listed_twice_takes_its_first_row(tmp_path):
    taxonomy = [tax_row(5, 'Gamma one', 'Gamma', 'F1'), tax_row(5, 'Delta one', 'Delta', 'F1')]
    write_inputs(tmp_path, {'5': 1.0}, taxonomy, [cn_row(5, 'Gamma one', 3.0), cn_row(5, 'Delta one', 9.0)])
    result = run_table(tmp_path)
    assert result.returncode == 0, result.stderr
    _, cell = read_cell(tmp_path / 's_cell.tsv')
    assert cell[['species', 'copy_number']].values.tolist() == [['Gamma one', 3.0]]


def test_gtdb_metadata_skips_blank_and_short_lines(tmp_path):
    path = tmp_path / 'm.tsv'
    path.write_text('accession\tgtdb_taxonomy\tncbi_assembly_level\tssu_count\n'
                    f'RS_1\t{LIN.format(g="Alpha", s="Alpha one")}\tComplete Genome\t5\n\nRS_2\tcut\n')
    assert gtdb_metadata_values([path], lambda lin: build_gtdb_db.parse_gtdb_lineage(lin)['species']) == \
        {'Alpha one': (5.0, 1)}


def test_tax_ids_without_abundance_do_not_count():
    rows = copy_numbers.row_cells({'1': 0.0, '2': 0.0, '3': 0.5}, {'1': 0, '2': 0, '3': 1},
                                  {'1': (4.0, 'species'), '2': (8.0, 'species'), '3': (2.0, 'genus')})
    assert rows == {0: None, 1: (0.25, 0.5, 2.0, 'genus')}


def test_a_species_name_shared_by_tax_ids_counts_once_in_a_median():
    taxonomy = [tax_row(1, 'A a', 'A', 'F'), tax_row(2, 'A a', 'A', 'F'), tax_row(3, 'A a', 'A', 'F'),
                tax_row(4, 'A b', 'A', 'F'), tax_row(5, 'A c', 'A', 'F')]
    rows = fill_copy_numbers(taxonomy, {'A a': (7.0, 1), 'A b': (2.0, 1)}, 'test')
    assert by_tax(rows)['5'] == (4.5, 'genus')


# ─── Second review ───────────────────────────────────────────────────────────

def test_short_table_lines_are_skipped_not_a_crash(tmp_path):
    path = tmp_path / 't.tsv'
    path.write_text('species\tcopy_number\nA a\t3\nB b\n')
    assert copy_numbers.table_values(path) == {'A a': (3.0, '')}


@pytest.mark.parametrize('level', ['', 'manual'])
def test_an_unknown_level_is_refused_when_read(tmp_path, level):
    path = tmp_path / 'copy_numbers.tsv'
    path.write_text(f'tax_id\tspecies\tcopy_number\tlevel\n1\tA a\t3\t{level}\n')
    with pytest.raises(ValueError, match='line 2: level'):
        read_copy_numbers(path)


def test_a_source_that_fails_to_parse_leaves_no_old_file(tmp_path):
    write_taxonomy(tmp_path / 'taxonomy.tsv', TAXONOMY[:2])
    (tmp_path / 'copy_numbers.tsv').write_text('old build\n')
    bad = tmp_path / 'bad.tsv'
    bad.write_text('species\tcopy_number\nAlpha one\t-1\n')
    with pytest.raises(ValueError):
        copy_numbers.build_copy_numbers(str(tmp_path), tmp_path / 'taxonomy.tsv', table=str(bad))
    assert not (tmp_path / 'copy_numbers.tsv').exists()


def test_a_group_of_one_has_exactly_its_rows_cell_share(tmp_path):
    # Abundances listed out of the species table's (alphabetical) order, with values whose sum
    # depends on the order of addition
    taxonomy = [tax_row(i, f'Sp{c} x', f'G{c}', f'F{c}') for i, c in enumerate('abcdef', 1)]
    freq = {'6': 0.1, '3': 0.2, '5': 0.3, '1': 0.15, '4': 0.05, '2': 0.2}
    groups = [{'tax_ids': [t], 'shared_reads': 0.0, 'max_shared': 0.0} for t in freq]
    cns = [cn_row(i, f'Sp{c} x', [3.0, 7.0, 11.0, 13.0, 17.0, 19.0][i - 1]) for i, c in enumerate('abcdef', 1)]
    write_inputs(tmp_path, freq, taxonomy, cns, groups=groups)
    result = run_table(tmp_path, '--groups-output', 's_groups')
    assert result.returncode == 0, result.stderr
    _, cell = read_cell(tmp_path / 's_cell.tsv')
    g = pd.read_csv(tmp_path / 's_groups.tsv', sep='\t').dropna(subset=['n_members'])
    by_species = dict(zip(cell['species'], cell['abundance']))
    assert all(v == by_species[m] for m, v in zip(g['members'], g['cell_abundance']))


# ─── Third review ────────────────────────────────────────────────────────────

@pytest.mark.parametrize('source', ['bad.tsv', 'missing.tsv'])
def test_build_database_reads_the_copy_number_source_before_building(tmp_path, source):
    build_database_inputs(tmp_path)
    (tmp_path / 'lineages.tsv').write_text(
        'tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n'
        '287\tAlpha one\tAlpha\tF1\tO\tC\tP\tBacteria\n1001\tAlpha two\tAlpha\tF1\tO\tC\tP\tBacteria\n')
    (tmp_path / 'bad.tsv').write_text('species\tcopy_number\nAlpha one\t0\n')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'main.py'), 'build-database', 'mydb', '--sequences', 'seqs.fa',
         '--seq2tax', 'seq2tax.tsv', '--taxonomy-list', 'lineages.tsv', '--copy-numbers', source],
        cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert not (tmp_path / 'mydb').exists()


def test_build_gtdb_db_reads_the_metadata_before_building(tmp_path):
    ssu, meta, out = tmp_path / 'ssu.fna', tmp_path / 'bac120_metadata_r0.tsv', tmp_path / 'db'
    ssu_fasta(ssu)
    meta.write_text('accession\tgtdb_taxonomy\n')
    with pytest.raises(ValueError):
        build_gtdb_db.build_database(str(ssu), 'db_gtdb', str(out), gtdb_metadata=[str(meta)])
    assert not out.exists()


def test_rrndb_refuses_a_repeated_species_taxid(tmp_path):
    path = tmp_path / 'pantaxa.tsv'
    path.write_text('taxid\trank\tname\tchildcount\tmedian\n'
                    '287\tspecies\tA a\t3\t4\n287\tspecies\tA b\t2\t6\n')
    with pytest.raises(ValueError, match="taxid '287' is listed more than once"):
        copy_numbers.rrndb_pantaxa_values(path)


def test_only_the_samples_tax_ids_are_read(tmp_path):
    path = tmp_path / 'copy_numbers.tsv'
    path.write_text('tax_id\tspecies\tcopy_number\tlevel\n1\tA a\t3\tspecies\n2\tB b\t4\tmanual\n')
    assert read_copy_numbers(path, {'1'}) == {'1': ('A a', 3.0, 'species')}
    with pytest.raises(ValueError, match='line 3: level'):
        read_copy_numbers(path)
