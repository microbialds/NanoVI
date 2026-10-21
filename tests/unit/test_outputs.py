"""Unit tests for the output tables and the GTDB database builder.

Covers freq_to_lineage_df and combine_outputs (bin/taxonomy.py),
build_consolidated_tables (bin/generate_phyloseq_tables.py) and the parsers
and builder in bin/build_gtdb_db.py. Tests marked xfail describe known bugs
(F11, F12) and are expected to pass once those are fixed.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
import build_gtdb_db
import taxonomy
from generate_phyloseq_tables import build_consolidated_tables, parse_sample_file
from taxonomy import combine_outputs, freq_to_lineage_df

RANKS = ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom']


@pytest.fixture
def taxonomy_df():
    """taxonomy.tsv as WRITE_OUTPUT reads it: string tax_id index, all columns str."""
    rows = [
        ('1', 'Escherichia coli', 'Escherichia', 'Enterobacteriaceae', 'Enterobacterales',
         'Gammaproteobacteria', 'Pseudomonadota', 'Bacteria'),
        ('2', 'Escherichia coli O157', 'Escherichia', 'Enterobacteriaceae', 'Enterobacterales',
         'Gammaproteobacteria', 'Pseudomonadota', 'Bacteria'),
        ('3', 'Staphylococcus aureus', 'Staphylococcus', 'Staphylococcaceae', 'Bacillales',
         'Bacilli', 'Bacillota', 'Bacteria'),
    ]
    df = pd.DataFrame(rows, columns=['tax_id', 'species'] + RANKS, dtype=str)
    return df.set_index('tax_id')


def sample(rows, counts=True):
    cols = ['species', 'abundance'] + RANKS + (['estimated_counts'] if counts else [])
    return pd.DataFrame(rows, columns=cols)


def write_sample(path, rows, counts=True):
    sample(rows, counts).to_csv(path, sep='\t', index=False)


def lineage(genus, family='F', order='O', cls='C', phylum='P', kingdom='Bacteria'):
    return [genus, family, order, cls, phylum, kingdom]


# ─── freq_to_lineage_df ──────────────────────────────────────────────────────

class TestFreqToLineageDf:
    def test_writes_tsv_with_expected_columns(self, tmp_path, taxonomy_df):
        out = tmp_path / 'S1_rel-abundance'
        df = freq_to_lineage_df({'1': 0.6, '3': 0.4}, str(out), taxonomy_df, 100, 5, counts=True)
        assert (tmp_path / 'S1_rel-abundance.tsv').exists()
        expected_cols = ['species', 'abundance'] + RANKS + ['estimated_counts']
        assert list(df.columns) == expected_cols
        on_disk = pd.read_csv(tmp_path / 'S1_rel-abundance.tsv', sep='\t')
        assert list(on_disk.columns) == expected_cols

    def test_no_counts_column_without_counts(self, tmp_path, taxonomy_df):
        df = freq_to_lineage_df({'1': 1.0}, str(tmp_path / 's'), taxonomy_df, 100, 0, counts=False)
        assert 'estimated_counts' not in df.columns

    def test_strains_merge_into_two_word_species(self, tmp_path, taxonomy_df):
        df = freq_to_lineage_df({'1': 0.3, '2': 0.2, '3': 0.5}, str(tmp_path / 's'), taxonomy_df, 100, 0)
        ecoli = df[df['species'] == 'Escherichia coli']
        assert len(ecoli) == 1
        assert ecoli['abundance'].iloc[0] == pytest.approx(0.5)

    def test_estimated_counts_from_assigned_count(self, tmp_path, taxonomy_df):
        df = freq_to_lineage_df({'1': 0.6, '3': 0.4}, str(tmp_path / 's'), taxonomy_df, 200, 0, counts=True)
        counts = dict(zip(df['species'], df['estimated_counts']))
        assert counts['Escherichia coli'] == 120
        assert counts['Staphylococcus aureus'] == 80

    def test_assigned_abundances_sum_to_one(self, tmp_path, taxonomy_df):
        df = freq_to_lineage_df({'1': 0.25, '2': 0.25, '3': 0.5}, str(tmp_path / 's'), taxonomy_df, 100, 0)
        assert df['abundance'].sum() == pytest.approx(1.0)

    def test_unassigned_row_kept_with_count(self, tmp_path, taxonomy_df):
        df = freq_to_lineage_df({'1': 0.6, '3': 0.4}, str(tmp_path / 's'), taxonomy_df, 100, 7, counts=True)
        unassigned = df[df['species'] == 'Unassigned']
        assert len(unassigned) == 1
        assert unassigned['estimated_counts'].iloc[0] == 7
        assert unassigned['abundance'].iloc[0] == 0
        assert df['species'].iloc[-1] == 'Unassigned'
        on_disk = pd.read_csv(tmp_path / 's.tsv', sep='\t')
        assert on_disk['species'].iloc[-1] == 'Unassigned'

    def test_unassigned_row_without_counts(self, tmp_path, taxonomy_df):
        df = freq_to_lineage_df({'1': 1.0}, str(tmp_path / 's'), taxonomy_df, 100, 7, counts=False)
        assert list(df['species']) == ['Escherichia coli', 'Unassigned']

    def test_taxid_missing_from_taxonomy_kept(self, tmp_path, taxonomy_df, capsys):
        df = freq_to_lineage_df({'1': 0.6, '99': 0.4}, str(tmp_path / 's'), taxonomy_df, 100, 0)
        assert df['abundance'].sum() == pytest.approx(1.0)
        unknown = df[df['species'] == 'Unknown taxid 99']
        assert len(unknown) == 1
        assert unknown['abundance'].iloc[0] == pytest.approx(0.4)
        assert 'tax_id 99 is missing from taxonomy.tsv' in capsys.readouterr().err

    def test_taxid_without_species_name_kept(self, tmp_path, taxonomy_df, capsys):
        taxonomy_df.loc['3', 'species'] = np.nan
        df = freq_to_lineage_df({'1': 0.6, '3': 0.4}, str(tmp_path / 's'), taxonomy_df, 100, 0)
        row = df[df['species'] == 'Unknown taxid 3']
        assert row['genus'].iloc[0] == 'Staphylococcus'
        assert 'tax_id 3 has no species name' in capsys.readouterr().err

    def test_row_with_empty_rank_kept(self, tmp_path, taxonomy_df):
        taxonomy_df.loc['3', 'family'] = np.nan
        df = freq_to_lineage_df({'1': 0.6, '3': 0.4}, str(tmp_path / 's'), taxonomy_df, 100, 0)
        assert 'Staphylococcus aureus' in set(df['species'])
        assert df.loc[df['species'] == 'Staphylococcus aureus', 'family'].iloc[0] == ''
        assert df['abundance'].sum() == pytest.approx(1.0)


# ─── build_consolidated_tables ───────────────────────────────────────────────

class TestBuildConsolidatedTables:
    def test_tables_shapes_and_values(self):
        s1 = sample([['Escherichia coli', 0.7] + lineage('Escherichia') + [70],
                     ['Staphylococcus aureus', 0.3] + lineage('Staphylococcus') + [30]])
        s2 = sample([['Escherichia coli', 1.0] + lineage('Escherichia') + [50]])
        otu_abund, otu_counts, tax, meta = build_consolidated_tables([('S1', s1), ('S2', s2)])

        assert list(otu_abund.columns) == ['S1', 'S2']
        assert list(otu_abund.index) == ['Escherichia coli', 'Staphylococcus aureus', 'Unassigned']
        assert otu_abund.loc['Staphylococcus aureus', 'S2'] == 0.0
        assert otu_counts.loc['Escherichia coli', 'S1'] == 70
        assert list(tax.columns) == RANKS
        assert tax.loc['Escherichia coli', 'genus'] == 'Escherichia'
        assert meta.loc['S1', 'n_taxa_detected'] == 2
        assert meta.loc['S1', 'assigned_counts'] == 100
        assert meta.loc['S2', 'total_abundance'] == pytest.approx(1.0)

    def test_nan_species_becomes_unassigned(self):
        s1 = sample([['Escherichia coli', 1.0] + lineage('Escherichia') + [90],
                     [np.nan, 0.0] + [np.nan] * 6 + [10]])
        otu_abund, otu_counts, _, meta = build_consolidated_tables([('S1', s1)])
        assert otu_counts.loc['Unassigned', 'S1'] == 10
        assert meta.loc['S1', 'unassigned_counts'] == 10

    def test_unassigned_counts_from_unassigned_row(self, tmp_path, taxonomy_df):
        freq_to_lineage_df({'1': 0.6, '3': 0.4}, str(tmp_path / 'S1_rel-abundance'),
                           taxonomy_df, 100, 30, counts=True)
        _, otu_counts, tax, meta = build_consolidated_tables(
            [parse_sample_file(tmp_path / 'S1_rel-abundance.tsv')])
        assert otu_counts.loc['Unassigned', 'S1'] == 30
        assert meta.loc['S1', 'unassigned_counts'] == 30
        assert meta.loc['S1', 'assigned_counts'] == 100
        assert meta.loc['S1', 'n_taxa_detected'] == 2
        assert list(otu_counts.index)[-1] == 'Unassigned'
        assert (tax.loc['Unassigned'] == '').all()

    def test_parse_sample_file(self, tmp_path):
        path = tmp_path / 'S9_rel-abundance.tsv'
        write_sample(path, [['Escherichia coli', 1.0] + lineage('Escherichia') + [12]])
        name, df = parse_sample_file(path)
        assert name == 'S9'
        assert df['abundance'].dtype == float

    @pytest.mark.xfail(strict=True, reason='F11: duplicate species in one sample overwrite instead of summing')
    def test_duplicate_species_are_summed(self):
        s1 = sample([['Bacillus cereus', 0.3] + lineage('Bacillus') + [30],
                     ['Bacillus cereus', 0.2] + lineage('Bacillus', family='Other') + [20],
                     ['Escherichia coli', 0.5] + lineage('Escherichia') + [50]])
        otu_abund, otu_counts, _, _ = build_consolidated_tables([('S1', s1)])
        assert otu_abund.loc['Bacillus cereus', 'S1'] == pytest.approx(0.5)
        assert otu_counts.loc['Bacillus cereus', 'S1'] == 50


# ─── combine_outputs ─────────────────────────────────────────────────────────

class TestCombineOutputs:
    @pytest.fixture
    def two_samples(self, tmp_path):
        write_sample(tmp_path / 'A_rel-abundance.tsv',
                     [['Escherichia coli', 0.6] + lineage('Escherichia') + [60],
                      ['Staphylococcus aureus', 0.4] + lineage('Staphylococcus') + [40]])
        write_sample(tmp_path / 'B_rel-abundance.tsv',
                     [['Escherichia coli', 1.0] + lineage('Escherichia') + [80]])
        return tmp_path

    def test_combines_samples_at_genus(self, two_samples):
        df = combine_outputs(str(two_samples), 'genus')
        assert (two_samples / 'vi-combined-genus.tsv').exists()
        assert set(df.columns) >= {'genus', 'A', 'B'}
        row = df.set_index('genus').loc['Escherichia']
        assert float(row['A']) == pytest.approx(0.6)
        assert float(row['B']) == pytest.approx(1.0)

    def test_count_table(self, two_samples):
        df = combine_outputs(str(two_samples), 'species', count_table=True)
        assert (two_samples / 'vi-combined-species-counts.tsv').exists()
        assert float(df.set_index('species').loc['Escherichia coli', 'B']) == 80

    def test_split_files(self, two_samples):
        combine_outputs(str(two_samples), 'genus', split_files=True)
        assert (two_samples / 'vi-combined-taxonomy-genus.tsv').exists()
        assert (two_samples / 'vi-combined-abundance-genus.tsv').exists()

    def test_ignores_other_files(self, two_samples):
        (two_samples / 'notes.tsv').write_text('x\n1\n')
        df = combine_outputs(str(two_samples), 'genus')
        assert 'notes' not in df.columns

    @pytest.mark.xfail(strict=True, reason='F12: --counts on tables without estimated_counts raises a bare KeyError')
    def test_missing_metric_column_names_the_file(self, tmp_path):
        write_sample(tmp_path / 'C_rel-abundance.tsv',
                     [['Escherichia coli', 1.0] + lineage('Escherichia')], counts=False)
        with pytest.raises(Exception) as excinfo:
            combine_outputs(str(tmp_path), 'species', count_table=True)
        assert 'C_rel-abundance.tsv' in str(excinfo.value)

    @pytest.mark.xfail(strict=True, reason='F12: taxa absent from a sample are left blank instead of 0')
    def test_absent_taxa_are_zero(self, two_samples):
        df = combine_outputs(str(two_samples), 'genus')
        value = df.set_index('genus').loc['Staphylococcus', 'B']
        assert not pd.isna(value) and float(value) == 0.0

    @pytest.mark.xfail(strict=True, reason='F12: files are read in os.listdir order, so sample columns are not sorted')
    def test_sample_columns_sorted(self, tmp_path, monkeypatch):
        for name in ['C', 'A', 'B']:
            write_sample(tmp_path / f'{name}_rel-abundance.tsv',
                         [['Escherichia coli', 1.0] + lineage('Escherichia') + [10]])
        real_listdir = os.listdir
        monkeypatch.setattr(taxonomy.os, 'listdir', lambda p: sorted(real_listdir(p), reverse=True))
        df = combine_outputs(str(tmp_path), 'genus')
        samples = [c for c in df.columns if c in {'A', 'B', 'C'}]
        assert samples == ['A', 'B', 'C']


# ─── build_gtdb_db ───────────────────────────────────────────────────────────

HEADER = ('>RS_GCF_000005845.2~NC_000913.3 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;'
          'o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli '
          '[location=4035531..4037072] [ssu_len=1542] [contig_len=4641652]')


class TestGtdbParsers:
    @pytest.mark.parametrize('name,expected', [
        ('Pseudomonas_E aeruginosa', 'Pseudomonas aeruginosa'),
        ('Listeria monocytogenes_B', 'Listeria monocytogenes'),
        ('Escherichia coli', 'Escherichia coli'),
    ])
    def test_normalize_gtdb_name(self, name, expected):
        assert build_gtdb_db.normalize_gtdb_name(name) == expected

    def test_parse_gtdb_lineage(self):
        lin = build_gtdb_db.parse_gtdb_lineage('d__Bacteria;p__Bacillota_A;g__Clostridium_P;s__Clostridium_P sp1')
        assert lin['superkingdom'] == 'Bacteria'
        assert lin['phylum'] == 'Bacillota'
        assert lin['genus'] == 'Clostridium'
        assert lin['species'] == 'Clostridium sp1'
        assert lin['family'] == ''

    def test_parse_ssu_header(self):
        genome_id, lin = build_gtdb_db.parse_ssu_header(HEADER)
        assert genome_id == 'RS_GCF_000005845.2'
        assert lin['species'] == 'Escherichia coli'
        assert lin['order'] == 'Enterobacterales'

    @pytest.mark.parametrize('header', ['>only_an_id', '>id some description without lineage'])
    def test_parse_ssu_header_unparseable(self, header):
        assert build_gtdb_db.parse_ssu_header(header) == (None, None)

    @pytest.mark.parametrize('gz', [False, True])
    def test_fasta_records(self, tmp_path, gz):
        text = '>a desc\nACGT\nAC\n>b\nGG\n'
        path = tmp_path / ('x.fna.gz' if gz else 'x.fna')
        if gz:
            with gzip.open(path, 'wt') as fh:
                fh.write(text)
        else:
            path.write_text(text)
        assert list(build_gtdb_db.fasta_records(str(path))) == [('>a desc', 'ACGTAC'), ('>b', 'GG')]


class TestBuildDatabase:
    def test_builds_fasta_and_taxonomy(self, tmp_path):
        lin = 'd__Bacteria;p__P;c__C;o__O;f__F;g__{g};s__{s}'
        seq_a, seq_b = 'ACGT' * 300, 'TTGCA' * 240
        records = [
            (f'>g1~c1 {lin.format(g="Zeta", s="Zeta alpha")} [ssu_len=1200]', seq_a),
            (f'>g2~c2 {lin.format(g="Alpha_A", s="Alpha_A beta")} [ssu_len=1200]', seq_b),
            (f'>g3~c3 {lin.format(g="Zeta", s="Zeta alpha")} [ssu_len=1200]', seq_a.lower()),  # duplicate
            (f'>g4~c4 {lin.format(g="Zeta", s="Zeta alpha")} [ssu_len=40]', 'ACGT' * 10),  # too short
            ('>g5~c5 no lineage here', seq_a),
        ]
        ssu = tmp_path / 'ssu.fna'
        ssu.write_text(''.join(f'{h}\n{s}\n' for h, s in records))
        out = tmp_path / 'db'
        build_gtdb_db.build_database(str(ssu), 'testdb', str(out), min_length=900)

        tax = pd.read_csv(out / 'taxonomy.tsv', sep='\t', dtype=str)
        assert list(tax.columns) == build_gtdb_db.RANKS_OUT
        assert list(tax['species']) == ['Alpha beta', 'Zeta alpha']
        assert list(tax['tax_id']) == ['1', '2']

        headers = [h for h, _ in build_gtdb_db.fasta_records(str(out / 'species_taxid.fasta'))]
        assert headers == ['>2:testdb:1', '>1:testdb:2']
