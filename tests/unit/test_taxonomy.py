"""Unit tests for taxonomy utility functions (bin/taxonomy.py).

Covers lineage resolution, species-name unification, and rank collapse.
File I/O functions that require NCBI dump files are tested with temporary
files or skipped in favour of integration tests.
"""
import os
import sys
import tempfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
from taxonomy import (
    TAXONOMY_RANKS,
    collapse_rank,
    get_species_tid,
    freq_to_lineage_df,
    lineage_dict_from_tid,
    species_row_keys,
    unify_species_abundance,
)


class TestLineageDictFromTid:
    def test_returns_tuple(self, mock_nodes_dict, mock_names_dict):
        result = lineage_dict_from_tid('287', mock_nodes_dict, mock_names_dict)
        assert isinstance(result, tuple)

    def test_length_equals_ranks_plus_taxid(self, mock_nodes_dict, mock_names_dict):
        result = lineage_dict_from_tid('287', mock_nodes_dict, mock_names_dict)
        assert len(result) == len(TAXONOMY_RANKS) + 1

    def test_species_name_in_result(self, mock_nodes_dict, mock_names_dict):
        result = lineage_dict_from_tid('287', mock_nodes_dict, mock_names_dict)
        assert 'Pseudomonas aeruginosa' in result

    def test_superkingdom_in_result(self, mock_nodes_dict, mock_names_dict):
        result = lineage_dict_from_tid('287', mock_nodes_dict, mock_names_dict)
        assert 'Bacteria' in result

    def test_genus_in_result(self, mock_nodes_dict, mock_names_dict):
        result = lineage_dict_from_tid('287', mock_nodes_dict, mock_names_dict)
        assert 'Pseudomonas' in result

    def test_taxid_is_first_element(self, mock_nodes_dict, mock_names_dict):
        result = lineage_dict_from_tid('287', mock_nodes_dict, mock_names_dict)
        assert result[0] == '287'


class TestGetSpeciesTid:
    def test_species_taxid_returns_itself(self, mock_nodes_dict):
        result = get_species_tid('287', mock_nodes_dict)
        assert result == '287'

    def test_genus_taxid_returns_genus(self, mock_nodes_dict):
        result = get_species_tid('286', mock_nodes_dict)
        assert result == '286'

    def test_unknown_taxid_raises_value_error(self, mock_nodes_dict):
        with pytest.raises(ValueError, match="not found"):
            get_species_tid('99999', mock_nodes_dict)


class TestUnifySpeciesAbundance:
    def test_full_names_kept_by_default(self, abundance_df):
        # Species names are no longer cut to two words (F17)
        result = unify_species_abundance(abundance_df.copy())
        assert sorted(result['species']) == sorted(abundance_df['species'].unique())

    def test_strain_names_merged_at_two_words(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy(), binomial=True)
        ecoli = result[result['species'] == 'Escherichia coli']
        assert len(ecoli) == 1

    def test_merged_abundance_is_sum_of_strains(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy(), binomial=True)
        ecoli_abundance = result.loc[result['species'] == 'Escherichia coli', 'abundance'].iloc[0]
        assert abs(ecoli_abundance - 0.5) < 1e-9

    @pytest.mark.parametrize('name,expected', [
        ('Escherichia coli O157', 'Escherichia coli'),
        ('Escherichia coli', 'Escherichia coli'),
        ('Candidatus Pelagibacter ubique', 'Candidatus Pelagibacter ubique'),
        ('Candidatus Pelagibacter ubique HTCC1062', 'Candidatus Pelagibacter ubique'),
        ('Pseudomonas_E aeruginosa', 'Pseudomonas_E aeruginosa'),
    ])
    def test_binomial_name_keeps_candidatus_whole(self, name, expected):
        from taxonomy import binomial_name
        assert binomial_name(name) == expected

    def test_candidatus_species_stay_apart(self):
        df = pd.DataFrame({'species': ['Candidatus Pelagibacter ubique', 'Candidatus Pelagibacter bermudensis'],
                           'abundance': [0.6, 0.4]})
        result = unify_species_abundance(df, binomial=True)
        assert sorted(result['species']) == ['Candidatus Pelagibacter bermudensis', 'Candidatus Pelagibacter ubique']

    def test_unambiguous_species_unchanged(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy())
        pa = result[result['species'] == 'Pseudomonas aeruginosa']
        assert len(pa) == 1
        assert abs(pa.iloc[0]['abundance'] - 0.5) < 1e-9

    def test_all_species_names_at_most_two_words(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy(), binomial=True)
        for name in result['species']:
            assert len(str(name).split()) <= 2

    def test_raises_without_species_column(self):
        df = pd.DataFrame({'genus': ['Escherichia'], 'abundance': [1.0]})
        with pytest.raises(ValueError, match="species"):
            unify_species_abundance(df)

    def test_missing_higher_rank_filled_with_unknown(self):
        df = pd.DataFrame({'species': ['Escherichia coli'], 'abundance': [1.0]})
        result = unify_species_abundance(df)
        assert 'genus' in result.columns


class TestCollapseRank:
    def _write_tmp_tsv(self, df):
        f = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        df.to_csv(f.name, sep='\t', index=False)
        f.close()
        return f.name

    def test_collapse_to_genus_produces_file(self, abundance_df):
        tmp = self._write_tmp_tsv(abundance_df)
        out = tmp.replace('.tsv', '-genus.tsv')
        try:
            collapse_rank(tmp, 'genus')
            assert os.path.exists(out)
        finally:
            for p in [tmp, out]:
                if os.path.exists(p):
                    os.remove(p)

    def test_collapsed_genus_has_abundance_column(self, abundance_df):
        tmp = self._write_tmp_tsv(abundance_df)
        out = tmp.replace('.tsv', '-genus.tsv')
        try:
            collapse_rank(tmp, 'genus')
            result = pd.read_csv(out, sep='\t')
            assert 'abundance' in result.columns
        finally:
            for p in [tmp, out]:
                if os.path.exists(p):
                    os.remove(p)

    def test_invalid_rank_raises_value_error(self, abundance_df):
        tmp = self._write_tmp_tsv(abundance_df)
        try:
            with pytest.raises(ValueError, match="rank"):
                collapse_rank(tmp, 'domain')
        finally:
            os.remove(tmp)

    def test_collapse_reduces_row_count(self, abundance_df):
        # Two E.coli strains → after collapsing to genus, only one Escherichia row
        tmp = self._write_tmp_tsv(abundance_df)
        out = tmp.replace('.tsv', '-genus.tsv')
        try:
            collapse_rank(tmp, 'genus')
            result = pd.read_csv(out, sep='\t')
            n_escherichia = (result.index.get_level_values('genus') == 'Escherichia').sum() \
                if 'genus' in result.index.names \
                else (result.get('genus', pd.Series()) == 'Escherichia').sum()
            assert n_escherichia <= 1
        finally:
            for p in [tmp, out]:
                if os.path.exists(p):
                    os.remove(p)


# ─── species_row_keys (F46) ─────────────────────────────────────────────────

@pytest.mark.parametrize('binomial', [False, True])
def test_species_row_keys_partition_tax_ids_as_the_table_does(tmp_path, binomial):
    # 1 and 2 share name and lineage; 3 has 1's name in another genus; 4 is a strain
    # of 1 (merged only with binomial names); 5 has no species name; 6 is missing
    taxonomy = pd.DataFrame(
        [['1', 'Escherichia coli', 'Escherichia', 'Enterobacteriaceae'],
         ['2', 'Escherichia coli', 'Escherichia', 'Enterobacteriaceae'],
         ['3', 'Escherichia coli', 'Shigella', 'Enterobacteriaceae'],
         ['4', 'Escherichia coli O157', 'Escherichia', 'Enterobacteriaceae'],
         ['5', '', 'Escherichia', 'Enterobacteriaceae']],
        columns=['tax_id', 'species', 'genus', 'family']).set_index('tax_id')
    tax_ids = [1, 2, 3, 4, 5, 6]
    keys = species_row_keys(tax_ids, taxonomy, binomial)
    # Each tax ID's abundance is a distinct power of two, so a row's abundance names its tax IDs
    freq = {t: 2.0 ** -t for t in tax_ids}
    table = freq_to_lineage_df(freq, str(tmp_path / 'out'), taxonomy, 10, 0, binomial=binomial)
    rows = [frozenset(t for t in tax_ids if int(a * 2 ** 6) >> (6 - t) & 1)
            for a in table['abundance'][:-2]]
    by_key = {}
    for t, k in zip(tax_ids, keys):
        by_key.setdefault(k, set()).add(t)
    assert sorted(map(sorted, rows)) == sorted(map(sorted, by_key.values()))
    assert (keys[0] == keys[3]) == binomial and keys[0] == keys[1] and keys[0] != keys[2]
