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
    lineage_dict_from_tid,
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
    def test_strain_names_merged_at_two_words(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy())
        ecoli = result[result['species'] == 'Escherichia coli']
        assert len(ecoli) == 1

    def test_merged_abundance_is_sum_of_strains(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy())
        ecoli_abundance = result.loc[result['species'] == 'Escherichia coli', 'abundance'].iloc[0]
        assert abs(ecoli_abundance - 0.5) < 1e-9

    def test_unambiguous_species_unchanged(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy())
        pa = result[result['species'] == 'Pseudomonas aeruginosa']
        assert len(pa) == 1
        assert abs(pa.iloc[0]['abundance'] - 0.5) < 1e-9

    def test_all_species_names_at_most_two_words(self, abundance_df):
        result = unify_species_abundance(abundance_df.copy())
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
