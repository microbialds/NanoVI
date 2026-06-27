"""Shared pytest fixtures for NanoVI unit tests."""
import math
import pytest
import pandas as pd


@pytest.fixture
def mock_nodes_dict():
    """Minimal NCBI-like taxonomy tree rooted at taxid '1'.

    Lineage: root(1) -> Bacteria(2) -> Proteobacteria(1224)
             -> Gammaproteobacteria(1236) -> Pseudomonadales(72274)
             -> Pseudomonadaceae(135621) -> Pseudomonas(286)
             -> Pseudomonas aeruginosa(287)
    """
    return {
        '1':      ('1',      'no rank'),
        '2':      ('1',      'superkingdom'),
        '1224':   ('2',      'phylum'),
        '1236':   ('1224',   'class'),
        '72274':  ('1236',   'order'),
        '135621': ('72274',  'family'),
        '286':    ('135621', 'genus'),
        '287':    ('286',    'species'),
    }


@pytest.fixture
def mock_names_dict():
    return {
        '1':      'root',
        '2':      'Bacteria',
        '1224':   'Proteobacteria',
        '1236':   'Gammaproteobacteria',
        '72274':  'Pseudomonadales',
        '135621': 'Pseudomonadaceae',
        '286':    'Pseudomonas',
        '287':    'Pseudomonas aeruginosa',
    }


@pytest.fixture
def abundance_df():
    """DataFrame with two E. coli strains that should merge and one Pseudomonas."""
    return pd.DataFrame({
        'species':      ['Escherichia coli O157', 'Escherichia coli K12', 'Pseudomonas aeruginosa'],
        'genus':        ['Escherichia', 'Escherichia', 'Pseudomonas'],
        'family':       ['Enterobacteriaceae', 'Enterobacteriaceae', 'Pseudomonadaceae'],
        'order':        ['Enterobacterales', 'Enterobacterales', 'Pseudomonadales'],
        'class':        ['Gammaproteobacteria'] * 3,
        'phylum':       ['Proteobacteria'] * 3,
        'superkingdom': ['Bacteria'] * 3,
        'abundance':    [0.3, 0.2, 0.5],
    })


@pytest.fixture
def simple_log_p_rgs():
    """Two reads each mapping unambiguously to one species."""
    return {
        'read1': ([1], [math.log(0.9)]),
        'read2': ([2], [math.log(0.8)]),
    }


@pytest.fixture
def ambiguous_log_p_rgs():
    """One read mapping to two species with different log-scores."""
    return {
        'read1': ([1, 2], [math.log(0.9), math.log(0.1)]),
    }


@pytest.fixture
def uniform_freq():
    return {1: 0.5, 2: 0.5}
