"""Unit tests for bin/identical_16s.py: species that share an identical 16S (D23, D24)."""
import csv
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
import identical_16s

S1, S2, S3, S4, S5, S6 = ('ACGT' * 300, 'TTGCA' * 240, 'GGCAT' * 240, 'CATG' * 300, 'AAGGT' * 240, 'TGCAG' * 240)


def write_db(tmp_path, records, lineages):
    """records: (tax_id, sequence); lineages: tax_id -> (species, genus, family, phylum)."""
    fasta = tmp_path / 'species_taxid.fasta'
    fasta.write_text(''.join(f'>{t}:db:{i}\n{seq}\n' for i, (t, seq) in enumerate(records, 1)))
    tax = tmp_path / 'taxonomy.tsv'
    with open(tax, 'w') as fh:
        fh.write('tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n')
        for t, (sp, genus, family, phylum) in lineages.items():
            fh.write(f'{t}\t{sp}\t{genus}\t{family}\tO\tC\t{phylum}\tBacteria\n')
    return str(fasta), str(tax)


@pytest.fixture
def db(tmp_path):
    lineages = {
        '1': ('Alpha one', 'Alpha', 'F1', 'P1'),
        '2': ('Alpha two', 'Alpha', 'F1', 'P1'),
        '3': ('Alpha three', 'Alpha', 'F1', 'P1'),
        '4': ('Beta one', 'Beta', 'F2', 'P2'),
        '5': ('Delta one', 'Delta_E', 'F3', 'P1'),
        '6': ('Delta two', 'Delta', 'F3', 'P1'),
        '7': ('Nogenus one', '', 'F4', 'P1'),
        '8': ('Gamma one', 'Gamma', 'F5', 'P1'),
    }
    records = [
        ('1', S1), ('2', S1.lower()),  # Alpha one and two share S1 (case does not matter)
        ('2', S2), ('3', S2),          # Alpha two and three share S2: one group of three
        ('1', S3), ('4', S3),          # S3 is in two genera: cross-genus, not a link
        ('5', S4), ('6', S4),          # Delta_E is Delta once the GTDB suffix is stripped
        ('7', S5), ('8', S5),          # Nogenus one has no genus: left out
        ('8', S6),                     # unique sequence
    ]
    return write_db(tmp_path, records, lineages)


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def test_groups_link_species_of_one_genus(db):
    rows, _ = identical_16s.identical_reports(*db)
    by_species = {r['species']: r for r in rows}
    assert set(by_species) == {'Alpha one', 'Alpha two', 'Alpha three', 'Delta one', 'Delta two'}
    # The largest group comes first
    assert {by_species[s]['group'] for s in ('Alpha one', 'Alpha two', 'Alpha three')} == {'G1'}
    assert by_species['Alpha one']['group_size'] == 3
    assert {by_species[s]['group'] for s in ('Delta one', 'Delta two')} == {'G2'}
    # identical_16s_species names the species a species shares a sequence with directly
    assert by_species['Alpha one']['identical_16s_species'] == 'Alpha two'
    assert by_species['Alpha two']['identical_16s_species'] == 'Alpha one; Alpha three'
    assert by_species['Delta one']['identical_16s_species'] == 'Delta two'


def test_sequences_in_several_genera_are_reported_apart(db):
    _, cross = identical_16s.identical_reports(*db)
    assert cross == [{
        'md5': cross[0]['md5'], 'species_count': 2, 'sequence_ids': '1:db:5;4:db:6', 'tax_ids': '1;4',
        'species': 'Alpha one;Beta one', 'genera': 'Alpha;Beta', 'families': 'F1;F2', 'phyla': 'P1;P2',
    }]


def test_script_writes_both_files(db, tmp_path):
    identical_16s.write_reports(*db, str(tmp_path / 'identical_16s.tsv'), str(tmp_path / 'cross.tsv'))
    groups = read_tsv(tmp_path / 'identical_16s.tsv')
    assert list(groups[0]) == identical_16s.GROUP_COLUMNS
    assert len(groups) == 5
    assert [r['tax_ids'] for r in read_tsv(tmp_path / 'cross.tsv')] == ['1;4']


def test_no_shared_sequences_gives_header_only(tmp_path):
    fasta, tax = write_db(tmp_path, [('1', S1), ('2', S2)],
                          {'1': ('Alpha one', 'Alpha', 'F', 'P'), '2': ('Alpha two', 'Alpha', 'F', 'P')})
    identical_16s.write_reports(fasta, tax, str(tmp_path / 'g.tsv'), str(tmp_path / 'c.tsv'))
    assert (tmp_path / 'g.tsv').read_text() == '\t'.join(identical_16s.GROUP_COLUMNS) + '\n'
    assert read_tsv(tmp_path / 'c.tsv') == []


def test_a_contaminant_in_another_genus_does_not_hide_an_identity(tmp_path):
    # Alpha one and Alpha two share S1; a Beta species carries S1 too (a contaminant contig)
    fasta, tax = write_db(tmp_path, [('1', S1), ('2', S1), ('3', S1)], {
        '1': ('Alpha one', 'Alpha', 'F1', 'P1'),
        '2': ('Alpha two', 'Alpha_B', 'F1', 'P1'),
        '3': ('Beta one', 'Beta', 'F2', 'P2'),
    })
    rows, cross = identical_16s.identical_reports(fasta, tax)
    assert {r['species']: r['identical_16s_species'] for r in rows} == {
        'Alpha one': 'Alpha two', 'Alpha two': 'Alpha one'}
    assert [c['species'] for c in cross] == ['Alpha one;Alpha two;Beta one']


def test_strip_gtdb_suffixes_matches_the_builder():
    import build_gtdb_db
    for name in ('Pseudomonas_E aeruginosa', 'Listeria monocytogenes_B', 'Escherichia coli',
                 'Lachnospiraceae_NK4A136_group'):
        assert build_gtdb_db.normalize_gtdb_name(name) == identical_16s.strip_gtdb_suffixes(name)


def test_genera_of_other_databases_are_not_cut(tmp_path):
    # Two genera that differ only inside the name stay two genera (an unanchored
    # suffix pattern made both 'Lachnospiraceae4A136_group')
    fasta, tax = write_db(tmp_path, [('1', S1), ('2', S1)], {
        '1': ('Lachnospiraceae_NK4A136_group sp1', 'Lachnospiraceae_NK4A136_group', 'F', 'P'),
        '2': ('Lachnospiraceae_XY4A136_group sp2', 'Lachnospiraceae_XY4A136_group', 'F', 'P'),
    })
    rows, cross = identical_16s.identical_reports(fasta, tax)
    assert rows == []
    assert len(cross) == 1

