"""Unit tests for bin/identical_16s.py: species that share an identical 16S (D23, D24)."""
import csv
import gzip
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
                 'Lachnospiraceae_NK4A136_group', 'Clostridium_AQ', 'Family_XIII', 'Ruminococcaceae_UCG'):
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


def test_fasta_records_is_the_builders_reader():
    import build_gtdb_db
    assert build_gtdb_db.fasta_records is identical_16s.fasta_records


# ─── Cross-genus rules (F56) ──────────────────────────────────────────────────

X, Y, W = 'ACGTT' * 240, 'GATTACA' * 172, 'CCGTA' * 240
OWN1, OWN2, OWN3 = 'TTAGC' * 240, 'AGCTA' * 240, 'CGATG' * 240


@pytest.fixture
def cross_db(tmp_path):
    """X: Alpha one, Alpha two, Beta one. Y: Alpha one, Gamma one. W: Delta one, Epsilon one."""
    lineages = {
        '1': ('Alpha one', 'Alpha', 'F1', 'P1'),
        '2': ('Alpha two', 'Alpha_A', 'F1', 'P1'),
        '3': ('Beta one', 'Beta', 'F2', 'P2'),
        '5': ('Gamma one', 'Gamma', 'F3', 'P2'),
        '6': ('Delta one', 'Delta', 'F4', 'P3'),
        '7': ('Epsilon one', 'Epsilon', 'F5', 'P3'),
    }
    records = [
        ('1', X), ('2', X), ('3', X),   # 1:db:1, 2:db:2, 3:db:3
        ('1', Y), ('5', Y),             # 1:db:4, 5:db:5
        ('1', OWN1), ('2', OWN2), ('3', OWN3),
        ('6', W), ('7', W),             # 6:db:9, 7:db:10
    ]
    fasta, tax = write_db(tmp_path, records, lineages)
    return fasta, tax, tmp_path


def build_reports(cross_db, rule, copies=None):
    fasta, tax, out = cross_db
    identical_16s.write_reports(fasta, tax, str(out / 'identical_16s.tsv'), str(out / 'cross_genus_identical.tsv'),
                                cross_genus_rule=rule, copies=copies)
    removed = {r['sequence_id']: r['kept_in'] for r in read_tsv(out / 'cross_genus_removed.tsv')}
    emptied = [r['species'] for r in read_tsv(out / 'cross_genus_emptied.tsv')]
    ids = [h[1:] for h, _ in identical_16s.fasta_records(fasta)]
    return removed, emptied, ids


def test_remove_drops_every_copy(cross_db):
    removed, emptied, ids = build_reports(cross_db, 'remove')
    assert removed == dict.fromkeys(['1:db:1', '2:db:2', '3:db:3', '1:db:4', '5:db:5', '6:db:9', '7:db:10'], '')
    assert ids == ['1:db:6', '2:db:7', '3:db:8']
    assert emptied == ['Delta one', 'Epsilon one', 'Gamma one']


def test_best_genus_keeps_the_genus_with_most_entries_and_drops_ties(cross_db):
    removed, emptied, ids = build_reports(cross_db, 'best-genus')
    # X: two Alpha entries (Alpha_A is Alpha) against one Beta; Y and W tie, so every copy goes
    assert removed == {'3:db:3': 'Alpha', '1:db:4': '', '5:db:5': '', '6:db:9': '', '7:db:10': ''}
    assert ids == ['1:db:1', '2:db:2', '1:db:6', '2:db:7', '3:db:8']
    assert emptied == ['Delta one', 'Epsilon one', 'Gamma one']


def test_most_copies_counts_input_copies_then_genus_entries(cross_db):
    # Beta one's X stands for ten input sequences; every other entry for one
    removed, emptied, ids = build_reports(cross_db, 'most-copies', copies={'3:db:3': 10})
    # X goes to Beta; Y ties on copies and goes to Alpha (5 entries against Gamma's 1);
    # W ties on both, so both copies go
    assert removed == {'1:db:1': 'Beta', '2:db:2': 'Beta', '5:db:5': 'Alpha', '6:db:9': '', '7:db:10': ''}
    assert ids == ['3:db:3', '1:db:4', '1:db:6', '2:db:7', '3:db:8']
    assert emptied == ['Delta one', 'Epsilon one', 'Gamma one']


def test_reports_describe_the_cleaned_database(cross_db):
    fasta, tax, out = cross_db
    build_reports(cross_db, 'best-genus')
    # X still links Alpha one and Alpha two; cross_genus_identical.tsv lists the sequences as built
    assert [r['species'] for r in read_tsv(out / 'identical_16s.tsv')] == ['Alpha one', 'Alpha two']
    assert len(read_tsv(out / 'cross_genus_identical.tsv')) == 3
    removed_rows = read_tsv(out / 'cross_genus_removed.tsv')
    assert list(removed_rows[0]) == identical_16s.REMOVED_COLUMNS
    assert removed_rows[0] == {'sequence_id': '3:db:3', 'tax_id': '3', 'species': 'Beta one', 'genus': 'Beta',
                               'md5': removed_rows[0]['md5'], 'rule': 'best-genus', 'kept_in': 'Alpha'}
    emptied_rows = read_tsv(out / 'cross_genus_emptied.tsv')
    assert list(emptied_rows[0]) == identical_16s.EMPTIED_COLUMNS
    assert emptied_rows[0] == {'tax_id': '6', 'species': 'Delta one', 'genus': 'Delta', 'sequences_removed': '1'}


def test_removing_a_shared_sequence_removes_the_group(cross_db):
    fasta, tax, out = cross_db
    build_reports(cross_db, 'most-copies', copies={'3:db:3': 10})
    assert read_tsv(out / 'identical_16s.tsv') == []


def test_keep_leaves_the_database_as_built(cross_db):
    fasta, tax, out = cross_db
    before = open(fasta).read()
    identical_16s.write_reports(fasta, tax, str(out / 'identical_16s.tsv'), str(out / 'cross.tsv'))
    assert open(fasta).read() == before
    assert not (out / 'cross_genus_removed.tsv').exists()
    assert not (out / 'cross_genus_emptied.tsv').exists()


def test_most_copies_without_genus_entries_drops_ties(cross_db):
    fasta, tax, _ = cross_db
    _, cross = identical_16s.identical_reports(fasta, tax)
    removed = identical_16s.cross_genus_removals(cross, identical_16s.read_taxonomy(tax), 'most-copies',
                                                 {'3:db:3': 10})
    # Without genus_entries a tie on copies cannot be broken: Y and W lose every copy
    assert removed == {'1:db:1': 'Beta', '2:db:2': 'Beta', '1:db:4': '', '5:db:5': '', '6:db:9': '', '7:db:10': ''}


def test_gzipped_database_stays_gzipped(cross_db):
    fasta, tax, out = cross_db
    gz = fasta + '.gz'
    with open(fasta, 'rb') as src, gzip.open(gz, 'wb') as dst:
        dst.write(src.read())
    identical_16s.write_reports(gz, tax, str(out / 'identical_16s.tsv'), cross_genus_rule='remove')
    with gzip.open(gz, 'rt') as fh:
        ids = [line[1:].strip() for line in fh if line.startswith('>')]
    assert ids == ['1:db:6', '2:db:7', '3:db:8']


def test_unknown_rule_is_refused():
    with pytest.raises(ValueError, match='unknown cross-genus rule'):
        identical_16s.cross_genus_removals([], {}, 'majority')
