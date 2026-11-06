"""A small NanoVI results directory for the run report's tests (F60): three samples (S1, S2, S10, so natural order
matters), species tables as WRITE_OUTPUT writes them, the consolidated tables from the real
generate_phyloseq_tables.py, fastplong JSON reports, software_versions.yml and, on request, the group, cell-level
and per-read tables."""
import gzip
import json
import os
import subprocess
import sys

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
RANKS = ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
LINEAGE = {
    'Alpha one': ['Alpha', 'Fam1', 'Ord1', 'Cls1', 'Phy1', 'Bacteria'],
    'Alpha two': ['Alpha', 'Fam1', 'Ord1', 'Cls1', 'Phy1', 'Bacteria'],
    'Beta three': ['Beta', 'Fam2', 'Ord1', 'Cls1', 'Phy1', 'Bacteria'],
    'Gamma five': ['Gamma', 'Fam3', 'Ord2', 'Cls2', 'Phy2', 'Bacteria'],
}
# sample -> [(species, abundance, estimated_counts, identical_16s_species, shared_read_fraction)]
SPECIES = {
    'S1': [('Alpha one', 0.6, 60, '', 0.1), ('Alpha two', 0.3, 30, '', 0.7), ('Beta three', 0.1, 10, 'Beta four', 0.0)],
    'S10': [('Alpha one', 0.5, 50, '', 0.0), ('Gamma five', 0.5, 50, '', 0.2)],
    'S2': [],
}
PLACEHOLDER_COUNTS = {'S1': (5, 5), 'S10': (0, 10), 'S2': (0, 20)}  # Below detection, Unassigned
FASTPLONG = {'S1': (120, 110, 8, 2), 'S10': (130, 110, 15, 5), 'S2': (25, 20, 5, 0)}  # raw, passed, short, long
VERSIONS = """ALIGN_READS:
  minimap2: 2.31-r1302
  samtools: '1.24'
FILTER_READS:
  fastplong: 0.7.0
RUN_VI:
  python: 3.11.17
  numpy: 2.4.6
  scipy: 1.17.1
Workflow:
    microbialds/NanoVI: 1.1.1
    Nextflow: 26.04.4
"""


def species_table(sample, counts=True):
    header = ['species', 'abundance'] + RANKS + (['estimated_counts'] if counts else []) + \
        ['identical_16s_species', 'shared_read_fraction']
    lines = ['\t'.join(header)]
    for sp, ab, n, ident, shared in SPECIES[sample]:
        lines.append('\t'.join([sp, repr(ab)] + LINEAGE[sp] + ([str(n)] if counts else []) + [ident, repr(shared)]))
    below, unassigned = PLACEHOLDER_COUNTS[sample]
    for label, n in (('Below detection', below), ('Unassigned', unassigned)):
        lines.append('\t'.join([label, '0.0'] + [label] * 6 + ([str(n)] if counts else []) + ['', '']))
    return '\n'.join(lines) + '\n'


def make_results(root, counts=True, groups=False, cell=False, reads=False, versions=True):
    so = root / 'sample_outputs'
    so.mkdir(parents=True)
    for s in SPECIES:
        (so / f'{s}_rel-abundance.tsv').write_text(species_table(s, counts))
    co = root / 'consolidated_output'
    subprocess.run([sys.executable, os.path.join(BIN, 'generate_phyloseq_tables.py'),
                    *sorted(str(p) for p in so.glob('*_rel-abundance.tsv')), '--output-dir', str(co)],
                   check=True, capture_output=True)
    qc = root / 'qc' / 'fastplong'
    qc.mkdir(parents=True)
    for s, (raw, passed, short, long_) in FASTPLONG.items():
        report = {'summary': {'fastplong_version': '0.7.0', 'before_filtering': {'total_reads': raw},
                              'after_filtering': {'total_reads': passed}},
                  'filtering_result': {'passed_filter_reads': passed, 'low_quality_reads': 0, 'too_many_N_reads': 0,
                                       'too_short_reads': short, 'too_long_reads': long_, 'chimeric_reads': 0}}
        with gzip.open(qc / f'{s}_report.json.gz', 'wt') as fh:
            json.dump(report, fh)
    if versions:
        (root / 'pipeline_info').mkdir()
        (root / 'pipeline_info' / 'software_versions.yml').write_text(VERSIONS)
    if groups:
        header = ['group', 'members', 'n_members', 'abundance', 'lowest_common_rank', 'lowest_common_taxon',
                  'estimated_counts', 'identical_16s_relatives', 'max_shared', 'shared_read_fraction']
        extra = ['cell_abundance', 'copy_number', 'copy_number_level'] if cell else []
        rows = {'S1': [['G1', 'Alpha one; Alpha two', '2', '0.9', 'genus', 'Alpha', '90', '', '0.8', '0.05'],
                       ['G2', 'Beta three', '1', '0.1', 'species', 'Beta three', '10', 'Beta four', '', '0.0']],
                'S10': [['G1', 'Alpha one', '1', '0.5', 'species', 'Alpha one', '50', '', '', '0.0'],
                        ['G2', 'Gamma five', '1', '0.5', 'species', 'Gamma five', '50', '', '', '0.2']],
                'S2': []}
        cells = {'Alpha one; Alpha two': ['0.8', '6.0', 'species'], 'Beta three': ['0.2', '2.0', 'genus'],
                 'Alpha one': ['0.4', '6.0', 'species'], 'Gamma five': ['0.6', '4.0', 'family']}
        for s, rs in rows.items():
            lines = ['\t'.join(header + extra)]
            for r in rs:
                lines.append('\t'.join(r + (cells[r[1]] if cell else [])))
            for label in ('Below detection', 'Unassigned'):
                lines.append('\t'.join([label, label] + [''] * (len(header) + len(extra) - 2)))
            (so / f'{s}_groups.tsv').write_text('\n'.join(lines) + '\n')
    if cell:
        header = ['species', 'abundance', 'read_abundance', 'copy_number', 'copy_number_level'] + RANKS
        copies = {'Alpha one': (6.0, 'species'), 'Alpha two': (6.0, 'species'), 'Beta three': (2.0, 'genus'),
                  'Gamma five': (4.0, 'family')}
        for s, rs in SPECIES.items():
            cell_share = {sp: ab / copies[sp][0] for sp, ab, *_ in rs}
            total = sum(cell_share.values())
            below, unassigned = PLACEHOLDER_COUNTS[s]
            lines = ['# Cell-level abundance: read-level abundance / 16S copy number, renormalised.',
                     f'# Not included: Below detection ({below} reads) and Unassigned ({unassigned} reads).',
                     '\t'.join(header)]
            for sp, ab, *_ in rs:
                lines.append('\t'.join([sp, repr(cell_share[sp] / total), repr(ab), repr(copies[sp][0]),
                                        copies[sp][1]] + LINEAGE[sp]))
            (so / f'{s}_cell-abundance.tsv').write_text('\n'.join(lines) + '\n')
    if reads:
        assignments = {
            'S1': [('r1', 'Alpha one', 1.0), ('r2', 'Alpha one', 0.6), ('r2', 'Alpha two', 0.4),
                   ('r3', 'Alpha two', 0.3), ('r3', 'Alpha one', 0.3), ('r3', 'Beta three', 0.4),
                   ('r4', 'Below detection', 1.0), ('r5', 'Unassigned', 1.0)],
            'S10': [('q1', 'Gamma five', 0.995), ('q2', 'Alpha one', 1.0)],
            'S2': [('u1', 'Unassigned', 1.0)],
        }
        for s, rows in assignments.items():
            text = 'read_id\tspecies\tprobability\n' + ''.join(f'{r}\t{sp}\t{p!r}\n' for r, sp, p in rows)
            with gzip.open(so / f'{s}_read_assignments.tsv.gz', 'wt') as fh:
                fh.write(text)
    return root
