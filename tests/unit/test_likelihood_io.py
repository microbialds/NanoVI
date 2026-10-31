"""Unit tests for the SAM-reading functions in bin/likelihood.py.

get_cigar_op_log_probabilities, log_prob_rgs_dict and score_reads (LIKELIHOOD,
which does both in one pass) are run on small SAM and BAM files written with
pysam in tmp_path.
"""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pysam
import pytest

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, BIN)
from likelihood import get_cigar_op_log_probabilities, log_prob_rgs_dict, read_scores_to_arrays, score_reads

REFS = ('1:Species_one:1', '2:Species_two:1', '3:Species_three:1')
PRIMARY, UNMAPPED, SECONDARY, SUPPLEMENTARY = 0, 4, 256, 2048


def write_sam(path, records, mode='w'):
    """Write records (query_name, flag, reference, cigar, NM) as a SAM file ('wb': BAM).

    NM is the edit distance (mismatches + inserted + deleted bases); pysam's
    get_cigar_stats reports it at index 10, from which the code derives mismatches.
    """
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': r, 'LN': 5000} for r in REFS]}
    with pysam.AlignmentFile(str(path), mode, header=header) as out:
        for name, flag, ref, cigar, nm in records:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.flag = flag
            if flag & UNMAPPED:
                a.reference_id = -1
                a.reference_start = -1
                a.query_sequence = 'A' * 50
            else:
                a.reference_id = REFS.index(ref)
                a.reference_start = 0
                a.mapping_quality = 60
                a.cigarstring = cigar
                a.query_sequence = 'A' * a.query_length
                a.set_tag('NM', nm)
            out.write(a)
    return str(path)


def log_shares(counts):
    """Expected log probabilities for counts of I, D, S and X (zeros dropped)."""
    kept = [c for c in counts if c]
    return [math.log(c / sum(kept)) for c in kept]


class TestGetCigarOpLogProbabilities:
    def test_longest_alignment_is_the_maximum_over_all_alignments(self, tmp_path):
        # r1's secondary alignment is its longest. With two threads the chunked
        # version put it in the first chunk and kept the second chunk's 102 (F04).
        sam = write_sam(tmp_path / 'a.sam', [
            ('r1', SECONDARY, REFS[2], '150M', 4),
            ('r2', PRIMARY, REFS[1], '80M1D', 2),
            ('r1', SUPPLEMENTARY, REFS[1], '50M', 1),
            ('r1', PRIMARY, REFS[0], '100M2I', 3),
        ])
        _, _, longest = get_cigar_op_log_probabilities(sam)
        assert longest == {'r1': 150, 'r2': 81}

    def test_error_profile_counts_primary_alignments_only(self, tmp_path):
        # F03: secondary and supplementary alignments (and unmapped reads) are
        # left out of the error profile, as in Emu
        sam = write_sam(tmp_path / 'a.sam', [
            ('r1', PRIMARY, REFS[0], '4S100M2I', 5),       # I 2, S 4, X 3
            ('r1', SECONDARY, REFS[1], '100M3I3D', 20),    # left out
            ('r2', PRIMARY, REFS[1], '100M1D', 2),         # D 1, X 1
            ('r2', SUPPLEMENTARY, REFS[2], '10S50M', 9),   # left out
            ('r3', UNMAPPED, None, None, 0),
        ])
        log_probs, zero_locs, longest = get_cigar_op_log_probabilities(sam)
        assert zero_locs == []
        assert log_probs == pytest.approx(log_shares([2, 1, 4, 4]))
        assert longest == {'r1': 106, 'r2': 101, 'r3': 0}

    def test_operation_seen_only_in_secondary_alignments_is_a_zero_loc(self, tmp_path):
        sam = write_sam(tmp_path / 'a.sam', [
            ('r1', PRIMARY, REFS[0], '100M2I', 3),
            ('r1', SECONDARY, REFS[1], '100M2D', 2),
        ])
        log_probs, zero_locs, _ = get_cigar_op_log_probabilities(sam)
        # D and S never occur in a primary alignment
        assert zero_locs == [1, 2]
        assert log_probs == pytest.approx(log_shares([2, 0, 0, 1]))

    def test_single_alignment(self, tmp_path):
        # F04: the old chunked version divided by zero when there were fewer alignments than threads
        sam = write_sam(tmp_path / 'a.sam', [('r1', PRIMARY, REFS[0], '5S100M2I', 3)])
        log_probs, zero_locs, longest = get_cigar_op_log_probabilities(sam)
        # I = 2, D = 0, S = 5, X = 3 - 2 = 1
        assert zero_locs == [1]
        assert log_probs == pytest.approx(log_shares([2, 0, 5, 1]))
        assert longest == {'r1': 107}

    def test_empty_sam_raises(self, tmp_path):
        sam = write_sam(tmp_path / 'empty.sam', [])
        with pytest.raises(ValueError, match='No primary mapped alignments'):
            get_cigar_op_log_probabilities(sam)

    def test_only_unmapped_reads_raises(self, tmp_path):
        sam = write_sam(tmp_path / 'u.sam', [('r1', UNMAPPED, None, None, 0), ('r2', UNMAPPED, None, None, 0)])
        with pytest.raises(ValueError, match='No primary mapped alignments'):
            get_cigar_op_log_probabilities(sam)

    def test_no_error_operations_gives_an_empty_profile(self, tmp_path):
        # Error-free primary alignments no longer stop LIKELIHOOD (F45): every operation is
        # a zero loc, a perfect alignment scores 0 and one with any operation is skipped
        sam = write_sam(tmp_path / 'p.sam', [('r1', PRIMARY, REFS[0], '100M', 0),
                                            ('r1', SECONDARY, REFS[1], '95M5S', 0),
                                            ('r2', PRIMARY, REFS[1], '100M', 0)])
        log_probs, zero_locs, _ = get_cigar_op_log_probabilities(sam)
        assert log_probs == [] and zero_locs == [0, 1, 2, 3]
        assert score_reads(sam) == ({'r1': ([1], [0.0]), 'r2': ([2], [0.0])}, 0, 2)


class TestLogProbRgsDict:
    LOG_P = [math.log(0.2), math.log(0.3), math.log(0.1), math.log(0.4)]

    def score(self, stats, longest, align_len, log_p=None):
        log_p = self.LOG_P if log_p is None else log_p
        return sum(lp * n for lp, n in zip(log_p, stats)) * longest / align_len

    def test_scores_and_best_alignment_per_species(self, tmp_path):
        sam = write_sam(tmp_path / 'a.sam', [
            ('r1', PRIMARY, REFS[0], '100M2I', 3),       # I 2, X 1
            ('r1', SECONDARY, REFS[1], '100M1D', 4),     # D 1, X 3
            ('r1', SECONDARY, REFS[1], '100M', 6),       # X 6: worse, same species
            ('r2', UNMAPPED, None, None, 0),
        ])
        longest = {'r1': 102, 'r2': 0}
        log_p_rgs, unassigned, assigned = log_prob_rgs_dict(sam, self.LOG_P, longest, [])
        assert assigned == 1 and unassigned == 1
        sids, scores = log_p_rgs['r1']
        assert sids == [1, 2]
        assert scores == pytest.approx([
            self.score([2, 0, 0, 1], 102, 102),
            self.score([0, 1, 0, 3], 102, 101),
        ])

    def test_zero_locs_none_and_empty_agree(self, tmp_path):
        sam = write_sam(tmp_path / 'a.sam', [('r1', PRIMARY, REFS[0], '100M2I', 3)])
        longest = {'r1': 102}
        assert log_prob_rgs_dict(sam, self.LOG_P, longest) == log_prob_rgs_dict(sam, self.LOG_P, longest, [])

    def test_alignment_using_a_zero_count_operation(self, tmp_path):
        # F10: the profile has no deletions (zero_locs [1]). r1's only alignment
        # has a deletion, so r1 is unassigned; r2 keeps its other alignment.
        sam = write_sam(tmp_path / 'a.sam', [
            ('r1', PRIMARY, REFS[0], '100M1D', 2),
            ('r2', PRIMARY, REFS[1], '100M2I', 3),
            ('r2', SECONDARY, REFS[2], '100M1D', 1),
        ])
        log_p = [math.log(0.3), math.log(0.2), math.log(0.5)]  # I, S, X
        longest = {'r1': 101, 'r2': 102}
        log_p_rgs, unassigned, assigned = log_prob_rgs_dict(sam, log_p, longest, [1])
        assert set(log_p_rgs) == {'r2'}
        assert log_p_rgs['r2'][0] == [2]
        assert log_p_rgs['r2'][1] == pytest.approx([self.score([2, 0, 1], 102, 102, log_p)])
        assert assigned == 1 and unassigned == 1


# ─── score_reads and read_likelihood.py: one pass (F24) ─────────────────────

MIXED = [
    ('r1', SECONDARY, REFS[2], '150M', 4),
    ('r1', PRIMARY, REFS[0], '100M2I', 3),
    ('r2', PRIMARY, REFS[1], '80M1D', 2),
    ('r3', UNMAPPED, None, None, 0),
    ('r1', SUPPLEMENTARY, REFS[1], '50M', 1),
    ('r4', PRIMARY, REFS[1], '90M5S', 2),
    ('r4', SECONDARY, REFS[1], '95M', 7),
]


@pytest.mark.parametrize('mode,suffix', [('w', 'sam'), ('wb', 'bam')])
def test_score_reads_equals_the_two_passes(tmp_path, mode, suffix):
    path = write_sam(tmp_path / f'a.{suffix}', MIXED, mode)
    log_probs, zero_locs, longest = get_cigar_op_log_probabilities(path)
    assert score_reads(path) == log_prob_rgs_dict(path, log_probs, longest, zero_locs)


@pytest.mark.parametrize('records,unassigned', [
    ([], 0),                                                              # header only: no read passed the filter
    ([('r1', UNMAPPED, None, None, 0), ('r2', UNMAPPED, None, None, 0)], 2),  # nothing aligned
])
def test_score_reads_without_primary_alignments_leaves_every_read_unassigned(tmp_path, records, unassigned):
    # A blank, or a sample whose reads all failed the length filter, must not stop the run (F45)
    path = write_sam(tmp_path / 'a.bam', records, 'wb')
    assert score_reads(path) == ({}, unassigned, 0)


def test_read_scores_to_arrays():
    lp = {'b': ([3, 1], [-2.0, -5.0]), 'a': ([2], [-1.5])}
    offsets, tax_ids, scores = read_scores_to_arrays(lp)
    assert offsets.tolist() == [0, 2, 3]
    assert tax_ids.tolist() == [3, 1, 2]
    assert scores.tolist() == [-2.0, -5.0, -1.5]


def test_score_reads_lists_the_unassigned_reads(tmp_path):
    # In file order, and only on request: the three-value form stays (F55)
    records = MIXED + [('r5', UNMAPPED, None, None, 0)]
    path = write_sam(tmp_path / 'a.bam', records, 'wb')
    log_p_rgs, unassigned, assigned, names = score_reads(path, unassigned_names=True)
    assert names == ['r3', 'r5'] and unassigned == 2
    assert score_reads(path) == (log_p_rgs, unassigned, assigned)
    blank = write_sam(tmp_path / 'b.bam', [('u1', UNMAPPED, None, None, 0), ('u2', UNMAPPED, None, None, 0)], 'wb')
    assert score_reads(blank, unassigned_names=True) == ({}, 2, 0, ['u1', 'u2'])


def test_load_read_scores_with_and_without_names(tmp_path):
    from read_scores import load_read_scores, save_read_scores
    lp = {'b': ([3, 1], [-2.0, -5.0]), 'a': ([2], [-1.5])}
    offsets, tax_ids, scores = read_scores_to_arrays(lp)
    save_read_scores(tmp_path / 'new.npz', offsets, tax_ids, scores, read_names=list(lp), unassigned_names=['u'])
    # One long name does not widen the others: the names are stored joined, not fixed-width
    long = 'x' * 300
    save_read_scores(tmp_path / 'long.npz', offsets, tax_ids, scores, read_names=['b', long], unassigned_names=[])
    with np.load(tmp_path / 'long.npz') as data:
        assert data['read_names'].nbytes == 302 and data['unassigned_names'].size == 0
    assert load_read_scores(tmp_path / 'long.npz')[3:] == (['b', long], [])
    loaded = load_read_scores(tmp_path / 'new.npz')
    assert loaded[3:] == (['b', 'a'], ['u'])
    assert load_read_scores(tmp_path / 'new.npz', names=False)[3:] == (None, None)
    for got, expected in zip(loaded[:3], (offsets, tax_ids, scores)):
        np.testing.assert_array_equal(got, expected)
    # A file written before the names were kept still loads
    np.savez_compressed(tmp_path / 'old.npz', offsets=offsets, tax_ids=tax_ids, scores=scores)
    assert load_read_scores(tmp_path / 'old.npz')[3:] == (None, None)
    with pytest.raises(ValueError, match='1 read names for 2 reads'):
        save_read_scores(tmp_path / 'bad.npz', offsets, tax_ids, scores, read_names=['b'])


def test_read_likelihood_script(tmp_path):
    path = write_sam(tmp_path / 'a.bam', MIXED, 'wb')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
         '--scores', 'scores.npz', '--counts', 'counts.json'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    log_p_rgs, unassigned, assigned = score_reads(path)
    with np.load(tmp_path / 'scores.npz') as data:
        assert set(data.files) == {'offsets', 'tax_ids', 'scores'}
        expected = read_scores_to_arrays(log_p_rgs)
        for key, value in zip(('offsets', 'tax_ids', 'scores'), expected):
            np.testing.assert_array_equal(data[key], value)
    assert json.loads((tmp_path / 'counts.json').read_text()) == {
        'assigned_count': assigned, 'unassigned_count': unassigned}
    assert (assigned, unassigned) == (3, 1)


def test_read_likelihood_script_stores_read_names_on_request(tmp_path):
    # The name of each CSR row, and the unassigned reads, in file order (F55)
    from read_scores import load_read_scores
    path = write_sam(tmp_path / 'a.bam', MIXED, 'wb')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
         '--scores', 'scores.npz', '--counts', 'counts.json', '--read-names'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    with np.load(tmp_path / 'scores.npz') as data:
        assert set(data.files) == {'offsets', 'tax_ids', 'scores', 'read_names', 'unassigned_names'}
        assert data['read_names'].tobytes() == b'r1\nr2\nr4'
    log_p_rgs = score_reads(path)[0]
    loaded = load_read_scores(tmp_path / 'scores.npz')
    for got, expected in zip(loaded[:3], read_scores_to_arrays(log_p_rgs)):
        np.testing.assert_array_equal(got, expected)
    assert loaded[3:] == (list(log_p_rgs), ['r3']) == (['r1', 'r2', 'r4'], ['r3'])


def test_read_likelihood_script_without_primary_alignments(tmp_path):
    # An empty read-score npz and the read counts, with a warning, instead of a failure (F45)
    path = write_sam(tmp_path / 'S1_vi_alignments.bam', [('r1', UNMAPPED, None, None, 0)], 'wb')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
         '--scores', 'scores.npz', '--counts', 'counts.json'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert 'S1_vi_alignments.bam' in result.stderr and 'no primary alignment' in result.stderr
    with np.load(tmp_path / 'scores.npz') as data:
        assert data['offsets'].tolist() == [0]
        assert data['tax_ids'].size == 0 and data['scores'].size == 0
    assert json.loads((tmp_path / 'counts.json').read_text()) == {'assigned_count': 0, 'unassigned_count': 1}


# ─── Options removed after the Phase 7a benchmark ──────────────────────────

def test_error_profile_option_is_removed(tmp_path):
    # The error profile comes from primary alignments only (D3); 'all' lost in Phase 7a (D31)
    path = write_sam(tmp_path / 'a.bam', MIXED, 'wb')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
         '--scores', 'scores.npz', '--counts', 'counts.json', '--error-profile', 'all'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert 'unrecognized arguments: --error-profile' in result.stderr


def test_per_base_score_is_removed(tmp_path):
    # Emu's score stays; the per-base score lost in the Phase 7a benchmark (D29)
    path = write_sam(tmp_path / 'a.bam', MIXED, 'wb')
    result = subprocess.run(
        [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
         '--scores', 'scores.npz', '--counts', 'counts.json', '--score', 'per-base'],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert 'unrecognized arguments: --score' in result.stderr


def test_identity_and_length_filters_are_removed(tmp_path):
    # Identity and aligned-length filters lost in the Phase 7a benchmark (D25)
    path = write_sam(tmp_path / 'a.bam', MIXED, 'wb')
    argv = [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
            '--scores', 'scores.npz', '--counts', 'counts.json']
    for option, value in (('--min-identity', '95'), ('--min-aligned-length', '1000')):
        result = subprocess.run(argv + [option, value], cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 2
        assert f'unrecognized arguments: {option}' in result.stderr


def test_read_likelihood_script_has_no_unknown_component(tmp_path):
    # The unknown-taxon prototype was dropped after the Phase 7a benchmark (D22)
    path = write_sam(tmp_path / 'a.bam', MIXED, 'wb')
    argv = [sys.executable, os.path.join(BIN, 'read_likelihood.py'), '--alignments', path,
            '--scores', 'scores.npz', '--counts', 'counts.json']
    result = subprocess.run(argv + ['--unknown-delta', '0.05'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 2
    assert 'unrecognized arguments: --unknown-delta' in result.stderr
    result = subprocess.run(argv, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with np.load(tmp_path / 'scores.npz') as data:
        assert set(data.files) == {'offsets', 'tax_ids', 'scores'}
