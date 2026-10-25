# Portions of this file are adapted from Emu (https://github.com/treangenlab/emu),
# Copyright (c) 2021 Kristen Curry, Qi Wang, Michael G. Nute, Alona Tyshaieva,
# Qinglong Wu, Tor Savidge, Alexander Dilthey, Todd J. Treangen. MIT License;
# see THIRD_PARTY_LICENSES for the full text.
# Adapted here: get_align_stats, get_align_len, cigar_profile
# (get_cigar_op_log_probabilities) and read_scores_from_records (log_prob_rgs_dict):
# Emu's CIGAR-based read score.

"""
Read likelihoods for NanoVI.

1. CIGAR operation log-probabilities estimated from the sample's primary alignments
2. Per-read log-likelihood of each candidate species

score_reads does both from one pass over the alignments (LIKELIHOOD, F24); the
estimators (em_estimator.py, variational_inference_bayesian.py) take the
per-read scores it returns.
"""

import math
from array import array

import numpy as np
import pysam
from operator import mul
from sys import stdout

CIGAR_OPS = [1, 2, 4, 10]
CIGAR_OPS_ALL = [0, 1, 2, 4]

def get_align_stats(alignment, cigar_stats=None):
    if cigar_stats is None:
        cigar_stats = alignment.get_cigar_stats()[0]
    n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]
    return [cigar_stats[1], cigar_stats[2], cigar_stats[4], n_mismatch]

def get_align_len(alignment, cigar_stats=None):
    if cigar_stats is None:
        cigar_stats = alignment.get_cigar_stats()[0]
    return sum(cigar_stats[cigar_op] for cigar_op in CIGAR_OPS_ALL)

class AlignmentRecords:
    """What the read score needs from each alignment, read in one pass (F24).

    Per alignment, in file order: the read (an index into names), the
    reference (an index into ref_names, -1 when unmapped), the aligned length
    (M + I + D + S), the counts of I, D, S and mismatches, whether it is a
    primary mapped alignment, and the counts of aligned columns (M, = and X
    operations) and hard-clipped bases, which the experimental scores and
    filters use. Integers are kept in array.array, about 80 bytes per
    alignment, so the alignments are not held as pysam objects.
    """

    def __init__(self, sam_path):
        self.sam_path = sam_path
        self.names = []
        self.read = array("q")
        self.ref = array("q")
        self.length = array("q")
        self.stats = array("q")
        self.primary = array("b")
        self.extra = array("q")
        read_index = {}
        with pysam.AlignmentFile(sam_path) as sam_pysam:
            self.ref_names = list(sam_pysam.references)
            for alignment in sam_pysam:
                query_name = alignment.query_name
                r = read_index.get(query_name)
                if r is None:
                    r = read_index[query_name] = len(self.names)
                    self.names.append(query_name)
                self.read.append(r)
                self.ref.append(alignment.reference_id if alignment.reference_name else -1)
                # The CIGAR is parsed once per alignment
                cigar_stats = alignment.get_cigar_stats()[0]
                self.length.append(get_align_len(alignment, cigar_stats))
                self.stats.extend(get_align_stats(alignment, cigar_stats))
                self.primary.append(
                    not alignment.is_unmapped and not alignment.is_secondary and not alignment.is_supplementary
                )
                self.extra.extend((cigar_stats[0] + cigar_stats[7] + cigar_stats[8], cigar_stats[5]))

    def __len__(self):
        return len(self.read)

    def align_stats(self, i):
        return list(self.stats[4 * i:4 * i + 4])

    def columns(self):
        """Per alignment: aligned columns (M, = and X) and hard-clipped bases, as an (n, 2) array."""
        return np.frombuffer(self.extra, dtype=np.int64).reshape(-1, 2)

    def profile_mask(self, error_profile="primary"):
        """The alignments an error profile is estimated from: primary mapped ones, or
        with error_profile 'all' (experimental, Phase 7a) every mapped alignment."""
        if error_profile == "all":
            return np.frombuffer(self.ref, dtype=np.int64) >= 0
        if error_profile != "primary":
            raise ValueError(f"error_profile must be 'primary' or 'all', not {error_profile!r}")
        return np.frombuffer(self.primary, dtype=np.int8).astype(bool)

    def op_counts(self):
        """Per alignment: matches, I, D, S (hard clips included) and mismatches, as an (n, 5) array.

        Matches are aligned columns minus mismatches (NM - I - D). Hard-clipped
        bases count as soft clips, so every alignment of a read covers the whole read.
        """
        stats = np.frombuffer(self.stats, dtype=np.int64).reshape(-1, len(CIGAR_OPS))
        cols = self.columns()
        return np.column_stack([cols[:, 0] - stats[:, 3], stats[:, 0], stats[:, 1],
                                stats[:, 2] + cols[:, 1], stats[:, 3]])

    def identity(self):
        """Per alignment: Emu's percent identity, 100 x (1 - NM / aligned query bases).

        Aligned query bases are M, =, X and I columns (clips excluded); 0 when there are none.
        """
        stats = np.frombuffer(self.stats, dtype=np.int64).reshape(-1, len(CIGAR_OPS))
        aligned_query = self.columns()[:, 0] + stats[:, 0]
        nm = stats[:, 0] + stats[:, 1] + stats[:, 3]
        return np.where(aligned_query > 0, 100.0 * (1 - nm / np.maximum(aligned_query, 1)), 0.0)

    def aligned_query_length(self):
        """Per alignment: aligned query bases (M, =, X and I; clips excluded), as Emu's --min-align-len."""
        stats = np.frombuffer(self.stats, dtype=np.int64).reshape(-1, len(CIGAR_OPS))
        return self.columns()[:, 0] + stats[:, 0]


def read_alignment_records(sam_path):
    """Read a SAM or BAM file once into AlignmentRecords."""
    return AlignmentRecords(sam_path)


def cigar_profile(records, error_profile="primary"):
    """Error profile and longest alignments from AlignmentRecords.

    Returns what get_cigar_op_log_probabilities returns. error_profile 'all'
    (experimental, Phase 7a) counts every mapped alignment, secondary and
    supplementary included, as NanoVI did before F03.
    """
    longest = np.full(len(records.names), -1, dtype=np.int64)
    np.maximum.at(longest, np.frombuffer(records.read, dtype=np.int64),
                  np.frombuffer(records.length, dtype=np.int64))
    # The error profile comes from primary mapped alignments only, as in Emu (F03, D3)
    n_mapped = int(records.profile_mask("primary").sum())
    primary = records.profile_mask(error_profile)
    stats = np.frombuffer(records.stats, dtype=np.int64).reshape(-1, len(CIGAR_OPS))
    cigar_stats_primary = [int(x) for x in stats[primary].sum(axis=0)]

    if n_mapped == 0:
        raise ValueError(
            f"No primary mapped alignments in {records.sam_path}. Check that the reads overlap with the "
            "reference database and that filtering (--min_length, --max_length) left reads to align."
        )

    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]
    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]

    n_char = sum(cigar_stats_primary)
    if n_char == 0:
        raise ValueError(
            f"No insertions, deletions, soft clips or mismatches in the primary alignments of {records.sam_path}; "
            "the error profile cannot be estimated."
        )
    dict_longest_align = dict(zip(records.names, longest.tolist()))
    return [math.log(x) for x in np.array(cigar_stats_primary) / n_char], zero_locs, dict_longest_align


def get_cigar_op_log_probabilities(sam_path):
    """Compute global log-probabilities for each CIGAR operation from the alignments.

    Reads the SAM or BAM once into AlignmentRecords (memory grows with the number
    of alignments, about 80 bytes each). Returns (log_probs, zero_locs, dict_longest_align):
    the log share of each error operation (I, D, S, X) over the primary mapped
    alignments, with the operations of zero count removed (F03); the indices of
    those operations in CIGAR_OPS; and each read's longest alignment over all of
    its alignments, secondary and supplementary included (F04).

    Raises ValueError when the SAM holds no primary mapped alignment or when those
    alignments have no error operation.
    """
    return cigar_profile(read_alignment_records(sam_path))

def per_base_log_probabilities(records, error_profile="primary", pseudocount=0.5):
    """Per-base log-probabilities of match, I, D, S and X (experimental, M4).

    Each operation's share of all aligned and clipped bases over the primary
    mapped alignments (every mapped alignment with error_profile 'all'), with a
    pseudocount so that an operation unseen there still has a finite cost.
    """
    mask = records.profile_mask(error_profile)
    counts = records.op_counts()[mask].sum(axis=0) + pseudocount
    return np.log(counts / counts.sum())


def read_scores_from_records(records, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None,
                             per_base_log_p=None, keep=None):
    """log_prob_rgs_dict on AlignmentRecords: the same scores, in the same order.

    Experimental (Phase 7a): with per_base_log_p (from per_base_log_probabilities)
    the score is Σ count_op × ln p_op over match, I, D, S and X, with no length
    scaling (clips count as S, so every alignment covers the whole read); keep
    is a boolean array over the alignments, and an alignment outside it is
    skipped like one that uses a zero-count operation.
    """
    zero_locs = sorted(p_cigar_op_zero_locs or [], reverse=True)
    log_p_rgs, unassigned_set = {}, set()
    species_tids = {}
    if per_base_log_p is not None:
        per_base_scores = records.op_counts() @ per_base_log_p

    for i in range(len(records)):
        query_name = records.names[records.read[i]]
        align_len = records.length[i]
        ref = records.ref[i]
        if not (ref >= 0 and align_len):
            unassigned_set.add(query_name)
            continue
        if keep is not None and not keep[i]:
            unassigned_set.add(query_name)
            continue
        if per_base_log_p is not None:
            log_score = float(per_base_scores[i])
        else:
            cigar_stats = records.align_stats(i)
            if any(cigar_stats[j] for j in zero_locs):
                unassigned_set.add(query_name)
                continue
            for j in zero_locs:
                del cigar_stats[j]
            log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats))) * (dict_longest_align[query_name] / align_len)
        species_tid = species_tids.get(ref)
        if species_tid is None:
            species_tid = species_tids[ref] = int(records.ref_names[ref].split(":")[0])
        if query_name not in log_p_rgs:
            log_p_rgs[query_name] = ([species_tid], [log_score])
        elif species_tid not in log_p_rgs[query_name][0]:
            log_p_rgs[query_name][0].append(species_tid)
            log_p_rgs[query_name][1].append(log_score)
        else:
            logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
            if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                log_p_rgs[query_name][1][logprgs_idx] = log_score

    assigned_reads = set(log_p_rgs.keys())
    unassigned_reads = unassigned_set - assigned_reads
    unassigned_count = len(unassigned_reads)
    stdout.write(f"Unassigned read count: {unassigned_count}\n")

    return log_p_rgs, unassigned_count, len(assigned_reads)

def log_prob_rgs_dict(sam_path, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None):
    """Build dict mapping read_name -> (species_ids, log_scores) from SAM alignments.

    Each species keeps the best score among the read's alignments to it. An
    alignment that uses an operation with no count in the error profile
    (p_cigar_op_zero_locs) has no score and is skipped; its read counts as
    unassigned unless another alignment of the read is scored (F10).
    """
    return read_scores_from_records(
        read_alignment_records(sam_path), log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs
    )

def score_reads(sam_path, error_profile="primary", score="emu", min_identity=None, min_aligned_length=None,
                unknown_delta=None, records=None):
    """Error profile and read scores from one pass over a SAM or BAM file (F24).

    Returns (log_p_rgs, unassigned_count, assigned_count), as log_prob_rgs_dict
    does after get_cigar_op_log_probabilities, which read the file twice.

    The other arguments are experimental (Phase 7a benchmark; defaults give the
    pipeline's scores): error_profile 'all' estimates the error mix from every
    mapped alignment; score 'per-base' uses per_base_log_probabilities;
    min_identity (percent) and min_aligned_length skip alignments below them,
    as Emu's --min-pid and --min-align-len; unknown_delta adds a fourth value,
    the per-read floor scores of the unknown-taxon component (unknown_floor_scores),
    in the order of log_p_rgs.
    """
    if records is None:
        records = read_alignment_records(sam_path)
    log_probs, zero_locs, dict_longest_align = cigar_profile(records, error_profile)
    per_base = per_base_log_probabilities(records, error_profile) if score == "per-base" else None
    if score not in ("emu", "per-base"):
        raise ValueError(f"score must be 'emu' or 'per-base', not {score!r}")
    keep = None
    if min_identity is not None or min_aligned_length is not None:
        keep = np.ones(len(records), dtype=bool)
        if min_identity is not None:
            keep &= records.identity() >= min_identity
        if min_aligned_length is not None:
            keep &= records.aligned_query_length() >= min_aligned_length
    result = read_scores_from_records(records, log_probs, dict_longest_align, zero_locs,
                                      per_base_log_p=per_base, keep=keep)
    if unknown_delta is None:
        return result
    floors = unknown_floor_scores(records, result[0], dict_longest_align, log_probs,
                                  per_base, unknown_delta, error_profile)
    return (*result, floors)


def _median_rate(errors, covered, delta):
    ok = covered > 0
    return min(float(np.median(errors[ok] / covered[ok])) + delta, 0.99)


def unknown_floor_scores(records, log_p_rgs, dict_longest_align, log_p_cigar_op,
                         per_base_log_p, delta, error_profile="primary"):
    """Per-read score of the unknown-taxon component (experimental, M5; docs/design/unknown_component.md).

    The score of a virtual alignment of the read, as long as its longest
    alignment, whose error rate is the sample's median per-read error rate plus
    delta, with the sample's mix of error operations. A read from a species
    missing from the database aligns to its nearest relative with roughly its
    own errors plus the divergence between the two, so its best score falls
    below this floor once that divergence exceeds delta.
    """
    mask = records.profile_mask(error_profile)
    if per_base_log_p is not None:
        ops = records.op_counts()[mask]
        errors, covered = ops[:, 1:].sum(axis=1), ops.sum(axis=1)
        rate = _median_rate(errors, covered, delta)
        mix = ops[:, 1:].sum(axis=0) / ops[:, 1:].sum()
        per_base_cost = (1 - rate) * per_base_log_p[0] + rate * float(mix @ per_base_log_p[1:])
    else:
        # Emu's score charges error operations only (I, D, S, X over M + I + D + S
        # columns), each at the log of its share, so the virtual alignment costs
        # rate x Σ p_op ln p_op per column
        stats = np.frombuffer(records.stats, dtype=np.int64).reshape(-1, len(CIGAR_OPS))[mask]
        lengths = np.frombuffer(records.length, dtype=np.int64)[mask]
        rate = _median_rate(stats.sum(axis=1), lengths, delta)
        log_p = np.array(log_p_cigar_op)
        per_base_cost = rate * float(np.exp(log_p) @ log_p)
    return np.array([dict_longest_align[r] * per_base_cost for r in log_p_rgs], dtype=np.float64)


def read_scores_to_arrays(log_p_rgs):
    """Per-read lists as flat arrays, in read order (the read-scores npz, F24).

    Returns (offsets, tax_ids, scores): read r's candidates are
    tax_ids[offsets[r]:offsets[r + 1]], with their scores.
    """
    lengths = [len(sids) for sids, _ in log_p_rgs.values()]
    offsets = np.zeros(len(lengths) + 1, dtype=np.int64)
    np.cumsum(lengths, out=offsets[1:])
    tax_ids = np.fromiter((sid for sids, _ in log_p_rgs.values() for sid in sids),
                          dtype=np.int64, count=int(offsets[-1]))
    scores = np.fromiter((score for _, scores in log_p_rgs.values() for score in scores),
                         dtype=np.float64, count=int(offsets[-1]))
    return offsets, tax_ids, scores
