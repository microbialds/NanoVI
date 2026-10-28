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

# The read-scores layout, in a module without pysam that the estimators share (F52)
from read_scores import read_scores_to_arrays

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
    (M + I + D + S), the counts of I, D, S and mismatches, and whether it is a
    primary mapped alignment. Integers are kept in array.array, about 60 bytes
    per alignment, so the alignments are not held as pysam objects.
    """

    def __init__(self, sam_path):
        self.sam_path = sam_path
        self.names = []
        self.read = array("q")
        self.ref = array("q")
        self.length = array("q")
        self.stats = array("q")
        self.primary = array("b")
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

    def __len__(self):
        return len(self.read)

    def align_stats(self, i):
        return list(self.stats[4 * i:4 * i + 4])


def read_alignment_records(sam_path):
    """Read a SAM or BAM file once into AlignmentRecords."""
    return AlignmentRecords(sam_path)


def cigar_profile(records):
    """Error profile and longest alignments from AlignmentRecords.

    Returns what get_cigar_op_log_probabilities returns.
    """
    longest = np.full(len(records.names), -1, dtype=np.int64)
    np.maximum.at(longest, np.frombuffer(records.read, dtype=np.int64),
                  np.frombuffer(records.length, dtype=np.int64))
    # The error profile comes from primary mapped alignments only, as in Emu (F03, D3)
    primary = np.frombuffer(records.primary, dtype=np.int8).astype(bool)
    n_mapped = int(primary.sum())
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
    of alignments, about 60 bytes each). Returns (log_probs, zero_locs, dict_longest_align):
    the log share of each error operation (I, D, S, X) over the primary mapped
    alignments, with the operations of zero count removed (F03); the indices of
    those operations in CIGAR_OPS; and each read's longest alignment over all of
    its alignments, secondary and supplementary included (F04).

    Raises ValueError when the SAM holds no primary mapped alignment or when those
    alignments have no error operation.
    """
    return cigar_profile(read_alignment_records(sam_path))

def read_scores_from_records(records, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None):
    """log_prob_rgs_dict on AlignmentRecords: the same scores, in the same order."""
    zero_locs = sorted(p_cigar_op_zero_locs or [], reverse=True)
    log_p_rgs, unassigned_set = {}, set()
    species_tids = {}

    for i in range(len(records)):
        query_name = records.names[records.read[i]]
        align_len = records.length[i]
        ref = records.ref[i]
        if not (ref >= 0 and align_len):
            unassigned_set.add(query_name)
            continue
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

def score_reads(sam_path):
    """Error profile and read scores from one pass over a SAM or BAM file (F24).

    Returns (log_p_rgs, unassigned_count, assigned_count), as log_prob_rgs_dict
    does after get_cigar_op_log_probabilities, which read the file twice. With no
    primary mapped alignment (a blank, or a sample whose reads all failed the
    length filter) no read is scored and every read is unassigned (F45).
    """
    records = read_alignment_records(sam_path)
    if not any(records.primary):
        stdout.write(f"Unassigned read count: {len(records.names)}\n")
        return {}, len(records.names), 0
    log_probs, zero_locs, dict_longest_align = cigar_profile(records)
    return read_scores_from_records(records, log_probs, dict_longest_align, zero_locs)

