# Portions of this file are adapted from Emu (https://github.com/treangenlab/emu),
# Copyright (c) 2021 Kristen Curry, Qi Wang, Michael G. Nute, Alona Tyshaieva,
# Qinglong Wu, Tor Savidge, Alexander Dilthey, Todd J. Treangen. MIT License;
# see THIRD_PARTY_LICENSES for the full text.
# Adapted here: get_align_stats, get_align_len,
# get_cigar_op_log_probabilities, compute_log_prob_rgs and log_prob_rgs_dict (Emu's
# CIGAR-based read score).

"""
Read likelihoods for NanoVI.

1. CIGAR operation log-probabilities estimated from the sample's primary alignments
2. Per-read log-likelihood of each candidate species

The estimators (em_estimator.py, variational_inference_bayesian.py) take the
per-read scores that log_prob_rgs_dict returns.
"""

import math
import numpy as np
import pysam
from operator import add, mul
from sys import stdout

CIGAR_OPS = [1, 2, 4, 10]
CIGAR_OPS_ALL = [0, 1, 2, 4]

def get_align_stats(alignment):
    cigar_stats = alignment.get_cigar_stats()[0]
    n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]
    return [cigar_stats[1], cigar_stats[2], cigar_stats[4], n_mismatch]

def get_align_len(alignment):
    return sum(alignment.get_cigar_stats()[0][cigar_op] for cigar_op in CIGAR_OPS_ALL)

def get_cigar_op_log_probabilities(sam_path, threads=None):
    """Compute global log-probabilities for each CIGAR operation from the alignments.

    Reads the SAM once, streaming. Returns (log_probs, zero_locs, dict_longest_align):
    the log share of each error operation (I, D, S, X) over the primary mapped
    alignments, with the operations of zero count removed (F03); the indices of
    those operations in CIGAR_OPS; and each read's longest alignment over all of
    its alignments, secondary and supplementary included (F04).

    threads is ignored: the pass is serial (F25). It is kept so that existing
    callers keep working.

    Raises ValueError when the SAM holds no primary mapped alignment or when those
    alignments have no error operation.
    """
    cigar_stats_primary = [0] * len(CIGAR_OPS)
    dict_longest_align = {}
    n_mapped = 0

    with pysam.AlignmentFile(sam_path) as sam_pysam:
        for alignment in sam_pysam:
            align_len = get_align_len(alignment)
            query_name = alignment.query_name
            if dict_longest_align.get(query_name, -1) < align_len:
                dict_longest_align[query_name] = align_len
            # The error profile comes from primary mapped alignments only, as in Emu (F03, D3)
            if not alignment.is_unmapped and not alignment.is_secondary and not alignment.is_supplementary:
                n_mapped += 1
                cigar_stats_primary = list(map(add, cigar_stats_primary, get_align_stats(alignment)))

    if n_mapped == 0:
        raise ValueError(
            f"No primary mapped alignments in {sam_path}. Check that the reads overlap with the "
            "reference database and that filtering (--min_length, --max_length) left reads to align."
        )

    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]
    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]

    n_char = sum(cigar_stats_primary)
    if n_char == 0:
        raise ValueError(
            f"No insertions, deletions, soft clips or mismatches in the primary alignments of {sam_path}; "
            "the error profile cannot be estimated."
        )
    return [math.log(x) for x in np.array(cigar_stats_primary) / n_char], zero_locs, dict_longest_align

def compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len):
    ref_name, query_name = alignment.reference_name, alignment.query_name
    log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats))) * (dict_longest_align[query_name] / align_len)
    species_tid = int(ref_name.split(":")[0])
    return log_score, query_name, species_tid

def log_prob_rgs_dict(sam_path, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None):
    """Build dict mapping read_name -> (species_ids, log_scores) from SAM alignments.

    Each species keeps the best score among the read's alignments to it. An
    alignment that uses an operation with no count in the error profile
    (p_cigar_op_zero_locs) has no score and is skipped; its read counts as
    unassigned unless another alignment of the read is scored (F10).
    """
    zero_locs = sorted(p_cigar_op_zero_locs or [], reverse=True)
    log_p_rgs, unassigned_set = {}, set()

    with pysam.AlignmentFile(sam_path) as sam_pysam:
        for alignment in sam_pysam:
            align_len = get_align_len(alignment)
            if not (alignment.reference_name and align_len):
                unassigned_set.add(alignment.query_name)
                continue
            cigar_stats = get_align_stats(alignment)
            if any(cigar_stats[i] for i in zero_locs):
                unassigned_set.add(alignment.query_name)
                continue
            for i in zero_locs:
                del cigar_stats[i]
            log_score, query_name, species_tid = compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len)
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
