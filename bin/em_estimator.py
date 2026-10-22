# Portions of this file are adapted from Emu (https://github.com/treangenlab/emu),
# Copyright (c) 2021 Kristen Curry, Qi Wang, Michael G. Nute, Alona Tyshaieva,
# Qinglong Wu, Tor Savidge, Alexander Dilthey, Todd J. Treangen. MIT License;
# see THIRD_PARTY_LICENSES for the full text.
# Adapted here: variational_inference, update_q_distribution, maximize_elbo and
# variational_inference_iterations (Emu's EM loop, expectation_maximization and
# expectation_maximization_iterations).

"""
Emu-style expectation-maximization estimator for NanoVI (--estimator em).

A point estimate of the species frequencies, kept as a baseline for the
Bayesian variational inference in variational_inference_bayesian.py, which the
pipeline runs by default (decision D2).
"""

import math
import numpy as np
from sys import stdout


def variational_inference(log_p_rgs, freq, max_iterations=20, tolerance=5e-4):
    """Run inner VI/EM loop: update q(s|r) and frequencies until convergence."""
    for iteration in range(max_iterations):
        q_dist, elbo = update_q_distribution(log_p_rgs, freq)
        new_freq = maximize_elbo(q_dist)
        for key in freq:
            if key not in new_freq:
                new_freq[key] = 0.0
        diff = sum(abs(new_freq[key] - freq.get(key, 0.0)) for key in new_freq)
        if diff < tolerance:
            print(f"Converged at iteration {iteration + 1} with ELBO = {elbo}")
            break
        freq = new_freq
        print(f"Iteration {iteration + 1}: ELBO = {elbo}")
    return freq, elbo

def update_q_distribution(log_p_rgs, freq):
    q_dist = {}
    elbo = 0
    for read in log_p_rgs:
        valid_seqs, log_p_rns = [], []
        for seq in range(len(log_p_rgs[read][0])):
            s_val = log_p_rgs[read][0][seq]
            if s_val in freq and freq[s_val] != 0:
                logprns_val = log_p_rgs[read][1][seq] + math.log(freq[s_val])
                valid_seqs.append(s_val)
                log_p_rns.append(logprns_val)
        if len(valid_seqs) != 0:
            logc = -np.max(log_p_rns)
            prnsc = np.exp(log_p_rns + logc)
            prc = np.sum(prnsc)
            elbo += (np.log(prc) - logc)
            for seq in enumerate(valid_seqs):
                tax_id, prob = seq[1], prnsc[seq[0]] / prc
                if tax_id not in q_dist:
                    q_dist[tax_id] = {}
                q_dist[tax_id][read] = prob
    return q_dist, elbo

def maximize_elbo(q_dist):
    """M-step: re-estimate species frequencies from q(s|r) assignments."""
    freq = {}
    for tax_id, read_id in q_dist.items():
        freq[tax_id] = sum(read_id.values())
    total = sum(freq.values())
    for tax_id in freq:
        freq[tax_id] /= total
    return freq

def count_below_detection(log_p_rgs, reported):
    """Reads with a score none of whose species is reported with a positive abundance (F06)."""
    return sum(
        1 for sids, _ in log_p_rgs.values()
        if not any(reported.get(sid, 0) > 0 for sid in sids)
    )

def variational_inference_iterations(log_p_rgs, db_ids, lli_thresh, input_threshold, prune=True):
    """EM with pruning, as in Emu.

    Returns (freq_full, freq_set_thresh, below_detection_count): the frequencies
    after pruning species below 1/R (10/R above 1,000 reads); the frequencies
    refitted on the species at or above input_threshold, or None when that
    threshold is not above the pruning one; and the number of reads whose
    species are all missing from the reported set (the second if given, else
    the first), which are reported as 'Below detection' (F06).

    prune=False (experimental, Phase 7a) sets the pruning threshold to 0, so
    every species with data is kept until input_threshold applies.
    """
    n_db = len(db_ids)
    n_reads = len(log_p_rgs)
    stdout.write(f"Assigned read count: {n_reads}\n")
    if n_reads == 0:
        raise ValueError(
            "0 reads were assigned to any reference taxon. "
            "Check that input reads overlap with the reference database "
            "and that filtering parameters (--min_length, --max_length) are appropriate."
        )
    freq, counter = dict.fromkeys(db_ids, 1 / n_db), 1
    freq_thresh = 1 / n_reads
    if n_reads > 1000:
        freq_thresh = 10 / n_reads
    if not prune:
        freq_thresh = 0.0
    total_elbo = -math.inf
    while True:
        freq, _ = variational_inference(log_p_rgs, freq)
        _, new_elbo = update_q_distribution(log_p_rgs, freq)
        elbo_diff = new_elbo - total_elbo
        total_elbo = new_elbo
        if elbo_diff < 0:
            raise ValueError(
                f"ELBO decreased by {elbo_diff:.6f} from prior iteration, "
                "indicating numerical instability. Check input data quality."
            )
        if elbo_diff < lli_thresh:
            stdout.write(f"Number of VI iterations: {counter}\n")
            freq = {k: v for k, v in freq.items() if v >= freq_thresh}
            freq_full, _ = variational_inference(log_p_rgs, freq)
            freq_set_thresh = None
            if freq_thresh < input_threshold:
                freq = {k: v for k, v in freq_full.items() if v >= input_threshold}
                freq_set_thresh, _ = variational_inference(log_p_rgs, freq)
            reported = freq_set_thresh if freq_set_thresh else freq_full
            below_detection_count = count_below_detection(log_p_rgs, reported)
            stdout.write(f"Reads below detection (all their species pruned): {below_detection_count}\n")
            return freq_full, freq_set_thresh, below_detection_count
        counter += 1
