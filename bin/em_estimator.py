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

from variational_inference_bayesian import posterior_read_summary, row_groups


def variational_inference(log_p_rgs, freq, max_iterations=20, tolerance=5e-4):
    """Run inner VI/EM loop: update q(s|r) and frequencies until convergence.

    Returns the frequencies of the last update, which sum to 1 (F47; the input,
    for example a pruned set, need not), and the log-likelihood of the
    frequencies that update started from.
    """
    for iteration in range(max_iterations):
        q_dist, elbo = update_q_distribution(log_p_rgs, freq)
        new_freq = maximize_elbo(q_dist)
        for key in freq:
            if key not in new_freq:
                new_freq[key] = 0.0
        diff = sum(abs(new_freq[key] - freq.get(key, 0.0)) for key in new_freq)
        freq = new_freq
        if diff < tolerance:
            print(f"Converged at iteration {iteration + 1} with ELBO = {elbo}")
            break
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

def read_summary(idx, logL, freq, species_ids, reported_ids, species_rows=None, resp=None):
    """Expected and shared reads of each reported species under EM's responsibilities (D27).

    idx/logL are the padded likelihood over the species that survived pruning
    (columns in the order of species_ids, as sparse_likelihood_from_csr gives them)
    and freq their frequencies (freq_full). A read's responsibility is
    ∝ freq_s × p(read | s) over those species, as the Bayesian VI's φ comes from its
    final fit; expected and shared reads are then counted over the reported species
    (reported_ids) as in posterior_read_summary, without alpha; species_rows
    (species ID → table row) counts as shared only the reads split between rows (F46).
    resp: (idx_reported, phi) from responsibilities with the same arguments, when
    the caller has them already.
    """
    idx_reported, phi = resp if resp is not None else responsibilities(idx, logL, freq, species_ids, reported_ids)
    return posterior_read_summary(None, idx_reported, np.where(idx_reported >= 0, phi, 0.0), list(reported_ids),
                                  groups=row_groups(list(reported_ids), species_rows))


def responsibilities(idx, logL, freq, species_ids, reported_ids):
    """EM's responsibilities, ∝ freq_s × p(read | s), as read_summary uses them.

    Returns (idx_reported, phi): phi over the columns of idx (species_ids, with
    frequencies freq), each read's row summing to 1 (0 for a read none of whose
    species survived pruning), and each entry's column among reported_ids, -1
    for padding and for species that are not reported (the per-read
    assignments, F55).
    """
    finite = np.isfinite(logL)
    with np.errstate(divide="ignore"):
        log_w = np.log(np.asarray(freq, dtype=np.float64))
    score = np.where(finite, logL + log_w[idx], -np.inf)
    top = score.max(axis=1, keepdims=True)
    phi = np.exp(score - np.where(np.isfinite(top), top, 0.0))
    total = phi.sum(axis=1, keepdims=True)
    phi = np.divide(phi, total, out=np.zeros_like(phi), where=total > 0)
    column = {sid: i for i, sid in enumerate(reported_ids)}
    reported_col = np.array([column.get(sid, -1) for sid in species_ids], dtype=np.int64)
    idx_reported = np.where(np.isfinite(score), reported_col[idx], -1)
    return idx_reported, phi


def variational_inference_iterations(log_p_rgs, db_ids, lli_thresh, input_threshold):
    """EM with pruning, as in Emu.

    Returns (freq_full, freq_set_thresh, below_detection_count): the frequencies
    after pruning species below 1/R (10/R above 1,000 reads); the frequencies
    refitted on the species at or above input_threshold, or None when that
    threshold is not above the pruning one; and the number of reads whose
    species are all missing from the reported set (the second if given, else
    the first), which are reported as 'Below detection' (F06).
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
