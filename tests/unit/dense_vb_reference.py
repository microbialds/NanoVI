"""Dense reference for the Bayesian VI: dev's implementation before F23.

bin/variational_inference_bayesian.py stores the likelihood as padded
reads x candidates arrays (F23). This module keeps the original dense
reads x species version, unchanged apart from comments and progress output,
so the tests can check that the sparse version gives the same results.
"""
import itertools

import numpy as np
from scipy.special import digamma, gammaln

SHARED_READ_MIN_SPLIT = 0.01


def build_likelihood_matrix(log_p_rgs, species_ids):
    """Dense (R x S) log-likelihood matrix, -inf where a read has no alignment."""
    sid_to_col = {sid: j for j, sid in enumerate(species_ids)}
    read_names = list(log_p_rgs.keys())
    L = np.full((len(read_names), len(species_ids)), -np.inf, dtype=np.float64)
    for r, rname in enumerate(read_names):
        sids, scores = log_p_rgs[rname]
        for sid, score in zip(sids, scores):
            if sid in sid_to_col:
                col = sid_to_col[sid]
                L[r, col] = max(L[r, col], score)
    return L, read_names, sid_to_col


def compute_dirichlet_expectation(alpha):
    return digamma(alpha) - digamma(np.sum(alpha))


def update_responsibilities(L, E_ln_pi):
    log_phi = L + E_ln_pi[np.newaxis, :]
    log_phi_max = np.max(log_phi, axis=1, keepdims=True)
    finite_mask = np.isfinite(log_phi_max.ravel())
    phi = np.zeros_like(log_phi)
    if np.any(finite_mask):
        phi_unnorm = np.exp(log_phi[finite_mask] - log_phi_max[finite_mask])
        phi[finite_mask] = phi_unnorm / phi_unnorm.sum(axis=1, keepdims=True)
    return phi


def update_dirichlet_parameters(alpha_0, phi):
    return alpha_0 + phi.sum(axis=0)


def compute_elbo(L, phi, alpha, alpha_0, E_ln_pi):
    log_phi_safe = np.where(phi > 0, np.log(phi), 0.0)
    mask = phi > 0
    term1 = np.sum(
        phi[mask] * (L[mask] + np.broadcast_to(E_ln_pi, L.shape)[mask] - log_phi_safe[mask])
    )
    log_B_alpha0 = np.sum(gammaln(alpha_0)) - gammaln(np.sum(alpha_0))
    log_B_alpha = np.sum(gammaln(alpha)) - gammaln(np.sum(alpha))
    term2 = log_B_alpha - log_B_alpha0 + np.sum((alpha_0 - alpha) * E_ln_pi)
    return term1 + term2


def cavi_inner_loop(L, alpha_0, alpha_init, max_iterations=1000, tolerance=1e-7):
    alpha = alpha_init.copy()
    mean = alpha / np.sum(alpha)
    E_ln_pi = compute_dirichlet_expectation(alpha)
    for iteration in range(1, max_iterations + 1):
        phi = update_responsibilities(L, E_ln_pi)
        alpha = update_dirichlet_parameters(alpha_0, phi)
        E_ln_pi = compute_dirichlet_expectation(alpha)
        elbo = compute_elbo(L, phi, alpha, alpha_0, E_ln_pi)
        new_mean = alpha / np.sum(alpha)
        max_change = float(np.max(np.abs(new_mean - mean)))
        mean = new_mean
        if max_change < tolerance:
            return alpha, phi, E_ln_pi, elbo, iteration
    return alpha, phi, E_ln_pi, elbo, max_iterations


def posterior_read_summary(alpha, phi, species_ids):
    expected = phi.sum(axis=0)
    row_total = phi.sum(axis=1, keepdims=True)
    phi_reported = np.divide(phi, row_total, out=np.zeros_like(phi), where=row_total > 0)
    split = (phi_reported > 0) & (phi_reported <= 1.0 - SHARED_READ_MIN_SPLIT)
    shared = np.where(split, phi, 0.0).sum(axis=0)
    return {
        sid: {
            "alpha": float(alpha[i]),
            "expected_reads": float(expected[i]),
            "shared_reads": float(shared[i]),
        }
        for i, sid in enumerate(species_ids)
    }


def bayesian_vi_iterations(log_p_rgs, db_species_ids, abundance_threshold=1e-3, alpha_prior=1.0,
                           max_outer_iterations=10, max_inner_iterations=1000, inner_tolerance=1e-7,
                           posterior_summary=False):
    """Dev's outer loop on the dense matrix. Returns (freq_full, freq_thresholded,
    below_detection_count, inner_iterations[, summary]); the unused credible
    intervals are left out."""
    n_reads = len(log_p_rgs)
    if n_reads == 0:
        raise ValueError("0 reads were assigned to any reference taxon.")
    L, _, _ = build_likelihood_matrix(log_p_rgs, list(db_species_ids))
    current_species = list(db_species_ids)

    freq_thresh = max(1.0 / n_reads, abundance_threshold)
    if n_reads > 1000:
        freq_thresh = max(10.0 / n_reads, abundance_threshold)

    col_has_data = np.any(np.isfinite(L), axis=0)
    if not np.all(col_has_data):
        L = L[:, col_has_data]
        current_species = [s for s, keep in zip(current_species, col_has_data) if keep]
    scored = np.any(np.isfinite(L), axis=1)
    n_scored = int(np.sum(scored))

    expected_counts = None
    inner_iterations = []
    final_fit = False
    for outer_iter in itertools.count(1):
        S = len(current_species)
        alpha_0 = np.full(S, alpha_prior / S)
        alpha_init = alpha_0 + (n_scored / S) if expected_counts is None else alpha_0 + expected_counts
        alpha, phi, E_ln_pi, elbo, n_iters = cavi_inner_loop(
            L, alpha_0, alpha_init, max_iterations=max_inner_iterations, tolerance=inner_tolerance)
        inner_iterations.append(n_iters)
        posterior_mean = alpha / np.sum(alpha)
        keep_mask = posterior_mean >= freq_thresh
        if final_fit:
            break
        if not np.any(keep_mask):
            keep_mask[np.argmax(posterior_mean)] = True
        if int(S - np.sum(keep_mask)) == 0:
            break
        keep_indices = np.where(keep_mask)[0]
        current_species = [current_species[i] for i in keep_indices]
        L = L[:, keep_indices].copy()
        expected_counts = (alpha - alpha_0)[keep_indices]
        final_fit = outer_iter >= max_outer_iterations

    posterior_mean = alpha / np.sum(alpha)
    freq_full = {sid: float(posterior_mean[i]) for i, sid in enumerate(current_species)}
    total = sum(freq_full.values())
    if total > 0:
        freq_full = {k: v / total for k, v in freq_full.items()}

    freq_thresholded = None
    if abundance_threshold > 0:
        kept = {k: v for k, v in freq_full.items() if v >= abundance_threshold}
        if len(kept) < len(freq_full):
            total_t = sum(kept.values())
            if total_t > 0:
                freq_thresholded = {k: v / total_t for k, v in kept.items()}

    reported = freq_thresholded if freq_thresholded is not None else freq_full
    reported_cols = [i for i, sid in enumerate(current_species) if sid in reported]
    below_detection_count = int(np.sum(scored & (phi[:, reported_cols].sum(axis=1) == 0)))

    result = (freq_full, freq_thresholded, below_detection_count, inner_iterations)
    if posterior_summary:
        summary = posterior_read_summary(
            alpha[reported_cols], phi[:, reported_cols], [current_species[i] for i in reported_cols])
        return result + (summary,)
    return result
