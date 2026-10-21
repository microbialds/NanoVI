import math
import sys
import numpy as np
from scipy.special import digamma, gammaln
from scipy.stats import beta as beta_dist
from sys import stdout
from typing import Dict, Tuple, Optional, List

# Largest ELBO fall between CAVI iterations, relative to |ELBO|, taken as rounding
ELBO_DECREASE_TOLERANCE = 1e-9


# ============================================================================
#  DATA STRUCTURES
# ============================================================================

def build_likelihood_matrix(
    log_p_rgs: Dict[str, Tuple[List[int], List[float]]],
    species_ids: List[int],
) -> Tuple[np.ndarray, List[str], Dict[int, int]]:
    """
    Convert the sparse dict-of-lists log-likelihood representation into a dense
    (R × S) matrix for efficient vectorized CAVI updates.

    Parameters
    ----------
    log_p_rgs : dict
        Keys are read names. Values are (list_of_species_ids, list_of_log_scores).
        This is the same format produced by log_prob_rgs_dict().
    species_ids : list of int
        Full list of species taxonomy IDs in the reference database.

    Returns
    -------
    L : np.ndarray, shape (R, S)
        Log-likelihood matrix. L[r, s] = ln p(x_r | z_r = s).
        Entries with no alignment are set to -inf (zero probability).
    read_names : list of str
        Ordered read names corresponding to rows of L.
    sid_to_col : dict
        Mapping from species taxonomy ID to column index in L.
    """
    sid_to_col = {sid: j for j, sid in enumerate(species_ids)}
    read_names = list(log_p_rgs.keys())
    R = len(read_names)
    S = len(species_ids)

    # Initialize with -inf (log(0) = impossible alignment)
    L = np.full((R, S), -np.inf, dtype=np.float64)

    for r, rname in enumerate(read_names):
        sids, scores = log_p_rgs[rname]
        for sid, score in zip(sids, scores):
            if sid in sid_to_col:
                col = sid_to_col[sid]
                # Keep the best score if multiple alignments to same species
                L[r, col] = max(L[r, col], score)

    return L, read_names, sid_to_col


# ============================================================================
#  CORE CAVI ALGORITHM
# ============================================================================

def compute_dirichlet_expectation(alpha: np.ndarray) -> np.ndarray:
    """
    Compute E_q[ln π_s] under q(π) = Dirichlet(α).

    E[ln π_s] = ψ(α_s) − ψ(Σ_j α_j)

    where ψ(·) is the digamma function.

    This is the key quantity that distinguishes VI from EM: in EM, one uses
    ln(π_s) directly (the point estimate), whereas in VI, the digamma
    function introduces Bayesian shrinkage toward the prior.

    Parameters
    ----------
    alpha : np.ndarray, shape (S,)
        Dirichlet variational parameters.

    Returns
    -------
    E_ln_pi : np.ndarray, shape (S,)
        Expected log-frequencies under the variational Dirichlet.
    """
    return digamma(alpha) - digamma(np.sum(alpha))


def update_responsibilities(
    L: np.ndarray,
    E_ln_pi: np.ndarray,
) -> np.ndarray:
    """
    E-step: update variational categorical parameters φ_{r,s} for each read.

    φ_{r,s} ∝ exp( ln p(x_r | s) + E_q[ln π_s] )

    Normalization is performed in log-space for numerical stability using
    the log-sum-exp trick.

    Parameters
    ----------
    L : np.ndarray, shape (R, S)
        Log-likelihood matrix.
    E_ln_pi : np.ndarray, shape (S,)
        Expected log-frequencies from compute_dirichlet_expectation().

    Returns
    -------
    phi : np.ndarray, shape (R, S)
        Responsibility matrix. phi[r, s] = q(z_r = s).
        Each row sums to 1.
    """
    # Unnormalized log-responsibilities
    log_phi = L + E_ln_pi[np.newaxis, :]  # broadcasting: (R, S)

    # Log-sum-exp normalization per read (row-wise)
    log_phi_max = np.max(log_phi, axis=1, keepdims=True)

    # Guard against rows where all entries are -inf (read maps to nothing)
    finite_mask = np.isfinite(log_phi_max.ravel())
    phi = np.zeros_like(log_phi)

    if np.any(finite_mask):
        safe_rows = finite_mask
        log_phi_safe = log_phi[safe_rows]
        max_safe = log_phi_max[safe_rows]
        phi_unnorm = np.exp(log_phi_safe - max_safe)
        phi[safe_rows] = phi_unnorm / phi_unnorm.sum(axis=1, keepdims=True)

    return phi


def update_dirichlet_parameters(
    alpha_0: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    """
    M-step: update variational Dirichlet parameters α.

    α_s = α₀_s + Σ_r φ_{r,s}

    The sufficient statistic Σ_r φ_{r,s} is the expected number of reads
    assigned to species s under the current variational distribution.

    Parameters
    ----------
    alpha_0 : np.ndarray, shape (S,)
        Prior Dirichlet parameters.
    phi : np.ndarray, shape (R, S)
        Current responsibility matrix.

    Returns
    -------
    alpha : np.ndarray, shape (S,)
        Updated variational Dirichlet parameters.
    """
    return alpha_0 + phi.sum(axis=0)


def compute_elbo(
    L: np.ndarray,
    phi: np.ndarray,
    alpha: np.ndarray,
    alpha_0: np.ndarray,
    E_ln_pi: np.ndarray,
) -> float:
    """
    Compute the Evidence Lower BOund (ELBO).

    ELBO = E_q[ln p(x, z, π)] − E_q[ln q(z, π)]

    Decomposed into three terms:

    (1) Data + assignment term:
        Σ_r Σ_s φ_{r,s} [ ln p(x_r|s) + E[ln π_s] − ln φ_{r,s} ]

    (2) Prior − variational Dirichlet KL divergence:
        ln B(α) − ln B(α₀) + Σ_s (α₀_s − α_s) E[ln π_s]

        where ln B(α) = Σ_s ln Γ(α_s) − ln Γ(Σ_s α_s)

    The ELBO is guaranteed to be non-decreasing across CAVI iterations
    (up to numerical precision), providing a convergence diagnostic.

    Parameters
    ----------
    L : np.ndarray, shape (R, S)
        Log-likelihood matrix.
    phi : np.ndarray, shape (R, S)
        Responsibility matrix.
    alpha : np.ndarray, shape (S,)
        Current variational Dirichlet parameters.
    alpha_0 : np.ndarray, shape (S,)
        Prior Dirichlet parameters.
    E_ln_pi : np.ndarray, shape (S,)
        Expected log-frequencies.

    Returns
    -------
    elbo : float
        The evidence lower bound.
    """
    # Avoid log(0) in entropy term
    log_phi_safe = np.where(phi > 0, np.log(phi), 0.0)

    # Term 1: E_q[ln p(x, z)] + H[q(z)]
    # Only sum over entries where phi > 0 (finite contributions)
    mask = phi > 0
    term1 = np.sum(
        phi[mask] * (L[mask] + np.broadcast_to(E_ln_pi, L.shape)[mask]
                     - log_phi_safe[mask])
    )

    # Term 2: KL(q(π) || p(π))  [note: ELBO subtracts KL, so this is negative KL]
    # ln B(α) − ln B(α₀) + Σ_s (α₀_s − α_s) E[ln π_s]
    log_B_alpha0 = np.sum(gammaln(alpha_0)) - gammaln(np.sum(alpha_0))
    log_B_alpha = np.sum(gammaln(alpha)) - gammaln(np.sum(alpha))
    term2 = log_B_alpha - log_B_alpha0 + np.sum((alpha_0 - alpha) * E_ln_pi)

    return term1 + term2


def cavi_inner_loop(
    L: np.ndarray,
    alpha_0: np.ndarray,
    alpha_init: np.ndarray,
    max_iterations: int = 1000,
    tolerance: float = 1e-7,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, int]:
    """
    Run the inner CAVI loop until convergence.

    Alternates between:
      1. Update E[ln π] from current α  (Dirichlet expectation)
      2. Update φ (responsibilities)      (E-step analogue)
      3. Update α from φ                  (M-step analogue)
      4. Compute ELBO for monitoring

    The loop stops when the largest change in posterior mean abundance,
    max_s |Δ(α_s / Σ_j α_j)|, falls below `tolerance` (F05). A relative ELBO
    change is not used: the ELBO grows with the number of reads, so a relative
    rule stops large samples early. CAVI never lowers the ELBO; a fall of more
    than ELBO_DECREASE_TOLERANCE × |ELBO| is reported on stderr as a sign of
    numerical trouble.

    Parameters
    ----------
    L : np.ndarray, shape (R, S)
        Log-likelihood matrix.
    alpha_0 : np.ndarray, shape (S,)
        Prior Dirichlet parameters.
    alpha_init : np.ndarray, shape (S,)
        Initial variational Dirichlet parameters.
    max_iterations : int
        Maximum number of CAVI iterations.
    tolerance : float
        Convergence threshold on the largest change in posterior mean abundance.
    verbose : bool
        If True, print progress.

    Returns
    -------
    alpha : np.ndarray, shape (S,)
        Converged variational Dirichlet parameters.
    phi : np.ndarray, shape (R, S)
        Final responsibility matrix.
    E_ln_pi : np.ndarray, shape (S,)
        Final expected log-frequencies.
    elbo : float
        Final ELBO value.
    n_iters : int
        Number of iterations until convergence.
    """
    alpha = alpha_init.copy()
    mean = alpha / np.sum(alpha)
    prev_elbo = None
    E_ln_pi = compute_dirichlet_expectation(alpha)

    for iteration in range(1, max_iterations + 1):
        # Step 1: Dirichlet expectations of the current α (E_ln_pi, computed
        # before the loop or at the end of the previous iteration)

        # Step 2: Update responsibilities (E-step)
        phi = update_responsibilities(L, E_ln_pi)

        # Step 3: Update Dirichlet parameters (M-step)
        alpha = update_dirichlet_parameters(alpha_0, phi)

        # Step 4: Compute ELBO
        E_ln_pi = compute_dirichlet_expectation(alpha)  # Recompute with updated α
        elbo = compute_elbo(L, phi, alpha, alpha_0, E_ln_pi)

        new_mean = alpha / np.sum(alpha)
        max_change = float(np.max(np.abs(new_mean - mean)))
        mean = new_mean

        if prev_elbo is not None and elbo < prev_elbo - ELBO_DECREASE_TOLERANCE * abs(prev_elbo):
            sys.stderr.write(
                f"  WARNING: ELBO fell from {prev_elbo:.6f} to {elbo:.6f} "
                f"at CAVI iteration {iteration}; this may indicate numerical instability.\n"
            )
        prev_elbo = elbo

        if verbose and (iteration == 1 or iteration % 100 == 0):
            stdout.write(
                f"  CAVI iteration {iteration}: "
                f"ELBO = {elbo:.4f}, "
                f"max Δ(mean) = {max_change:.2e}\n"
            )
        if max_change < tolerance:
            if verbose:
                stdout.write(
                    f"  Converged at iteration {iteration} "
                    f"(ELBO = {elbo:.4f}, max Δ(mean) = {max_change:.2e})\n"
                )
            return alpha, phi, E_ln_pi, elbo, iteration

    if verbose:
        stdout.write(
            f"  Reached max iterations ({max_iterations}). "
            f"ELBO = {elbo:.4f}, max Δ(mean) = {max_change:.2e}\n"
        )
    return alpha, phi, E_ln_pi, elbo, max_iterations


# ============================================================================
#  OUTER OPTIMIZATION WITH PRUNING
# ============================================================================

def bayesian_vi_iterations(
    log_p_rgs: Dict[str, Tuple[List[int], List[float]]],
    db_species_ids: List[int],
    elbo_threshold: float = 1e-4,
    abundance_threshold: float = 1e-3,
    alpha_prior: float = 1.0,
    max_outer_iterations: int = 10,
    max_inner_iterations: int = 1000,
    inner_tolerance: float = 1e-7,
) -> Tuple[Dict[int, float], Optional[Dict[int, float]], Dict[int, Tuple[float, float]], int]:
    """
    Outer optimization loop with iterative species pruning.

    This function:
    1. Runs CAVI to convergence on the full species set.
    2. Prunes species with posterior mean abundance below a threshold.
    3. Re-runs CAVI on the reduced set, warm-started from the previous α.
    4. Repeats until no species is below the threshold (F07). After
       max_outer_iterations prunings it fits the pruned set once more and
       stops, with a warning if species are still below the threshold.
    5. Returns both full and thresholded abundance estimates,
       plus posterior credible intervals.

    Parameters
    ----------
    log_p_rgs : dict
        Read-to-species log-likelihood dictionary (from log_prob_rgs_dict).
    db_species_ids : list of int
        All species taxonomy IDs in the reference database.
    elbo_threshold : float
        Ignored. The outer loop no longer compares ELBOs across species sets
        (F07); kept so existing callers keep working.
    abundance_threshold : float
        Minimum posterior mean abundance to retain a species.
    alpha_prior : float
        Total Dirichlet prior concentration (shared across all species).
        Each species receives alpha_prior / S pseudo-counts, so the total
        prior contribution is always alpha_prior regardless of S.
        - alpha_prior = 1.0: weakly informative (1 pseudo-read total, default).
        - alpha_prior < 1.0: sparser prior (encourages fewer species).
        - alpha_prior > 1.0: smoothing prior (shrinks toward uniform).
    max_outer_iterations : int
        Maximum number of pruning iterations; one more fit follows the last.
    max_inner_iterations : int
        Maximum CAVI iterations per inner loop.
    inner_tolerance : float
        Convergence tolerance for the inner CAVI loop: the largest change in
        posterior mean abundance between iterations (F05).

    Returns
    -------
    freq_full : dict
        Species ID → posterior mean abundance (after pruning, renormalized).
    freq_thresholded : dict or None
        Species ID → posterior mean abundance after applying abundance_threshold.
        None if no thresholding needed.
    credible_intervals : dict
        Species ID → (lower_95, upper_95) credible intervals on abundance.
    below_detection_count : int
        Number of reads whose responsibilities are all zero at the end: every
        species they align to was pruned (F06). They are not part of the
        abundances and should not be counted as assigned to the reported species.

    Raises
    ------
    ValueError
        If no reads are assigned.
    """
    n_reads = len(log_p_rgs)
    stdout.write(f"Assigned read count: {n_reads}\n")

    if n_reads == 0:
        raise ValueError(
            "0 reads were assigned to any reference taxon. "
            "Check that input reads overlap with the reference database "
            "and that filtering parameters (--min_length, --max_length) "
            "are appropriate."
        )

    # --- Build dense log-likelihood matrix ---
    active_species = list(db_species_ids)
    L_full, read_names, sid_to_col_full = build_likelihood_matrix(
        log_p_rgs, active_species
    )

    # --- Determine adaptive pruning threshold ---
    freq_thresh = max(1.0 / n_reads, abundance_threshold)
    if n_reads > 1000:
        freq_thresh = max(10.0 / n_reads, abundance_threshold)

    stdout.write(
        f"Running Bayesian VI with Dirichlet prior "
        f"(total α₀ = {alpha_prior}, per-species = {alpha_prior}/{len(active_species)})\n"
        f"Initial species count: {len(active_species)}\n"
        f"Pruning threshold: {freq_thresh:.2e}\n"
    )

    # --- Pre-filter species with zero evidence ---
    L = L_full.copy()
    current_species = list(active_species)
    col_has_data = np.any(np.isfinite(L), axis=0)
    if not np.all(col_has_data):
        n_empty = int(np.sum(~col_has_data))
        L = L[:, col_has_data]
        current_species = [
            current_species[i]
            for i in range(len(current_species))
            if col_has_data[i]
        ]
        stdout.write(
            f"Removed {n_empty} species with zero aligned reads "
            f"(remaining: {len(current_species)})\n"
        )

    # Reads with at least one finite score; only these start the uniform α (F06)
    n_scored = int(np.sum(np.any(np.isfinite(L), axis=1)))

    # --- Outer loop: CAVI + prune ---
    # The loop stops only when no species is below the pruning threshold. ELBOs
    # of successive iterations are not compared: each belongs to a different
    # species set and prior (F07). α is warm-started from the previous solution.
    expected_counts = None
    inner_iterations = []
    final_fit = False

    for outer_iter in range(1, max_outer_iterations + 2):
        S = len(current_species)
        stdout.write(f"\n=== Outer iteration {outer_iter} (S = {S}) ===\n")

        # Symmetric Dirichlet prior scaled so total pseudo-count = alpha_prior
        alpha_0 = np.full(S, alpha_prior / S)

        if expected_counts is None:
            # Initialize α uniformly (prior + uniform pseudo-counts)
            alpha_init = alpha_0 + (n_scored / S)
        else:
            # Warm start: the surviving species keep their expected read counts
            alpha_init = alpha_0 + expected_counts

        # Run inner CAVI
        alpha, phi, E_ln_pi, elbo, n_iters = cavi_inner_loop(
            L, alpha_0, alpha_init,
            max_iterations=max_inner_iterations,
            tolerance=inner_tolerance,
        )
        inner_iterations.append(n_iters)
        stdout.write(f"  ELBO = {elbo:.4f} after {n_iters} CAVI iterations\n")

        # Compute posterior mean abundances: E[π_s] = α_s / Σ_j α_j
        posterior_mean = alpha / np.sum(alpha)
        keep_mask = posterior_mean >= freq_thresh

        if final_fit:
            n_below = int(S - np.sum(keep_mask))
            if n_below:
                sys.stderr.write(
                    f"WARNING: reached {max_outer_iterations} pruning iterations; "
                    f"{n_below} species remain below the pruning threshold "
                    f"({freq_thresh:.2e}).\n"
                )
            break

        if not np.any(keep_mask):
            # Every species is below the threshold (more species than reads can
            # support); keep the most abundant one rather than none
            keep_mask[np.argmax(posterior_mean)] = True
        n_pruned = int(S - np.sum(keep_mask))

        if n_pruned == 0:
            stdout.write("  Converged: no species below the pruning threshold.\n")
            break

        stdout.write(f"  Pruning {n_pruned} species below threshold.\n")

        # Update active species, likelihood matrix and warm-start counts
        keep_indices = np.where(keep_mask)[0]
        current_species = [current_species[i] for i in keep_indices]
        L = L[:, keep_indices].copy()
        expected_counts = (alpha - alpha_0)[keep_indices]

        # Remove columns that are all -inf (no reads map to them)
        col_has_data = np.any(np.isfinite(L), axis=0)
        if not np.all(col_has_data):
            L = L[:, col_has_data]
            expected_counts = expected_counts[col_has_data]
            current_species = [
                current_species[i]
                for i in range(len(current_species))
                if col_has_data[i]
            ]

        # The last allowed pruning is followed by one fit on the pruned set
        final_fit = outer_iter == max_outer_iterations

    # --- Final results ---
    # Reads left with no surviving species (F06)
    below_detection_count = int(np.sum(phi.sum(axis=1) == 0))

    # Posterior mean and credible intervals from the Dirichlet
    alpha_sum = np.sum(alpha)
    posterior_mean = alpha / alpha_sum

    # 95% credible intervals from exact Beta marginals:
    # π_s | α ~ Beta(α_s, α_sum - α_s)
    credible_intervals = {}
    freq_full = {}

    for i, sid in enumerate(current_species):
        mean_s = posterior_mean[i]
        a_s = alpha[i]
        b_s = alpha_sum - a_s
        lower, upper = beta_dist.ppf([0.025, 0.975], a_s, b_s)
        freq_full[sid] = float(mean_s)
        credible_intervals[sid] = (float(lower), float(upper))

    # Renormalize (should already sum to 1, but ensures consistency)
    total = sum(freq_full.values())
    if total > 0:
        freq_full = {k: v / total for k, v in freq_full.items()}

    # Apply user-facing abundance threshold
    freq_thresholded = None
    if abundance_threshold > 0:
        freq_thresh_dict = {
            k: v for k, v in freq_full.items() if v >= abundance_threshold
        }
        if len(freq_thresh_dict) < len(freq_full):
            total_t = sum(freq_thresh_dict.values())
            if total_t > 0:
                freq_thresholded = {
                    k: v / total_t for k, v in freq_thresh_dict.items()
                }

    stdout.write(
        f"\nFinal species count: {len(freq_full)}\n"
        f"Reads below detection (all their species pruned): {below_detection_count}\n"
        f"Number of VI outer iterations: {outer_iter} "
        f"(CAVI iterations: {', '.join(map(str, inner_iterations))})\n"
    )

    return freq_full, freq_thresholded, credible_intervals, below_detection_count


# ============================================================================
#  POSTERIOR ANALYSIS UTILITIES
# ============================================================================

def posterior_entropy(alpha: np.ndarray) -> float:
    """
    Compute the differential entropy of the posterior Dirichlet.

    H[Dir(α)] = ln B(α) + (α₀ − S) ψ(α₀) − Σ_s (α_s − 1) ψ(α_s)

    where α₀ = Σ_s α_s and B(α) = ∏ Γ(α_s) / Γ(α₀).

    Higher entropy indicates more uncertainty in the abundance estimates.

    Parameters
    ----------
    alpha : np.ndarray
        Dirichlet parameters.

    Returns
    -------
    entropy : float
        Differential entropy in nats.
    """
    alpha_0 = np.sum(alpha)
    S = len(alpha)
    log_B = np.sum(gammaln(alpha)) - gammaln(alpha_0)
    return (
        log_B
        + (alpha_0 - S) * digamma(alpha_0)
        - np.sum((alpha - 1) * digamma(alpha))
    )


def kl_divergence_dirichlet(alpha: np.ndarray, beta: np.ndarray) -> float:
    """
    KL divergence between two Dirichlet distributions.

    KL(Dir(α) || Dir(β)) = ln B(β)/B(α) + Σ_s (α_s − β_s)[ψ(α_s) − ψ(α₀)]

    Useful for measuring how far the posterior has moved from the prior.

    Parameters
    ----------
    alpha : np.ndarray
        Parameters of the first (approximate posterior) Dirichlet.
    beta : np.ndarray
        Parameters of the second (prior) Dirichlet.

    Returns
    -------
    kl : float
        KL divergence in nats.
    """
    alpha_0 = np.sum(alpha)
    beta_0 = np.sum(beta)
    log_B_beta = np.sum(gammaln(beta)) - gammaln(beta_0)
    log_B_alpha = np.sum(gammaln(alpha)) - gammaln(alpha_0)
    return (
        log_B_beta - log_B_alpha
        + np.sum((alpha - beta) * (digamma(alpha) - digamma(alpha_0)))
    )


def expected_effective_species(alpha: np.ndarray) -> float:
    """
    Compute the expected effective number of species (exp of Shannon entropy)
    under the posterior Dirichlet.

    This uses the first-order approximation:
        E[H(π)] ≈ −Σ_s E[π_s] ln E[π_s]

    and returns exp(E[H(π)]).

    Parameters
    ----------
    alpha : np.ndarray
        Posterior Dirichlet parameters.

    Returns
    -------
    n_eff : float
        Expected effective number of species.
    """
    pi_mean = alpha / np.sum(alpha)
    pi_mean = pi_mean[pi_mean > 0]
    H = -np.sum(pi_mean * np.log(pi_mean))
    return np.exp(H)


# ============================================================================
#  DROP-IN REPLACEMENT WRAPPER
# ============================================================================

def variational_inference_iterations(
    log_p_rgs: Dict[str, Tuple[List[int], List[float]]],
    db_ids: List[int],
    lli_thresh: float = 1e-4,
    input_threshold: float = 1e-3,
) -> Tuple[Dict[int, float], Optional[Dict[int, float]], int]:
    """
    Drop-in replacement for the original variational_inference_iterations().

    This function keeps the signature of the original EM-based implementation
    and returns a 3-tuple like it, but the third value is the number of reads
    below detection (F06) instead of the EM version's unused None.

    Parameters
    ----------
    log_p_rgs : dict
        Read-to-species log-likelihood dictionary.
    db_ids : list of int
        Species taxonomy IDs in the reference database.
    lli_thresh : float
        Ignored by the Bayesian VI (maps to elbo_threshold, which the outer
        loop no longer uses; F07).
    input_threshold : float
        Minimum abundance for output (maps to abundance_threshold).

    Returns
    -------
    freq_full : dict
        Species ID → posterior mean abundance.
    freq_thresholded : dict or None
        Abundance estimates after thresholding, or None.
    below_detection_count : int
        Reads whose species were all pruned (F06). This position held an
        unused read_dist placeholder (None).
    """
    freq_full, freq_thresholded, _, below_detection_count = bayesian_vi_iterations(
        log_p_rgs=log_p_rgs,
        db_species_ids=db_ids,
        elbo_threshold=lli_thresh,
        abundance_threshold=input_threshold,
        alpha_prior=1.0,  # Uniform Dirichlet prior
        max_outer_iterations=10,
        max_inner_iterations=1000,
        inner_tolerance=1e-7,
    )
    return freq_full, freq_thresholded, below_detection_count
