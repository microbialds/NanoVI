"""Unit tests for the Bayesian VI used by RUN_VI (bin/variational_inference_bayesian.py).

Covers every public function: the likelihood matrix, the CAVI updates, the
ELBO (checked against an independent calculation), the inner and outer loops,
the posterior utilities and the drop-in wrapper. The tests named after F05,
F06 and F07 were strict xfails until Phase 3 fixed those bugs.
"""
import io
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.special import digamma
from scipy.stats import dirichlet

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
import variational_inference_bayesian as vb
from taxonomy import freq_to_lineage_df

from .synthetic import synthetic_reads


@pytest.fixture(autouse=True)
def quiet(monkeypatch):
    """The module writes progress to sys.stdout bound at import; silence it."""
    monkeypatch.setattr(vb, 'stdout', io.StringIO())


# ─── build_likelihood_matrix ─────────────────────────────────────────────────

class TestBuildLikelihoodMatrix:
    def test_shape_and_order(self):
        lp = {'r1': ([1], [-1.0]), 'r2': ([2, 1], [-2.0, -3.0])}
        L, names, sid_to_col = vb.build_likelihood_matrix(lp, [1, 2, 3])
        assert L.shape == (2, 3)
        assert names == ['r1', 'r2']
        assert sid_to_col == {1: 0, 2: 1, 3: 2}

    def test_values_and_missing_are_neg_inf(self):
        lp = {'r1': ([1], [-1.0]), 'r2': ([2, 1], [-2.0, -3.0])}
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2, 3])
        expected = np.array([[-1.0, -np.inf, -np.inf], [-3.0, -2.0, -np.inf]])
        np.testing.assert_array_equal(L, expected)

    def test_duplicate_species_keeps_best_score(self):
        lp = {'r1': ([1, 1], [-5.0, -2.0])}
        L, _, _ = vb.build_likelihood_matrix(lp, [1])
        assert L[0, 0] == -2.0

    def test_species_not_in_list_ignored(self):
        lp = {'r1': ([7], [-1.0])}
        L, _, _ = vb.build_likelihood_matrix(lp, [1])
        assert np.isneginf(L[0, 0])


# ─── CAVI updates ────────────────────────────────────────────────────────────

class TestComputeDirichletExpectation:
    def test_matches_digamma_formula(self):
        alpha = np.array([1.0, 2.5, 10.0])
        expected = digamma(alpha) - digamma(alpha.sum())
        np.testing.assert_allclose(vb.compute_dirichlet_expectation(alpha), expected)

    def test_matches_monte_carlo(self):
        alpha = np.array([2.0, 5.0, 3.0])
        rng = np.random.default_rng(1)
        draws = rng.dirichlet(alpha, size=200_000)
        np.testing.assert_allclose(
            vb.compute_dirichlet_expectation(alpha), np.log(draws).mean(axis=0), atol=0.01
        )


class TestUpdateResponsibilities:
    def test_softmax_by_hand(self):
        L = np.array([[math.log(0.6), math.log(0.2)]])
        E_ln_pi = np.array([math.log(0.5), math.log(0.5)])
        phi = vb.update_responsibilities(L, E_ln_pi)
        np.testing.assert_allclose(phi, [[0.75, 0.25]])

    def test_rows_sum_to_one_and_unaligned_row_is_zero(self):
        L = np.array([[-1.0, -2.0, -np.inf], [-np.inf, -np.inf, -np.inf], [-3.0, -3.0, -3.0]])
        phi = vb.update_responsibilities(L, np.zeros(3))
        np.testing.assert_allclose(phi.sum(axis=1), [1.0, 0.0, 1.0])
        assert phi[0, 2] == 0.0

    def test_stable_for_very_negative_scores(self):
        L = np.array([[-5000.0, -5001.0]])
        phi = vb.update_responsibilities(L, np.zeros(2))
        np.testing.assert_allclose(phi, [[1 / (1 + math.exp(-1)), math.exp(-1) / (1 + math.exp(-1))]])


class TestUpdateDirichletParameters:
    def test_prior_plus_expected_counts(self):
        alpha_0 = np.array([0.5, 0.5])
        phi = np.array([[0.75, 0.25], [1.0, 0.0], [0.0, 0.0]])
        np.testing.assert_allclose(vb.update_dirichlet_parameters(alpha_0, phi), [2.25, 0.75])


# ─── ELBO ────────────────────────────────────────────────────────────────────

def elbo_by_hand(L, phi, alpha, alpha_0):
    """ELBO written term by term, with scipy's Dirichlet entropy for H[q(pi)]."""
    E = digamma(alpha) - digamma(alpha.sum())
    total = 0.0
    for r in range(L.shape[0]):
        for s in range(L.shape[1]):
            if phi[r, s] > 0:
                # E[ln p(x_r | z_r = s)] + E[ln p(z_r = s | pi)] - E[ln q(z_r = s)]
                total += phi[r, s] * (L[r, s] + E[s] - math.log(phi[r, s]))
    log_B_alpha_0 = sum(math.lgamma(a) for a in alpha_0) - math.lgamma(alpha_0.sum())
    E_log_prior = -log_B_alpha_0 + float(np.sum((alpha_0 - 1) * E))
    return total + E_log_prior + dirichlet(alpha).entropy()


class TestComputeElbo:
    def test_matches_hand_calculation(self):
        L = np.array([[-1.0, -2.0], [-3.0, -np.inf], [-0.5, -0.7]])
        alpha_0 = np.array([0.5, 0.5])
        alpha = np.array([2.0, 1.5])
        E = vb.compute_dirichlet_expectation(alpha)
        phi = vb.update_responsibilities(L, E)
        alpha = vb.update_dirichlet_parameters(alpha_0, phi)
        E = vb.compute_dirichlet_expectation(alpha)
        assert vb.compute_elbo(L, phi, alpha, alpha_0, E) == pytest.approx(
            elbo_by_hand(L, phi, alpha, alpha_0), rel=1e-12
        )

    def test_elbo_never_decreases(self):
        lp = synthetic_reads({1: 300, 2: 150, 3: 40, 4: 10}, seed=3)
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2, 3, 4])
        alpha_0 = np.full(4, 0.25)
        alpha = alpha_0 + L.shape[0] / 4
        elbos = []
        for _ in range(50):
            phi = vb.update_responsibilities(L, vb.compute_dirichlet_expectation(alpha))
            alpha = vb.update_dirichlet_parameters(alpha_0, phi)
            elbos.append(vb.compute_elbo(L, phi, alpha, alpha_0, vb.compute_dirichlet_expectation(alpha)))
        diffs = np.diff(elbos)
        assert np.all(diffs >= -1e-9 * abs(elbos[-1]))


# ─── Inner loop ──────────────────────────────────────────────────────────────

class TestCaviInnerLoop:
    def test_returns_converged_state(self):
        lp = synthetic_reads({1: 200, 2: 100}, seed=4)
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2])
        alpha_0 = np.full(2, 0.5)
        alpha, phi, E_ln_pi, elbo, n_iters = vb.cavi_inner_loop(
            L, alpha_0, alpha_0 + 150, max_iterations=100, tolerance=1e-10, verbose=False
        )
        assert alpha.shape == (2,)
        assert phi.shape == L.shape
        np.testing.assert_allclose(alpha, alpha_0 + phi.sum(axis=0))
        np.testing.assert_allclose(E_ln_pi, vb.compute_dirichlet_expectation(alpha))
        assert np.isfinite(elbo)
        assert 1 <= n_iters <= 100

    def test_stops_at_max_iterations(self):
        lp = synthetic_reads({1: 200, 2: 100}, seed=4)
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2])
        alpha_0 = np.full(2, 0.5)
        *_, n_iters = vb.cavi_inner_loop(L, alpha_0, alpha_0 + 150, max_iterations=3, tolerance=0.0, verbose=False)
        assert n_iters == 3


    def test_stops_on_change_in_posterior_mean(self):
        lp = synthetic_reads({1: 200, 2: 100, 3: 30}, seed=4)
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2, 3])
        alpha_0 = np.full(3, 1 / 3)
        alpha, *_, n_iters = vb.cavi_inner_loop(L, alpha_0, alpha_0 + 110, tolerance=1e-7, verbose=False)
        assert n_iters < 1000
        # one more update moves the posterior mean by less than the tolerance
        phi = vb.update_responsibilities(L, vb.compute_dirichlet_expectation(alpha))
        nxt = vb.update_dirichlet_parameters(alpha_0, phi)
        assert np.max(np.abs(nxt / nxt.sum() - alpha / alpha.sum())) < 1e-7

    def test_converged_start_stops_after_one_iteration(self):
        lp = synthetic_reads({1: 200, 2: 100}, seed=4)
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2])
        alpha_0 = np.full(2, 0.5)
        alpha, *_ = vb.cavi_inner_loop(L, alpha_0, alpha_0 + 150, tolerance=1e-12, verbose=False)
        *_, n_iters = vb.cavi_inner_loop(L, alpha_0, alpha, tolerance=1e-7, verbose=False)
        assert n_iters == 1

    def test_warns_when_elbo_falls(self, monkeypatch, capsys):
        lp = synthetic_reads({1: 200, 2: 100}, seed=4)
        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2])
        alpha_0 = np.full(2, 0.5)
        elbos = iter([-1000.0, -1000.0 - 1e-7, -1001.0])
        monkeypatch.setattr(vb, 'compute_elbo', lambda *args: next(elbos))
        vb.cavi_inner_loop(L, alpha_0, alpha_0 + 150, max_iterations=3, tolerance=0.0, verbose=False)
        err = capsys.readouterr().err
        # a fall of 1e-7 (1e-10 of |ELBO|) is rounding; a fall of 1 is reported once
        assert err.count('WARNING: ELBO fell') == 1
        assert 'iteration 3' in err


# ─── Outer loop and wrapper ──────────────────────────────────────────────────

class TestBayesianViIterations:
    def test_recovers_known_proportions(self):
        truth = {1: 0.40, 2: 0.30, 3: 0.15, 4: 0.10, 5: 0.05}
        lp = synthetic_reads({s: int(p * 2000) for s, p in truth.items()}, seed=5)
        freq_full, _, ci, _ = vb.bayesian_vi_iterations(lp, list(truth))
        assert sum(freq_full.values()) == pytest.approx(1.0)
        for sid, p in truth.items():
            assert freq_full[sid] == pytest.approx(p, abs=0.02)
            lower, upper = ci[sid]
            assert lower <= freq_full[sid] <= upper

    def test_prunes_species_without_support(self):
        lp = synthetic_reads({1: 1500, 2: 500}, seed=6)
        lp['stray'] = ([3, 1], [-90.0, -40.0])
        freq_full, _, _, _ = vb.bayesian_vi_iterations(lp, [1, 2, 3, 4])
        assert set(freq_full) == {1, 2}

    def test_empty_input_raises(self):
        with pytest.raises(ValueError, match='0 reads'):
            vb.bayesian_vi_iterations({}, [1, 2])

    def test_species_below_input_threshold_not_reported(self):
        lp = synthetic_reads({1: 1000, 2: 4}, seed=7, max_decoys=0)
        freq_full, freq_thresh, _, _ = vb.bayesian_vi_iterations(lp, [1, 2], abundance_threshold=0.01)
        assert set(freq_full) == {1}
        assert freq_thresh is None or set(freq_thresh) == {1}

    def test_result_is_a_fixed_point(self):
        # F05: 20,000 reads, two species, each read favours its source by 2 nats.
        # The ELBO is about -2e6, so the old relative tolerance of 5e-4 let the
        # loop stop after two iterations, far from the fixed point.
        rng = np.random.default_rng(8)
        z = rng.choice(2, size=20_000, p=[0.7, 0.3])
        lp = {
            f'r{i}': ([1, 2], [-100.0 if zi == 0 else -102.0, -102.0 if zi == 0 else -100.0])
            for i, zi in enumerate(z)
        }
        freq_full, _, _ = vb.variational_inference_iterations(lp, [1, 2])

        L, _, _ = vb.build_likelihood_matrix(lp, [1, 2])
        alpha_0 = np.full(2, 0.5)
        alpha = alpha_0 + L.shape[0] / 2
        for _ in range(5000):
            phi = vb.update_responsibilities(L, vb.compute_dirichlet_expectation(alpha))
            new_alpha = vb.update_dirichlet_parameters(alpha_0, phi)
            if np.max(np.abs(new_alpha / new_alpha.sum() - alpha / alpha.sum())) < 1e-12:
                break
            alpha = new_alpha
        fixed_point = alpha / alpha.sum()
        assert freq_full[1] == pytest.approx(fixed_point[0], abs=1e-4)

    def test_outer_loop_prunes_before_stopping(self, monkeypatch):
        # F07: the inner loop is scripted to isolate the outer loop's control flow.
        # After the first pruning, species 2 falls below the threshold and the
        # ELBO is lower (it is a different model); the loop must still prune it.
        # The fake takes the arguments the outer loop passes to cavi_inner_loop
        # (max_iterations and tolerance by keyword) and returns its 5-tuple.
        # With 2,000 reads the pruning threshold is 10 / 2000 = 0.005; 8 reads is 0.004.
        scripted = {3: np.array([1800.0, 192.0, 8.0]), 2: np.array([1992.0, 8.0]), 1: np.array([2000.0])}
        elbos = iter([-100.0, -200.0, -300.0, -400.0])

        def fake_inner_loop(L, alpha_0, alpha_init, max_iterations, tolerance):
            alpha = alpha_0 + scripted[L.shape[1]]
            E = vb.compute_dirichlet_expectation(alpha)
            return alpha, vb.update_responsibilities(L, E), E, next(elbos), 2

        monkeypatch.setattr(vb, 'cavi_inner_loop', fake_inner_loop)
        lp = {f'r{i}': ([1, 2, 3], [-1.0, -1.0, -1.0]) for i in range(2000)}
        freq_full, _, _, _ = vb.bayesian_vi_iterations(lp, [1, 2, 3])
        assert set(freq_full) == {1}

    def test_outer_loop_warm_starts(self, monkeypatch):
        calls = []
        real_inner_loop = vb.cavi_inner_loop

        def spy(L, alpha_0, alpha_init, **kwargs):
            calls.append(np.array(alpha_init))
            return real_inner_loop(L, alpha_0, alpha_init, **kwargs)

        monkeypatch.setattr(vb, 'cavi_inner_loop', spy)
        lp = synthetic_reads({1: 1500, 2: 450, 3: 5}, seed=9, max_decoys=1)
        vb.bayesian_vi_iterations(lp, [1, 2, 3])
        assert len(calls) >= 2
        second = calls[1]
        assert not np.allclose(second, second[0])
        # the surviving species start from their expected counts
        np.testing.assert_allclose(second.sum(), calls[0].sum() - 5, atol=5)

    def test_outer_loop_stops_after_max_iterations_with_a_final_fit(self, monkeypatch, capsys):
        # Every fit leaves its last species at 8 of 2,000 reads, below 10 / 2000,
        # so each pass prunes one. With max_outer_iterations=2 the loop prunes
        # twice, fits the pruned set once more and warns that one is still below.
        sizes = []

        def fake_inner_loop(L, alpha_0, alpha_init, max_iterations, tolerance):
            S = L.shape[1]
            sizes.append(S)
            alpha = alpha_0 + np.array([2000.0 - 8 - 40 * (S - 2)] + [40.0] * (S - 2) + [8.0])
            E = vb.compute_dirichlet_expectation(alpha)
            return alpha, vb.update_responsibilities(L, E), E, -1.0, 1

        monkeypatch.setattr(vb, 'cavi_inner_loop', fake_inner_loop)
        lp = {f'r{i}': ([1, 2, 3, 4, 5], [-1.0] * 5) for i in range(2000)}
        freq_full, _, _, _ = vb.bayesian_vi_iterations(lp, [1, 2, 3, 4, 5], max_outer_iterations=2)
        assert sizes == [5, 4, 3]
        assert set(freq_full) == {1, 2, 3}
        assert '1 species remain below the pruning threshold' in capsys.readouterr().err

    def test_all_species_below_threshold_keeps_the_most_abundant(self):
        # One read split evenly over three species: each has 1/3 < 1/n = 1
        lp = {'r1': ([1, 2, 3], [-1.0, -1.0, -1.1])}
        freq_full, _, _, _ = vb.bayesian_vi_iterations(lp, [1, 2, 3])
        assert list(freq_full) in ([1], [2])
        assert freq_full[list(freq_full)[0]] == pytest.approx(1.0)

    def test_reads_of_pruned_species_are_not_counted_as_assigned(self, tmp_path):
        # F06: the wrapper returns the reads whose species were all pruned, and
        # freq_to_lineage_df reports them in a 'Below detection' row
        lp = synthetic_reads({1: 1500, 2: 495}, seed=10, max_decoys=1)
        for i in range(5):  # five reads that align only to species 3, which is pruned
            lp[f'orphan{i}'] = ([3], [-40.0])
        freq_full, freq_thresh, below = vb.variational_inference_iterations(lp, [1, 2, 3])
        assert 3 not in freq_full
        assert below == 5
        taxonomy = pd.DataFrame(
            {'species': ['Species one', 'Species two', 'Species three'], 'genus': ['G1', 'G2', 'G3']},
            index=pd.Index(['1', '2', '3'], name='tax_id'),
        )
        freq = {str(k): v for k, v in (freq_thresh or freq_full).items()}
        df = freq_to_lineage_df(
            freq, str(tmp_path / 's'), taxonomy, assigned_count=len(lp), unassigned_count=0, counts=True,
            below_detection_count=below,
        )
        assigned = df[df['species'].isin(['Species one', 'Species two'])]
        assert assigned['estimated_counts'].sum() == len(lp) - 5
        below_row = df[df['species'] == 'Below detection']
        assert below_row['estimated_counts'].tolist() == [5]
        assert list(df['species'][-2:]) == ['Below detection', 'Unassigned']

    def test_no_reads_below_detection_when_every_read_keeps_a_species(self):
        lp = synthetic_reads({1: 1500, 2: 450}, seed=9, max_decoys=1)
        for i in range(5):  # species 3's reads also align to species 1
            lp[f'r3_{i}'] = ([3, 1], [-40.0, -45.0])
        freq_full, _, _, below = vb.bayesian_vi_iterations(lp, [1, 2, 3])
        assert 3 not in freq_full
        assert below == 0

    def test_initial_alpha_counts_only_scored_reads(self, monkeypatch):
        # 'elsewhere' aligns only to a species outside the list, so it has no
        # finite score; the uniform start spreads 100 reads, not 101, and the
        # read is not below detection either (it was never part of the fit)
        calls = []
        real_inner_loop = vb.cavi_inner_loop

        def spy(L, alpha_0, alpha_init, **kwargs):
            calls.append((alpha_0.copy(), alpha_init.copy()))
            return real_inner_loop(L, alpha_0, alpha_init, **kwargs)

        monkeypatch.setattr(vb, 'cavi_inner_loop', spy)
        lp = synthetic_reads({1: 60, 2: 40}, seed=13, max_decoys=1)
        lp['elsewhere'] = ([9], [-40.0])
        *_, below = vb.bayesian_vi_iterations(lp, [1, 2])
        alpha_0, alpha_init = calls[0]
        assert (alpha_init - alpha_0).sum() == pytest.approx(100.0)
        assert below == 0

    def test_reads_of_species_dropped_at_the_cap_are_below_detection(self, monkeypatch):
        # One pruning allowed: species 3 is pruned, and in the final fit species 2
        # is still below the threshold, so the abundance threshold drops it. Its
        # 10 reads of its own must not be spread over species 1.
        scripted = {3: np.array([1980.0, 12.0, 8.0]), 2: np.array([1994.0, 6.0])}

        def fake_inner_loop(L, alpha_0, alpha_init, max_iterations, tolerance):
            alpha = alpha_0 + scripted[L.shape[1]]
            E = vb.compute_dirichlet_expectation(alpha)
            return alpha, vb.update_responsibilities(L, E), E, -1.0, 1

        monkeypatch.setattr(vb, 'cavi_inner_loop', fake_inner_loop)
        lp = {f'r{i}': ([1, 2, 3], [-1.0, -1.0, -1.0]) for i in range(1990)}
        lp.update({f'own2_{i}': ([2], [-1.0]) for i in range(10)})
        freq_full, freq_thresh, _, below = vb.bayesian_vi_iterations(
            lp, [1, 2, 3], abundance_threshold=0.004, max_outer_iterations=1)
        assert set(freq_full) == {1, 2}
        assert set(freq_thresh) == {1}
        assert below == 10

    def test_zero_max_outer_iterations_still_refits_after_pruning(self):
        lp = synthetic_reads({1: 1500, 2: 495}, seed=10, max_decoys=1)
        lp.update({f'orphan{i}': ([3], [-40.0]) for i in range(5)})
        freq_full, _, _, below = vb.bayesian_vi_iterations(lp, [1, 2, 3], max_outer_iterations=0)
        reference, _, _, ref_below = vb.bayesian_vi_iterations(lp, [1, 2, 3])
        assert set(freq_full) == {1, 2}
        assert freq_full == pytest.approx(reference, abs=1e-6)
        assert below == ref_below == 5


class TestVariationalInferenceIterationsWrapper:
    def test_return_shape(self):
        lp = synthetic_reads({1: 800, 2: 200}, seed=11)
        result = vb.variational_inference_iterations(lp, [1, 2])
        assert isinstance(result, tuple) and len(result) == 3
        freq_full, freq_thresh, below = result
        assert below == 0
        assert isinstance(freq_full, dict)
        assert set(freq_full) <= {1, 2}
        assert all(isinstance(v, float) for v in freq_full.values())
        assert sum(freq_full.values()) == pytest.approx(1.0)
        assert freq_thresh is None or isinstance(freq_thresh, dict)

    def test_empty_input_raises(self):
        with pytest.raises(ValueError, match='0 reads'):
            vb.variational_inference_iterations({}, [])


# ─── Posterior utilities ─────────────────────────────────────────────────────

class TestPosteriorUtilities:
    def test_posterior_entropy_matches_scipy(self):
        alpha = np.array([2.0, 3.5, 1.2])
        assert vb.posterior_entropy(alpha) == pytest.approx(dirichlet(alpha).entropy(), rel=1e-10)

    def test_kl_divergence_zero_for_identical(self):
        alpha = np.array([2.0, 3.0])
        assert vb.kl_divergence_dirichlet(alpha, alpha) == pytest.approx(0.0, abs=1e-12)

    def test_kl_divergence_matches_monte_carlo(self):
        a, b = np.array([3.0, 2.0, 5.0]), np.array([1.0, 1.0, 1.0])
        draws = np.random.default_rng(12).dirichlet(a, size=100_000)
        mc = np.mean(dirichlet(a).logpdf(draws.T) - dirichlet(b).logpdf(draws.T))
        assert vb.kl_divergence_dirichlet(a, b) == pytest.approx(mc, abs=0.02)

    def test_expected_effective_species(self):
        assert vb.expected_effective_species(np.full(4, 10.0)) == pytest.approx(4.0)
        assert vb.expected_effective_species(np.array([1000.0, 1e-9])) == pytest.approx(1.0, abs=1e-6)
