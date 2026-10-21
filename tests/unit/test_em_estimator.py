"""Unit tests for the Emu-style EM estimator (bin/em_estimator.py, --estimator em).

Tests cover update_q_distribution, maximize_elbo, variational_inference and
variational_inference_iterations, and compare the EM and Bayesian VI estimates
on synthetic reads. The SAM-reading functions are tested in test_likelihood_io.py.
"""
import math
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
from em_estimator import (
    update_q_distribution,
    maximize_elbo,
    variational_inference,
    variational_inference_iterations,
)
import variational_inference_bayesian as vb

from .synthetic import synthetic_reads


class TestUpdateQDistribution:
    def test_returns_non_empty_dict(self, simple_log_p_rgs, uniform_freq):
        q_dist, _ = update_q_distribution(simple_log_p_rgs, uniform_freq)
        assert isinstance(q_dist, dict)
        assert len(q_dist) > 0

    def test_per_read_probabilities_sum_to_one(self, ambiguous_log_p_rgs, uniform_freq):
        q_dist, _ = update_q_distribution(ambiguous_log_p_rgs, uniform_freq)
        total = sum(
            q_dist[sp]['read1']
            for sp in q_dist
            if 'read1' in q_dist[sp]
        )
        assert abs(total - 1.0) < 1e-9

    def test_elbo_is_finite(self, simple_log_p_rgs, uniform_freq):
        _, elbo = update_q_distribution(simple_log_p_rgs, uniform_freq)
        assert math.isfinite(elbo)

    def test_species_absent_from_freq_excluded(self, simple_log_p_rgs):
        # Only species 1 in freq → read2 (mapped to species 2) should be excluded
        freq = {1: 1.0}
        q_dist, _ = update_q_distribution(simple_log_p_rgs, freq)
        assert 2 not in q_dist

    def test_zero_frequency_species_excluded(self, ambiguous_log_p_rgs):
        freq = {1: 1.0, 2: 0.0}
        q_dist, _ = update_q_distribution(ambiguous_log_p_rgs, freq)
        assert 2 not in q_dist


class TestMaximizeElbo:
    def test_frequencies_sum_to_one(self, simple_log_p_rgs, uniform_freq):
        q_dist, _ = update_q_distribution(simple_log_p_rgs, uniform_freq)
        freq = maximize_elbo(q_dist)
        assert abs(sum(freq.values()) - 1.0) < 1e-9

    def test_returns_dict(self, simple_log_p_rgs, uniform_freq):
        q_dist, _ = update_q_distribution(simple_log_p_rgs, uniform_freq)
        assert isinstance(maximize_elbo(q_dist), dict)

    def test_equal_assignment_gives_equal_freqs(self):
        q_dist = {1: {'read1': 0.5, 'read2': 0.5}, 2: {'read3': 0.5, 'read4': 0.5}}
        freq = maximize_elbo(q_dist)
        assert abs(freq[1] - freq[2]) < 1e-9

    def test_all_reads_to_one_species(self):
        q_dist = {1: {'read1': 1.0, 'read2': 1.0}}
        freq = maximize_elbo(q_dist)
        assert abs(freq[1] - 1.0) < 1e-9


class TestVariationalInference:
    def test_converges_without_error(self, simple_log_p_rgs, uniform_freq):
        freq, elbo = variational_inference(simple_log_p_rgs, uniform_freq)
        assert isinstance(freq, dict)
        assert math.isfinite(elbo)

    def test_output_frequencies_sum_to_one(self, simple_log_p_rgs, uniform_freq):
        freq, _ = variational_inference(simple_log_p_rgs, uniform_freq)
        assert abs(sum(freq.values()) - 1.0) < 1e-9

    def test_dominant_species_gets_higher_frequency(self, ambiguous_log_p_rgs):
        # read1 maps to species 1 with 9× higher probability → species 1 should dominate
        freq_init = {1: 0.5, 2: 0.5}
        freq, _ = variational_inference(ambiguous_log_p_rgs, freq_init)
        assert freq[1] > freq[2]

    def test_deterministic_with_fixed_input(self, simple_log_p_rgs, uniform_freq):
        freq_a, elbo_a = variational_inference(simple_log_p_rgs, uniform_freq)
        freq_b, elbo_b = variational_inference(simple_log_p_rgs, dict(uniform_freq))
        assert freq_a == freq_b
        assert elbo_a == elbo_b


class TestVariationalInferenceIterations:
    def test_raises_on_empty_read_dict(self):
        with pytest.raises(ValueError, match="0 reads"):
            variational_inference_iterations(
                {}, db_ids=[1, 2], lli_thresh=1e-3, input_threshold=0.01
            )

    def test_returns_three_values(self, simple_log_p_rgs):
        result = variational_inference_iterations(
            simple_log_p_rgs, db_ids=[1, 2], lli_thresh=1e-3, input_threshold=0.001
        )
        assert len(result) == 3

    def test_freq_full_sums_to_one(self, simple_log_p_rgs):
        freq_full, _, _ = variational_inference_iterations(
            simple_log_p_rgs, db_ids=[1, 2], lli_thresh=1e-3, input_threshold=0.001
        )
        assert abs(sum(freq_full.values()) - 1.0) < 1e-6

    def test_threshold_removes_rare_species(self):
        # 200 reads for species 1, 1 read for species 2 → species 2 below threshold
        log_p = {f'read{i}': ([1], [math.log(0.95)]) for i in range(200)}
        log_p['read_rare'] = ([2], [math.log(0.05)])
        freq_full, freq_thresh, _ = variational_inference_iterations(
            log_p, db_ids=[1, 2], lli_thresh=1e-3, input_threshold=0.5
        )
        assert 1 in freq_full
        if freq_thresh is not None:
            assert 1 in freq_thresh

    def test_reads_of_pruned_species_are_below_detection(self):
        # 2,000 reads of species 1 and 5 reads that align only to species 2:
        # species 2 (5/2005) falls below the pruning threshold 10/2005
        log_p = {f'read{i}': ([1], [math.log(0.9)]) for i in range(2000)}
        log_p.update({f'rare{i}': ([2], [math.log(0.9)]) for i in range(5)})
        freq_full, freq_thresh, below = variational_inference_iterations(
            log_p, db_ids=[1, 2], lli_thresh=1e-4, input_threshold=0.001
        )
        assert 2 not in (freq_thresh or freq_full)
        assert below == 5

    def test_no_reads_below_detection_when_every_read_keeps_a_species(self, simple_log_p_rgs):
        _, _, below = variational_inference_iterations(
            simple_log_p_rgs, db_ids=[1, 2], lli_thresh=1e-3, input_threshold=0.001
        )
        assert below == 0


class TestEmAgainstBayesianVi:
    """The two estimators differ only for species with a handful of expected
    reads (M1); on well-supported species they agree."""

    TRUTH = {1: 0.40, 2: 0.30, 3: 0.15, 4: 0.10, 5: 0.05}

    def test_estimates_agree_on_synthetic_reads(self, monkeypatch, capsys):
        monkeypatch.setattr(vb, 'stdout', sys.stdout)
        lp = synthetic_reads({s: int(p * 2000) for s, p in self.TRUTH.items()}, seed=5)
        em_full, em_thresh, em_below = variational_inference_iterations(
            lp, list(self.TRUTH), lli_thresh=1e-4, input_threshold=0.001)
        vb_full, vb_thresh, _, vb_below = vb.bayesian_vi_iterations(lp, list(self.TRUTH))
        em_freq = em_thresh or em_full
        vb_freq = vb_thresh or vb_full
        assert set(em_freq) == set(vb_freq) == set(self.TRUTH)
        for sid, p in self.TRUTH.items():
            assert em_freq[sid] == pytest.approx(p, abs=0.02)
            assert em_freq[sid] == pytest.approx(vb_freq[sid], abs=0.001)
        assert em_below == vb_below == 0
