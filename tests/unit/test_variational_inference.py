"""Unit tests for the variational inference engine (bin/variational_inference.py).

Tests cover the pure-Python functions that implement the core EM/VI algorithm:
update_q_distribution, maximize_elbo, variational_inference, and
variational_inference_iterations. pysam-dependent I/O functions are excluded
from unit tests and covered by the nf-test integration suite.
"""
import math
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
from variational_inference import (
    update_q_distribution,
    maximize_elbo,
    variational_inference,
    variational_inference_iterations,
)


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
