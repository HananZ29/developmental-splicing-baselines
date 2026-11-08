"""Hand-checkable cases for the E15 local-pattern EM (PREREGISTRATION_E15 section 3).

Model: P(read on intron i) is proportional to (A w)_i, where A[i, p] = 1 if intron i is in splicing
pattern p. Junction usage is u_i = (A w)_i / (1^T A w). The weights are identifiable only when
rank(A) equals the number of patterns.

Run: python -m pytest tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "publication_extension" / "scripts"))
import E15_baseline_diagnostics as D  # noqa: E402
import E15_predict as P  # noqa: E402

# Cassette exon: the inclusion pattern uses introns 0 and 1; the skipping pattern uses intron 2.
CASSETTE = [(0, 1), (2,)]
# Three overlapping patterns over four introns (an alternative 3' site plus a skip), full rank.
OVERLAP = [(0, 1), (0, 2), (3,)]


def roundtrip(pats, n, w_true):
    A = P.pattern_matrix(pats, n)
    u = P.usage_from_weights(np.asarray(w_true, float), A)
    return A, u, P.fit_weights(u, A)


def test_pattern_matrix_layout():
    A = P.pattern_matrix(CASSETTE, 3)
    assert A.tolist() == [[1, 0], [1, 0], [0, 1]]


def test_usage_is_a_distribution():
    A = P.pattern_matrix(OVERLAP, 4)
    u = P.usage_from_weights(np.array([0.5, 0.3, 0.2]), A)
    assert u.min() >= 0
    assert u.sum() == pytest.approx(1.0)


def test_cassette_usage_by_hand():
    # Inclusion weight 0.6 contributes 2 reads per molecule and skipping 0.4 contributes 1:
    # u = (0.6, 0.6, 0.4) / 1.6.
    A = P.pattern_matrix(CASSETTE, 3)
    u = P.usage_from_weights(np.array([0.6, 0.4]), A)
    assert np.allclose(u, np.array([0.6, 0.6, 0.4]) / 1.6)


@pytest.mark.parametrize("pats,n,w", [
    (CASSETTE, 3, [0.6, 0.4]),
    (CASSETTE, 3, [0.05, 0.95]),
    (OVERLAP, 4, [0.5, 0.3, 0.2]),
    (OVERLAP, 4, [0.1, 0.1, 0.8]),
])
def test_em_recovers_true_weights(pats, n, w):
    _, _, w_hat = roundtrip(pats, n, w)
    assert np.allclose(w_hat, w, atol=1e-5)


def test_weights_are_a_distribution():
    _, _, w_hat = roundtrip(OVERLAP, 4, [0.2, 0.5, 0.3])
    assert w_hat.min() >= 0
    assert w_hat.sum() == pytest.approx(1.0)


def test_em_matches_noisy_usage_better_than_its_start():
    # With usage that no weight vector reproduces exactly, EM must still return a distribution
    # whose implied usage is close to the observed usage.
    A = P.pattern_matrix(OVERLAP, 4)
    u = np.array([0.40, 0.22, 0.18, 0.20])
    w = P.fit_weights(u, A)
    assert w.sum() == pytest.approx(1.0)
    assert np.abs(P.usage_from_weights(w, A) - u).max() < 0.05


def test_missing_usage_returns_none():
    A = P.pattern_matrix(CASSETTE, 3)
    assert P.fit_weights(None, A) is None
    assert P.fit_weights(np.full(3, np.nan), A) is None


def test_identifiability_check():
    assert D.full_rank(CASSETTE, 3)
    assert D.full_rank(OVERLAP, 4)
    # A single pattern carries no usage information.
    assert not D.full_rank([(0, 1)], 2)
    # Rank deficient: pattern {0, 1} is the sum of patterns {0} and {1}.
    assert not D.full_rank([(0,), (1,), (0, 1)], 2)


def test_diagnostics_em_agrees_with_prediction_em():
    # E15_baseline_diagnostics.em_weights and E15_predict.fit_weights implement the same model.
    pats, n, w = OVERLAP, 4, np.array([0.5, 0.3, 0.2])
    A = P.pattern_matrix(pats, n)
    u = P.usage_from_weights(w, A)
    w_d = np.asarray(D.em_weights(u, pats, n), float)
    assert np.allclose(w_d / w_d.sum(), P.fit_weights(u, A), atol=1e-4)
