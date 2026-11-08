"""Hand-checkable cases for the propagation model. Run: python -m pytest tests -q"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from propagation import intron_membership, intron_usage, propagate, smooth_baseline, transcript_factors  # noqa: E402
from splice_graph import Transcript  # noqa: E402


def cassette():
    """+ strand cassette exon: inclusion uses donor 100 / acceptor 200 / donor 300 / acceptor 400;
    skipping uses donor 100 / acceptor 400."""
    inc = Transcript("INC", "G", "chr1", "+", [(1, 100), (200, 300), (400, 500)])
    skp = Transcript("SKP", "G", "chr1", "+", [(1, 100), (400, 500)])
    return [inc, skp]


def tracks(window, changes):
    t = {k: np.full(2 * window + 1, 0.9, dtype=float) for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    for (key, offset), (ref, alt) in changes.items():
        t[f"{key}_ref"][window + offset], t[f"{key}_alt"][window + offset] = ref, alt
    return t


def test_sites_follow_strand_convention():
    inc, _ = cassette()
    assert inc.sites() == {(100, "donor"), (200, "acceptor"), (300, "donor"), (400, "acceptor")}
    minus = Transcript("M", "G", "chr1", "-", [(1, 100), (200, 300)])
    assert minus.sites() == {(200, "donor"), (100, "acceptor")}


def test_exon_acceptor_loss_shifts_to_skipping():
    ts, W, v = cassette(), 150, 200
    f = transcript_factors(ts, v, tracks(W, {("acc", 0): (0.9, 0.09)}), W, eps=0.0)
    assert np.allclose(f, [0.1, 1.0])
    pi = np.array([0.5, 0.5])
    new = propagate(pi, f)
    assert np.allclose(new, [0.05 / 0.55, 0.5 / 0.55])
    introns, index, m, comp = intron_membership(ts)
    du = intron_usage(new, m, comp) - intron_usage(pi, m, comp)
    assert du[index[(101, 399)]] > 0 and du[index[(101, 199)]] < 0


def test_baseline_changes_magnitude_not_factor():
    ts, W = cassette(), 150
    f = transcript_factors(ts, 200, tracks(W, {("acc", 0): (0.9, 0.09)}), W, eps=0.0)
    d_low = propagate(np.array([0.1, 0.9]), f) - np.array([0.1, 0.9])
    d_mid = propagate(np.array([0.5, 0.5]), f) - np.array([0.5, 0.5])
    assert abs(d_mid[0]) > abs(d_low[0])  # largest effect when the exon is at intermediate inclusion


def test_sites_outside_window_ignored():
    ts, W = cassette(), 50
    f = transcript_factors(ts, 200, tracks(W, {}), W)
    assert np.allclose(f, 1.0)


def test_smoothing_is_identical_and_normalized():
    p = smooth_baseline({"A": 1.0}, ["A", "B"], alpha=0.001)
    assert np.isclose(p.sum(), 1) and p[1] > 0
    assert np.allclose(smooth_baseline({}, ["A", "B"]), [0.5, 0.5])


def test_component_clusters_are_transitive():
    """A=(101,199) and B=(101,299) share a start; B and C=(151,299) share an end. A and C share
    nothing directly, but LeafCutter would place all three in one cluster."""
    t1 = Transcript("T1", "G", "chr1", "+", [(1, 100), (200, 400)])
    t2 = Transcript("T2", "G", "chr1", "+", [(1, 100), (300, 400)])
    t3 = Transcript("T3", "G", "chr1", "+", [(1, 150), (300, 400)])
    introns, index, m, comp_c = intron_membership([t1, t2, t3], "component")
    _, _, _, comp_d = intron_membership([t1, t2, t3], "direct")
    a, c = index[(101, 199)], index[(151, 299)]
    assert c in comp_c[a] and c not in comp_d[a]
    assert len(comp_c[a]) == 3
