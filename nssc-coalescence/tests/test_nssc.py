"""Consistency tests. Run from the repository root: python -m pytest tests"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nssc.config import group_t
from nssc.models import NSSC, Pnisland, create_glue_matrix, create_q_matrix
from nssc.ms import empirical_distribution
from nssc.stats import clean_density, coalescence_stats
from nssc.theory import model_iicr, stationary_nisland_iicr


def island_matrix(n, rate):
    return rate * (np.ones((n, n)) - np.eye(n))


def test_q_matrix_rows_sum_to_zero():
    for dist in (False, True):
        Q = create_q_matrix(island_matrix(3, 1.0), [1, 2, 0.5], dist)
        assert np.allclose(Q.sum(axis=1), 0)


@pytest.mark.parametrize("n_before, n_after", [(4, 2), (3, 2), (2, 3), (5, 1)])
def test_glue_matrix_rows_sum_to_one(n_before, n_after):
    G = create_glue_matrix(island_matrix(n_before, 1.0), n_after)
    assert np.allclose(G.sum(axis=1), 1)


def test_glue_two_vanishing_demes():
    # 4 -> 2 demes: the lineages in demes 2 and 3 each move to 0 or 1 (probability 1/2)
    G = create_glue_matrix(island_matrix(4, 1.0), 2)
    row_23 = 8  # state (2, 3)
    assert np.allclose(G[row_23], [0.25, 0.5, 0.25, 0])


@pytest.mark.parametrize("same_island", [True, False])
def test_pnisland_matches_exact_formula(same_island):
    t = np.linspace(0.01, 5, 200)
    model = Pnisland(same_island, [0], [10], [1.0], [1])
    exact = stationary_nisland_iicr(t, 10, 1.0, same_island)
    assert np.allclose(model_iicr(model, t), exact, rtol=1e-8)


def test_nssc_matches_pnisland():
    # Same n-island model, described by a migration matrix
    t = np.linspace(0.01, 10, 500)
    n, M = 5, 2.0
    nssc = NSSC([2, 0, 0, 0, 0], [0, 1], [island_matrix(n, M / (n - 1))] * 2, [[1] * n, [2] * n])
    pn = Pnisland(True, [0, 1], [n, n], [M, M], [1, 2])
    assert np.allclose(nssc.pdf(t), pn.pdf(t))
    assert np.allclose(nssc.cdf(t), pn.cdf(t))


def test_distinguishable_lineages_give_same_law():
    t = np.linspace(0.01, 10, 500)
    args = ([1, 1, 0], [0], [[[0, 1, 2], [1, 0, 0.5], [2, 0.5, 0]]], [[1, 2, 1]])
    assert np.allclose(NSSC(*args).pdf(t), NSSC(*args, lineages_are_dist=True).pdf(t))


def test_two_island_mean():
    # Two lineages in the same island: E[T2] = number of islands
    t = np.linspace(0, 100, 100_000)
    model = NSSC([2, 0], [0], [island_matrix(2, 1.0)], [[1, 1]])
    assert coalescence_stats(model, t)["mean"] == pytest.approx(2, rel=1e-4)


def test_clean_density():
    f = clean_density([np.nan, 0.5, -1, np.inf, 0.2])
    assert np.allclose(f, [1e-14, 0.5, 0.5, 0.5, 0.2])


def test_group_t():
    assert np.allclose(group_t(np.arange(10.0), "2*1+3+2"), [0, 1, 2, 5, 7])


def test_empirical_distribution():
    cdf, pdf = empirical_distribution([0.5, 1.5, 1.6, 2.5], [0, 1, 2, 3])
    assert np.allclose(cdf, [0, 1, 3, 4])
    assert np.allclose(pdf, [0, 1, 2, 1])  # cells [0, .5), [.5, 1.5), [1.5, 2.5), [2.5, 3.5]
