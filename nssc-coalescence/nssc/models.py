"""
Structured coalescent models for a pair of lineages.

Each model is a continuous-time Markov chain that is piecewise constant:
on each period [times[k], times[k+1]) the chain has a rate matrix Q_k. The
last state is absorbing (the two lineages have coalesced), so the
coalescence time T2 is the absorption time.

- `PiecewiseCoalescent`: shared machinery (density and CDF of T2);
- `NSSC`: non-stationary structured coalescent, with any migration matrix
  and a number of demes that may change from one period to the next;
- `Pnisland`: piecewise n-island model (only 3 states).
"""
import numpy as np


# --------------------------------------------------------------------------
# State numbering
# --------------------------------------------------------------------------

def pair_indices(n_demes, lineages_are_dist=False):
    """
    Number the states "lineage 1 in deme i, lineage 2 in deme j".
    Indistinguishable lineages: pairs (i, j) with i <= j.
    Distinguishable lineages:   all pairs (i, j).
    """
    if lineages_are_dist:
        pairs = [(i, j) for i in range(n_demes) for j in range(n_demes)]
    else:
        pairs = [(i, j) for i in range(n_demes) for j in range(i, n_demes)]
    return {pair: k for k, pair in enumerate(pairs)}


def ordered(a, b):
    return (a, b) if a <= b else (b, a)


def initial_state_index(sampling_vector, lineages_are_dist=False):
    """Index of the initial state: 2 lineages in one deme, or 1 in each of two demes."""
    if 2 in sampling_vector:
        a = b = sampling_vector.index(2)
    else:
        a, b = [i for i, x in enumerate(sampling_vector) if x == 1][:2]
    return pair_indices(len(sampling_vector), lineages_are_dist)[(a, b)]


# --------------------------------------------------------------------------
# Matrices
# --------------------------------------------------------------------------

def create_q_matrix(mig_matrix, deme_sizes, lineages_are_dist=False):
    """
    Rate matrix Q (the last row/column is the absorbing state).
    Each lineage in deme i migrates to deme l at rate M[i, l] / 2;
    two lineages in deme i coalesce at rate 1 / c[i].
    """
    mig = 0.5 * np.asarray(mig_matrix, dtype=float)
    n = mig.shape[0]
    idx = pair_indices(n, lineages_are_dist)
    key = (lambda a, b: (a, b)) if lineages_are_dist else ordered
    Q = np.zeros((len(idx) + 1, len(idx) + 1))
    for (i, j), row in idx.items():
        for l in range(n):
            if l != i:  # the lineage in deme i migrates to l
                Q[row, idx[key(l, j)]] += mig[i, l]
            if l != j:  # the lineage in deme j migrates to l
                Q[row, idx[key(i, l)]] += mig[j, l]
        if i == j:      # coalescence in deme i
            Q[row, -1] = 1.0 / deme_sizes[i]
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def destination_probs(mig_matrix, deme, n_kept):
    """
    Probabilities that a lineage in deme `deme` (which vanishes) is sent to
    each of the remaining demes (indices < n_kept), proportionally to the
    migration rates.
    """
    rates = np.asarray(mig_matrix, dtype=float)[deme, :n_kept]
    neighbours = np.flatnonzero(rates > 0)
    total = rates[neighbours].sum()
    return {int(k): rates[k] / total for k in neighbours}


def create_glue_matrix(mig_before, n_demes_after):
    """
    Transition ("glue") matrix between the states of one period and those of
    the next one when the number of demes changes (indistinguishable lineages).

    Demes with index >= n_demes_after vanish: each of their lineages moves,
    independently, to a remaining deme with probability proportional to the
    migration rates of the period that ends. If the number of demes
    increases, lineages stay where they are.
    """
    idx_before = pair_indices(len(mig_before))
    idx_after = pair_indices(n_demes_after)
    G = np.zeros((len(idx_before) + 1, len(idx_after) + 1))
    G[-1, -1] = 1

    def probs(deme):
        return destination_probs(mig_before, deme, n_demes_after)

    for (a, b), row in idx_before.items():
        if b < n_demes_after:                 # both demes remain
            G[row, idx_after[(a, b)]] = 1
        elif a < n_demes_after:               # only the lineage in deme b has to move
            for l, p in probs(b).items():
                G[row, idx_after[ordered(a, l)]] += p
        else:                                 # both lineages have to move
            p_a, p_b = probs(a), probs(b)
            for u, pu in p_a.items():
                for v, pv in p_b.items():
                    G[row, idx_after[ordered(u, v)]] += pu * pv
    return G


# --------------------------------------------------------------------------
# Shared machinery
# --------------------------------------------------------------------------

class PiecewiseCoalescent:
    """
    Piecewise-constant Markov chain.

    times          : start time of each period (the first one is 0)
    q_matrices     : one rate matrix per period
    initial_state  : index of the initial state
    glue_matrices  : one transition matrix per period change; it is only
                     used when the number of states changes.
    """

    def __init__(self, times, q_matrices, initial_state, glue_matrices=None):
        if len(times) != len(q_matrices):
            raise ValueError("There must be as many times as Q matrices.")
        if times[0] != 0:
            raise ValueError("The first period must start at t = 0.")
        self.times = np.asarray(times, dtype=float)
        self.q_matrices = [np.asarray(Q, dtype=float) for Q in q_matrices]
        glue_matrices = glue_matrices or [None] * (len(times) - 1)

        # For each period k, the distribution of T2 reads, with dt = t - times[k]:
        #   P(T2 <= t) = Re( sum_m cdf_w[m] * exp(eigvals[m] * dt) )
        #   f(t)       = Re( sum_m pdf_w[m] * exp(eigvals[m] * dt) )
        self._periods = []
        state = np.zeros(len(self.q_matrices[0]))  # state distribution at the start of the period
        state[initial_state] = 1
        for k, Q in enumerate(self.q_matrices):
            if k > 0:
                dt = self.times[k] - self.times[k - 1]
                state = ((state @ A) * np.exp(eigvals * dt)) @ A_inv
                if len(state) != len(Q):
                    G = glue_matrices[k - 1]
                    if G is None:
                        raise ValueError(f"Missing glue matrix at time {times[k]}.")
                    state = state @ np.asarray(G, dtype=float)
            eigvals, A = np.linalg.eig(Q)
            A_inv = np.linalg.inv(A)
            cdf_w = (state @ A) * A_inv[:, -1]
            self._periods.append((eigvals, cdf_w, eigvals * cdf_w))

    def _evaluate(self, t, which):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        period = np.searchsorted(self.times, t, side="right") - 1
        out = np.empty_like(t)
        for k, (eigvals, *weights) in enumerate(self._periods):
            mask = period == k
            dt = t[mask] - self.times[k]
            acc = np.zeros(dt.shape, dtype=complex)
            for lam, w in zip(eigvals, weights[which]):
                acc += w * np.exp(lam * dt)
            out[mask] = acc.real
        return out

    def cdf(self, t):
        """Cumulative distribution function of T2 (t in units of 2N0 generations)."""
        return self._evaluate(t, 0)

    def pdf(self, t):
        """Density of T2 (t in units of 2N0 generations)."""
        return self._evaluate(t, 1)

    def plateau(self):
        """Limit of the IICR: -1 / (largest non-zero eigenvalue of the last Q)."""
        eigvals = np.linalg.eigvals(self.q_matrices[-1]).real
        return -1.0 / eigvals[eigvals < -1e-12].max()


# --------------------------------------------------------------------------
# Models
# --------------------------------------------------------------------------

class NSSC(PiecewiseCoalescent):
    """
    Non-stationary structured coalescent.

    sampling_vector   : number of lineages sampled in each deme (sums to 2)
    times             : start time of each period (the first one is 0)
    mig_matrices      : one migration matrix per period
    deme_sizes        : one vector of deme sizes per period
    glue_matrices     : transition matrices between periods; those that are
                        None (or all of them, if None) are computed automatically
    lineages_are_dist : distinguishable lineages (ordered states (i, j))
    """

    def __init__(self, sampling_vector, times, mig_matrices, deme_sizes,
                 glue_matrices=None, lineages_are_dist=False):
        q_matrices = [create_q_matrix(M, c, lineages_are_dist)
                      for M, c in zip(mig_matrices, deme_sizes)]
        glue_matrices = list(glue_matrices or [None] * (len(times) - 1))
        for k, G in enumerate(glue_matrices):
            n_before, n_after = len(mig_matrices[k]), len(mig_matrices[k + 1])
            if G is None and n_before != n_after:
                if lineages_are_dist:
                    raise NotImplementedError(
                        "No automatic glue matrix for distinguishable lineages: "
                        "please provide the glue matrices.")
                glue_matrices[k] = create_glue_matrix(mig_matrices[k], n_after)
        super().__init__(times, q_matrices,
                         initial_state_index(sampling_vector, lineages_are_dist),
                         glue_matrices)


class Pnisland(PiecewiseCoalescent):
    """
    Piecewise n-island model: all islands have the same size and exchange
    migrants at the same rate. Three states: same island, different
    islands, coalescence.

    same_island : both lineages are sampled in the same island
    n_demes, mig_rates, deme_sizes : one value per period
    """

    def __init__(self, same_island, times, n_demes, mig_rates, deme_sizes):
        q_matrices = [self.create_q_matrix(n, M, 1.0 / c)
                      for n, M, c in zip(n_demes, mig_rates, deme_sizes)]
        super().__init__(times, q_matrices, 0 if same_island else 1)

    @staticmethod
    def create_q_matrix(n, M, c):
        return np.array([[-M - c, M, c],
                         [M / (n - 1), -M / (n - 1), 0],
                         [0, 0, 0]], dtype=float)
