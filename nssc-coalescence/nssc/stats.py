"""
Statistics of the coalescence time T2 of a pair of lineages: mean, variance
and probability of observing at least one difference.
"""
import numpy as np
from scipy import integrate

PDF_FLOOR = 1e-14  # below this, the density is treated as numerical noise


def clean_density(f):
    """
    Replace invalid density values (nan, inf, < PDF_FLOOR) with the last
    valid value (PDF_FLOOR if it is the first one).
    """
    f = np.array(f, dtype=float)
    with np.errstate(invalid="ignore"):
        bad = ~(np.isfinite(f) & (f >= PDF_FLOOR))
    if bad[0]:
        f[0], bad[0] = PDF_FLOOR, False
    last_valid = np.maximum.accumulate(np.where(bad, 0, np.arange(len(f))))
    return f[last_valid]


def coalescence_stats(model, t, theta=1e-4, n_sites=500):
    """
    Mean and variance of T2, and probability that a pair of sequences of
    `n_sites` sites shows at least one difference (theta per site), by
    numerical integration of the model's density on the grid t.
    """
    t = np.asarray(t, dtype=float)
    f = clean_density(model.pdf(t))
    mean = integrate.simpson(t * f, x=t)
    variance = integrate.simpson(t**2 * f, x=t) - mean**2
    p_no_difference = integrate.simpson(np.exp(-n_sites * theta * t) * f, x=t)
    return {"mean": mean,
            "variance": variance,
            "p_at_least_one_diff": 1 - p_no_difference}
