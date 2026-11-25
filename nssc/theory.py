"""
Theoretical IICR: IICR(t) = (1 - F(t)) / f(t), where F and f are the
cumulative distribution function and the density of T2 (time in units of
2N0 generations).
"""
import numpy as np

from .scenarios import general_model, is_array_like, nisland_model

PDF_FLOOR = 1e-14
PLATEAU_TOL = 1e-3      # below this gap, the IICR is set to its plateau
PLATEAU_MIN_TIME = 10   # no forced plateau before this time


def model_iicr(model, t):
    """
    IICR of a model (NSSC, Pnisland, ...) on the grid t.

    At large times, after the last event, F tends to 1 and f to 0, so the
    ratio becomes numerically unstable. It is then replaced by its limit
    (the plateau) as soon as it is close to it, or as soon as the previous
    values have been set to it.
    """
    t = np.asarray(t, dtype=float)
    F = model.cdf(t)
    f = model.pdf(t)
    if not 0 <= F[0] <= 1:
        F[0] = 1
    if not (np.isfinite(f[0]) and f[0] >= PDF_FLOOR):
        f[0] = PDF_FLOOR
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        iicr = (1 - F) / f

    plateau = model.plateau()
    can_be_plateau = (t > model.times[-1]) & (t > PLATEAU_MIN_TIME)
    close = np.abs(iicr - plateau) < PLATEAU_TOL
    # "already at the plateau" in the sense of np.allclose (default tolerances)
    at_plateau_tol = 1e-8 + 1e-5 * abs(plateau)
    streak = np.zeros(len(t), dtype=int)  # number of consecutive plateau values ending at i
    for i in range(len(t)):
        if i > 0 and can_be_plateau[i] and (close[i] or (i >= 10 and streak[i - 2] >= 9)):
            iicr[i] = plateau
        if abs(iicr[i] - plateau) <= at_plateau_tol:
            streak[i] = streak[i - 1] + 1 if i > 0 else 1
    return iicr


def stationary_nisland_iicr(t, n, M, same_island=True):
    """Exact IICR of the stationary n-island model (closed form of f)."""
    t = np.asarray(t, dtype=float)
    gamma = M / (n - 1)
    delta = (1 + n * gamma) ** 2 - 4 * gamma
    alpha = 0.5 * (1 + n * gamma + np.sqrt(delta))
    beta = 0.5 * (1 + n * gamma - np.sqrt(delta))
    ea, eb = np.exp(-alpha * t), np.exp(-beta * t)
    if same_island:
        numerator = (1 - beta) * ea + (alpha - 1) * eb
        denominator = (alpha - gamma) * ea + (gamma - beta) * eb
    else:
        numerator = beta * ea - alpha * eb
        denominator = gamma * (ea - eb)
    return numerator / denominator


def general_iicr(entry, t):
    """IICR of a "theoretical_IICR_general" entry."""
    return model_iicr(general_model(entry), t)


def nisland_iicr(entry, t):
    """IICR of a "theoretical_IICR_nisland" entry (exact formula if stationary)."""
    if is_array_like(entry["M"]):
        return model_iicr(nisland_model(entry), t)
    return stationary_nisland_iicr(t, entry["n"], entry["M"], entry["sampling_same_island"])


def theoretical_curves(p, t):
    """
    All theoretical IICRs of the parameter file.
    Returns a list of (family, index, JSON entry, IICR) with
    family = "nisland" or "gen".
    """
    curves = []
    for i, entry in enumerate(p.get("theoretical_IICR_nisland", [])):
        curves.append(("nisland", i, entry, nisland_iicr(entry, t)))
    for i, entry in enumerate(p.get("theoretical_IICR_general", [])):
        curves.append(("gen", i, entry, general_iicr(entry, t)))
    return curves


def distance_to_psmc(model, x, y, n_ref):
    """
    Distance (sum of squares) between the model's IICR, scaled by n_ref, and
    a step-wise PSMC curve (x, y). The IICR is evaluated at the middle of
    each step, converted to units of 2 * n_ref generations.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    midpoints = (x[1:] + x[:-1]) * 0.25 / n_ref
    iicr = model_iicr(model, midpoints)
    return np.sum((y[:len(iicr)] - n_ref * iicr) ** 2)
