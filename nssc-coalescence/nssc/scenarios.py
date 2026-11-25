"""
Building models from the entries of the JSON file.

"theoretical_IICR_general" entry:
    M        : one migration matrix per period
    size     : one vector of deme sizes per period
    tau      : start time of each period (the first one is 0)
    sampling : number of lineages sampled per deme, e.g. [2, 0, 0]
    glue_matrix       (optional): one glue matrix, used at every change in
                                  the number of demes; or a list of matrices,
                                  one per period change
    lineages_are_dist (optional): distinguishable lineages (default: false)
    Without "tau", M and size describe a stationary model (a single period).

"theoretical_IICR_nisland" entry:
    n, M      : number of islands and migration rate (numbers, or lists per period)
    tau, size : start time and island size of each period (if M is a list)
    sampling_same_island : both lineages come from the same island
"""
import numpy as np

from .models import NSSC, Pnisland


def is_array_like(obj):
    return hasattr(obj, "__len__") and not isinstance(obj, str)


def _glue_list(glue, n_transitions):
    if glue is None:
        return None
    if np.ndim(glue[0][0]) == 0:   # a single matrix, reused everywhere
        return [glue] * n_transitions
    return list(glue)              # one matrix per period change


def general_model(entry):
    """NSSC described by a "theoretical_IICR_general" entry."""
    if "tau" in entry:
        tau, M, size = entry["tau"], entry["M"], entry["size"]
    else:
        tau, M, size = [0], [entry["M"]], [entry["size"]]
    return NSSC(entry["sampling"], tau, M, size,
                glue_matrices=_glue_list(entry.get("glue_matrix"), len(tau) - 1),
                lineages_are_dist=entry.get("lineages_are_dist", False))


def nisland_model(entry):
    """Pnisland described by a "theoretical_IICR_nisland" entry."""
    same_island = entry["sampling_same_island"]
    if not is_array_like(entry["M"]):
        return Pnisland(same_island, [0], [entry["n"]], [entry["M"]], [entry.get("size", 1)])
    check_nisland_entry(entry)
    return Pnisland(same_island, entry["tau"], entry["n"], entry["M"], entry["size"])


def check_nisland_entry(entry):
    M, tau, size = entry["M"], entry.get("tau"), entry.get("size")
    if not (is_array_like(tau) and is_array_like(size)):
        raise TypeError("'tau' and 'size' must be lists when 'M' is a list.")
    if not len(M) == len(tau) == len(size):
        raise ValueError("'M', 'tau' and 'size' must have the same length.")
    if tau[0] != 0:
        raise ValueError("The first time in 'tau' must be 0.")
