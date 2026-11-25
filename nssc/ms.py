"""
ms simulations: simulated T2 values, empirical distribution, empirical IICR
and density.

"scenarios" entry of the JSON file:
    ms_command : ms command, e.g. "ms 2 100000 -T -L -I 2 2 0 1"
                 (run from the "path2ms" folder)
"""
import os
import subprocess

import numpy as np
from scipy.special import comb


def run_ms(command):
    """Run ms and return the times read on the "time:" lines (ms units of 4N0)."""
    output = subprocess.run(command, shell=True, capture_output=True, text=True, check=True).stdout
    return [float(line.split("\t")[1]) for line in output.splitlines()
            if line.startswith("time:\t")]


def sample_size(ms_command):
    """Number of sampled sequences: the first argument of ms."""
    return int(ms_command.split()[1])


def simulate_t2(scenario, path2ms=""):
    """Simulated T2 values, in units of 2N0 generations (ms uses 4N0)."""
    command = os.path.join(path2ms, scenario["ms_command"])
    return 2 * np.array(run_ms(command))


def empirical_distribution(obs, x_vector):
    """
    Empirical distribution of the observations, evaluated at the points of
    x_vector (the first point is treated as 0). Returns (cdf, pdf) as counts:
      cdf[k] = number of observations <= x_vector[k]
      pdf[k] = number of observations in the cell centred on x_vector[k],
               divided by the width of that cell.
    Observations beyond the last grid point are ignored.
    """
    x = np.array(x_vector, dtype=float)
    x[0] = 0
    half_dx = np.diff(x) / 2
    centred_edges = np.concatenate([[0], x[:-1] + half_dx, [x[-1] + half_dx[-1]]])

    counts = np.histogram(obs, bins=x)[0]
    cdf = np.concatenate([[0], counts.cumsum()])
    pdf = np.histogram(obs, bins=centred_edges)[0] / np.diff(centred_edges)
    return cdf, pdf


def empirical_iicr(obs, x_vector, n_samples=2):
    """
    Empirical IICR (N - F) / f at the points of x_vector. If more than 2
    sequences are sampled, it is multiplied by C(n_samples, 2) to bring it
    back to the scale of a pair.
    """
    cdf, pdf = empirical_distribution(obs, x_vector)
    with np.errstate(divide="ignore", invalid="ignore"):
        return comb(n_samples, 2) * (len(obs) - cdf) / pdf


def empirical_density(obs, x_vector):
    """Empirical density of T2 at the points of x_vector (integrates to ~1)."""
    return empirical_distribution(obs, x_vector)[1] / len(obs)


def ms_iicr(scenario, x_vector, path2ms=""):
    """Simulate an ms scenario and return its empirical IICR at the points of x_vector."""
    obs = simulate_t2(scenario, path2ms)
    return empirical_iicr(obs, x_vector, sample_size(scenario["ms_command"]))
