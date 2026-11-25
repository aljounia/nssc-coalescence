"""
Reading the JSON parameter file and the quantities derived from it
(time scales, evaluation grids, prefix of the output files).

Units: the models work in units of 2N0 generations. On the figures, time is
multiplied by 2 * N0 * generation_time and the IICR by N0.
"""
import json
from pathlib import Path

import numpy as np

N_THEORY_POINTS = 100_000


def load_params(params_file):
    with open(params_file) as f:
        return json.load(f)


def output_prefix(params_file):
    """Path of the parameter file without its extension (base name of the output files)."""
    return str(Path(params_file).with_suffix(""))


def time_scale(p):
    """Number of figure time units (e.g. years) per unit of 2N0 generations."""
    return 2 * p["scale_params"]["N0"] * p["scale_params"]["generation_time"]


def size_scale(p):
    return p["scale_params"]["N0"]


def theory_time_grid(p, n_points=N_THEORY_POINTS):
    """
    Grid (in units of 2N0 generations) on which theoretical quantities are
    evaluated: it covers the x axis of the figure, up to plot_limits[1].
    """
    t_max = p["plot_params"]["plot_limits"][1]
    if p["plot_params"]["time_log"]:
        grid = np.logspace(-2, np.log10(t_max), n_points)
    else:
        grid = np.linspace(0, t_max, n_points)
    return grid / time_scale(p)


def ms_time_grid(p):
    """Grid (in units of 2N0 generations) of the histograms of ms simulations."""
    custom = p.get("custom_x_vector", {})
    if custom.get("set_custom_xvector"):
        return np.array(custom["x_vector"], dtype=float)
    cp = p["computation_parameters"]
    t = compute_t_vector(cp["start"], cp["end"], cp["number_of_values"], cp["x_vector_type"])
    return group_t(t, cp["pattern"])


def compute_t_vector(start, end, number_of_values, vector_type):
    """Linear grid, or PSMC-like "log" grid (number_of_values + 1 points)."""
    if vector_type == "log":
        i = np.arange(number_of_values + 1)
        t = 0.1 * (np.exp(i * np.log(1 + 10 * end) / number_of_values) - 1)
        t[0] += start
        return t
    return np.linspace(start, end, number_of_values)


def group_t(time_interval, pattern):
    """
    Group time intervals following a PSMC-like pattern, e.g. "4*1+25*2+4+6":
    4 blocks of 1 interval, 25 blocks of 2, 1 block of 4, 1 block of 6.
    """
    t = list(time_interval) + [time_interval[-1]]
    grouped = [t[0]]
    pos = 0
    for block in pattern.split("+"):
        n_blocks, size = map(int, block.split("*")) if "*" in block else (1, int(block))
        for _ in range(n_blocks):
            pos += size
            grouped.append(t[pos])
    return np.array(grouped)
