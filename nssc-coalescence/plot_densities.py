#!/usr/bin/env python3
"""
Plot the densities of T2: empirical densities from ms simulations
("scenarios") and theoretical densities ("theoretical_IICR_nisland",
"theoretical_IICR_general", if plot_params.plot_theor_IICR is 1).

Time is in units of 2N0 generations (no rescaling).

Optional JSON keys, in "plot_densities":
    densities_to_plot : indices of the ms scenarios to plot (empty = all)
    x_lim, y_lim      : axis limits

Output: <file>_densities.png if plot_params.save_figure.
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from nssc import config
from nssc.ms import empirical_density, simulate_t2
from nssc.plotting import finish_figure, plot_curve
from nssc.scenarios import general_model, nisland_model

N_POINTS = 2000


def main():
    parser = argparse.ArgumentParser(description="Plot theoretical and/or simulated densities of T2.")
    parser.add_argument("params_file", help="JSON parameter file")
    args = parser.parse_args()

    p = config.load_params(args.params_file)
    options = p.get("plot_densities", {})
    fig, ax = plt.subplots()

    # Empirical densities (ms)
    scenarios = p.get("scenarios", [])
    if scenarios:
        x = config.ms_time_grid(p)
        for i in options.get("densities_to_plot") or range(len(scenarios)):
            obs = simulate_t2(scenarios[i], p.get("path2ms", ""))
            plot_curve(ax, x, empirical_density(obs, x), scenarios[i], drawstyle="steps-mid")

    # Theoretical densities
    if p["plot_params"].get("plot_theor_IICR"):
        x_max = options.get("x_lim", [0, config.theory_time_grid(p)[-1]])[1]
        t = np.linspace(0, x_max, N_POINTS)
        models = ([(e, nisland_model(e)) for e in p.get("theoretical_IICR_nisland", [])]
                  + [(e, general_model(e)) for e in p.get("theoretical_IICR_general", [])])
        for entry, model in models:
            plot_curve(ax, t, model.pdf(t), entry)

    plot_params = {
        "plot_xlabel": "time (2N0 generations)",
        "plot_ylabel": "density of T2",
        "plot_title": p["plot_params"].get("plot_title", ""),
        "save_figure": p["plot_params"].get("save_figure", 0),
        "show_plot": p["plot_params"].get("show_plot", 0),
    }
    if "x_lim" in options and "y_lim" in options:
        plot_params["plot_limits"] = [*options["x_lim"], *options["y_lim"]]
    finish_figure(fig, ax, plot_params, config.output_prefix(args.params_file) + "_densities.png")


if __name__ == "__main__":
    main()
