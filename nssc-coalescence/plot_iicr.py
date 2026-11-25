#!/usr/bin/env python3
"""
Plot on the same figure the IICRs described in the JSON file: empirical
IICRs from ms simulations ("scenarios"), theoretical IICRs
("theoretical_IICR_nisland", "theoretical_IICR_general", if
plot_params.plot_theor_IICR is 1), piecewise constant functions
("piecewise_constant_functions") and vertical lines ("vertical_lines").

Outputs (depending on the JSON file):
    <file>.png                              if plot_params.save_figure
    <file>_IICR_<family>_<i>_text_file.txt  if save_theor_IICR_as_file
"""
import argparse

import matplotlib.pyplot as plt

from nssc import config
from nssc.ms import ms_iicr
from nssc.output import write_curve
from nssc.plotting import add_piecewise_functions, add_vertical_lines, finish_figure, plot_curve
from nssc.theory import theoretical_curves


def main():
    parser = argparse.ArgumentParser(description="Plot theoretical and/or ms-simulated IICRs.")
    parser.add_argument("params_file", help="JSON parameter file")
    args = parser.parse_args()

    p = config.load_params(args.params_file)
    prefix = config.output_prefix(args.params_file)
    t_scale, n_scale = config.time_scale(p), config.size_scale(p)
    fig, ax = plt.subplots()

    # Empirical IICRs (ms)
    if p.get("scenarios"):
        x = config.ms_time_grid(p)
        x_plot = x.copy()
        if x_plot[0] == 0:          # avoid t = 0 on a logarithmic axis
            x_plot[0] = x_plot[1] / 100
        for scenario in p["scenarios"]:
            iicr = ms_iicr(scenario, x, p.get("path2ms", ""))
            plot_curve(ax, t_scale * x_plot, n_scale * iicr, scenario, drawstyle="steps-post")

    add_vertical_lines(ax, p.get("vertical_lines", []), t_scale)

    # Theoretical IICRs
    if p["plot_params"].get("plot_theor_IICR"):
        t = config.theory_time_grid(p)
        for family, i, entry, iicr in theoretical_curves(p, t):
            plot_curve(ax, t_scale * t, n_scale * iicr, entry)
            if p.get("save_theor_IICR_as_file"):
                write_curve(f"{prefix}_IICR_{family}_{i}_text_file.txt", t_scale * t, n_scale * iicr)

    add_piecewise_functions(ax, p.get("piecewise_constant_functions", []))
    finish_figure(fig, ax, p["plot_params"], f"{prefix}.png")


if __name__ == "__main__":
    main()
