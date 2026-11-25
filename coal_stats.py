#!/usr/bin/env python3
"""
Statistics of the coalescence time T2 of a pair of lineages for each
"theoretical_IICR_general" entry of the JSON file: mean, variance and
probability of observing at least one difference.

JSON keys used:
    scale_params, plot_params.time_log, plot_params.plot_limits[1] (integration grid)
    theoretical_IICR_general
    theta, L (optional, defaults 1e-4 and 500)

Output: <file>_stats.txt (one line per scenario).
"""
import argparse

from nssc import config
from nssc.output import write_table
from nssc.scenarios import general_model
from nssc.stats import coalescence_stats


def main():
    parser = argparse.ArgumentParser(description="Mean and variance of T2, and P(at least one difference).")
    parser.add_argument("params_file", help="JSON parameter file")
    args = parser.parse_args()

    p = config.load_params(args.params_file)
    t = config.theory_time_grid(p)
    theta, n_sites = p.get("theta", 1e-4), p.get("L", 500)

    header = ["tau_last", "mean", "variance", "p_at_least_one_diff"]
    rows = []
    for entry in p["theoretical_IICR_general"]:
        stats = coalescence_stats(general_model(entry), t, theta, n_sites)
        rows.append([entry.get("tau", [0])[-1], *stats.values()])

    out_file = config.output_prefix(args.params_file) + "_stats.txt"
    write_table(out_file, header, rows)
    print(" ".join(header))
    for row in rows:
        print(" ".join(f"{x:.6g}" for x in row))
    print(f"\nResults written to {out_file}")


if __name__ == "__main__":
    main()
