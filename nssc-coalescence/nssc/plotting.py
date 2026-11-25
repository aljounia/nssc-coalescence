"""
Shared plotting helpers. Each curve takes its style from its JSON entry
(keys color, linestyle, linewidth, alpha, label).
"""
import matplotlib.pyplot as plt

STYLE_KEYS = {"color": "color", "linestyle": "ls", "linewidth": "linewidth",
              "alpha": "alpha", "label": "label"}


def style(entry):
    """matplotlib arguments taken from a JSON entry (missing keys are ignored)."""
    return {mpl: entry[key] for key, mpl in STYLE_KEYS.items() if key in entry}


def plot_curve(ax, x, y, entry, **kwargs):
    ax.plot(x, y, **style(entry), **kwargs)


def add_vertical_lines(ax, times, scale=1.0):
    for t in times:
        ax.axvline(scale * t, color="k", ls="--", alpha=0.3)


def add_piecewise_functions(ax, functions):
    """Piecewise constant functions (already at the figure's scale)."""
    for f in functions:
        ax.step(f["x"], f["y"], where="post", **style(f))


def finish_figure(fig, ax, plot_params, output_file):
    """Legend, axes, title; saves and/or shows the figure according to plot_params."""
    ax.set_xlabel(plot_params.get("plot_xlabel", ""))
    ax.set_ylabel(plot_params.get("plot_ylabel", ""))
    if plot_params.get("time_log"):
        ax.set_xscale("log")
    if plot_params.get("y_scale") == "log":
        ax.set_yscale("log")
    if "plot_limits" in plot_params:
        x_a, x_b, y_a, y_b = plot_params["plot_limits"]
        ax.set_xlim(x_a, x_b)
        ax.set_ylim(y_a, y_b)
    if plot_params.get("plot_title"):
        ax.set_title(plot_params["plot_title"])
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="best")
    if plot_params.get("save_figure"):
        fig.savefig(output_file, format="png", dpi=300)
        print(f"Figure saved: {output_file}")
    if plot_params.get("show_plot"):
        plt.show()
