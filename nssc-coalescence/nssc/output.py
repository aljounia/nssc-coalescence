"""Writing results to text files."""
import numpy as np


def write_curve(path, x, y):
    """Two lines: the x values separated by spaces, then the y values."""
    with open(path, "w") as f:
        f.write(" ".join(str(v) for v in x) + "\n")
        f.write(" ".join(str(v) for v in y) + "\n")


def write_table(path, header, rows):
    """Table with a header line, columns separated by spaces."""
    np.savetxt(path, rows, delimiter=" ", header=" ".join(header), comments="")
