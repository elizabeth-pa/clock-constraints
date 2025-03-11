"""
    plot_options.py

    Various settings for matplotlib.
    IMPORTANT: for the final versions of the plots, turn up the
    resolution to something like 10k
"""

import matplotlib.pyplot as plt

import numpy as np

# Number of grid points in the plots.
# Around 1k looks ok and is fast.
# Set this to a high value (~10k) for the final versions
# of the plots.
resolution = 10000

# Color scheme
colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Most standard printers are 300dpi, so the figsize is the actual
# size on a page.
figsize = (5, 5)
dpi_setting = 300

title_box_properties = dict(boxstyle='round', facecolor = "white", alpha=1.0)

def shade_below(x, y, color, alpha=0.3, boundary=True):
    """ Plot a curve, and shade below it. """

    # Add points to the ends so the outline covers the sides
    # This looks better than plt.vlines(),
    # as it takes care of the corners
    x = np.concatenate( ([x[0]], x, [x[-1]]) )
    y = np.concatenate( ([0], y, [0]) )

    if boundary: plt.plot(x, y, color=color)
    plt.fill_between(x, -100, y, color=color, alpha=alpha)
    return

def shade_above(x, y, color, alpha=0.3, boundary=True):
    """ Plot a curve, and shade above it. """

    # Add points to the ends so the outline covers the sides
    # This looks better than plt.vlines(),
    # as it takes care of the corners
    x = np.concatenate( ([x[0]], x, [x[-1]]) )
    y = np.concatenate( ([0], y, [0]) )

    if boundary: plt.plot(x, y, color=color)
    plt.fill_between(x, y, 100, color=color, alpha=alpha)
    return

