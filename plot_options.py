"""
    plot_options.py

    Various settings for matplotlib.
    IMPORTANT: for the final versions of the plots, turn up the
    resolution to something like 10k
"""

import matplotlib.pyplot as plt

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
