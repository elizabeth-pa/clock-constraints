"""
    plot_DE.py

    Creates the two dark energy plots in the paper.
"""
import matplotlib.pyplot as plt
import numpy as np

from theory.physical_constants import *
import theory.bounds as bounds

import mu_constraints as mc
import plot_options as po

def clock_bound(M, M_e, clock_pair, w=-0.95):
    # Amplitude in s^-1.  Convert to eV
    A = mc.DE_max_amplitude(clock_pair) * hbar / c
    M_eff = (M**-1 - M_e**-1)**-1

    X = rho_DE * (1. + w) / (1. - w)

    M_eff = np.abs(M_eff)
    return A < np.sqrt(X) / M_eff   # False = 0, True = 1


def draw_clock_bound(X, Y, clocks, col, ls='solid', w=-0.95):
    M   = Mpl * np.pow(10, X)
    M_e = Mpl * np.pow(10, Y)

    Z = clock_bound(M, M_e, clocks, w)

    plt.contour(X, Y, Z,  levels=[0.5], colors = col, linestyles=ls)
    #plt.contourf(X, Y, Z, levels=[-1, 0], colors = col, alpha = 0.3)
    return



### First plot ###
# M_e vs M
plt.figure(1, figsize=po.figsize )

# logarithmic limits
# i.e. M = Mpl * 10^x
xmin, xmax = 4, 8
ymin, ymax = 3, 8

X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(ymin, ymax, po.resolution))

M  = Mpl * np.pow(10, X)
M_e = Mpl * np.pow(10, Y)
# Curves

# Microscope #
Z = bounds.microscope_massless(M, M_e)
plt.contour(X, Y, Z, levels=[0.5], colors='gray', linewidths = 1)
plt.contourf(X, Y, Z, levels=[0.5, 1.], colors='gray', alpha = 0.5)

# Clocks
draw_clock_bound(X, Y, 'CaF/Sr', po.colorcycle[3], 'dotted', w = -0.999)
draw_clock_bound(X, Y, 'CaF/Sr', po.colorcycle[1], 'dashed', w= -0.99)
draw_clock_bound(X, Y, 'CaF/Sr', po.colorcycle[0], w = -0.95)


# Style
plt.tick_params( direction = 'in' )

plt.xticks([4,5,6,7,8])

plt.xlabel(r"$\log_{10} M~/~M_\mathrm{Pl}$", fontsize=12)
plt.ylabel(r"$\log_{10} M_\mathrm{e}~/~M_\mathrm{Pl}$", fontsize=12)

x_pos = 7.1
pi = 3.14
plt.text(x_pos, 3.21, "Microscope", ha = "left")
plt.text(x_pos, 6.54, "$w = -0.95$", c=po.colorcycle[0], ha = "left")
plt.text(x_pos, 6.21, "$w = -0.99$", c=po.colorcycle[1], ha = "left")
plt.text(x_pos, 5.73, "$w = -0.999$", c=po.colorcycle[3], ha = "left")


props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
plt.text(6, 7.7, "Dark Energy",
         fontsize=12, ha='center', bbox=po.title_box_properties)

plt.savefig("plots/DE_Me_vs_M.png", dpi=po.dpi_setting)

### Second plot ###
# M_eff vs w

plt.figure(2, figsize=po.figsize)
plt.clf()

# semi log axes
# w on x axis
# M_eff on y axis
xmin, xmax = -1.005, -0.94
ymin, ymax = 3, 7

X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(ymin, ymax, po.resolution))

# Microscope #

w = X
Meff = Mpl * np.pow(10, Y)

Z = bounds.microscope_massless(Meff, np.inf)
plt.contour(X, Y, Z, levels=[0.5], colors='gray')
plt.contourf(X, Y, Z, levels=[0.5, 1], colors='gray', alpha=0.5)

Z = bounds.microscope_massless(np.inf, Meff)
plt.contour(X, Y, Z, levels=[0.5], colors='gray')
plt.contourf(X, Y, Z, levels=[0.5, 1], colors='gray', alpha=0.5)

plt.text(-0.999, 5.05, r"Microscope ($M_\mathrm{e} \to \infty$)", ha = "left", fontsize=6)
plt.text(-0.999, 3.225, r"Microscope ($M \to \infty$)", ha = "left", fontsize=6)


# Clocks #
Z = clock_bound(Meff, np.inf, 'CaF/Sr', w)
plt.contour(X, Y, Z,  levels=[0.5], colors=po.colorcycle[0])
plt.contourf(X, Y, Z,  levels=[0.5,1], colors=po.colorcycle[0], alpha=0.3)

Z = clock_bound(Meff, np.inf, 'Cs/Sr', w) 
plt.contour(X, Y, Z,  levels=[0.5], colors=po.colorcycle[1])
plt.contourf(X, Y, Z,  levels=[0.5,1], colors=po.colorcycle[1], alpha=0.3)

plt.text(-0.98, 6.15, "CaF/Sr clocks", color='black', rotation=7, fontsize=12)
plt.text(-0.98, 4.9, "Cs/Sr clocks", color='black', rotation=7, fontsize=12)


# Planck
Z = bounds.planck(w)
plt.contour(X, Y, Z,  levels=[0.5], colors = 'brown')
plt.contourf(X, Y, Z,  levels=[0.5, 1], colors = 'brown', alpha=0.3)
plt.text(-0.9495, 3.75, "Planck", rotation = "vertical", fontsize=8)

# Phantom dark energy
plt.axvline(x = -1, color = 'black', linestyle = 'dashed')
plt.text(-1.002, 6, "Phantom DE", rotation = 'vertical', fontsize=8)

# Style and annotations #

plt.xticks([-1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94])
plt.yticks([3,4,5,6,7])


x_pos = np.mean([xmin, xmax])
plt.text(x_pos, 6.76, "Dark Energy",
         fontsize=12, ha='center', bbox=po.title_box_properties)
plt.tick_params( direction = 'in' )

plt.xlabel(r"DE equation of state $w$", fontsize=12)
plt.ylabel(r"$\log_{10} M_{\rm eff} ~/~ M_\mathrm{Pl}$", fontsize=12)

plt.savefig("plots/DE_Meff_vs_w.png", dpi=po.dpi_setting)

#plt.show()
