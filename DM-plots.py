"""
    DM-plots.py

    Creates the two dark matter plots in the paper.
    One is for an ultralight dark matter theory,
    and the other is for the raw signal
        delta mu = A / w cos(w t)
"""
import matplotlib.pyplot as plt
import numpy as np

from theory.physical_constants import *
import theory.bounds as bounds

import mu_constraints as mc
import plot_options as po

### Ultralight dark matter bounds ###
plt.figure(1, figsize=po.figsize)

xmin, xmax = -25, -16
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(3, 10, po.resolution))

m = np.pow(10, X)
Meff = Mpl * np.pow(10, Y)

def clock_bound_raw(clocks):
    # Angular frequencies, and amplitudes, for the clock signals
    # Signals are ~ A / w * cos(w t), and both A, w are
    # in units of s^-1
    if clocks == 'Sherrill':
        w_vals, A_vals = mc.DM_max_amplitude_Sherrill()
    else:
        w_vals, A_vals = mc.DM_max_amplitude(clocks)

    return w_vals, A_vals

def clock_bound_ULDM(clocks, color):

    w_vals, A_vals = clock_bound_raw(clocks)

    # Convert s^-1 to eV
    m = w_vals * hbar / c
    A_vals = A_vals * hbar / c

    rho = rho_DM_local

    # Solve for the M that saturates the bound.
    # Note that the w is already accounted for:
    #   A / w = sqrt(2 rho_DM) / (Meff m)
    # and w = m, so we only need
    #   Meff = sqrt(2 rho) / A
    M = np.sqrt(2*rho) / (A_vals)

    x,y = np.log10(m), np.log10(M / Mpl)

    # Add points to the ends so the outline covers the sides
    # This looks better than plt.vlines(),
    # as it takes care of the corners
    x = np.concatenate( ([x[0]], x, [x[-1]]) )
    y = np.concatenate( ([0], y, [0]) )

    plt.plot(x, y, color=color)
    plt.fill_between(x, 0, y, alpha=0.3, color=color)
    return

# Clocks
clock_bound_ULDM('CaF/Sr', po.colorcycle[0])
clock_bound_ULDM('Cs/Sr', po.colorcycle[1])

clock_bound_ULDM('Sherrill', 'purple')

# CMB
col = 'brown'
Z = bounds.DM_CMB(m)
plt.contour(X, Y, Z, levels = [0.5], colors=col, linestyles='dashed')
plt.contourf(X, Y, Z, levels = [0.5, 1], alpha = 0.3, colors=col)

plt.text(-24.25, 6, r'CMB & LSS', rotation='vertical', fontsize=8)

# Text labels for clocks
plt.text(-22, 7.5, "CaF/Sr clocks", fontsize=10, rotation=-34)
plt.text(-22, 6.35, "Cs/Sr clocks", fontsize=10, rotation=-34)
plt.text(-18.65, 4.25, 'Yb/Cs\nclocks', ha='center', fontsize=8)

# Microscope
col = 'gray'
Z = bounds.microscope_massless(Meff, np.inf)
plt.contour(X, Y, Z, levels = [0.5], colors=col)
plt.contourf(X, Y, Z, levels = [0.5, 1], alpha = 0.3, colors=col)

plt.text(-22.25, 4.96, r'Microscope ($M_\mathrm{e} \to \infty$)',
         fontsize=8)
plt.text(-22.25, 3.15, r'Microscope ($M \to \infty$)',
         fontsize=8)

Z = bounds.microscope_massless(np.inf, Meff)
plt.contour(X, Y, Z, levels = [0.5], colors=col)
plt.contourf(X, Y, Z, levels = [0.5, 1], alpha = 0.3, colors=col)

plt.xlim(xmin, xmax)
plt.ylim(3, 10)

plt.tick_params(direction='in')

x_pos = np.mean([xmin, xmax])
props = po.title_box_properties
plt.text(x_pos, 9.55, "Ultralight Dark Matter",
         ha='center', fontsize=12, bbox=props)


# Second axis
ax = plt.gca()
# Functions for converting to/from the second axis units
def m_to_Hz(m):
    # m = w = 2 pi f
    m = np.pow(10, m)
    f = m / (2*PI)
    f *= c / hbar   # convert to s^-1
    return np.log10(f)

def Hz_to_m(f):
    f = np.pow(10, f)
    m = 2*PI*f
    m *= hbar / c
    return np.log10(m)

ax2 = ax.secondary_xaxis('top', functions=(m_to_Hz, Hz_to_m))
ax2.tick_params(direction='in')
ax2.set_xlabel(r'Frequency $\log_{10} \frac{m}{2\pi} ~/~ \mathrm{Hz}$',
               fontsize=10)

plt.text(-17.1, 3.5, r"$m = 2 \pi (10~\mathrm{min})^{-1}$",
         fontsize = 6, rotation='vertical', ha='left')
plt.text(-22.65, 3.5, r"$m = 2 \pi (3~\mathrm{yr})^{-1}$",
         fontsize = 6, rotation='vertical', ha='left')


plt.xlabel(r"$\log_{10} m ~/~ \mathrm{eV}$", fontsize=12)
plt.ylabel(r"$\log_{10} M_\mathrm{eff} ~/~  M_\mathrm{Pl}$", fontsize=12)

plt.savefig("plots/DM-M-vs-m.png", dpi=po.dpi_setting)

### Raw signal constraint plot ###

plt.figure(2, figsize=po.figsize)

xmin, xmax = -9, -2
ymin, ymax = -26, -15
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(ymin, ymax, po.resolution))

def draw_raw_clock_bound(clocks, color):
    w_vals, A_vals = clock_bound_raw(clocks)

    # Convert from angular frequency to Hz
    # w = 2 pi f
    f_vals = w_vals / (2*PI)

    x = np.log10(f_vals)
    y = np.log10(A_vals)

    x = np.concatenate( ([x[0]], x, [x[-1]]) )
    y = np.concatenate( ([0], y, [0]) )

    plt.plot(x, y, color=color)
    plt.fill_between(x, y, 0, alpha=0.3, color=color)
    return

draw_raw_clock_bound('CaF/Sr', po.colorcycle[0])
draw_raw_clock_bound('Cs/Sr', po.colorcycle[1])
draw_raw_clock_bound('Sherrill', 'purple')

plt.text(-7.5, -23.35, "CaF/Sr clocks", fontsize=10, rotation=17)
plt.text(-7.5, -22.1, "Cs/Sr clocks", fontsize=10, rotation=17)
plt.text(-4, -19, 'Yb/Cs\nclocks', ha='center', fontsize=8)

plt.text(-2.72, -18, r"$f = (10~\mathrm{min})^{-1}$",
         fontsize = 6, rotation='vertical', ha='left')
plt.text(-8.2, -18, r"$f = (3~\mathrm{yr})^{-1}$",
         fontsize = 6, rotation='vertical', ha='left')

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

plt.tick_params(direction='in')
plt.xlabel(r"Frequency $\log_{10} f ~/~ \mathrm{Hz}$", fontsize=12)
plt.ylabel(r"Amplitude $\log_{10} A ~/~  \mathrm{s}^{-1}$", fontsize=12)

props = po.title_box_properties
x_pos = np.mean([xmin, xmax])
plt.text(x_pos, -15.65, "Dark Matter Signal",
         ha='center', fontsize=12, bbox=props)

plt.savefig("plots/DM-A-vs-f.png", dpi=po.dpi_setting)

#plt.show()
