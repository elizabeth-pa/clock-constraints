import matplotlib.pyplot as plt
import numpy as np

from physical_constants import *
import mu_constraints as mc
import plot_options as po
#import field_solutions as dphi
import bounds


PI = np.pi

### Ultralight dark matter bounds ###
plt.figure(1, figsize=po.figsize)

xmin, xmax = -25, -16
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(3, 10, po.resolution))

m = np.pow(10, X)
Meff = Mpl * np.pow(10, Y)

def clock_bound_ULDM(clocks, color):
    # Angular frequencies, and amplitudes, for the clock signals
    # Signals are ~ A / w * cos(w t), and both A, w are
    # in units of s^-1
    w_vals, A_vals = mc.DM_max_amplitude(clocks)

    # Convert s^-1 to eV
    m = w_vals * hbar / c
    A_vals = A_vals * hbar / c

    rho = rho_DM_local

    # Solve for the M that saturates the bound.
    # Note that the w is already accounted for, as the signal
    # is = A / w * cos(w t)
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

# CMB
col = 'brown'
Z = bounds.DM_CMB(m)
plt.contour(X, Y, Z, levels = [0.5], colors=col, linestyles='dashed')
plt.contourf(X, Y, Z, levels = [0.5, 1], alpha = 0.3, colors=col)

plt.text(-24.25, 6, r'CMB & LSS', rotation='vertical', fontsize=8)

# Text labels for clocks
plt.text(-22, 7.5, "CaF/Sr clocks", fontsize=10, rotation=-34)
plt.text(-22, 6.35, "Cs/Sr clocks", fontsize=10, rotation=-34)

# Microscope
col = 'gray'
Z = bounds.microscope_massless(Meff, np.inf)
plt.contour(X, Y, Z, levels = [0.5], colors=col)
plt.contourf(X, Y, Z, levels = [0.5, 1], alpha = 0.3, colors=col)

plt.text(-22, 4.96, r'Microscope ($M_\mathrm{e} \to \infty$)', fontsize=8)
plt.text(-22, 3.15, r'Microscope ($M \to \infty$)', fontsize=8)

Z = bounds.microscope_massless(np.inf, Meff)
plt.contour(X, Y, Z, levels = [0.5], colors=col)
plt.contourf(X, Y, Z, levels = [0.5, 1], alpha = 0.3, colors=col)

plt.xlim(xmin, xmax)
plt.ylim(3, 10)

plt.tick_params(direction='in')

x_pos = np.mean([xmin, xmax])
props = po.title_box_properties
plt.text(x_pos, 9.55, "Dark matter",
         ha='center', fontsize=12, bbox=props)

def m_to_period(m):
    m = np.pow(10, m)
    T = hbar / (c * m)
    return np.log10(T)

def period_to_m(T):
    T = np.pow(10, T)
    m = hbar / (c * T)
    return np.log10(m)

ax = plt.gca()
ax2 = ax.secondary_xaxis('top', functions=(m_to_period, period_to_m))
ax2.tick_params(direction='in')

plt.xlabel(r"$\log_{10} m ~/~ \mathrm{eV}$", fontsize=12)
plt.ylabel(r"$\log_{10} M_\mathrm{eff} ~/~  M_\mathrm{Pl}$", fontsize=12)

plt.savefig("plots/DM-M-vs-m.png", dpi=po.dpi_setting)

plt.show()

