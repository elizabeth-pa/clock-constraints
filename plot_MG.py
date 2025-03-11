"""
    plot_MG.py

    Creates the two modified gravity plots in the paper.
    One is for the galileon (cubic and quartic), the second
    is for the generalized interaction model that is described
    in the paper.
"""
import matplotlib.pyplot as plt
import numpy as np

from theory.physical_constants import *
import theory.bounds as bounds
import theory.field_solutions as dphi

import mu_constraints as mc
import plot_options as po

def clock_bound(M, dphi):
    """ Modified gravity clock signal.
    Uses the effect of the Earth going around the Sun.
    Parameters:
        M       The scalar-matter coupling M.  This really is M_eff but
                we generally set M_e -> infinity for the modified
                gravity signal.
        dphi    dphi(m, r), which gives the field gradient around
                a spherical object of mass m as a function of distance r.
    """
    A = mc.MG_max_amplitude('CaF/Sr')   # dimensionless

    M_obj = solar_mass
    r = AU
    e = earth_eccentricity

    dp = dphi(M_obj, r)

    return A < (dp * r * e) / M
    

### Galileon ###
plt.figure(1, figsize=po.figsize)

xmin, xmax = -16, 0
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(-6, 6, po.resolution))

Lambda = np.pow(10, X)
M = Mpl * np.pow(10, Y)

# Clock

dp = lambda m, r: dphi.massless(m, r, M)
Z2 = clock_bound(M, dp)

dp = lambda m, r: dphi.galileon_3(m, r, M, Lambda)
Z3 = clock_bound(M, dp)

c4 = 1e-12
dp = lambda m, r:  dphi.galileon_4(m, r, M, Lambda, c4)
Z4 = clock_bound(M, dp)

Z = Z2 & Z3 & Z4

col = po.colorcycle[0]
plt.contour(X, Y, Z,  levels=[0.5], colors=col)
plt.contourf(X, Y, Z,  levels=[0.5, 1], colors=col, alpha=0.3)
plt.text(np.log10(2e-5), np.log10(3e3),
         "CaF/Sr clocks", va = "center")

# Microscope
dp = lambda m, r: dphi.massless(m, r, M)
Z2 = bounds.microscope(M, np.inf, dp)

dp = lambda m, r: dphi.galileon_3(m, r, M, Lambda)
Z3 = bounds.microscope(M, np.inf, dp)

dp = lambda m, r: dphi.galileon_4(m, r, M, Lambda, c4)
Z4 = bounds.microscope(M, np.inf, dp)

Z = Z2 & Z3 & Z4

col = 'gray'
#plt.contour(X, Y, Z, levels=[0.5], colors=col, linestyles='dashed')
plt.contourf(X, Y, Z, levels=[0.5, 1], colors=col, alpha=0.3)
plt.text(np.log10(5e-4), np.log10(8e4),
         "Microscope", va = "center", fontsize = 8)

# LLR
Z1 = bounds.LLR(M, Lambda)
rV3 = bounds.rV3(M, Lambda, earth_mass)
rV4 = bounds.rV4(M, Lambda, earth_mass, c4)

Z2 = earth_moon_distance < rV3
Z3 = earth_moon_distance > rV4**4 / rV3**3

Z = Z1 & Z2 & Z3

col = 'gray'
#plt.contour(X, Y, Z, levels=[0.5], colors=col, linestyles='dashed')
plt.contourf(X, Y, Z, levels=[0.5, 1], colors=col, alpha=0.3)

plt.text(-6.85, 4.5,
         "LLR", va = "center", fontsize = 8)

# Styling
plt.tick_params( direction = 'in' )

plt.xticks( range(-16, 1, 2) )

x_pos = np.mean([xmin, xmax])
plt.text(x_pos, 5.23, "Galileon",
         fontsize=12, ha='center', bbox=po.title_box_properties)

# Text box for c_4
#props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
props = po.title_box_properties
plt.text(-14, 4, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")

plt.axhline(0, linestyle='dashed', color='black')
plt.text(-15.5, -0.4, r"$M = M_\mathrm{Pl}$", fontsize=8)

plt.axvline(np.log10(1.3e-13), linestyle='dashed', color='black')
plt.text(-13.45, -2, r"$m_\mathrm{g} = H_0$", fontsize=8, rotation='vertical')

# Labels of Vainshtein regions
plt.text(np.log10(1e-2), np.log10(8.5e3),
         r"$R_{\rm V3} < {\rm AU}$", fontsize = 6)
plt.text(-9.925, 1,
         r"$R_{\rm V4}^4 / R_{\rm V3}^3 < {\rm AU} < R_{\rm V3}$",
         fontsize = 6, rotation = 53)
plt.text(-14.92, -5, r"${\rm AU} < R_{\rm V4}^4 / R_{\rm V3}^3$",
         fontsize = 6, rotation = 66)

plt.xlabel(r'$\log_{10} \Lambda ~/~ \mathrm{eV}$', fontsize=12)
plt.ylabel(r'$\log_{10} M ~/~ M_\mathrm{Pl}$', fontsize=12)

# Second set of axes
def Lambda_to_mg(L):
    L = np.pow(10, L)
    out = np.sqrt(L**3 / Mpl)
    out = np.log10(out)
    return out

def mg_to_Lambda(mg):
    mg = np.pow(10, mg)
    out = np.pow(mg**2 * Mpl, 1/3.)
    out = np.log10(out)
    return out

ax = plt.gca()
# For this to work automatically, two functions converting to/from
# the two axis coordinates are needed.
ax2 = ax.secondary_xaxis('top', functions=(Lambda_to_mg, mg_to_Lambda))

ax2.tick_params(direction='in')
ax2.set_xlabel(r'$\log_{10} m_g ~/~ \mathrm{eV}$', fontsize=12)

plt.savefig("plots/MG_gal_M_vs_Lambda.png", dpi=po.dpi_setting)


### Generalized model ###

plt.figure(2, figsize=po.figsize)


xmin, xmax = -0.05, 2.05
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(-0.05, 1.05, po.resolution))

betas  = X
alphas = Y

# Model parameters
M = Mpl
Lambda = 1e-10

## Microscope
dp = lambda m, r: dphi.generalized(m, r, M, Lambda, alphas, betas)
Z = bounds.microscope(M, np.inf, dp)
col = 'gray'
plt.contour(X, Y, Z, levels=[0.5], colors=col, linestyles='dashed')
plt.contourf(X, Y, Z, levels=[0.5, 1], colors=col, alpha=0.3)

plt.text(1.7, 0.665, "Microscope", rotation=12.5, fontsize = 8)

## Clocks

# We can reuse the same field gradient function dp(m, r) that was used for 
# microscope
Z = clock_bound(M, dp)

plt.text(1, 0.61, "CaF/Sr clocks", rotation=23, fontsize = 12)

col = po.colorcycle[0]
plt.contour(X, Y, Z,  levels=[0.5], colors=col)
plt.contourf(X, Y, Z,  levels=[0.5, 1], colors=col, alpha=0.3)

# Label theories in \alpha, \beta space
col = 'black'
fs = 8
plt.plot(0.5, 0.5, 'o', label='Cubic Galileon', color=col)
plt.text(0.515, 0.425, "Cubic\nGalileon", fontsize=fs, ha='left')

plt.plot(0., 0.4, 'o', label='Quartic Galileon', color=col)
plt.text(0.015, 0.32, "Quartic\nGalileon", fontsize=fs, ha='left')

plt.plot(0., 0.0, 'o', label='DBIon', color=col)
plt.text(0.02, 0.02, "DBIon", fontsize=fs, ha='left')

plt.plot(2., 1., 'o', label='Free scalar', color=col)
plt.text(1.625, 0.965, "Free scalar", fontsize=fs, ha='left')

# Text box with model parameter values
props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
s = r"$\Lambda = 10^{-10}~\mathrm{eV}$" + "\n" + r"$M = M_\mathrm{Pl}$"
plt.text(1.5, 0.025, s, bbox = props, ha = "left", fontsize = 10)

plt.tick_params(direction='in')

props = po.title_box_properties
plt.text(1.0, 0.975, 'Generalized interaction',
         fontsize=12, ha='center', bbox=props)

plt.xlabel(r"$\beta$", fontsize=12)
plt.ylabel(r"$\alpha$", fontsize=12)

plt.savefig("plots/MG_gen_alpha_vs_beta.png", dpi=po.dpi_setting)

#plt.show()
