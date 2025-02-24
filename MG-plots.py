import matplotlib.pyplot as plt
import numpy as np

from physical_constants import *
import mu_constraints as mc
import plot_options as po
import field_solutions as dphi
import bounds


### Galileon ###
plt.figure(1, figsize=po.figsize)

xmin, xmax = -16, 0
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(-6, 6, po.resolution))

Lambda = np.pow(10, X)
M = Mpl * np.pow(10, Y)

# Clock
A = mc.MG_max_amplitude('CaF/Sr')   # dimensionless

M_obj = solar_mass
r = AU
e = earth_eccentricity
c4 = 1e-12

dp = dphi.massless(M_obj, r, M)
Z2 = A < (dp * r * e) / M

dp = dphi.galileon_3(M_obj, r, M, Lambda)
Z3 = A < (dp * r * e) / M

dp = dphi.galileon_4(M_obj, r, M, Lambda, c4)
Z4 = A < (dp * r * e) / M

Z = Z2 & Z3 & Z4

col = po.colorcycle[0]
plt.contour(X, Y, Z,  levels=[0.5], colors=col)
plt.contourf(X, Y, Z,  levels=[0.5, 1], colors=col, alpha=0.3)
plt.text(np.log10(2e-5), np.log10(3e3),
         "CaF/Sr clocks", va = "center")

# Microscope
dp = lambda m, r : dphi.massless(m, r, M)
Z2 = bounds.microscope(M, np.inf, dp)

dp = lambda m, r : dphi.galileon_3(m, r, M, Lambda)
Z3 = bounds.microscope(M, np.inf, dp)

dp = lambda m, r : dphi.galileon_4(m, r, M, Lambda, c4)
Z4 = bounds.microscope(M, np.inf, dp)

Z = Z2 & Z3 & Z4

col = 'gray'
plt.contour(X, Y, Z,  levels=[0.5], colors=col, linestyles='dashed')
plt.contourf(X, Y, Z,  levels=[0.5, 1], colors=col, alpha=0.3)
plt.text(np.log10(5e-4), np.log10(8e4),
         "Microscope", va = "center", fontsize = 8)


# Styling
plt.tick_params( direction = 'in' )

plt.xticks( range(-16, 1, 2) )

x_pos = np.mean([xmin, xmax])
plt.text(x_pos, 5.25, "Galileon",
         fontsize=12, ha='center', bbox=po.title_box_properties)

# Text box for c_4
props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
plt.text(-14, 4, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")

plt.axhline(0, linestyle='dashed', color='black')
plt.text(-15.5, -0.4, r"$M = M_\mathrm{Pl}$", fontsize=8)

plt.axvline(np.log10(1.3e-13), linestyle='dashed', color='black')
plt.text(-13.45, -2, r"$m_\mathrm{g} = H_0$", fontsize=8, rotation='vertical')

# Labels of Vainshtein regions
plt.text(np.log10(1e-2), np.log10(8.5e3),
         r"$R_{\rm V3} < {\rm AU}$", fontsize = 6)
plt.text(np.log10(1.5e-10), np.log10(1e1),
         r"$R_{\rm V4} < {\rm AU} < R_{\rm V3}$",
         fontsize = 6, rotation = 53)
plt.text(-14.9, -5, r"${\rm AU} < R_{\rm V4}$",
         fontsize = 6, rotation = 65)

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

#
plt.savefig("plots/MG-gal-M-vs-Lambda.png", dpi=po.dpi_setting)

plt.show()

