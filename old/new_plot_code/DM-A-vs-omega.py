import matplotlib.pyplot as plt
import numpy as np
import DM_plotter

pi = np.pi

# Some times, in seconds
year_in_seconds = np.pi * 1e7
minute_in_seconds = 60

# Maximum and minimum times that we trust our projected measurements over
# TO DO: check these values and update if necessary
T_min = 5 * minute_in_seconds
T_max = 3 * year_in_seconds

omega_min = 2*pi / T_max
omega_max = 2*pi / T_min

# TO DO: replace these with the actual values,
# as these are simply placeholders!
cprime = 1e-5
h0 = 4.5e-30
hinv = 4.06e-35

def A_constraint(omega):
    A = cprime * omega
    A *= np.sqrt(h0 + 2 * np.pi * hinv / omega)
    return A


xs = np.logspace(np.log10(omega_min), np.log10(omega_max))
#ys = A_constraint(xs)

couples = np.genfromtxt("sigmas.csv", delimiter=",", dtype=float)[:,2]
ys=DM_plotter.A_thr(xs,2,sigmas=couples) # that's how to replace the curve The number 2 corresponds to the couple CaF/Sr

def add_edges(xs, ys, val = 1):
    """ Add points so that when plotted there is an
    outline along the edges """
    xs = np.insert(xs, 0, xs[0])
    ys = np.insert(ys, 0, val)

    xs = np.append(xs, xs[-1])
    ys = np.append(ys, val)
    return xs, ys


### Plot setup ###
plt.figure(figsize=(6, 6))
middle = np.power(10, (np.log10(omega_max) + np.log10(omega_min))/2)

## Main results
xs, ys = add_edges(xs, ys)
plt.loglog(xs, ys)
plt.fill_between(xs, ys, 1, alpha = 0.3)

plt.text(middle, 2.5e-24, "CaF/Sr\nclocks", ha = "center")

## Sherrill et al
# New J. Phys. 25 (2023) 093012
# Values drawn from Table 3
fmin = 8.7e-7
fmax = 8.3e-4

omega_min = 2 * pi * fmin
omega_max = 2 * pi * fmax

sigma_r = 1.6e-13 # which comes from 2302.04565, below equation (15)
h0_r = 2 * sigma_r**2 # the correct convertion to frequency domain white noise
sigma_A_over_om = np.sqrt(2*h0_r/T_max) # the Fisher forecast error for white noise over a total observation time T_max

xs_sherrill = np.logspace(np.log10(omega_min), np.log10(omega_max))
ys_sherrill = [sigma_A_over_om * omega for omega in xs_sherrill]    

xs_sherrill, ys_sherrill = add_edges(xs_sherrill, ys_sherrill)

plt.loglog(xs_sherrill, ys_sherrill)
plt.fill_between(xs_sherrill, ys_sherrill, 1, alpha = 0.3)

plt.text(1e-4, 5e-19, "Sr/Yb/Cs\nclocks", fontsize = 8, ha = "center")

### Plot formatting ###
ymin = ys[1] * 1e-1
#ymax = ys[-2] * 1e1
ymax = 1e-15

plt.ylim( ymin, ymax )
plt.xlabel(r"Frequency $\omega ~[\mathrm{rad} / s]$")
plt.ylabel(r"Amplitude $A ~[s^{-1}]$")

props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
plt.text(middle, ymax * 0.3, "Dark Matter Signal", ha = "center",
         fontsize = 12, bbox = props)

plt.tick_params(which="both", direction="in")
plt.savefig("DM_A_vs_omega.png", dpi = 300)
plt.show()
