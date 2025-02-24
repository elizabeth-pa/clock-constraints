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

xmin, xmax = -25, -17
X, Y = np.meshgrid(
        np.linspace(xmin, xmax, po.resolution),
        np.linspace(3, 10, po.resolution))

w_vals, A_vals = mc.DM_max_amplitude('CaF/Sr')

rho = rho_DM_local

m = w_vals * hbar / c

# TO DO: this step shouldn't be necessary
m /= 2*PI

A_vals = A_vals * hbar / c

# Solve for the M that saturates the bound.
# Note that the w is already accounted for, as the signal
# is = A / w * cos(w t)
M = np.sqrt(2*rho) / (A_vals)

plt.plot( np.log10(m), np.log10(M / Mpl))

plt.show()

