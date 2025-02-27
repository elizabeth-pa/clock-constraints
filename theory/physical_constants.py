'''
    physical_constants.py

    Used for the theory plots.

    All quantities are converted to eV, that is,
    natural units in which c = hbar = 1.

    Here Mpl = (8 pi G)^{-1/2} is the *reduced* Planck mass.
'''

from math import pi as PI

### Fundamental constants ###
hbar = 1.973e-5     # cm eV
c    = 2.998e10     # cm s^-1
Mpl  = 2.435e27     # eV

### Conversion factors ###
g_to_eV = 1e9 / 1.8e-24

### Time ###
second = c / hbar
minute = 60 * second
hour = 60 * minute
day = 24 * hour
year = 365 * day

### Solar System ###
AU = 1.5e13 / hbar  # cm to eV^-1
solar_mass = 2e33 * g_to_eV
earth_mass = 1.7e28 * g_to_eV
earth_eccentricity = 0.0167

### Cosmology ###
Lambda_DE = 2.4e-3  # eV
rho_DE = Lambda_DE**4.  # all densities are in eV^4
rho_DM = 1e-11
rho_DM_local = 2.6e-6
