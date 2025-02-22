import numpy as np

from physical_constants import *

PI = np.pi

microscope_eta = np.sqrt(2.3**2 + 1.5**2) * 1e-15
microscope_r = 7000 * 1e5 / hbar

# Epsilon values, summarized in our table
eps_Pt = 2e-4
eps_Ti = 2.3e-4
eps_earth = 2.33e-4 # Iron

def microscope_massless(M, M_e):
    """ A simplified form of the MICROSCOPE bound, for
    a massless scalar
    """
    eta = microscope_eta

    out = 2 * Mpl**2 * (M**-1 - M_e**-1)
    out *= (M**-1 + eps_earth / M_e)
    out *= (eps_Pt - eps_Ti)
    out = np.abs(out)

    return eta < out

def microscope(M, M_e, dphi):
    """ MICROSCOPE satellite.
    dphi(m, r) is a function for the field gradent as a function
    of distance from a spherical object for the theory at hand.
    """
    eta = microscope_eta
    R = microscope_r
    m = earth_mass

    out = 8 * PI * Mpl**2 * R**2 / m
    out *= (M**-1 - M_e**-1)
    out *= (eps_Pt - eps_Ti)
    out *= dphi(m, R)

    out = np.abs(out)

    return eta < out

def planck(w):
    return -0.95 < w
