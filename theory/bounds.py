"""
    bounds.py

    A number of existing bounds on the various theories
    under consideration.
    Particularly MICROSCOPE and CMB bounds.
"""
import numpy as np

from theory.physical_constants import *

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

    # TO DO: m should be fixed up to only be sourced by the 
    # appropriate combination of m, M, M_e
    # Not currently needed since we almost always  set M_e \to \infty.
    # The single time this isn't the case is handled with
    # microscope_massless().

    # The dphi functions generally set M_e \to \infty
    #Q = m * (M**-1 + eps_earth / M_e**-1) * 1 / (1 + eps_earth)

    out = 8 * PI * Mpl**2 * R**2 / m
    out *= (M**-1 - M_e**-1)
    out *= (eps_Pt - eps_Ti)
    out *= dphi(m, R)

    out = np.abs(out)

    return eta < out

def planck(w):
    return -0.95 < w

def DM_CMB(m):
    return m < 1e-24
