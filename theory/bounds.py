"""
    bounds.py

    A number of existing bounds on the various theories
    under consideration.
    Particularly MICROSCOPE and CMB bounds.
"""

import numpy as np

# This expects to be run from the parent folder, i.e. in plot_DE.py etc
from theory.physical_constants import *
import theory.field_solutions as solutions

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


def LLR(M, L):
    """ Lunar laser ranging bound on the cubic galileon
    As described in 1407.0059 """

    rV = rV3(M, L, earth_mass)

    g = Mpl / M
    R = earth_moon_distance
    return 1e-11 < g**2 * (R / rV)**(3/2.)

def rV3(M, L, M_obj):
    """ Cubic Vainshtein radius """
    return (M_obj / (2*PI*M))**(1/3.) / L


def rV4(M, L, M_obj, c4):
    """ Quartic Vainshtein radius """
    return (c4/2)**(1/6.) * rV3(M, L, M_obj)

delta_phi_LLR = 2.4e-11
def LLR3(M, L):
    """ Lunar laser ranging bound for cubic galileon.
    Described in 1106.1573 """ 
    m = earth_mass
    R = earth_moon_distance
    dphi = delta_phi_LLR

    rhs = 3*PI/2 * Mpl**2 * (8*PI * L**3 * R**3 /(m * M**3))**0.5
    return dphi < rhs

def LLR4(M, L, c4):
    """ LLR quartic galileon """
    m = earth_mass
    R = earth_moon_distance
    dphi = delta_phi_LLR

    rhs = 2*PI * Mpl**2 * L**2
    rhs *= ((8*PI)**2 / (c4 * M**4 * m**2))**(1/3.) * R**2
    return dphi < rhs

def LLR_generalized(M, L, a, b):
    """ LLR bound, generalized interaction model """
    m = earth_mass
    R = earth_moon_distance
    dphi = delta_phi_LLR

    rhs = PI * (2 - b) * (Mpl / M)**2
    rhs *= (m / (8*PI*M))**(a-1) * (L*R)**(2-b)
    return dphi < rhs


def planck(w):
    return -0.95 < w

def DM_CMB(m):
    return m < 1e-24


def load_amplitude_helper(fname):
    """
    Load a .csv file of m, d_m_e values.
    Convert to M_eff.
    """

    # m is in eV, d_m_e is dimensionless
    # both are the log10 of their real values
    import csv
    m, d_m_e = [], []
    with open(fname) as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            m.append(row[0])
            d_m_e.append(row[1])

    # These are log10 their actual values
    m = np.pow(10, m)
    d_m_e = np.pow(10, d_m_e)

    # Convert to our M_eff
    M_eff = np.sqrt(2) * Mpl / d_m_e
    return m, M_eff

def Meff_to_amplitude(m, M_eff):
    """ Generalized amplitude measurement bounds 
    Converts from the m, Meff bound to our generalized signal
        delta mu / mu = A / w cos(w t)
    Returns the frequency f and corresponding A values.
    Both are in units of s^-1
    """
    
    f = m / (2*PI)
    A = np.sqrt(2 * rho_DM_local) / M_eff

    # Convert to s^-1
    f *= c / hbar
    A *= c / hbar

    return f, A

def HSi_ULDM():
    """ H/Si clock bounds
    https://arxiv.org/pdf/2008.08773
    Converts to our M_eff coupling.  Returns an array of m values,
    and an array of the corresponding M_eff values that are the bound.
    """
    return load_amplitude_helper('theory/H_Si.csv')

def HSi_amplitude():
    """
    Converts H/Si bound to frequency and amplitude bounds.
    """
    m, M_eff = HSi_ULDM()
    return Meff_to_amplitude(m, M_eff)

def YbCs_ULDM():
    """ https://arxiv.org/pdf/2212.05721 """
    return load_amplitude_helper('theory/Yb_Cs.csv')

def YbCs_amplitude():
    m, M_eff = YbCs_ULDM()
    return Meff_to_amplitude(m, M_eff)


def NANOGrav_ULDM():
    """ https://arxiv.org/abs/2306.16219 """
    return load_amplitude_helper('theory/NANOGrav.csv')

