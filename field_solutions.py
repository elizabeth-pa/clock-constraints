"""
    field_solutions.py

    Contains solutions for the gradient of the scalar field,
    vs distance to a spherical object, in various theories.

    Used for the modified gravity clock signals, as well as for
    some of the existing bounds like MICROSCOPE.
"""

import numpy as np

PI = np.pi

def massless(M_obj, r, M):
    """ Massless scalar """
    return M_obj / (4*PI * M * r**2)

def yukawa(M_obj, r, M, m):
    """ Massive scalar """
    return dphi_massless(M_obj, r, M) * np.exp(- m * r)

def galileon_3(M_obj, r, M, Lambda):
    """ Cubic galileon """
    L = Lambda
    rhs = M_obj / (4*PI * M * r**3)
    return r * np.sqrt(L**3 / 2. * rhs)

def galileon_4(M_obj, r, M, Lambda, c4):
    """ Quartic galileon """
    L = Lambda
    rhs = M_obj / (4*PI * M * r**3)
    return r * np.pow(L**6 / (2*c4) * rhs, 1/3.)

def generalized(M_obj, r, M, Lambda, alpha, beta):
    """ Generalized force model """
    out = Lambda**2
    out *= np.pow( M_obj / (8*PI*M), alpha )
    out *= np.pow( Lambda * r, -beta)
    return out
