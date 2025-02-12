"""
mu_constraints.py

Implements three main functions to be exported as a module:
    Function name       Signal in mu = m_p / m_e

    DE_max_amplitude    Linear drift
    DM_max_amplitude    Oscillating on arbitrary frequency w
    MG_max_amplitude    Oscillating on a 1-year period

Each of these takes in a clock pair and returns a maximum signal that
is detectable by that clock pair.

"""

import math
import numpy as np

PI = np.pi

def DE_max_amplitude(clock_pair = "CaF/Sr"):
    """Maximum dark energy signal that could have escaped detection.
        signal = A * t

    The signal is dimensionless, A has units of sec^-1
    Returns A for the given clock pair.
    """
    match clock_pair:
        case "CaF/Sr":
            return 1.7e-25
        case "Sr/Cs":
            return 3.1e-24

    raise Exception("Unknown clock pair.")
    return -1


def DM_max_amplitude(clock_pair = "CaF/Sr", nPoints = 100):
    """Maximum dark matter signal that could have escaped detection.
        signal = A/w * cos(w t)
    This signal in delta mu / mu is dimensionless.

    Inputs:
        clock_pair: Name of a particular clock pair
        nPoints:    The number of A and omega values to return

    Returns a tuple of two arrays of length nPoints:
        w_vals: Angular frequencies w = 2*pi*freq that correspond to each 
                value in the A_vals.  Covers the whole range of omega
                that the clock pair is able to constrain.
                These should be logarithmically distributed within the
                range that is constrained.
                Units of sec^-1
        A_vals: Corresponding signal amplitudes for each angular
                frequency in w_vals.
                Units of sec^-1
    """

    # Minimum and maximum angular frequencies that are constrained
    # These in general depend on the clock pair in question.
    w_min = 2*PI / 1e7     # year
    w_max = 2*PI / 100     # 100 sec
        
    # Values logarithmically spaced between w_min and w_max
    w_vals = np.logspace(math.log10(w_min),
                         math.log10(w_max), nPoints)

    # Fill in A_vals depending on the clock pair
    A_vals = []
    for w in w_vals:
        A = DM_helper(w, clock_pair)
        A_vals.append(A)
    return w_vals, A_vals

def MG_max_amplitude(clock_pair = "CaF/Sr"):
    """Maximum modified gravity signal that could have escaped detection.
        signal = A * cos(2 pi t / year)

    The signal and A are dimensionless.

    Returns A for the given clock pair.
    """
    match clock_pair:
        case "CaF/Sr":
            return 7.8e-18
        case "Sr/Cs":
            return 1.4e-16

    raise Exception("Unknown clock pair.")
    return -1


if __name__ == "__main__":
    print("Dark energy constraint", DE_max_amplitude())
    print("Modified gravity constraint:", MG_max_amplitude())
    w, A = DM_max_amplitude()
    print("Dark matter omega, A constraints:")
    print(w)
    print(A)

