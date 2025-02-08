"""
mu_constraints.py
"""

def DE_max_amplitude(clock_pair = "CaF/Sr"):
    match clock_pair:
        case "CaF/Sr":
            return 1.7e-25
        case "Sr/Cs":
            return 3.1e-24

    raise Exception("Unknown clock pair.")
    return -1

def DM_max_amplitude(w = 1, clock_pair = "CaF/Sr"):
    """
    Takes in a clock pair, as well as an angular frequency w = 2*pi*freq
    Returns a maximum amplitude for the signal at that frequency.
    """
    match clock_pair:
        case "CaF/Sr":
            return 1
        case "Sr/Cs":
            return 0

    raise Exception("Unknown clock pair.")
    return -1

def MG_max_amplitude(clock_pair = "CaF/Sr"):
    match clock_pair:
        case "CaF/Sr":
            return 7.8e-18
        case "Sr/Cs":
            return 1.4e-16

    raise Exception("Unknown clock pair.")
    return -1
