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

import pandas as pd
import numpy as np
from stats import utils

PI = np.pi

def DE_max_amplitude(clock_pair):
    """Maximum dark energy signal that could have escaped detection.
        signal = A * t

    The signal is dimensionless, A has units of sec^-1
    Returns A for the given clock pair.
    Input:
        clock_pair: Name of a particular clock pair
    """
    file_path = "stats/sigmas.csv"
    
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")
    
    if clock_pair not in df.iloc[:, 0].values:
        raise ValueError(f"Clock name '{clock_pair}' not found in the file.")
    
    return np.array(df[df.iloc[:, 0] == clock_pair]['sigma_A_DE'])[0]

def DM_max_amplitude(clock_pair, nPoints = 100):
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
    T_year=3600*24*365 # year in seconds
    w_min = 2*PI / (3*T_year)     # 3 years
    w_max = 2*PI / 600     # 10 min

    w_ref = 2*PI / (3*T_year)     # year, reference value for sigmas.csv file
        
    # Values logarithmically spaced between w_min and w_max
    w_vals = np.logspace(np.log10(w_min),
                         np.log10(w_max), nPoints)
    
    file_path = "stats/sigmas.csv"
    
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")
    
    if clock_pair not in df.iloc[:, 0].values:
        raise ValueError(f"Clock name '{clock_pair}' not found in the file.")
    
    sigma_A_0=np.array(df[df.iloc[:, 0] == clock_pair]['sigma_A_DM'])[0]

    clock_couple=clock_pair.split('/')

    noise_pars_list=pd.read_csv('stats/clocks_pars.csv')

    row1 = noise_pars_list.loc[noise_pars_list['Clock_name'] == clock_couple[0]].iloc[0, 1:].tolist()
    row2 = noise_pars_list.loc[noise_pars_list['Clock_name'] == clock_couple[1]].iloc[0, 1:].tolist()
    noise_pars = row1 + row2

    # Fill in A_vals depending on the clock pair
    den=np.sqrt(utils.noise_PSD(w_ref/(2*np.pi),clocks_pars=noise_pars))
    PSD_vals=utils.noise_PSD(w_vals/(2*np.pi),clocks_pars=noise_pars)
    A_vals=sigma_A_0*(w_vals/w_ref)*np.sqrt(PSD_vals)/den
    
    return w_vals, A_vals

def MG_max_amplitude(clock_pair):
    """Maximum modified gravity signal that could have escaped detection.
        signal = A * cos(2 pi t / year)
    The signal and A are dimensionless.
    Input:
        clock_pair: Name of a particular clock pair

    Returns A for the given clock pair.
    """

    file_path = "stats/sigmas.csv"
    
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")
    
    if clock_pair not in df.iloc[:, 0].values:
        raise ValueError(f"Clock name '{clock_pair}' not found in the file.")
    
    return np.array(df[df.iloc[:, 0] == clock_pair]['sigma_A_MOD'])[0]


if __name__ == "__main__":
    print("Dark energy constraint:", DE_max_amplitude('N2+/Sr'))
    print("Modified gravity constraint:", MG_max_amplitude('N2+/Sr'))
    from matplotlib import pyplot as plt
    ws,As=DM_max_amplitude('CaF/Sr', nPoints = 100)
    ws2,As2=DM_max_amplitude('Cs/Sr', nPoints = 100)
    plt.loglog(ws,As,label='CaF/Sr')
    plt.loglog(ws2,As2,label='Cs/Sr')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'A')
    plt.legend()
    plt.show()
