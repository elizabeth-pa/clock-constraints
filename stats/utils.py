###### Constants ######

Ts = 3600*24*365 # year in seconds
T_obs_yr = 3.01 # default number of years of observation
dt = 1 # default sampling time in seconds

'''
Default noise clock couple values for h0, h-1 and K
'''

h0a = 3e-33
hm1a = 7.2e-37
Ka=0
h0b = 8e-28
hm1b = 2.8e-32
Kb=1

import numpy as np
import emcee

from scipy.fft import fftfreq

########## Frequencies array for freq-domain analysis ##########

def freqs_used(Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    T= Ts * T_obs_yr
    N = int(T / dt)
    freqs = fftfreq(N, dt)[0:N//2]
    fredd=freqs[freqs<5e-6] # low frequency for data compression and faster analysis
                            # change it if data has high frequency information
    return fredd[1:] # avoids f=0, which is irrelevant

########## MCMC ##########

def emcee_sampler(log_post, x0, nburn=300, steps=10000, **kwargs):
    
    '''
    Wrapper for emcee MCMC
    Inputs:
        - log_post: one parameter function (see the usage in big_sampler below)
        - x0: the initial proposal for the MCMC chain
    '''

    mnwalkers, mndim = x0.shape

    sampler = emcee.EnsembleSampler(mnwalkers,
                                    mndim,
                                    log_post)
    sampler.run_mcmc(x0, nburn, progress=True,tune=True)
    state = sampler.get_chain()[-1, :, :]
    sampler.reset()
    print('Finished initial run, burn-in dropped and starting real run')
    sampler.run_mcmc(state, steps, progress=True,tune=True)

    return sampler.get_log_prob(), sampler.get_chain()

def big_sampler(x,y,sigma_y,model, theta0, log_prior=lambda _:0.0, nburn=300, steps=3000):

    '''
    Wrapper for emcee MCMC to a posterior distribution with gaussian likelihood
    Inputs:
        - log_post: one parameter function (see the usage in big_sampler below)
        - x: array of x sampled values (e.g. frequencies)
        - y: array data samples values, datastream
        - sigma_y: array of the errors on the measured datastream y
        - model: one parameter function, the model where data come from
            i.e. y = model + noise
            where noise follows gaussian statistics with errors sigma_y
        - theta0: initial proposal of the parameters sampled by MCMC
    '''

    def log_posterior(theta): # builds the gaussian likelihood function in the noise
        res=np.abs(y-model(x,theta))/sigma_y
        return -2*np.sum(res**2)/2 + log_prior(theta)

    ndim=len(theta0)
    nwalkers=8*ndim
    p0=np.zeros((nwalkers,ndim))
    for i in range(nwalkers):
        p0[i]=theta0*np.random.normal(loc=1.0, scale=0.001, size=ndim)
        # sets a close enough initial guess for MCMC for convenience 

    return emcee_sampler(log_posterior, p0, nburn=nburn, steps=steps)


########## Clocks Noise ##########

'''
Being the noise in the two clocks independent, the overall noise Power Spectral Density
    for the datastream is just the sum of the ones of the two individual clocks
'''

def h0_comb(h0a=h0a, h0b=h0b):
    return h0a + h0b
def hm1_comb(hm1a=hm1a, hm1b=hm1b):
    return hm1a + hm1b
def Ktot_comb(Ka=Ka, Kb=Kb):
    return abs(Ka - Kb)

clocks_pars_0=[h0a,hm1a,Ka,h0b,hm1b,Kb]

def noise_PSD(f,clocks_pars=clocks_pars_0):
    '''
    Noise PSD for a clock pair.
    Inputs:
        - f: array of frequencies
        - clocks_pars array of 3 noise parameters for each of the two clocks
    '''
    h0a=clocks_pars[0]
    hm1a=clocks_pars[1]
    Ka=clocks_pars[2]
    h0b=clocks_pars[3]
    hm1b=clocks_pars[4]
    Kb=clocks_pars[5]

    return (h0_comb(h0a=h0a, h0b=h0b) + hm1_comb(hm1a=hm1a, hm1b=hm1b) / f) / Ktot_comb(Ka=Ka, Kb=Kb)


def noise_generator(clocks_pars=clocks_pars_0,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    Returns the array of frequencies and a realization of the noise
        given a pair of clocks with noise parameters clocks_pars
    '''
    fred=freqs_used(Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
    N = int(Ts*T_obs_yr / dt)

    return fred,np.sqrt(
        noise_PSD(fred,clocks_pars=clocks_pars) 
        * N) * np.exp(-2j * np.pi * np.random.uniform(size = len(fred)))


################## Phenomenological models ##################

# Modified gravity

MOD_log_prior=lambda _:0.0

def MOD_model(f, pars = np.array([5e-16]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    Modified-Gravity-inspired sinusoidal model d=A*sin(omega*t+phi) in the frequency domain
    pars is the only free parameter, A
    '''
    A = pars[0]
    omega = 2*np.pi/Ts
    phi = 0.
    T=Ts*T_obs_yr

    numerator = np.exp(-2j * np.pi * f * T) * (omega * np.cos(omega * T + phi) - np.exp(2j * np.pi * f * T) * (omega * np.cos(phi) +
           2j * np.pi * f * np.sin(phi)) + 2j * np.pi * f * np.sin(omega * T + phi))
    denominator = (4 * np.pi * np.pi * f**2 - omega * omega)

    return (1/dt) * A * numerator / denominator

def MOD_model_d(f, pars = np.array([5e-16]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    First derivative of the Modified-Gravity-inspired sinusoidal model MOD_model
        with respect to A, in the frequency domain
        pars is the only free parameter, A
    '''
    omega = 2*np.pi/Ts
    phi = 0.
    T=Ts*T_obs_yr

    numerator = np.exp(-2j * np.pi * f * T) * (omega * np.cos(omega * T + phi) - np.exp(2j * np.pi * f * T) * (omega * np.cos(phi) +
           2j * np.pi * f * np.sin(phi)) + 2j * np.pi * f * np.sin(omega * T + phi))
    denominator = (4 * np.pi * np.pi * f**2 - omega * omega)

    return (1/dt) * numerator / denominator

# DE

DE_log_prior=lambda _:0.0

def DE_model(f,pars=np.array([1e-33]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    Dark-Energy-inspired linear model d=a*t in the frequency domain
    pars is the only free parameter, a
    '''
    T=Ts*T_obs_yr
    a=pars[0]

    numerator = np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)
    denominator = (2*np.pi*f)**2

    return -(1/dt)*a*numerator/denominator

def DE_model_d(f,pars=np.array([1e-33]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    First derivative of the Dark-Energy-inspired linear model DE_model
        with respect to a, in the frequency domain
        pars is the only free parameter, a
    '''
    T=Ts*T_obs_yr
    numerator = np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)
    denominator = (2*np.pi*f)**2

    return -(1/dt) * numerator / denominator

# DM

def DM_log_prior(pars):
    if pars[1]<0 or pars[2]<0 or pars[2]>2*np.pi:
        return -np.inf
    return 0.

def DM_model(f,pars=np.array([2e-21,2*np.pi/(Ts),1.]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    Dark-Matter-inspired sinusoidal model d=(A/omega)*sin(omega*t+phi) in the frequency domain
    pars is the array of the three free parameters, [A,omega,phi]
    '''

    A=pars[0]
    om=pars[1]
    phi=pars[2]
    T=Ts*T_obs_yr

    numerator = np.exp(-2j*np.pi*f*T)*(
        om*np.cos(om*T+phi)-np.exp(2j*np.pi*f*T)*(om*np.cos(phi)+2j*np.pi*f*np.sin(phi))+2j*np.pi*f*np.sin(om*T+phi))
    denominator = (4*np.pi*np.pi*f*f-om*om)

    return (1/dt)*(A/om) * numerator / denominator

def DM_model_dA(f, pars = np.array([2e-21,2*np.pi/Ts,1.]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    First derivative of the Dark-Matter-inspired sinusoidal model MOD_model
        with respect to A, in the frequency domain
        pars is the array of the three free parameters, [A,omega,phi] at fiducial value
    '''

    A = pars[0]
    omega = pars[1]
    phi = pars[2]
    T=Ts*T_obs_yr

    numerator = np.exp(-2j*np.pi*f*T)*(
        omega*np.cos(omega*T+phi)-np.exp(2j*np.pi*f*T)*(omega*np.cos(phi)+2j*np.pi*f*np.sin(phi))+2j*np.pi*f*np.sin(omega*T+phi))
    denominator = (4*np.pi*np.pi*f*f-omega*omega)
    
    return (1/dt)*(1/omega)* numerator / denominator

def DM_model_dom(f, pars = np.array([2e-21,2*np.pi/Ts,1.]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    First derivative of the Dark-Matter-inspired sinusoidal model MOD_model
        with respect to omega, in the frequency domain
        pars is the array of the three free parameters, [A,omega,phi] at fiducial value
    '''

    A = pars[0]
    omega = pars[1]
    phi = pars[2]
    T=Ts*T_obs_yr
    
    pi = np.pi
    I = 1j

    exp1 = -4 * f**2 * pi**2 * omega + omega**3
    exp2 = -4 * f**2 * pi**2 + 3 * omega**2

    term1 = -exp2 * (-omega * np.cos(phi + T * omega) +
                    np.exp(2 * I * f * pi * T) * (omega * np.cos(phi) + 2 * I * f * pi * np.sin(phi)) -
                    2 * I * f * pi * np.sin(phi + T * omega))

    term2 = exp1 * (np.exp(2 * I * f * pi * T) * np.cos(phi) +
                    (-1 - 2 * I * f * pi * T) * np.cos(phi + T * omega) +
                    T * omega * np.sin(phi + T * omega))

    numerator = term1 + term2
    denominator = (exp1)**2

    expression = np.exp(-2 * I * f * pi * T) * numerator / denominator
    
    return (1/dt)*A*expression

def DM_model_dphi(f, pars = np.array([2e-21,2*np.pi/Ts,1.]),Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
    '''
    First derivative of the Dark-Matter-inspired sinusoidal model MOD_model
        with respect to phi, in the frequency domain
        pars is the array of the three free parameters, [A,omega,phi] at fiducial value
    '''

    A = pars[0]
    omega = pars[1]
    phi = pars[2]
    T=Ts*T_obs_yr
    
    pi = np.pi
    I = 1j

    numerator = np.exp(-2 * I * f * pi * T) *(-2 * I * f * pi * np.cos(phi + T * omega) + 
             np.exp(2 * I * f * pi * T) * (2 * I * f * pi * np.cos(phi) - omega * np.sin(phi)) + 
             omega * np.sin(phi + T * omega))

    denominator = (-4 * f**2 * pi**2 * omega + omega**3)

    expression = numerator / denominator
    
    return (1/dt)*A*expression