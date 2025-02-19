###### Constants ######
Ts=3600*24*365

h0a = 3e-33
hm1a = 7.2e-37
Ka=0
h0b = 8e-28
hm1b = 2.8e-32
Kb=1

import numpy as np
import emcee

from scipy.fft import fft, fftfreq

########## Frequency domain ##########

def freqs_used(Ts=Ts,nyr=3.01,dt=1):
    T= Ts * nyr
    N = int(T / dt)
    freqs = fftfreq(N, dt)[0:N//2]
    fredd=freqs[freqs<5e-6]
    return fredd[1:]

########## MCMC ##########

def emcee_sampler(log_post, x0, nburn=300, steps=10000, **kwargs):

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
    def log_posterior(theta):
        res=np.abs(y-model(x,theta))/sigma_y
        return -2*np.sum(res**2)/2 + log_prior(theta)

    ndim=len(theta0)
    nwalkers=8*ndim
    p0=np.zeros((nwalkers,ndim))
    for i in range(nwalkers):
        p0[i]=theta0*np.random.normal(loc=1.0, scale=0.001, size=ndim)

    return emcee_sampler(log_posterior, p0, nburn=nburn, steps=steps)


########## Clocks Noise ##########

def h0_comb(h0a=h0a, h0b=h0b):
    return h0a + h0b
def hm1_comb(hm1a=hm1a, hm1b=hm1b):
    return hm1a + hm1b
def Ktot_comb(Ka=Ka, Kb=Kb):
    return abs(Ka - Kb)

clocks_pars_0=[h0a,hm1a,Ka,h0b,hm1b,Kb]

def noise_PSD(f,clocks_pars=clocks_pars_0):
    h0a=clocks_pars[0]
    hm1a=clocks_pars[1]
    Ka=clocks_pars[2]
    h0b=clocks_pars[3]
    hm1b=clocks_pars[4]
    Kb=clocks_pars[5]

    return (h0_comb(h0a=h0a, h0b=h0b) + hm1_comb(hm1a=hm1a, hm1b=hm1b) / f) / Ktot_comb(Ka=Ka, Kb=Kb)


def noise_generator(clocks_pars=clocks_pars_0,Ts= Ts,nyr=3.01,dt=1):
    
    fred=freqs_used(Ts=Ts,nyr=nyr,dt=dt)
    N = int(Ts*nyr / dt)

    return fred,np.sqrt(
        noise_PSD(fred,clocks_pars=clocks_pars) 
        * N) * np.exp(-2j * np.pi * np.random.uniform(size = len(fred)))


################## Theory models ##################

# Modified gravity

MOD_log_prior=lambda _:0.0

def MOD_model(f, pars = np.array([5e-16]),Ts=Ts,nyr=3.01,dt=1):
    A = pars[0]
    omega = 2*np.pi/Ts
    phi = 0.
    T=Ts*nyr

    return (1/dt) * A * np.exp(-2j * np.pi * f * T) * (omega * np.cos(omega * T + phi) - np.exp(2j * np.pi * f * T) * (omega * np.cos(phi) +
           2j * np.pi * f * np.sin(phi)) + 2j * np.pi * f * np.sin(omega * T + phi)) / (4 * np.pi * np.pi * f**2 - omega * omega)

def MOD_model_d(f, pars = np.array([5e-16]),Ts=Ts,nyr=3.01,dt=1):
    omega = 2*np.pi/Ts
    phi = 0.
    T=Ts*nyr

    return (1/dt) * np.exp(-2j * np.pi * f * T) * (omega * np.cos(omega * T + phi) - np.exp(2j * np.pi * f * T) * (omega * np.cos(phi) +
           2j * np.pi * f * np.sin(phi)) + 2j * np.pi * f * np.sin(omega * T + phi)) / (4 * np.pi * np.pi * f**2 - omega * omega)

# DE

DE_log_prior=lambda _:0.0

def DE_model(f,pars=np.array([1e-33]),Ts=Ts,nyr=3.01,dt=1):
    T=Ts*nyr
    a=pars[0]
    return -(1/dt)*a*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

def DE_model_d(f,pars=np.array([1e-33]),Ts=Ts,nyr=3.01,dt=1):
    T=Ts*nyr
    return -(1/dt)*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

# DM

def DM_log_prior(pars):
    if pars[1]<0 or pars[2]<0 or pars[2]>2*np.pi:
        return -np.inf
    return 0.

def DM_model(f,pars=np.array([2e-21,2*np.pi/(Ts),1.]),Ts=Ts,nyr=3.01,dt=1):
    A=pars[0]
    om=pars[1]
    phi=pars[2]
    T=Ts*nyr
    
    return (1/dt)*(A/om)*np.exp(-2j*np.pi*f*T)*(
        om*np.cos(om*T+phi)-np.exp(2j*np.pi*f*T)*(om*np.cos(phi)+2j*np.pi*f*np.sin(phi))+2j*np.pi*f*np.sin(om*T+phi)
        )/(4*np.pi*np.pi*f*f-om*om)

def DM_model_dA(f, pars = np.array([2e-21,2*np.pi/Ts,1.]),Ts=Ts,nyr=3.01,dt=1):
    A = pars[0]
    omega = pars[1]
    phi = pars[2]
    T=Ts*nyr
    
    return (1/dt)*(1/omega)*np.exp(-2j*np.pi*f*T)*(
        omega*np.cos(omega*T+phi)-np.exp(2j*np.pi*f*T)*(omega*np.cos(phi)+2j*np.pi*f*np.sin(phi))+2j*np.pi*f*np.sin(omega*T+phi)
        )/(4*np.pi*np.pi*f*f-omega*omega)

def DM_model_dom(f, pars = np.array([2e-21,2*np.pi/Ts,1.]),Ts=Ts,nyr=3.01,dt=1):
    A = pars[0]
    omega = pars[1]
    phi = pars[2]
    T=Ts*nyr
    
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

def DM_model_dphi(f, pars = np.array([2e-21,2*np.pi/Ts,1.]),Ts=Ts,nyr=3.01,dt=1):
    A = pars[0]
    omega = pars[1]
    phi = pars[2]
    T=Ts*nyr
    
    pi = np.pi
    I = 1j

    numerator = np.exp(-2 * I * f * pi * T) *(-2 * I * f * pi * np.cos(phi + T * omega) + 
             np.exp(2 * I * f * pi * T) * (2 * I * f * pi * np.cos(phi) - omega * np.sin(phi)) + 
             omega * np.sin(phi + T * omega))

    denominator = (-4 * f**2 * pi**2 * omega + omega**3)

    expression = numerator / denominator
    
    return (1/dt)*A*expression