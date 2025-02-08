import constants
import utils

import numpy as np
from scipy.fft import fft, fftfreq
from scipy.integrate import simps

# d=a*t

################## MCMC ##################

def hFa(f,pars=np.array([1e-33]),Ts=constants.Ts,nyr=3.01,dt=1):
    T=Ts*nyr
    a=pars[0]
    return -(1/dt)*a*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

def hFaMCMC(f,pars=np.array([1e-33])):
    Ts=constants.Ts,nyr=3.01,dt=1
    T=Ts*nyr
    a=pars[0]
    return -(1/dt)*a*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

def df_gen(pars=np.array([1e-33]),h0a=constants.h0a,hm1a=constants.hm1a,Ka=constants.Ka,
              h0b=constants.h0b,hm1b=constants.hm1b, Kb=constants.Kb,Ts= constants.Ts,nyr=3.01,dt=1):
    freqs,noise=utils.noise_generator(h0a=h0a,hm1a=hm1a,Ka=Ka,
              h0b=h0b,hm1b=hm1b, Kb=Kb,Ts= Ts,nyr=nyr,dt=dt)
    signal=hFa(freqs,pars=pars,Ts=Ts,nyr=nyr,dt=dt)
    sigmas=np.abs(noise)
    return freqs,signal+noise,sigmas

def DE_chain(pars=np.array([1e-33]),h0a=constants.h0a,hm1a=constants.hm1a,Ka=constants.Ka,
              h0b=constants.h0b,hm1b=constants.hm1b, Kb=constants.Kb,Ts= constants.Ts,nyr=3.01,dt=1):
    
    freqs,data,sigmas=df_gen(pars=pars,h0a=h0a,hm1a=hm1a,Ka=Ka,
              h0b=h0b,hm1b=hm1b, Kb=Kb,Ts= Ts,nyr=nyr,dt=dt)
    
    chains_fa = utils.big_sampler(freqs, data, sigmas, hFaMCMC, pars * 1.1)[1]
    samples = chains_fa.reshape(-1,chains_fa.shape[-1]).flatten()
    return samples


################## Fisher ##################

def dhFa(f,Ts=constants.Ts,nyr=3.01,dt=1):
    T=Ts*nyr
    return -(1/dt)*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

def Fisher(h0a=constants.h0a,hm1a=constants.hm1a,Ka=constants.Ka,
              h0b=constants.h0b,hm1b=constants.hm1b, Kb=constants.Kb,Ts=constants.Ts,nyr=3.01,dt=1):

    fred=utils.freqs_used(Ts=Ts,nyr=nyr,dt=dt)

    integrand=2*(np.abs(dhFa(fred,Ts=Ts,nyr=nyr,dt=dt))**2)/utils.noise_PSD(fred,h0a=h0a,hm1a=hm1a,Ka=Ka,h0b=h0b,hm1b=hm1b, Kb=Kb)

    return simps(integrand,fred)

print(1/np.sqrt(Fisher()))
#print(np.sqrt(np.var(DE_chain(pars=np.array([1e-33])))))
print("hello")