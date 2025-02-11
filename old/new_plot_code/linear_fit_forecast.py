import constants
import utils

import numpy as np
from scipy.fft import fft, fftfreq
from scipy.integrate import simps

# d=a*t

################## MCMC ##################

def hFa(f,pars=np.array([1e-23]),Ts=constants.Ts,nyr=3.01,dt=1):
    T=Ts*nyr
    a=pars[0]
    return -(1/dt)*a*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

def df_gen(pars=np.array([1e-23]),clocks_pars=utils.clocks_pars_0,Ts= constants.Ts,nyr=3.01,dt=1):
    freqs,noise=utils.noise_generator(clocks_pars=clocks_pars,Ts= Ts,nyr=nyr,dt=dt)
    signal=hFa(freqs,pars=pars,Ts=Ts,nyr=nyr,dt=dt)
    sigmas=np.abs(noise)
    return freqs,signal+noise,sigmas

def DE_chain(pars=np.array([1e-23]),clocks_pars=utils.clocks_pars_0,Ts= constants.Ts,nyr=3.01,dt=1):
    
    freqs,data,sigmas=df_gen(pars=pars,clocks_pars=clocks_pars,Ts= Ts,nyr=nyr,dt=dt)
    
    def model(f,pars=np.array([1e-23])):
        return hFa(f,pars=pars,Ts=Ts,nyr=nyr,dt=dt)
    
    chains_fa = utils.big_sampler(freqs, data, sigmas, model, pars * 1.1)[1]
    samples = chains_fa.reshape(-1,chains_fa.shape[-1]).flatten()
    return samples


################## Fisher ##################

def dhFa(f,Ts=constants.Ts,nyr=3.01,dt=1):
    T=Ts*nyr
    return -(1/dt)*np.exp(-2j*np.pi*f*T)*(np.exp(2j*np.pi*f*T)-1-2j*np.pi*f*T)/(2*np.pi*f)**2

def Fisher(clocks_pars=utils.clocks_pars_0,Ts=constants.Ts,nyr=3.01,dt=1):

    fred=utils.freqs_used(Ts=Ts,nyr=nyr,dt=dt)

    integrand=2*(np.abs(dhFa(fred,Ts=Ts,nyr=nyr,dt=dt))**2)/utils.noise_PSD(fred,clocks_pars=clocks_pars)

    return simps(integrand,fred)


'''
print(np.sqrt(np.var(DE_chain(pars=np.array([1e-23])))))
print(1/np.sqrt(Fisher()))

import matplotlib.pyplot as plt
plt.loglog(utils.freqs_used(),np.abs(df_gen(pars=np.array([1e-23]))[1]))
plt.show()
print(np.abs(df_gen(pars=np.array([1e-23]))).shape)

plt.loglog(utils.freqs_used(),utils.noise_PSD(utils.freqs_used()))
plt.show()
'''

print(np.sqrt(np.var(DE_chain(pars=np.array([1e-23])))))
print(1/np.sqrt(Fisher()))