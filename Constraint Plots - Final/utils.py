import constants
import numpy as np
import emcee

from scipy.fft import fft, fftfreq

########## Frequency domain ##########

def freqs_used(Ts=constants.Ts,nyr=3.01,dt=1):
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

def h0_comb(h0a=constants.h0a, h0b=constants.h0b):
    return h0a + h0b
def hm1_comb(hm1a=constants.hm1a, hm1b=constants.hm1b):
    return hm1a + hm1b
def Ktot_comb(Ka=constants.Ka, Kb=constants.Kb):
    return abs(Ka - Kb)


def noise_PSD(f,h0a=constants.h0a,hm1a=constants.hm1a,Ka=constants.Ka,
              h0b=constants.h0b,hm1b=constants.hm1b, Kb=constants.Kb):

    return (h0_comb(h0a=h0a, h0b=h0b) + hm1_comb(hm1a=hm1a, hm1b=hm1b) / f) / Ktot_comb(Ka=Ka, Kb=Kb)


def noise_generator(h0a=constants.h0a,hm1a=constants.hm1a,Ka=constants.Ka,
              h0b=constants.h0b,hm1b=constants.hm1b, Kb=constants.Kb,Ts= constants.Ts,nyr=3.01,dt=1):
    
    fred=freqs_used(Ts=Ts,nyr=nyr,dt=dt)

    return fred,np.sqrt(
        noise_PSD(fred,h0a=h0a,hm1a=hm1a,Ka=Ka,h0b=h0b,hm1b=hm1b, Kb=Kb) 
        * np.exp(-2j * np.pi * np.random.uniform(size = len(fred)))
        )

########## Fisher ##########

def samples_from_Fisher(fid,Fish,N): #needs to be tested
    return np.random.multivariate_normal(fid, Fish, N).T