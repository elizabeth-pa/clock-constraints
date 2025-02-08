import constants
import utils

import numpy as np
from scipy.fft import fft, fftfreq
from scipy.integrate import simps

import matplotlib.pyplot as plt
import corner

noise_pars_list=np.array([
    [2.88e-30,5.26e-36,0.], # Yb+    0
    [4.61e-33,2.89e-36,0.], # Sr     1
    [6.48e-28,2.89e-32,1.], # Cs     2
    [4.5e-30,4.06e-35,0.5], # CaF    3
    [2.88e-28,1.1e-35,0.5], # N2+    4
    [1.8e-29,1.62e-38,0.], # Cf15+   5
    [1.25e-29,1.62e-38,0.], # Cf17+  6
])

names_list=[
    'Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs'
]

idxs_list=[
    [2,1], [4,1], [3,1], [3,2], [0,2], [3,0], [5,2]
]

def noise_pars_wrap(couple):
    return list(noise_pars_list[idxs_list[couple][0]])+list(noise_pars_list[idxs_list[couple][1]])

################## MCMC ##################

def data_gen(wich_model=utils.MOD_model,pars=np.array([5e-16]),clocks_pars=utils.clocks_pars_0,Ts= constants.Ts,nyr=3.01,dt=1):
    freqs,noise=utils.noise_generator(clocks_pars=clocks_pars,Ts=Ts,nyr=nyr,dt=dt)
    signal=wich_model(freqs,pars=pars,Ts=Ts,nyr=nyr,dt=dt)
    sigmas=np.abs(noise)
    return freqs,signal+noise,sigmas

def chain_gen(wich_model=utils.MOD_model,pars=np.array([5e-16]),clocks_pars=utils.clocks_pars_0,Ts= constants.Ts,nyr=3.01,dt=1,log_prior=lambda _:0.0):
    
    freqs,data,sigmas=data_gen(wich_model=wich_model,pars=pars,clocks_pars=clocks_pars,Ts= Ts,nyr=nyr,dt=dt)
    
    def model(f,parsmod):
        return wich_model(f,pars=parsmod,Ts=Ts,nyr=nyr,dt=dt)
    
    chains_fa = utils.big_sampler(freqs, data, sigmas, model, pars*1. ,log_prior=log_prior)[1]
    samples = chains_fa.reshape(-1,chains_fa.shape[-1])
    return samples

########## Fisher ##########

def Fisher(wich_dmodel=[utils.MOD_model_d],pars = np.array([5e-16]),clocks_pars=utils.clocks_pars_0,Ts=constants.Ts,nyr=3.01,dt=1):

    fred=utils.freqs_used(Ts=Ts,nyr=nyr,dt=dt)
    if len(pars)==1:
        integrand=2*(np.abs(wich_dmodel[0](fred,pars=pars,Ts=Ts,nyr=nyr,dt=dt))**2
                     /utils.noise_PSD(fred,clocks_pars=clocks_pars))

        return np.array([simps(integrand,fred)])
    else:
        def diag_F(h1,pars=pars,Ts=Ts,nyr=nyr,dt=dt):
            return simps(2*(np.abs(h1)**2)/(utils.noise_PSD(fred,clocks_pars=clocks_pars)),fred)
        
        def offdiag_F(h1,h2,pars=pars,Ts=Ts,nyr=nyr,dt=dt):
            return simps(2*(np.real(h1*np.conjugate(h2)))/(utils.noise_PSD(fred,clocks_pars=clocks_pars)),fred)
        
        FFop=np.zeros((len(pars),len(pars)))
        
        for i in range(len(pars)):
            FFop[i,i]=diag_F(wich_dmodel[i](fred,pars=pars,Ts=Ts,nyr=nyr,dt=dt),pars=pars,Ts=Ts,nyr=nyr,dt=dt)
            for j in range(i+1,len(pars)):
                FFop[i,j]=offdiag_F(wich_dmodel[i](fred,pars=pars,Ts=Ts,nyr=nyr,dt=dt),wich_dmodel[j](fred,pars=pars,Ts=Ts,nyr=nyr,dt=dt),pars=pars,Ts=Ts,nyr=nyr,dt=dt)
                FFop[j,i]=FFop[i,j]

        return FFop


def samps_from_Fish(fish,mean,size=10000):
    resc=1/mean[0]
    if len(mean)==1:
        return np.random.multivariate_normal(mean, 1/fish[None], size=size)
    else:
        mean_resc=mean*1.
        mean_resc[0]=mean_resc[0]*resc
        fish_resc=fish*1.
        fish_resc[0,:]=fish_resc[0,:]/resc
        fish_resc[:,0]=fish_resc[:,0]/resc
        samps=np.random.multivariate_normal(mean_resc, np.linalg.inv(fish_resc), size=size)
        samps[:,0]=samps[:,0]/resc
        return samps

########## Plotters ##########

def plotter_posterior(wich_model=utils.MOD_model,pars=np.array([5e-16]),
            wich_dmodel=[utils.MOD_model_d],
            clocks_pars=utils.clocks_pars_0,Ts=constants.Ts,nyr=3.01,dt=1,log_prior=lambda _:0.0, labels=['A']):
    
    if len(wich_dmodel)!=len(pars) or len(wich_dmodel)!=len(labels):
        raise Exception("Be careful with the dimensions of pars, derivatives, labels")
    
    samples_MCMC=chain_gen(wich_model=wich_model,pars=pars,
                         clocks_pars=clocks_pars,Ts=Ts,nyr=nyr,dt=dt,log_prior=log_prior)
    
    samples_rescaled=samples_MCMC-np.mean(samples_MCMC,axis=0)+pars

    samples_Fisher=samps_from_Fish(Fisher(wich_dmodel=wich_dmodel, pars = pars,clocks_pars=clocks_pars,Ts=Ts,nyr=nyr,dt=dt),pars,size=len(samples_rescaled[:,0]))
    

    if len(pars)==1:
        plt.hist(samples_rescaled,alpha=0.5,bins=50,label='MCMC')
        plt.hist(samples_Fisher,color='r',alpha=0.5,bins=50,label='Fisher')
        plt.axvline(x=pars[0],c='k')
        plt.xlabel(labels[0])
        plt.legend()
        plt.show()
        print('MCMC stdev on '+labels+':',np.sqrt(np.var(samples_MCMC)),'s^-1')
        print('Fisher stdev on '+labels+':',np.sqrt(np.var(samples_Fisher)),'s^-1')
    else:
        my_figureomphi = corner.corner(samples_Fisher, labels=labels, show_titles=True,truths=pars,color='b')
        corner.corner(samples_rescaled,fig=my_figureomphi)
        plt.show()
    
    return 


def plotter_wrapper(which='MOD'):
    if which=='MOD':
        pars_ex=np.array([5e-16])
        plotter_posterior(wich_model=utils.MOD_model,wich_dmodel=[utils.MOD_model_d],pars = pars_ex, labels=['A'])
    elif which=='DE':
        pars_ex=np.array([1e-23])
        plotter_posterior(wich_model=utils.DE_model,wich_dmodel=[utils.DE_model_d],pars = pars_ex, labels=['a'])
    elif which=='DM':
        pars_ex=np.array([2e-21,2*np.pi/constants.Ts,1.])
        plotter_posterior(wich_model=utils.DM_model,wich_dmodel=[utils.DM_model_dA,utils.DM_model_dom,utils.DM_model_dphi],pars = pars_ex, labels=['A','om','phi'])
    return

########## Table producer ##########

def sigmas_table(which_clocks=names_list,save=False):
    tab_out=np.zeros((len(which_clocks),3))

    for i in range(len(which_clocks)):
        samps_MOD=chain_gen(wich_model=utils.MOD_model,pars=np.array([5e-16]),clocks_pars=noise_pars_wrap(i))
        samps_DE=chain_gen(wich_model=utils.DE_model,pars=np.array([1e-23]),clocks_pars=noise_pars_wrap(i))
        samps_DM=chain_gen(wich_model=utils.DM_model,pars=np.array([2e-21,2*np.pi/constants.Ts,1.]),clocks_pars=noise_pars_wrap(i),log_prior=utils.DM_log_prior)
        tab_out[i,0]=np.sqrt(np.var(samps_MOD[:,0]))
        tab_out[i,1]=np.sqrt(np.var(samps_DE[:,0]))
        tab_out[i,2]=np.sqrt(np.var(samps_DM[:,0]))
    if save==True:
        np.savetxt('sigmas.csv',tab_out,delimiter=",")
    return tab_out