import utils

import numpy as np
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.integrate import simpson as simps

import matplotlib.pyplot as plt
import corner

noise_pars_list=pd.read_csv('stats/clocks_pars.csv')

couples_list=np.array([ # default clock couples considered
    # IMPORTANT: these names must be in agreement with the values used in the clocks_pars.csv file 
    ['Cs','Sr'],['N2+','Sr'],['CaF','Sr'],['CaF','Cs'],['Yb+','Cs'],['CaF','Yb+']
])


def noise_pars_wrap(idx,which_clocks=couples_list):
    '''
    Wrapper that given the idx term in the list of couples returns the 6 parameters for the noise model
    i.e. h0,h-1,K for each of the two clocks, concatenated in a list
    '''
    clock1=which_clocks[idx,0]
    clock2=which_clocks[idx,1]
    row1 = noise_pars_list.loc[noise_pars_list['Clock_name'] == clock1].iloc[0, 1:].tolist()
    row2 = noise_pars_list.loc[noise_pars_list['Clock_name'] == clock2].iloc[0, 1:].tolist()
    return row1 + row2

################## MCMC ##################

def data_gen(which_model=utils.MOD_model,pars=np.array([5e-16]),clocks_pars=utils.clocks_pars_0,Ts=utils.Ts,T_obs_yr=utils.T_obs_yr,dt=utils.dt):
    '''
    Generates datastream for a chosen phenomenological model and some fiducial value of the parameters
        in the frequency domain for a chosen pair of clocks
    Inputs:
        - which_model: the model considered, as a function of frequency and parameters
        - pars: array of the fiducial parameters for the model forecast
        - clocks_pars: array of the clock noise parameters
    '''

    freqs,noise=utils.noise_generator(clocks_pars=clocks_pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
    signal=which_model(freqs,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
    sigmas=np.abs(noise)
    return freqs,signal+noise,sigmas

def chain_gen(which_model=utils.MOD_model,pars=np.array([5e-16]),clocks_pars=utils.clocks_pars_0,Ts=utils.Ts,T_obs_yr=utils.T_obs_yr,dt=utils.dt,log_prior=lambda _:0.0):
    '''
    Generates posterior chain samples with MCMC on the phenomenological parameters for a given model
        given the fiducial value of the parameters for a chosen pair of clocks
    Inputs:
        - which_model: the model considered, as a function of frequency and parameters
        - pars: array of the fiducial parameters for the model forecast
        - clocks_pars: array of the clock noise parameters
    '''

    
    freqs,data,sigmas=data_gen(which_model=which_model,pars=pars,clocks_pars=clocks_pars,Ts= Ts,T_obs_yr=T_obs_yr,dt=dt)
    
    def model(f,parsmod):
        return which_model(f,pars=parsmod,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
    
    chains_fa = utils.big_sampler(freqs, data, sigmas, model, pars*1. ,log_prior=log_prior)[1]
    samples = chains_fa.reshape(-1,chains_fa.shape[-1])
    return samples

########## Fisher ##########

def Fisher(which_dmodel=[utils.MOD_model_d],pars = np.array([5e-16]),clocks_pars=utils.clocks_pars_0,Ts=utils.Ts,T_obs_yr=utils.T_obs_yr,dt=utils.dt):
    '''
    Computes the Fisher matrix of the phenomenological parameters for a given model
        given the fiducial value of the parameters for a chosen pair of clocks
    Inputs:
        - which_dmodel: the derivatives with respect to the phenomenological parameters
            as a LIST of the partial derivatives of the model
        - pars: array of the fiducial parameters for the model forecast
        - clocks_pars: array of the clock noise parameters
    '''
    fred=utils.freqs_used(Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
    if len(pars)==1:
        integrand=2*(np.abs(which_dmodel[0](fred,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt))**2
                     /utils.noise_PSD(fred,clocks_pars=clocks_pars))

        return np.array([simps(integrand,fred)])
    else:
        def diag_F(h1,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
            return simps(2*(np.abs(h1)**2)/(utils.noise_PSD(fred,clocks_pars=clocks_pars)),fred)
        
        def offdiag_F(h1,h2,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt):
            return simps(2*(np.real(h1*np.conjugate(h2)))/(utils.noise_PSD(fred,clocks_pars=clocks_pars)),fred)
        
        FFop=np.zeros((len(pars),len(pars)))
        
        for i in range(len(pars)):
            FFop[i,i]=diag_F(which_dmodel[i](fred,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt),pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
            for j in range(i+1,len(pars)):
                FFop[i,j]=offdiag_F(which_dmodel[i](fred,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt),which_dmodel[j](fred,pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt),pars=pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt)
                FFop[j,i]=FFop[i,j]

        return FFop


def samps_from_Fish(fish,mean,size=10000):
    '''
    Wrapper that produces samples from a multivariate Gaussian distribution with a mean value mean
    and a Covariance matrix C as the inverse of the Fisher matrix as input: C=fish^-1
    '''

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

def plotter_posterior(which_model=utils.MOD_model,pars=np.array([5e-16]),
            which_dmodel=[utils.MOD_model_d],
            clocks_pars=utils.clocks_pars_0,Ts=utils.Ts,T_obs_yr=utils.T_obs_yr,dt=utils.dt,log_prior=lambda _:0.0, labels=['A']):
    '''
    Plotter of the posterior distribution with MCMC and with Fisher matrix approach

    Inputs:
        - which_model: the model considered, as a function of frequency and parameters
        - pars: array of the fiducial parameters for the model forecast
        - clocks_pars: array of the clock noise parameters
    '''

    if len(which_dmodel)!=len(pars) or len(which_dmodel)!=len(labels):
        raise Exception("Be careful with the dimensions of pars, derivatives, labels")
    
    samples_MCMC=chain_gen(which_model=which_model,pars=pars,
                         clocks_pars=clocks_pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt,log_prior=log_prior)
    
    samples_rescaled=samples_MCMC-np.mean(samples_MCMC,axis=0)+pars

    samples_Fisher=samps_from_Fish(Fisher(which_dmodel=which_dmodel, pars = pars,clocks_pars=clocks_pars,Ts=Ts,T_obs_yr=T_obs_yr,dt=dt),pars,size=len(samples_rescaled[:,0]))
    

    if len(pars)==1:
        plt.hist(samples_rescaled,alpha=0.5,bins=50,label='MCMC')
        plt.hist(samples_Fisher,color='r',alpha=0.5,bins=50,label='Fisher')
        plt.axvline(x=pars[0],c='k')
        plt.xlabel(labels[0])
        plt.legend()
        plt.show()
        print('MCMC stdev on '+labels[0]+':',np.sqrt(np.var(samples_MCMC)),'s^-1')
        print('Fisher stdev on '+labels[0]+':',np.sqrt(np.var(samples_Fisher)),'s^-1')
    else:
        my_figureomphi = corner.corner(samples_Fisher, labels=labels, show_titles=True,truths=pars,color='b')
        corner.corner(samples_rescaled,fig=my_figureomphi)
        plt.show()
    
    return 


def plotter_wrapper(which='MOD'):
    '''
    Wrapper for the plotter function above plotter_posterior, for the three default models considered
    '''

    if which=='MOD':
        pars_ex=np.array([5e-16])
        plotter_posterior(which_model=utils.MOD_model,which_dmodel=[utils.MOD_model_d],pars = pars_ex, labels=['A'])
    elif which=='DE':
        pars_ex=np.array([1e-23])
        plotter_posterior(which_model=utils.DE_model,which_dmodel=[utils.DE_model_d],pars = pars_ex, labels=['a'])
    elif which=='DM':
        pars_ex=np.array([2e-21,2*np.pi/utils.Ts,1.])
        plotter_posterior(which_model=utils.DM_model,which_dmodel=[utils.DM_model_dA,utils.DM_model_dom,utils.DM_model_dphi],pars = pars_ex, labels=['A','om','phi'])
    return

########## Forecast errors ##########

def sigmas_table(which_clocks=couples_list):
    '''
    Produces the table of the forecast errors sigma_A, for the three models considered on the parameter A
        - Mod: d=A*sin(omega*t)
        - DE: d=A*t
        - DM: d=(A/omega)*sin(omega*t+phi)
    These values are obtained by taking the variance of the samples of A (which corresponds
        to the standard deviation of the marginalized posterior distribution of that parameter)
    A total observation time of T_obs=3yr and a sampling rate of dt=1s are assumed
    '''

    tab_out=np.zeros((len(which_clocks),3))

    for i in range(len(which_clocks)):
        samps_MOD=chain_gen(which_model=utils.MOD_model,pars=np.array([5e-16]),clocks_pars=noise_pars_wrap(i,which_clocks=which_clocks))
        samps_DE=chain_gen(which_model=utils.DE_model,pars=np.array([1e-23]),clocks_pars=noise_pars_wrap(i,which_clocks=which_clocks))
        samps_DM=chain_gen(which_model=utils.DM_model,pars=np.array([2e-21,2*np.pi/utils.Ts,1.]),clocks_pars=noise_pars_wrap(i,which_clocks=which_clocks),log_prior=utils.DM_log_prior)
        tab_out[i,0]=np.sqrt(np.var(samps_MOD[:,0]))
        tab_out[i,1]=np.sqrt(np.var(samps_DE[:,0]))
        tab_out[i,2]=np.sqrt(np.var(samps_DM[:,0]))
    return tab_out

def save_sigmas_to_csv(col_names=['sigma_A_MOD','sigma_A_DE','sigma_A_DM'], which_clocks=couples_list):
    '''
    Saves the value computed by sigmas_table of the forecast errors in the table sigmas.csv
    '''

    row_names=[]
    for i in range(len(which_clocks)):
        row_names.append(which_clocks[i,0]+'/'+which_clocks[i,1])

    tab_out=sigmas_table(which_clocks=which_clocks)

    df = pd.DataFrame(tab_out, index=row_names, columns=col_names)
    df.index.name = 'Clock_pair'
    
    df.to_csv('sigmas.csv')
    return 

if __name__ == "__main__":
    # save_sigmas_to_csv(which_clocks=couples_list)
    plotter_wrapper(which='DM')