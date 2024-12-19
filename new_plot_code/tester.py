import constants
import utils
import Data_to_model as dtm
import matplotlib.pyplot as plt

import numpy as np


couples = np.genfromtxt("sigmas.csv", delimiter=",", dtype=float)[:,2]

om0=2*np.pi/(3*constants.Ts)
om1=2*np.pi/(10)
om_range=np.logspace(np.log10(om0),np.log10(om1))

def A_thr(om,idx,sigmas=couples):
    omref=2*np.pi/constants.Ts
    den=np.sqrt(utils.noise_PSD(omref/(2*np.pi),clocks_pars=dtm.noise_pars_wrap(idx)))
    return sigmas[idx]*(om/omref)*np.sqrt(utils.noise_PSD(om/(2*np.pi),clocks_pars=dtm.noise_pars_wrap(idx)))/den

plt.loglog(om_range,A_thr(om_range,0,sigmas=couples),label=dtm.names_list[0])
plt.loglog(om_range,A_thr(om_range,2,sigmas=couples),label=dtm.names_list[2])
plt.axvline(x=2*np.pi/(constants.Ts),label=r'$\text{yr}^{-1}$',color='k')
plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$A_{\rm thr}$ [Hz]')
plt.legend()
plt.grid()
plt.show()