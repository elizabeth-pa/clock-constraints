import constants
import utils
import Data_to_model as dtm

from constants import hbar, rhoDM, Mpl

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

couples = np.genfromtxt("sigmas.csv", delimiter=",", dtype=float)[:,2]
# Computed for the pairs of clocks in Data_to_model and therefore for T_obs=10yr

colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


idx_pair=2
idx_pair_bad=0

def A_thr(om,idx,sigmas=couples):
    om0=2*np.pi/constants.Ts
    den=np.sqrt(utils.noise_PSD(om0/(2*np.pi),clocks_pars=dtm.noise_pars_wrap(idx)))
    return sigmas[idx]*(om/om0)*np.sqrt(utils.noise_PSD(om/(2*np.pi),clocks_pars=dtm.noise_pars_wrap(idx)))/den


def Mmass(a): # in units of Mpl
    return np.sqrt(2*constants.rhoDM)/(constants.hbar*a*constants.Mpl)

def Mmass_of_m(m,idx,sigmas=couples):
    return Mmass(A_thr(m/constants.hbar,idx,sigmas=sigmas))


''' MICROSCOPE '''

def Microscope_limits():
    def MicroscopeConstraint(M, Me):
        QE=constants.Mearth * (1 / M + constants.epsilon_Fe / Me) / (1 + constants.epsilon_Fe)
        dPhiMassless=QE/ (4 * np.pi * constants.Rmic**2)
        return (8*np.pi * constants.Mpl**2 * constants.Rmic**2 / constants.Mearth
                )  * (1/M - 1/Me) * (constants.epsilon_Pt - constants.epsilon_Ti)*dPhiMassless

    def lowcurve(M,Me):
        return MicroscopeConstraint(M, Me)-constants.eta_microscope
    def upcurve(M,Me):
        return MicroscopeConstraint(M, Me)+constants.eta_microscope

    Mrange = np.logspace(1,10,100)*constants.Mpl
    MeMrange = []
    MMrange = []

    for Mx in Mrange:
        def func0(Me):
            return lowcurve(Mx,Me)
        def func1(M):
            return upcurve(M,Mx)
        
        MeMrange.append(fsolve(func0, 1e30)[0])
        MMrange.append(fsolve(func1, 1e30)[0])
        
    line1x=np.log10(Mrange/constants.Mpl)
    line1y=np.log10(np.array(MeMrange)/constants.Mpl)
    line2x=np.log10(np.array(MMrange)/constants.Mpl)
    line2y=np.log10(Mrange/constants.Mpl)

    Mmicr_lim=line2x[-1]
    Mmicr_lim2=line1y[-1]

    return Mmicr_lim, Mmicr_lim2


''' Sherrill '''


def Sherrill_curve():

    sigma_r = 1.6e-13 # which comes from 2302.04565, below equation (15)
    h0_r = 2 * sigma_r**2 # the correct convertion to frequency domain white noise
    sigma_A_over_om = np.sqrt(2*h0_r/(3.01*constants.Ts)) # the Fisher forecast error for white noise over a total observation time T_max

    mminS=constants.hbar/(80000)
    mmaxS=constants.hbar/(600)

    mplot_range=np.logspace(np.log10(mminS),np.log10(mmaxS))

    xvals = mplot_range
    yvals = Mmass(sigma_A_over_om*mplot_range/constants.hbar)
    # Make it downturn at each end
    xvals = np.concatenate(([xvals[0]], xvals, [xvals[-1]]))
    yvals = np.concatenate(([0], yvals, [0]))

    return xvals, yvals


''' Clocks '''


def Clocks_curves(idx_pair=2):
    m3yr=constants.hbar/(3.01*constants.Ts)
    m10min=constants.hbar/(60*10)

    mplot_range=np.logspace(np.log10(m3yr),np.log10(m10min))

    ## CaF/Sr clock 
    xvals = mplot_range
    yvals = Mmass_of_m(mplot_range,idx_pair)
    # Make it downturn at each end
    xvals = np.concatenate(([xvals[0]], xvals, [xvals[-1]]))
    yvals = np.concatenate(([0], yvals, [0]))

    return xvals, yvals

def DM_plotter(idx_pair=2, idx_pair_bad=0, colorss = [colorcycle[i] for i in range(3)], savefig=False, showfig=True):

    # Square aspect ratio
    plt.figure(figsize=(5,5))

    ## Clocks
    xvals, yvals, = Clocks_curves(idx_pair=idx_pair)
    xvals2, yvals2 = Clocks_curves(idx_pair=idx_pair_bad)
    xvalsS, yvalsS = Sherrill_curve()

    plt.loglog(xvals, yvals, label=dtm.names_list[idx_pair], color=colorss[0])
    plt.fill_between(xvals, yvals, 0, color=colorss[0], alpha=0.25)

    plt.loglog(xvals2, yvals2, label=dtm.names_list[idx_pair_bad], color=colorss[1])
    plt.fill_between(xvals2, yvals2, 0, color=colorss[1], alpha=0.25)

    plt.loglog(xvalsS, yvalsS, label='Sherrill', color=colorss[2])
    plt.fill_between(xvalsS, yvalsS, 0, color=colorss[2], alpha=0.25)

    ## Microscope
    Mmicr_lim, Mmicr_lim2 = Microscope_limits()

    microscope_range = [1e-25, 1e-17]
    plt.axhline(10**Mmicr_lim, color='gray')
    plt.axhline(10**Mmicr_lim2, color='gray')

    plt.fill_between(microscope_range, np.array(microscope_range)*0 + 10**Mmicr_lim, color='gray', alpha=0.25)
    plt.fill_between(microscope_range, np.array(microscope_range)*0 + 10**Mmicr_lim2, color='gray', alpha=0.25)

    plt.xlabel(r"$\log_{10} m ~/~ \mathrm{eV}$", fontsize=12)
    plt.ylabel(r"$\log_{10} M_\mathrm{eff} ~/~  M_\mathrm{Pl}$", fontsize=12)

    ## Tick labels
    ax1 = plt.gca()

    xticks = range(-25, -16)
    xticklabels = ["%i" % x for x in xticks]
    xticks = [10**x for x in xticks]
    ax1.set_xticks( xticks )
    ax1.set_xticklabels(xticklabels)

    yticks = range(3, 11)
    yticklabels = ["%i" % y for y in yticks]
    yticks = [10**y for y in yticks]
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticklabels)

    '''
    ## Top ticks (Disabled for now)
    ax2 = ax1.twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(ax1.get_xlim())
    ax2.tick_params(direction = "in", which='both')
    ax2.minorticks_off()

    # Year, month, day, hour
    year = 3.15e7
    month = 2.6e6
    day = 86400
    hour = 3600
    minute = 60

    tick_times = [10*year, year, month, day, hour]
    tick_times = [constants.hbar / t for t in tick_times]

    time_labels = ['decade', 'year', 'month', 'day', 'hour']
    time_labels = [r"$\mathrm{%s}^{-1}$" % t for t in time_labels]

    ax2.set_xticks(tick_times)
    ax2.set_xticklabels(time_labels, fontsize=6, ha='center')
    '''
    ax2 = ax1.twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(ax1.get_xlim())
    ax2.tick_params(direction = "in")
    ax2.axvline(-18)

    # Define powers of 10 for distance ticks (in eV)
    yrplot=86400*365.
    monthplot=86400*30.
    dayplot=86400.
    hourplot=3600.
    minplot=70.
    ts_top_ticks = np.array([yrplot,monthplot,dayplot,hourplot,minplot])
    ts_names=[r'${\rm yr}^{-1}$',r'${\rm mth}^{-1}$',r'${\rm day}^{-1}$',r'${\rm hr}^{-1}$',r'${\rm min}^{-1}$']

    # Calculate corresponding Theta values for those distances
    m_top_ticks = constants.hbar/ts_top_ticks

    # Set ticks and labels for the top axis (only powers of 10)
    ax2.set_xticks(m_top_ticks)
    #ax2.set_xticklabels([f'$10^{{{int(np.log10(dist))}}}$' for dist in mgr_top_ticks])
    ax2.set_xticklabels(ts_names[inam] for inam in range(len(ts_names)))

    # Label for the secondary axis
    #ax2.set_xlabel(r'Times', fontsize=12)

    # Text labels for Microscope
    plt.text(2e-24, 9e4, r'Microscope ($M_\mathrm{e} \to \infty$)', fontsize=8)
    plt.text(2e-24, 1.4e3, r'Microscope ($M \to \infty$)', fontsize=8)


    plt.text(1e-22, 4.0e7, "CaF/Sr clocks", fontsize=10, rotation=-35)
    plt.text(1e-22, 2.5e6, "Cs/Sr clocks", fontsize=10, rotation=-35)

    # A vertical label for the cutoff at 10 mins
    plt.text(1.1e-18, 5e3, r"$m = (10~\mathrm{min})^{-1}$", fontsize = 6, rotation='vertical')

    plt.axvline(1e-24, color='brown', linestyle='dashed')
    plt.fill_betweenx([0, 1e11], 1e-24, color='brown',alpha=0.25)

    plt.text(6e-25, 5e5, r'CMB & LSS', rotation='vertical', fontsize=8)

    
    #plt.xlim(xlims)
    plt.ylim(1e3, 1e10)
    plt.tick_params(direction='in', which='both')
    #xp = np.power(10, np.mean(np.log10(xlims)))
    ax1.set_xlim(xmin=1e-25, xmax=1e-17)
    ax2.set_xlim(xmin=1e-25, xmax=1e-17)
    #plt.text(xp, 4e9, "Dark matter", ha='center', fontsize=12)

    #plt.axvline(x=2*np.pi*3.15e7**(-1)*hbar/3,c='grey',label='m=1/(3yr)',linestyle='dashed')
    #plt.title('DM constraints for the model '+r"$d=\frac{\sqrt{2\rho}}{m*M}\cos(m*t+\phi)$")
    #plt.legend(loc='lower right')
    #plt.grid()

    if savefig:
        plt.savefig("DM-M-vs-m.png", dpi = 300)
    if showfig:
        plt.show()
    #plt.clf()

if __name__ == "__main__":
    DM_plotter(savefig=True, showfig=True)
