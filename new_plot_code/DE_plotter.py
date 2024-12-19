import constants
import Data_to_model as dtm

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

couples = np.genfromtxt("sigmas.csv", delimiter=",", dtype=float)[:,1]
# Computed for the pairs of clocks in Data_to_model and therefore for T_obs=10yr

wss = [-0.95,-0.99,-0.999]
Mrange = np.logspace(1,10,100)*constants.Mpl
colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


def MicroscopeConstraint(M, Me):
    QE=constants.Mearth * (1 / M + constants.epsilon_Fe / Me) / (1 + constants.epsilon_Fe)
    dPhiMassless=QE/ (4 * np.pi * constants.Rmic**2)
    return (8*np.pi * constants.Mpl**2 * constants.Rmic**2 / constants.Mearth
            )  * (1/M - 1/Me) * (constants.epsilon_Pt - constants.epsilon_Ti)*dPhiMassless

def lowcurveMic(M,Me):
    return MicroscopeConstraint(M, Me)-constants.eta_microscope
def upcurveMic(M,Me):
    return MicroscopeConstraint(M, Me)+constants.eta_microscope

def lowcurveClocks(M,Me,w,sigmas=couples,idx_pair=2):
    return np.sqrt(((1.+w)/(1-w))*constants.rhoDE) * (1/M - 1/Me) /(constants.hbar*np.pi*2) +sigmas[idx_pair]

def upcurveClocks(M,Me,w,sigmas=couples,idx_pair=2):
    return np.sqrt(((1.+w)/(1-w))*constants.rhoDE) * (1/M - 1/Me) /(constants.hbar*np.pi*2) -sigmas[idx_pair]

def Meff_of_w(w,i,sigmas=couples):
    return (np.sqrt(constants.rhoDE*(1+w)/(1-w))/(sigmas[i]*constants.hbar*np.pi*2))/constants.Mpl


def Me_M_curve_Mic(lowcurv,upcurv,range_M=Mrange):

    MeMrange = []
    MMrange = []

    for Mx in range_M:
        def func0(Me):
            return lowcurv(Mx,Me)
        def func1(M):
            return upcurv(M,Mx)
        
        MeMrange.append(fsolve(func0, 1e30)[0])
        MMrange.append(fsolve(func1, 1e30)[0])
        
    line1x = np.log10(range_M / constants.Mpl)
    line1y = np.log10(np.array(MeMrange) / constants.Mpl)
    line2x = np.log10(np.array(MMrange) / constants.Mpl)
    line2y = np.log10(range_M / constants.Mpl)
    return line1x, line1y, line2x, line2y

def Me_M_curve_clocks(w,sigmas=couples,range_M=Mrange,idx_pair=2):
    MeMrange = []
    MMrange = []
    for Mx in range_M:
        def func0(Me):
            return lowcurveClocks(Mx,Me,w,sigmas=sigmas,idx_pair=idx_pair)
        def func1(M):
            return upcurveClocks(M,Mx,w,sigmas=sigmas,idx_pair=idx_pair)
        MeMrange.append(fsolve(func0, 1e30)[0])
        MMrange.append(fsolve(func1, 1e30)[0])
    line1x = np.log10(range_M / constants.Mpl)
    line1y = np.log10(np.array(MeMrange) / constants.Mpl)
    line2x = np.log10(np.array(MMrange) / constants.Mpl)
    line2y = np.log10(range_M / constants.Mpl)
    return line1x, line1y, line2x, line2y


''' Plotters '''

def Microscope_plotter(range_M=Mrange):
    line1xMic, line1yMic, line2xMic, line2yMic = Me_M_curve_Mic(lowcurveMic,upcurveMic,range_M=range_M)
    plt.plot(line1xMic,line1yMic, label='Microscope',c='gray',alpha=0.5)
    plt.plot(line2xMic,line2yMic, c='gray',alpha=0.5)
    plt.fill_between(line1xMic,line1yMic,color='gray',alpha=0.5)
    plt.fill_betweenx(line2yMic,line2xMic,color='gray',alpha=0.5)
    return

def Clocks_plotter(idx_pair=2,sigmas=couples,w=wss,range_M=Mrange,
                   linestyles = ["solid", "dashed", "dotted"], colorss = [colorcycle[i] for i in range(3)]):
    for i in range(len(w)):
        line1x, line1y, line2x, line2y =Me_M_curve_clocks(w[i],sigmas=sigmas,range_M=range_M,idx_pair=idx_pair)
        plt.plot(line1x,line1y,label='w = ' + str(w[i]), color=colorss[i], linestyle = linestyles[i])
        plt.plot(line2x,line2y,color=colorss[i], linestyle = linestyles[i])
    return


def DE1_plotter(idx_pair=2,sigmas=couples,w=wss,range_M=Mrange,
                   linestyles = ["solid", "dashed", "dotted"], colorss = [colorcycle[i] for i in range(3)],
                   savefig=False):

    plt.figure(figsize=(5, 5))
    Microscope_plotter(range_M=range_M)
    Clocks_plotter(idx_pair=idx_pair,sigmas=sigmas,w=w,range_M=range_M,linestyles=linestyles,colorss=colorss)
    x_pos = 7.1

    plt.text(x_pos, 5.775, "$w = -0.95$", color = colorss[0], ha = "left")
    plt.text(x_pos, 5.425, "$w = -0.99$", color = colorss[1], ha = "left")
    plt.text(x_pos, 4.925, "$w = -0.999$", color = colorss[2], ha = "left")


    props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
    plt.text(6, 7.7, "Dark energy", fontsize=12, ha='center', bbox=props)
    plt.xlabel(r'$\log_{10} M ~/~ M_{\rm Pl}$', fontsize=12)
    plt.ylabel(r'$\log_{10} M_\mathrm{e} ~/~ M_{\rm Pl} $', fontsize=12)
    plt.xlim([4,8])
    plt.ylim([3,8])
    plt.tick_params(direction = "in")

    plt.xticks(range(4,9))
    if savefig==True:
        plt.savefig("DE-Me-vs-M.png", dpi = 300)
    plt.show()
    #plt.clf()
    return


def DE2_plotter(idx_pair=2,idx_pair_bad=0,sigmas=couples,ws=np.linspace(-0.999999,-0.94,1000),
                range_M=Mrange, colorss = [colorcycle[i] for i in range(3)], savefig=False):
    line1xMic, line1yMic, line2xMic, line2yMic =Me_M_curve_Mic(lowcurveMic,upcurveMic,range_M=range_M)
    Mmicr_lim = line2xMic[-1]
    Mmicr_lim2 = line1yMic[-1]

    # Clocks bounds
    col = colorss[0]
    xvals = ws
    yvals = np.log10(Meff_of_w(ws,idx_pair,sigmas=sigmas))
    plt.plot(xvals, yvals, label=dtm.names_list[idx_pair], color=col)
    plt.fill_between(xvals, yvals, color=col, alpha = 0.25)
    plt.text(-0.98, 5.375, "CaF/Sr clocks", color='black', rotation=7, fontsize=12)

    col = colorss[1]
    yvals = np.log10(Meff_of_w(ws,idx_pair_bad,sigmas=sigmas))
    plt.plot(xvals, yvals, label=dtm.names_list[idx_pair_bad], color=col)
    plt.fill_between(xvals, yvals, color=col, alpha = 0.25)
    plt.text(-0.98, 4.1, "Cs/Sr clocks", color='black', rotation=7, fontsize=12)

    # Microscope bounds
    w_range_full = np.linspace(-1.1, -0.94, 100)
    plt.plot(w_range_full, w_range_full * 0 + Mmicr_lim, label='Microscope M', color='gray', alpha=0.5)
    plt.plot(w_range_full, w_range_full * 0 + Mmicr_lim2, label='Microscope Me', color='gray', alpha=0.9)

    plt.text(-0.97, 5.05, r"Microscope ($M_\mathrm{e} \to \infty$)", ha = "left", fontsize=8)
    plt.text(-0.97, 3.225, r"Microscope ($M \to \infty$)", ha = "left", fontsize=8)

    # Planck
    plt.fill_between(w_range_full, w_range_full * 0 + Mmicr_lim, color='gray', alpha=0.5)
    plt.axvline(x=-0.95, color='brown', label='Planck')
    plt.fill_betweenx(np.linspace(0,7), -0.95 + 0 * np.linspace(0,7), color='brown',alpha=0.5)
    plt.text(-0.9495, 3.75, "Planck", rotation = "vertical", fontsize=8)

    # Phantom DE indication
    plt.axvline(x = -1, color = 'black', linestyle = 'dashed')
    plt.text(-1.002, 6, "Phantom DE", rotation = 'vertical', fontsize=8)

    xlims = [-1.005, -0.94]

    xp = np.mean(xlims)
    plt.text(xp, 6.8, "Dark energy", fontsize=12, ha='center')

    plt.xlabel(r"DE equation of state $w$", fontsize=12)
    plt.ylabel(r"$\log_{10} M_{\rm eff} ~/~ M_\mathrm{Pl}$", fontsize=12)
    plt.xlim(xlims)
    plt.ylim([3,7])
    #plt.legend(loc='lower right')
    plt.tick_params(direction = "in")
    #plt.grid()
    if savefig==True:
        plt.savefig("DE-Meff-vs-w.png", dpi = 300)
    plt.show()
    #plt.clf()
    return

DE1_plotter(idx_pair = 2)
DE2_plotter(idx_pair = 2,idx_pair_bad = 0)