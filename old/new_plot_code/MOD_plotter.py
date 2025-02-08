import constants
import Data_to_model as dtm

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

couples = np.genfromtxt("sigmas.csv", delimiter=",", dtype=float)[:,0]
# Computed for the pairs of clocks in Data_to_model and therefore for T_obs=10yr

colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

def Lambda_of_M(M,i,sigmas=couples,alp=0.5,bet=0.5):
    return (M*(constants.Msun/(8*np.pi*M))**(-alp) * constants.AUev**(-1+bet) * sigmas[i]/constants.eps)**(1/(2-bet))

def phipInt(Lam,M,alp=0.5,bet=0.5):
    return Lam**2 * (constants.Mearth/(8*np.pi*M))**alp * (Lam*constants.Rmic)**(-bet)

def resMic(Lam, M,alp=0.5,bet=0.5):
    return (8*np.pi * constants.Mpl**2 * constants.Rmic**2 / constants.Mearth)  * (1/M) * (constants.epsilon_Pt - constants.epsilon_Ti)*phipInt(Lam,M,alp=alp,bet=bet)+constants.eta_microscope

def Gen_Int1_plotter(idx_pair=2, idx_pair_bad=0, M_range=np.logspace(-2,4),
                     sigmas=couples,alp=0.5,bet=0.5, savefig=False):
    Lamrange=[]

    for M in M_range:
        def func10(Lam):
            return resMic(Lam,M*constants.Mpl,alp=alp,bet=bet)
        
        Lamrange.append(fsolve(func10, 1e-8)[0])
    plt.plot(np.log10(Lambda_of_M(M_range*constants.Mpl,idx_pair_bad,sigmas=sigmas,alp=alp,bet=bet)),np.log10(M_range),label=dtm.names_list[idx_pair_bad])
    plt.plot(np.log10(Lambda_of_M(M_range*constants.Mpl,idx_pair,sigmas=sigmas,alp=alp,bet=bet)),np.log10(M_range),label=dtm.names_list[idx_pair])
    plt.fill_betweenx(np.log10(M_range),np.log10(Lambda_of_M(M_range*constants.Mpl,idx_pair,sigmas=sigmas,alp=alp,bet=bet)),1.,alpha=0.25,color='orange')
    plt.plot(np.log10(np.array(Lamrange)),np.log10(M_range),label='Microscope',color='gray')

    plt.title('Generalized interaction '+r"$d=\frac{\Lambda^{2-\beta}}{M}\left(\frac{M_{\rm sun}}{8\pi M}\right)^{\alpha}a^{1-\beta}\epsilon\cos(\omega t)$")
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlim([-12.1,-4])
    plt.xlabel(r'$\log_{10}(\Lambda/{\rm eV})$', fontsize=12)
    plt.ylabel(r'$\log_{10}(M/M_{\rm pl})$', fontsize=12)
    plt.legend(loc='lower right')
    plt.tick_params(direction = "in")
    plt.grid()
    if savefig==True:
        plt.savefig("MG-gen-Meff-vs-Lambda.png", dpi = 300)
    plt.show()
    #plt.clf()
    return 

def Gen_Int2_plotter(idx_pair=2, bet_range=np.linspace(-0.2,2.5),Mval = constants.Mpl, Lval = 1e-10,
                     sigmas=couples,colorss = [colorcycle[i] for i in range(3)], savefig=False):

    def alp_of_bet(bet,i,M=Mval,Lam=Lval):
        return np.log(M*constants.AUev**(-1+bet)*Lam**(-2+bet)*sigmas[i]/constants.eps)/np.log(constants.Msun/(8*np.pi*M))

    def fpmic(M=Mval):
        return constants.eta_microscope*M*constants.Mearth/((constants.epsilon_Ti-constants.epsilon_Pt)*8*np.pi*constants.Mpl**2 *constants.Rmic**2)

    def alp_of_betMic(bet,M=Mval,Lam=Lval):
        return (np.log(fpmic(M=M)/(Lam**2))+bet*np.log(Lam*constants.Rmic))/np.log(constants.Mearth/(8*np.pi*M))


    plt.figure(figsize = (5, 5))
    


    # CaF/Sr clock
    xvals = bet_range
    yvals = alp_of_bet(bet_range, idx_pair, M=Mval, Lam=Lval)
    col = colorss[0]
    plt.plot(xvals, yvals, label=dtm.names_list[idx_pair], color = col)
    plt.fill_between(xvals, yvals, 10, alpha=0.25, color=col)
    plt.text(1, 0.61, "CaF/Sr clocks", rotation=23, fontsize = 12)

    # Microscope
    xvals = bet_range
    yvals = alp_of_betMic(bet_range, M=Mval, Lam=Lval)
    col = 'grey'
    plt.plot(xvals, yvals, label='Microscope', color=col, linestyle='dashed')
    plt.fill_between(xvals, yvals, 10, alpha=0.25, color=col)
    plt.text(1.7, 0.66, "Microscope", rotation=12.5, fontsize = 8)


    # Label theories in \alpha, \beta space
    col = 'black'
    fs = 8
    plt.plot(0.5, 0.5, 'o', label='Cubic Galileon', color=col)
    plt.text(0.515, 0.425, "Cubic\nGalileon", fontsize=fs, ha='left')

    plt.plot(0., 0.4, 'o', label='Quartic Galileon', color=col)
    plt.text(0.015, 0.32, "Quartic\nGalileon", fontsize=fs, ha='left')

    plt.plot(0., 0.0, 'o', label='DBIon', color=col)
    plt.text(0.02, 0.02, "DBIon", fontsize=fs, ha='left')

    plt.plot(2., 1., 'o', label='Free scalar', color=col)
    plt.text(1.625, 0.965, "Free scalar", fontsize=fs, ha='left')
    #plt.title('Generalized interaction '+r"$d=\frac{\Lambda^{2-\beta}}{M}\left(\frac{M_{\rm sun}}{8\pi M}\right)^{\alpha}a^{1-\beta}\epsilon\cos(\omega t)$")
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.title('Generalized interaction')
    props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
    plt.text(1.0, 0.975, 'Generalized interaction', fontsize=12, ha='center', bbox=props)

    # Text box with values
    props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
    s = r"$\Lambda = 10^{-10}~\mathrm{eV}$" + "\n" + "$M = M_\mathrm{Pl}$"
    plt.text(1.5, 0.025, s, bbox = props, ha = "left", fontsize = 10)

    plt.xlim(-0.1, 2.1)
    plt.ylim(-0.05, 1.05)
    plt.xlabel(r"$\beta$", fontsize=12)
    plt.ylabel(r"$\alpha$", fontsize=12)
    #plt.legend(loc='lower right')
    plt.tick_params(direction = "in")
    #plt.grid()
    if savefig==True:
        plt.savefig("MG-gen-alpha-vs-beta.png", dpi=300)
    plt.show()
    #plt.clf()
    return

def Gal_plotter(idx_pair=2,Lamrange = np.logspace(-20, 0, 100),c4=1e-10,
                sigmas=couples, colorss = [colorcycle[i] for i in range(3)], savefig=False):

    ''' Clocks '''

    def fp(M,sig):
        return sig*M/(constants.eps*constants.AUev)

    ''' Microscope '''

    def fpmic(M):
        return constants.eta_microscope*M*constants.Mearth/((constants.epsilon_Ti-constants.epsilon_Pt)*8*np.pi*constants.Mpl**2*constants.Rmic**2)

    ''' Forecast '''

    def myf(Lam,M,sig,c4=c4):
        return (fp(M,sig)/constants.AUev) + (2/Lam**3)*(fp(M,sig)/constants.AUev)**2 + (2*c4/Lam**6)*(fp(M,sig)/constants.AUev)**3 - constants.Msun/(4*np.pi*M*constants.AUev**3)

    def myfm(Lam,M,c4=c4):
        return (fpmic(M)/constants.Rmic) + (2/Lam**3)*(fpmic(M)/constants.Rmic)**2 + (2*c4/Lam**6)*(fpmic(M)/constants.Rmic)**3 - constants.Mearth/(4*np.pi*M*constants.Rmic**3)

    Mrange10 = []
    MrangeMic10 = []

    for Lam in Lamrange:
        def func10(M):
            return myf(Lam,M,sigmas[idx_pair], c4=c4)
        
        def mfunc10(M):
            return myfm(Lam, M, c4=c4)
        
        Mrange10.append(fsolve(func10, 1e20)[0])
        MrangeMic10.append(fsolve(mfunc10, 1e20)[0])
        

    fig, ax1 = plt.subplots(figsize=(5,5))


    # Clocks line
    #print("AU is %e" % AUev**-1)
    # To do: check Earth mass, Sun mass, etc
    #print("%e" % Mearth)
    Lvals = Lamrange[Lamrange > 1 / constants.AUev]
    Mvals = np.array(Mrange10)[Lamrange > 1 / constants.AUev] / constants.Mpl

    # Add a point to give a vertical strong-coupling line
    Lvals = np.insert(Lvals, 0, [Lvals[0]])
    Mvals = np.insert(Mvals, 0, [0])

    col = colorss[0]
    ax1.loglog(Lvals, Mvals, label=dtm.names_list[idx_pair], color = col)
    #ax1.axvline(x = Lvals[0], ymax = line_height, color = col)#, label='Clocks strong coupling',linestyle='dotted')
    ax1.fill_between(Lvals, Mvals, 0, alpha = 0.25, color = col)
    #ax1.text(2e-17, 1e-10, "CaF/Sr\nclocks", va = "center")
    ax1.text(2e-5, 3e3, "CaF/Sr clocks", va = "center")

    # Labels of Vainshtein regions
    ax1.text(1e-2, 8e3, r"$R_{\rm V3} < {\rm AU}$", fontsize = 6)
    ax1.text(1.85e-10, 1e1, r"$R_{\rm V4} < {\rm AU} < R_{\rm V3}$", fontsize = 6, rotation = 53)
    ax1.text(3e-15, 9e-6, r"${\rm AU} < R_{\rm V4}$", fontsize = 6, rotation = 65)

    # Repeat the above for microscope
    # All is the same except Mrange10 --> MrangeMic10 and AUev --> Rmic
    #Lvals = Lamrange[Lamrange > 1 / Rmic]
    #Mvals = np.array(MrangeMic10)[Lamrange > 1 / Rmic] / Mpl

    # Use these lines to disable the vertical line at Rmic
    Lvals = Lamrange[Lamrange > 1 / constants.AUev]
    Mvals = np.array(MrangeMic10)[Lamrange > 1 / constants.AUev] / constants.Mpl

    Lvals = np.insert(Lvals, 0, [Lvals[0]])
    Mvals = np.insert(Mvals, 0, [0])

    #col = colorss[1]
    col = 'gray'
    ax1.loglog(Lvals, Mvals,  label = "Microscope", color = col, linestyle = "dashed")
    ax1.fill_between(Lvals, Mvals, 0, alpha = 0.25, color = col)
    #ax1.text(1e-8, 1e-10, "Microscope", va = "center")
    ax1.text(5e-4, 8e4, "Microscope", va = "center", fontsize = 8)

    '''
    ax1.loglog(Lamrange[Lamrange>1/Rmic], np.array(MrangeMic10)[Lamrange>1/Rmic]/Mpl, label='Microscope')

    #ax1.fill_between(Lamrange[Lamrange>1/AUev],np.array(Mrange10)[Lamrange>1/AUev]/Mpl,0,alpha=0.25,color='blue')
    ax1.fill_between(Lamrange[Lamrange>1/Rmic],np.array(MrangeMic10)[Lamrange>1/Rmic]/Mpl,0,alpha=0.25,color='orange')

    # Plot the lines where the models become strongly coupled.
    #strong_coupling_Mmax = Mrange10[1/AUev] / Mpl
    #ax1.axvline(x=1/AUev, ymax = strong_coupling_Mmax, color='k',label='Clocks strong coupling',linestyle='dotted')
    ax1.axvline(x=1/Rmic, color='brown',label='Microscope strong coupling',linestyle='dotted')
    '''

    ax1.set_xlabel(r'$\log_{10} \Lambda ~/~ \mathrm{eV}$', fontsize=12)
    ax1.set_ylabel(r'$\log_{10} M ~/~ M_\mathrm{Pl}$', fontsize=12)

    xticks = np.logspace(-20, 0, 11)
    ax1.set_xticks(xticks)
    yticks = np.logspace(-14, 6, 11)
    ax1.set_yticks(yticks)
    #ax1.set_xticks(np.logspace(-20, 0, 10))
    ax1.tick_params(direction = "in")

    #ax1.set_xlim(1e-20, 1e0)
    #ax1.set_ylim(1e-14, 1e6)
    ax1.set_xlim(1e-16, 1e0)
    ax1.set_ylim(1e-6, 1e6)


    ax1.set_xticklabels(str(int(x)) for x in np.log10(xticks))
    ax1.set_yticklabels(str(int(y)) for y in np.log10(yticks))

    ax1.axhline(1, linestyle='dashed', color='black')
    plt.text(5e-16, 4e-1, r"$M = M_\mathrm{Pl}$", fontsize=8)

    ax1.axvline(1.3e-13, linestyle='dashed', color='black')
    plt.text(3.5e-14, 1.5e-2, r"$m_\mathrm{g} = H_0$", fontsize=8, rotation='vertical')


    # Title, inset
    # Placed in the center, requires a logarithmic mean of x axis limits
    x_placement = np.power(10, np.mean( np.log10(ax1.get_xlim())))
    #ax1.text(x_placement, 8e4, "Galileon constraints", fontsize = 12, ha = "center")
    ax1.text(x_placement, 2e5, "Galileon constraints", fontsize = 12, ha = "center")

    # Text box for c_4
    props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
    #ax1.text(1e-17, 1e2, r"$c_4 = 10^{-12}$", bbox = props)
    #ax1.text(x_placement, 1e3, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")
    #ax1.text(4e-18, 1e3, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")
    ax1.text(4e-14, 1e4, r"$c_4 = $"+str(c4), bbox = props, ha = "center")

    # Labels for the strong-coupling scales
    #ax1.text(2.25e-18, 2e-14, r"$\Lambda = \mathrm{AU}^{-1}$",
            #fontsize = 8, rotation = "vertical", ha = "left") 
    ax1.text(2.6e-18, 2e-14, r"$\Lambda = \mathrm{AU}^{-1}$",
            fontsize = 8, rotation = "vertical", ha = "left") 
    ax1.text(4.6e-14, 3.1e-14, r"$\Lambda = R_{\mu\mathrm{scope}}^{-1}$",
            fontsize = 8, rotation = "vertical", ha = "left") 

    #ax1.legend()
    #ax1.grid(True, which="both", ls="--")

    ax2 = ax1.twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(ax1.get_xlim())
    ax2.tick_params(direction = "in")

    ax2.axvline(-33) # Hubble rate

    # Define powers of 10 for distance ticks (in eV)
    mgr_top_ticks = np.array([1e-35, 1e-30, 1e-25, 1e-20, 1e-15])

    # Calculate corresponding Theta values for those distances
    Lam_top_ticks = (mgr_top_ticks**2 * constants.Mpl)**(1/3)

    # Set ticks and labels for the top axis (only powers of 10)
    ax2.set_xticks(Lam_top_ticks)
    #ax2.set_xticklabels([f'$10^{{{int(np.log10(dist))}}}$' for dist in mgr_top_ticks])
    ax2.set_xticklabels(str(int(mg)) for mg in np.log10(mgr_top_ticks))

    # Label for the secondary axis
    ax2.set_xlabel(r'$\log_{10} m_{\rm g} ~/~ \mathrm{eV}$', fontsize=12)
    if savefig==True:
        plt.savefig("MG-gal-M-vs-Lambda.png", dpi = 300)
    plt.show()
    #plt.clf()



Gen_Int1_plotter()
Gen_Int2_plotter()
Gal_plotter()