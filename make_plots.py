"""
To do:
    - Split this large script into several smaller ones, one file per plot.  This will avoid namespace collisions,
      and make it easier to re-run for just a single plot at a time.
    - Follow common style practices, particularly add a single space around = signs and most operators
        - A good baseline is:
        - https://peps.python.org/pep-0008/
        - Particularly important is the usage of spaces around = signs etc.
    - Only define constants in one place (Mpl, Msun, etc). Normally this is done up at the top, or
      in a separate file ("constants.py") which is then imported, such as
        import constants as c
      or
        from constants import *
      Same for model parameters.
    - Don't rely on implicit conversions, e.g. the definition of AUev already had been
      converted to eV.  Instead do the conversions in the code.  For example, rather than
        AUev = 7.59e17
      it's better to do
        AUev = 1.5e13 * c / hbar
      where the first number is the value in cm.  This way it's very easy to double-check
      the numbers by googling "AU in cm"
    - Suggest factoring such that each plot is made in its own function, and probably its  to avoid namespace
      collisions and unintended side effects
    - Avoid duplicating calculations, as changes to the code then have to be repeated.  E.g. when calculating
      Lambda values for plotting, they were repeated on multiple lines.  Since these are reused, it's
      better to store the values in a variable, so that if we change which values are used at some point
      it only needs to be changed in a single place.  Multiple definitions lead to bugs very easily.
    - Use multi-step calculations in lines sparingly -- optimize for readability, not efficiency/terseness

Other comments:
    - I believe AUev was off by a factor of 2 \pi.  It was set as 1.21e17, but I'm getting 7.59e17
    - Galileon plot: is c4 = 1e-10 or 1e-12?  I see both values in the code.
    - I changed the DM plot y axis label from M to M_eff
    - For the DM plot, I get (3 yr)^-1 = 7e-24, not 4e-23.  These values differ by ~ 2\pi
"""

# Plot output resolution
dpi_setting = 300

####### PLOTTER DE1 ###########
'''

The model is d = [sqrt(X)/(M_eff)] * t = sqrt((1+w)*0.5*rhoDE) / Meff * t

where we used:
w=1-2beta
beta=phiprime^2/(2V)=X/V -> X=beta*V=(1-w)*V/2 -> sqrt(X)=sqrt((1+w)*0.5*rhoDE)

parametrized d = a * t
'''

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

couples = np.genfromtxt("1par_file.csv", delimiter=",", dtype=float)[:,1]
couples_list = ['Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs']

# Constants, in natural units (eV)
# Assumed c = hbar = 1
Msun = 1.12e66
#AUev  = 1.21e17
AUev = 7.59e17
rhoDE = (2.4e-3)**4
hbar = 4.135667696e-15/(2*np.pi)
Mpl = 1.22e28/np.sqrt(8*np.pi)
A0f = 2e-18
om0f = 1.99e-7

def Clowcurve(M,Me,i,w):
    return np.sqrt(((1.+w)/(1-w))*rhoDE) * (1/M - 1/Me) /(hbar*np.pi*2) +couples[i]

def Cupcurve(M,Me,i,w):
    return np.sqrt(((1.+w)/(1-w))*rhoDE) * (1/M - 1/Me) /(hbar*np.pi*2) -couples[i]

Merange = np.logspace(2,10)*Mpl
wss = [-0.95,-0.99,-0.999]
#colorss = ['blue','orange','red']

# Use the default cycle, skipping one because red looks nicer
colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
colorss = [colorcycle[i] for i in [0, 1, 3]]

idx_pair0 = 0 # Change it to change the clocks pair
idx_pair1 = 2


''' MICROSCOPE '''

Rmic = AUev / (1.496e11 / 7e6)
Mearth = Msun * (5.972 / 1988000)
eta_microscope = np.sqrt(2.3**2 + 1.5**2) * 1e-15

#print("Rmicroscope is %e" % Rmic)

gToeV = 1e9 / 1.8e-24
Mearth = 6e27 * gToeV

# Atomic Properties

me = 0.5e6
mp = 1e9
ZH, AH = 1, 1
mH = AH * mp

ZTi, ATi = 22, 47.9
mTi = ATi * mp
epsilon_Ti = ZTi / ATi * me / mp

ZPt, APt = 78, 195
mPt = APt * mp
epsilon_Pt = ZPt / APt * me / mp

ZFe, AFe = 26, 55.8
mFe = 55.845 * mp
epsilon_Fe = ZFe / AFe * me / mp

def QE(M,Me):
    return Mearth * (1 / M + epsilon_Fe / Me) / (1 + epsilon_Fe)

def dPhiMassless(M,Me):
    return QE(M,Me) / (4 * np.pi * Rmic**2)

def MicroscopeConstraint(M, Me):
    return (8*np.pi * Mpl**2 * Rmic**2 / Mearth)  * (1/M - 1/Me) * (epsilon_Pt - epsilon_Ti)*dPhiMassless(M,Me)

def lowcurve(M,Me):
    return MicroscopeConstraint(M, Me)-eta_microscope
def upcurve(M,Me):
    return MicroscopeConstraint(M, Me)+eta_microscope

Mrange = np.logspace(1,10,100)*Mpl
MeMrange = []
MMrange = []

for Mx in Mrange:
    def func0(Me):
        return lowcurve(Mx,Me)
    def func1(M):
        return upcurve(M,Mx)
    
    MeMrange.append(fsolve(func0, 1e30)[0])
    MMrange.append(fsolve(func1, 1e30)[0])
    
line1x = np.log10(Mrange / Mpl)
line1y = np.log10(np.array(MeMrange) / Mpl)
line2x = np.log10(np.array(MMrange) / Mpl)
line2y = np.log10(Mrange / Mpl)

Mmicr_lim = line2x[-1]
Mmicr_lim2 = line1y[-1]


''' PLOTS '''

plt.figure(figsize=(5, 5))

plt.plot(line1x,line1y, label='Microscope',c='gray',alpha=0.5)
plt.plot(line2x,line2y, c='gray',alpha=0.5)
plt.fill_between(line1x,line1y,color='gray',alpha=0.5)
plt.fill_betweenx(line2y,line2x,color='gray',alpha=0.5)

#linestyles = ["-", "--", "."]
linestyles = ["solid", "dashed", "dotted"]

for i in range(len(wss)):
    MeMrange=[]
    MMrange=[]

    for Mx in Mrange:
        def func0(Me):
            return Clowcurve(Mx,Me,idx_pair0,wss[i])
        def func1(M):
            return Cupcurve(M,Mx,idx_pair0,wss[i])

        MeMrange.append(fsolve(func0, 1e30)[0])
        MMrange.append(fsolve(func1, 1e30)[0])
    
    line1x=np.log10(Mrange/Mpl)
    line1y=np.log10(np.array(MeMrange)/Mpl)
    line2x=np.log10(np.array(MMrange)/Mpl)
    line2y=np.log10(Mrange/Mpl)
    
    #plt.plot(line1x,line1y,label=couples_list[idx_pair0]+' w='+str(wss[i]),color=colorss[i],linestyle='dashed')
    #plt.plot(line2x,line2y,color=colorss[i],linestyle='dashed')
    
    MeMrange=[]
    MMrange=[]

    for Mx in Mrange:
        def func0(Me):
            return Clowcurve(Mx,Me,idx_pair1,wss[i])
        def func1(M):
            return Cupcurve(M,Mx,idx_pair1,wss[i])

        MeMrange.append(fsolve(func0, 1e30)[0])
        MMrange.append(fsolve(func1, 1e30)[0])
    
    line1x=np.log10(Mrange/Mpl)
    line1y=np.log10(np.array(MeMrange)/Mpl)
    line2x=np.log10(np.array(MMrange)/Mpl)
    line2y=np.log10(Mrange/Mpl)
    
    #plt.plot(line1x,line1y,label=couples_list[idx_pair1]+' w='+str(wss[i]),color=colorss[i])
    plt.plot(line1x,line1y,label='w = ' + str(wss[i]), color=colorss[i], linestyle = linestyles[i])
    plt.plot(line2x,line2y,color=colorss[i], linestyle = linestyles[i])

x_pos = 7.1
plt.text(x_pos, 3.225, "Microscope", ha = "left")
plt.text(x_pos, 5.775, "$w = -0.95$", color = colorss[0], ha = "left")
plt.text(x_pos, 5.425, "$w = -0.99$", color = colorss[1], ha = "left")
plt.text(x_pos, 4.925, "$w = -0.999$", color = colorss[2], ha = "left")


props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
plt.text(6, 7.7, "Dark energy", fontsize=12, ha='center', bbox=props)

#plt.title('DE constraints for the model '+r"$d=\left(\frac{1}{M}-\frac{1}{M_e}\right)\sqrt{X}\,t$")
#plt.title("Dark energy coupling strengths vs. equation of state for a CaF/Sr clock pair")
#plt.title("Dark energy coupling strengths vs. equation of state for a CaF/Sr clock pair")
plt.xlabel(r'$\log_{10} M ~/~ M_{\rm Pl}$', fontsize=12)
plt.ylabel(r'$\log_{10} M_\mathrm{e} ~/~ M_{\rm Pl} $', fontsize=12)
plt.xlim([4,8])
plt.ylim([3,8])
#plt.legend(loc='upper right')
plt.tick_params(direction = "in")

plt.xticks( range(4, 9) )
#plt.grid()
plt.savefig("plots/DE-Me-vs-M.png", dpi = dpi_setting)
#plt.show()
plt.clf()

 #%%
 
####### PLOTTER DE2 ###########

'''

The model is d = [sqrt(X)/(M_eff)] * t = sqrt((1+w)*0.5*rhoDE) / Meff * t

where we used:
w=1-2beta
beta=phiprime^2/(2V)=X/V -> X=beta*V=(1-w)*V/2 -> sqrt(X)=sqrt((1+w)*0.5*rhoDE)

parametrized d = a * t


'''

couples=np.genfromtxt("1par_file.csv", delimiter=",", dtype=float)[:,1]
couples_list=['Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs']

rhoDE=(2.4e-3)**4
hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7

def Meff_of_w(w,i):
    return (np.sqrt(rhoDE*(1+w)/(1-w))/(couples[i]*hbar*np.pi*2))/Mpl

w_range=np.linspace(-0.999999,-0.94,1000)

# Clocks bounds
col = colorss[0]
xvals = w_range
yvals = np.log10(Meff_of_w(w_range,2))
plt.plot(xvals, yvals, label=couples_list[2], color=col)
plt.fill_between(xvals, yvals, color=col, alpha = 0.25)
plt.text(-0.98, 5.375, "CaF/Sr clocks", color='black', rotation=7, fontsize=12)

col = colorss[1]
yvals = np.log10(Meff_of_w(w_range,0))
plt.plot(xvals, yvals, label=couples_list[0], color=col)
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

#plt.title('DE constraints for the model '+r"$d=\frac{\sqrt{X}}{M_{eff}}t$")
#plt.title("Dark energy constraints for multiple clock combinations")
plt.xlabel(r"DE equation of state $w$", fontsize=12)
plt.ylabel(r"$\log_{10} M_{\rm eff} ~/~ M_\mathrm{Pl}$", fontsize=12)
plt.xlim(xlims)
plt.ylim([3,7])
#plt.legend(loc='lower right')
plt.tick_params(direction = "in")
#plt.grid()
plt.savefig("plots/DE-Meff-vs-w.png", dpi = dpi_setting)
#plt.show()
plt.clf()

 #%%
 
####### PLOTTER general interaction ###########

'''
The model is d = Lambda^(2+b) * (r^(b+1))/(b+1) * (Msun/4pi)^a * M^(-a-1)
'''

import numpy as np
from matplotlib import pyplot as plt

couples=np.genfromtxt("1par_file.csv", delimiter=",", dtype=float)[:,0]
couples_list=['Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs']

rhoDE=(2.4e-3)**4
rhoDM=2.6e-6
Msun=1.12e66
#AUev=1.21e17
eps=0.0167
hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7

def Lambda_of_M(M,i,alp=0.5,bet=0.5):
    return (M*(Msun/(8*np.pi*M))**(-alp) * AUev**(-1+bet) * couples[i]/eps)**(1/(2-bet))

def phipInt(Lam,M,alp,bet):
    return Lam**2 * (Mearth/(8*np.pi*M))**alp * (Lam*Rmic)**(-bet)

def resMic(Lam, M):
    return (8*np.pi * Mpl**2 * Rmic**2 / Mearth)  * (1/M) * (epsilon_Pt - epsilon_Ti)*phipInt(Lam,M,1/2,1/2)+eta_microscope

M_range=np.logspace(-2,4)
Lamrange=[]

for M in M_range:
    def func10(Lam):
        return resMic(Lam,M*Mpl)
    
    Lamrange.append(fsolve(func10, 1e-8)[0])


plt.plot(np.log10(Lambda_of_M(M_range*Mpl,0)),np.log10(M_range),label=couples_list[0])
plt.plot(np.log10(Lambda_of_M(M_range*Mpl,2)),np.log10(M_range),label=couples_list[2])
plt.fill_betweenx(np.log10(M_range),np.log10(Lambda_of_M(M_range*Mpl,1)),1.,alpha=0.25,color='orange')
plt.plot(np.log10(np.array(Lamrange)),np.log10(M_range),label='Microscope',color='gray')

plt.title('Generalized interaction '+r"$d=\frac{\Lambda^{2-\beta}}{M}\left(\frac{M_{\rm sun}}{8\pi M}\right)^{\alpha}a^{1-\beta}\epsilon\cos(\omega t)$")
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlim([-12.1,-4])
plt.xlabel(r'$\log_{10}(\Lambda/{\rm eV})$', fontsize=12)
plt.ylabel(r'$\log_{10}(M/M_{\rm pl})$', fontsize=12)
plt.legend(loc='lower right')
plt.tick_params(direction = "in")
plt.grid()
plt.savefig("plots/MG-gen-Meff-vs-Lambda.png", dpi = dpi_setting)
#plt.show()
plt.clf()

bet_range=np.linspace(-0.2,2.5)

def alp_of_bet(bet,i,M=Mpl,Lam=5e-10):
    return np.log(M*AUev**(-1+bet)*Lam**(-2+bet)*couples[i]/eps)/np.log(Msun/(8*np.pi*M))

def fpmic(M=Mpl):
    return eta_microscope*M*Mearth/((epsilon_Ti-epsilon_Pt)*8*np.pi*Mpl*Mpl*Rmic*Rmic)

def alp_of_betMic(bet,M=Mpl,Lam=5e-10):
    return (np.log(fpmic(M=M)/(Lam**2))+bet*np.log(Lam*Rmic))/np.log(Mearth/(8*np.pi*M))


plt.figure(figsize = (5, 5))
Mval = Mpl
Lval = 1e-10


# CaF/Sr clock
xvals = bet_range
yvals = alp_of_bet(bet_range, 2, M=Mval, Lam=Lval)
col = colorss[0]
plt.plot(xvals, yvals, label=couples_list[2], color = col)
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
plt.savefig("plots/MG-gen-alpha-vs-beta.png", dpi=dpi_setting)
#plt.show()
plt.clf()

 #%%
 
rhoDE=(2.4e-3)**4
rhoDM=2.6e-6
Msun=1.12e66
#AUev=1.21e17
eps=0.0167
hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7


def fp(M,sig):
    return sig*M/(eps*AUev)

''' Microscope '''

Rmic = AUev/(1.496e11/7e6)
Mearth = Msun*(5.972/1988000)

def fpmic(M):
    return eta_microscope*M*Mearth/((epsilon_Ti-epsilon_Pt)*8*np.pi*Mpl*Mpl*Rmic*Rmic)

''' Forecast '''

myc4=1e-10

def myf(Lam,M,sig,c4=myc4):
    return (fp(M,sig)/AUev) + (2/Lam**3)*(fp(M,sig)/AUev)**2 + (2*c4/Lam**6)*(fp(M,sig)/AUev)**3 - Msun/(4*np.pi*M*AUev**3)

def myfm(Lam,M,c4=myc4):
    return (fpmic(M)/Rmic) + (2/Lam**3)*(fpmic(M)/Rmic)**2 + (2*c4/Lam**6)*(fpmic(M)/Rmic)**3 - Mearth/(4*np.pi*M*Rmic**3)

Lamrange = np.logspace(-20, 0, 100)
Mrange10 = []
MrangeMic10 = []

for Lam in Lamrange:
    def func10(M):
        return myf(Lam,M,couples[2], c4=1e-10)
    
    def mfunc10(M):
        return myfm(Lam, M, c4=1e-10)
    
    Mrange10.append(fsolve(func10, 1e20)[0])
    MrangeMic10.append(fsolve(mfunc10, 1e20)[0])
    

def mgraviton(Lam):
    return np.sqrt(Lam**3 / Mpl)

mgr=mgraviton(Lamrange)

fig, ax1 = plt.subplots(figsize=(5,5))


# Clocks line
#print("AU is %e" % AUev**-1)
# To do: check Earth mass, Sun mass, etc
#print("%e" % Mearth)
Lvals = Lamrange[Lamrange > 1 / AUev]
Mvals = np.array(Mrange10)[Lamrange > 1 / AUev] / Mpl

# Add a point to give a vertical strong-coupling line
Lvals = np.insert(Lvals, 0, [Lvals[0]])
Mvals = np.insert(Mvals, 0, [0])

col = colorss[0]
ax1.loglog(Lvals, Mvals, label=couples_list[2], color = col)
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
Lvals = Lamrange[Lamrange > 1 / AUev]
Mvals = np.array(MrangeMic10)[Lamrange > 1 / AUev] / Mpl

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
x_placement = np.pow(10, np.mean( np.log10(ax1.get_xlim())))
#ax1.text(x_placement, 8e4, "Galileon constraints", fontsize = 12, ha = "center")
ax1.text(x_placement, 2e5, "Galileon constraints", fontsize = 12, ha = "center")

# Text box for c_4
props = dict(boxstyle='round', facecolor = "white", alpha=1.0)
#ax1.text(1e-17, 1e2, r"$c_4 = 10^{-12}$", bbox = props)
#ax1.text(x_placement, 1e3, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")
#ax1.text(4e-18, 1e3, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")
ax1.text(4e-14, 1e4, r"$c_4 = 10^{-12}$", bbox = props, ha = "center")

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

ax2.axvline(-33)

# Define powers of 10 for distance ticks (in kpc)
mgr_top_ticks = np.array([1e-35, 1e-30, 1e-25, 1e-20, 1e-15])

# Calculate corresponding Theta values for those distances
Lam_top_ticks = (mgr_top_ticks**2 * Mpl)**(1/3)

# Set ticks and labels for the top axis (only powers of 10)
ax2.set_xticks(Lam_top_ticks)
#ax2.set_xticklabels([f'$10^{{{int(np.log10(dist))}}}$' for dist in mgr_top_ticks])
ax2.set_xticklabels(str(int(mg)) for mg in np.log10(mgr_top_ticks))

# Label for the secondary axis
ax2.set_xlabel(r'$\log_{10} m_{\rm g} ~/~ \mathrm{eV}$', fontsize=12)

plt.savefig("plots/MG-gal-M-vs-Lambda.png", dpi = dpi_setting)
#plt.show()
plt.clf()

 #%%
 
couples=np.genfromtxt("3par_file.csv", delimiter=",", dtype=float)
couples_list=['Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs']

hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7

# BE: I can't make heads or tails of this.  Presumably for the model signal A \cos( \omega t ) then there
# is a different constraints for each A given a particular frequency \omega.
# What I suggest doing is convert A and \omega to eV, so that then the calculation of M for each m can
# be done in natural units.  Also please add comments & clearly named variables so that it's easier
# for someone to figure out how this works.
# The arrangement is spaghetti-ish: Mmass_of_m calls Mmass, which calls A_thr.  Such an arrangement
# might be ok if Mmass and A_thr were used elsewhere but they only appear in one place so they
# should probably not be their own functions.
def A_thr(idx,omf):
    return np.sqrt( couples[idx,0]**2 * (omf/om0f) + couples[idx,1]**2 * (A0f/omf)**2 * (omf/om0f) )

def Mmass(a,rhoDM=2.6e-6): # in units of Mpl
    return np.sqrt(2*rhoDM)/(hbar*a*Mpl)

def Mmass_of_m(m,idx,rhoDM=2.6e-6):
    return Mmass(A_thr(idx,m/hbar),rhoDM=rhoDM)

mplot_range=np.logspace(-24,-18)

''' MICROSCOPE '''

eta_microscope = np.sqrt(2.3**2 + 1.5**2) * 1e-15

gToeV = 1e9 / (1.8e-24)
Mearth = 6e27 * gToeV

# Atomic Properties

me = 0.5e6
mp = 1e9
ZH, AH = 1, 1
mH = AH * mp

ZTi, ATi = 22, 47.9
mTi = ATi * mp
epsilon_Ti = ZTi / ATi * me / mp

ZPt, APt = 78, 195
mPt = APt * mp
epsilon_Pt = ZPt / APt * me / mp

ZFe, AFe = 26, 55.8
mFe = 55.845 * mp
epsilon_Fe = ZFe / AFe * me / mp

def QE(M,Me):
    return Mearth * (1 / M + epsilon_Fe / Me) / (1 + epsilon_Fe)

def dPhiMassless(M,Me):
    return QE(M,Me) / (4 * np.pi * Rmic**2)

def MicroscopeConstraint(M, Me):
    return (8*np.pi * Mpl**2 * Rmic**2 / Mearth)  * (1/M - 1/Me) * (epsilon_Pt - epsilon_Ti)*dPhiMassless(M,Me)

def lowcurve(M,Me):
    return MicroscopeConstraint(M, Me)-eta_microscope
def upcurve(M,Me):
    return MicroscopeConstraint(M, Me)+eta_microscope

Mrange = np.logspace(1,10,100)*Mpl
MeMrange = []
MMrange = []

for Mx in Mrange:
    def func0(Me):
        return lowcurve(Mx,Me)
    def func1(M):
        return upcurve(M,Mx)
    
    MeMrange.append(fsolve(func0, 1e30)[0])
    MMrange.append(fsolve(func1, 1e30)[0])
    
line1x=np.log10(Mrange/Mpl)
line1y=np.log10(np.array(MeMrange)/Mpl)
line2x=np.log10(np.array(MMrange)/Mpl)
line2y=np.log10(Mrange/Mpl)

Mmicr_lim=line2x[-1]
Mmicr_lim2=line1y[-1]

xlims = [1e-25, 1e-17]
plt.xlim(xlims)
plt.ylim(1e3, 1e10)

plt.tick_params(direction='in', which='both')


xp = pow(10, np.mean(np.log10(xlims)))
plt.text(xp, 4e9, "Dark matter", ha='center', fontsize=12)


## CaF/Sr clock 
xvals = mplot_range
yvals = Mmass_of_m(mplot_range,2)
col = colorss[0]

# Make it downturn at each end
xvals = np.concatenate(([xvals[0]], xvals, [xvals[-1]]))
yvals = np.concatenate(([0], yvals, [0]))

plt.loglog(xvals, yvals, label=couples_list[2], color=col)
plt.fill_between(xvals, yvals, 0, color=col, alpha=0.25)

## Cs/Sr clock
xvals = mplot_range
yvals = Mmass_of_m(mplot_range,0)
col = colorss[1]

# Make it downturn at each end
xvals = np.concatenate(([xvals[0]], xvals, [xvals[-1]]))
yvals = np.concatenate(([0], yvals, [0]))

plt.loglog(xvals, yvals, label=couples_list[0], color=col)
plt.fill_between(xvals, yvals, 0, color=col, alpha=0.25)


# Text labels for clocks
plt.text(1e-21, 2.5e7, "CaF/Sr clocks", fontsize=10, rotation=-30)
plt.text(1e-21, 1.6e6, "Cs/Sr clocks", fontsize=10, rotation=-30)

#plt.loglog(mplot_range,Mmass_of_m(mplot_range,0),label=couples_list[0])



## Microscope
microscope_range = [1e-26, 1e-17]
plt.axhline(10**Mmicr_lim, color='gray')
plt.axhline(10**Mmicr_lim2, color='gray')

plt.fill_between(microscope_range, np.array(microscope_range)*0 + 10**Mmicr_lim, color='gray', alpha=0.25)
plt.fill_between(microscope_range, np.array(microscope_range)*0 + 10**Mmicr_lim2, color='gray', alpha=0.25)

plt.text(2e-24, 9e4, r'Microscope ($M_\mathrm{e} \to \infty$)', fontsize=8)
plt.text(2e-24, 1.4e3, r'Microscope ($M \to \infty$)', fontsize=8)

#plt.loglog(microscope_range, microscope_range*0+10**Mmicr_lim, label='Microscope M', color='gray',alpha=0.5)
#plt.loglog(microscope_range, microscope_range*0+10**Mmicr_lim2, label='Microscope Me', color='gray',alpha=0.9)

## Clocks 2302.04565
# The stated bound is sigma_mu < 1.6e-13 / sqrt(tau / s)
# for 600 s < tau < 80000
tmin = 600
tmax = 80000

# Convert these times to eV
m_min = hbar / tmax
m_max = hbar / tmin

# Convert their bound to use m in eV
sigma = 1.6e-13 / np.sqrt(hbar)

m_sherrill = np.logspace(np.log10(m_min), np.log10(m_max), 10)

M_sherrill = np.sqrt(2 * rhoDM) / (sigma * np.pow(m_sherrill, 3/2.)) / Mpl

# Add a point at the front to get a boundary on the left

m_sherrill = np.insert(m_sherrill, 0, [m_sherrill[0]])
M_sherrill = np.insert(M_sherrill, 0, [0])


#col = 'olivedrab'
col = 'slateblue'
plt.loglog(m_sherrill, M_sherrill, color=col)
plt.fill_between(m_sherrill, M_sherrill, alpha=0.25, color=col)

plt.text(2.25e-20, 5e3, 'Yb/Cs\nclocks', ha='center', fontsize=8)

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

## Top ticks
ax2 = ax1.twiny()
ax2.set_xscale('log')
ax2.set_xlim(ax1.get_xlim())
ax2.tick_params(direction = "in", which='both')
ax2.minorticks_off()

# Define powers of 10 for distance ticks (in kpc)
# Year, month, day, hour
year = 3.15e7
month = 2.6e6
day = 86400
hour = 3600
minute = 60

tick_times = [10*year, year, month, day, hour]
tick_times = [hbar / t for t in tick_times]

time_labels = ['decade', 'year', 'month', 'day', 'hour']
time_labels = [r"$\mathrm{%s}^{-1}$" % t for t in time_labels]

ax2.set_xticks(tick_times)
ax2.set_xticklabels(time_labels, fontsize=6, ha='center')

# A vertical label for the cutoff at 10 mins
plt.text(1.1e-18, 5e3, r"$m = (10~\mathrm{min})^{-1}$", fontsize = 6, rotation='vertical')


col = 'brown'
plt.axvline(1e-24, color=col, linestyle='dashed')
plt.fill_betweenx([0, 1e11], 1e-24, color='brown',alpha=0.25)

plt.text(6e-25, 5e5, r'CMB & LSS', rotation='vertical', fontsize=8)

#plt.axvline(x=2*np.pi*3.15e7**(-1)*hbar/3,c='grey',label='m=1/(3yr)',linestyle='dashed')
#plt.title('DM constraints for the model '+r"$d=\frac{\sqrt{2\rho}}{m*M}\cos(m*t+\phi)$")
#plt.legend(loc='lower right')
#plt.grid()
plt.savefig("plots/DM-M-vs-m.png", dpi = dpi_setting)
#plt.show()
#plt.clf()
