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

couples=np.genfromtxt("1par_file.csv", delimiter=",", dtype=float)[:,1]
couples_list=['Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs']

Msun=1.12e66
AUev=1.21e17
rhoDE=(2.4e-3)**4
hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7

def Clowcurve(M,Me,i,w):
    return np.sqrt(((1.+w)/(1-w))*rhoDE) * (1/M - 1/Me) /(hbar*np.pi*2) +couples[i]

def Cupcurve(M,Me,i,w):
    return np.sqrt(((1.+w)/(1-w))*rhoDE) * (1/M - 1/Me) /(hbar*np.pi*2) -couples[i]

Merange=np.logspace(2,10)*Mpl
wss=[-0.95,-0.99,-0.999]
colorss=['blue','orange','red']

idx_pair0=0 # Change it to change the clocks pair
idx_pair1=2


''' MICROSCOPE '''

Rmic=AUev/(1.496e11/7e6)
Mearth=Msun*(5.972/1988000)
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

Mrange=np.logspace(1,10,100)*Mpl
MeMrange=[]
MMrange=[]

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


''' PLOTS '''

plt.figure(figsize=(10, 6))

plt.plot(line1x,line1y, label='Microscope',c='gray',alpha=0.5)
plt.plot(line2x,line2y, c='gray',alpha=0.5)
plt.fill_between(line1x,line1y,color='gray',alpha=0.5)
plt.fill_betweenx(line2y,line2x,color='gray',alpha=0.5)

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
    
    plt.plot(line1x,line1y,label=couples_list[idx_pair1]+' w='+str(wss[i]),color=colorss[i])
    plt.plot(line2x,line2y,color=colorss[i])

plt.title('DE constraints for the model '+r"$d=\left(\frac{1}{M}-\frac{1}{M_e}\right)\sqrt{X}\,t$")
plt.xlabel(r'$\log_{10}(M/M_{\rm pl})$')
plt.ylabel(r'$\log_{10}(Me/M_{\rm pl})$')
plt.xlim([3,10])
plt.ylim([3,10])
plt.legend(loc='upper right')
plt.tick_params(direction = "in")
plt.grid()
plt.show()

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

plt.plot(w_range,np.log10(Meff_of_w(w_range,0)),label=couples_list[0])
plt.plot(w_range,np.log10(Meff_of_w(w_range,2)),label=couples_list[2])

plt.plot(w_range,w_range*0+Mmicr_lim,label='Microscope M',color='gray',alpha=0.5)
plt.plot(w_range,w_range*0+Mmicr_lim2,label='Microscope Me',color='gray',alpha=0.9)
plt.fill_between(w_range,w_range*0+Mmicr_lim,color='gray',alpha=0.5)
plt.axvline(x=-0.95,color='brown',label='Planck')
plt.fill_betweenx(np.linspace(0,6),-0.95+0*np.linspace(0,6),color='brown',alpha=0.5)



plt.title('DE constraints for the model '+r"$d=\frac{\sqrt{X}}{M_{eff}}t$")
plt.xlabel(r"w")
plt.ylabel(r"$M_{\rm eff}$ [Mpl]")
plt.xlim([-1,-0.94])
plt.ylim([3,6])
plt.legend(loc='lower right')
plt.tick_params(direction = "in")
plt.grid()
plt.show()

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
AUev=1.21e17
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
plt.xlabel(r'$\log_{10}(\Lambda/{\rm eV})$')
plt.ylabel(r'$\log_{10}(M/M_{\rm pl})$')
plt.legend(loc='lower right')
plt.tick_params(direction = "in")
plt.grid()
plt.show()

bet_range=np.linspace(-0.2,2.5)

def alp_of_bet(bet,i,M=Mpl,Lam=2e-9):
    return np.log(M*AUev**(-1+bet)*Lam**(-2+bet)*couples[i]/eps)/np.log(Msun/(8*np.pi*M))

def fpmic(M=Mpl):
    return eta_microscope*M*Mearth/((epsilon_Ti-epsilon_Pt)*8*np.pi*Mpl*Mpl*Rmic*Rmic)

def alp_of_betMic(bet,M=Mpl,Lam=2e-9):
    return (np.log(fpmic(M=M)/(Lam**2))+bet*np.log(Lam*Rmic))/np.log(Mearth/(8*np.pi*M))


plt.plot(bet_range,alp_of_bet(bet_range,0),label=couples_list[0])
plt.plot(bet_range,alp_of_bet(bet_range,2),label=couples_list[2])
plt.plot(bet_range,alp_of_betMic(bet_range),label='Microscope',color='grey')
plt.plot(0.5, 0.5, 'o', label='Cubic Galileon')
plt.plot(0., 0.4, 'o', label='Quartic Galileon')
plt.plot(0., 0.0, 'o', label='DBIon')
plt.plot(2., 1., 'o', label='Free scalar')
plt.fill_between(bet_range,alp_of_bet(bet_range,1),10,alpha=0.25,color='orange')
plt.title('Generalized interaction '+r"$d=\frac{\Lambda^{2-\beta}}{M}\left(\frac{M_{\rm sun}}{8\pi M}\right)^{\alpha}a^{1-\beta}\epsilon\cos(\omega t)$")
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlim([-0.2,2.5])
plt.ylim([-0.2,1.5])
plt.xlabel(r"$\beta$")
plt.ylabel(r"$\alpha$")
plt.legend(loc='lower right')
plt.tick_params(direction = "in")
plt.grid()
plt.show()

 #%%
 
rhoDE=(2.4e-3)**4
rhoDM=2.6e-6
Msun=1.12e66
AUev=1.21e17
eps=0.0167
hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7

def fp(M,sig):
    return sig*M/(eps*AUev)

''' Microscope '''

Rmic=AUev/(1.496e11/7e6)
Mearth=Msun*(5.972/1988000)

def fpmic(M):
    return eta_microscope*M*Mearth/((epsilon_Ti-epsilon_Pt)*8*np.pi*Mpl*Mpl*Rmic*Rmic)

''' Forecast '''

myc4=1e-10

def myf(Lam,M,sig,c4=myc4):
    return (fp(M,sig)/AUev) + (2/Lam**3)*(fp(M,sig)/AUev)**2 + (2*c4/Lam**6)*(fp(M,sig)/AUev)**3 - Msun/(4*np.pi*M*AUev**3)

def myfm(Lam,M,c4=myc4):
    return (fpmic(M)/Rmic) + (2/Lam**3)*(fpmic(M)/Rmic)**2 + (2*c4/Lam**6)*(fpmic(M)/Rmic)**3 - Mearth/(4*np.pi*M*Rmic**3)

Lamrange=np.logspace(-18,-2,100)
Mrange10=[]
MrangeMic10=[]

for Lam in Lamrange:
    def func10(M):
        return myf(Lam,M,couples[2],c4=1e-12)
    
    def mfunc10(M):
        return myfm(Lam,M,c4=1e-12)
    
    Mrange10.append(fsolve(func10, 1e20)[0])
    MrangeMic10.append(fsolve(mfunc10, 1e20)[0])
    

def mgraviton(Lam):
    return np.sqrt(Lam**3 / Mpl)

mgr=mgraviton(Lamrange)


fig, ax1 = plt.subplots()

ax1.loglog(Lamrange[Lamrange>1/AUev], np.array(Mrange10)[Lamrange>1/AUev]/Mpl, label=couples_list[2])
ax1.loglog(Lamrange[Lamrange>1/Rmic], np.array(MrangeMic10)[Lamrange>1/Rmic]/Mpl, label='Microscope')
ax1.fill_between(Lamrange[Lamrange>1/AUev],np.array(Mrange10)[Lamrange>1/AUev]/Mpl,0,alpha=0.25,color='blue')
ax1.fill_between(Lamrange[Lamrange>1/Rmic],np.array(MrangeMic10)[Lamrange>1/Rmic]/Mpl,0,alpha=0.25,color='orange')
ax1.axvline(x=1/AUev,color='k',label='Clocks strong coupling',linestyle='dotted')
ax1.axvline(x=1/Rmic,color='brown',label='Microscope strong coupling',linestyle='dotted')
ax1.set_xlabel(r'$\Lambda$ [eV]')
ax1.set_ylabel(r'$M$ [Mpl]')
ax1.legend()
ax1.grid(True, which="both", ls="--")
#ax1.tick_params(direction = "in")

ax2 = ax1.twiny()
ax2.set_xscale('log')
ax2.set_xlim(ax1.get_xlim())

# Define powers of 10 for distance ticks (in kpc)
mgr_top_ticks = np.array([1e-40, 1e-35, 1e-30, 1e-25, 1e-20, 1e-15])

# Calculate corresponding Theta values for those distances
Lam_top_ticks = (mgr_top_ticks**2 * Mpl)**(1/3)

# Set ticks and labels for the top axis (only powers of 10)
ax2.set_xticks(Lam_top_ticks)
ax2.set_xticklabels([f'$10^{{{int(np.log10(dist))}}}$' for dist in mgr_top_ticks])

# Label for the secondary axis
ax2.set_xlabel(r'$m_{\rm g}\;[eV]$')

plt.show()

 #%%
 
couples=np.genfromtxt("3par_file.csv", delimiter=",", dtype=float)
couples_list=['Cs/Sr','N2+/Sr','CaF/Sr','CaF/Cs','Yb+/Cs','CaF/Yb+','Cf15+/Cs']

hbar=4.135667696e-15/(2*np.pi)
Mpl=1.22e28/np.sqrt(8*np.pi)
A0f=2e-18
om0f=1.99e-7

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

Mrange=np.logspace(1,10,100)*Mpl
MeMrange=[]
MMrange=[]

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



plt.loglog(mplot_range,Mmass_of_m(mplot_range,0),label=couples_list[0])
plt.loglog(mplot_range,Mmass_of_m(mplot_range,2),label=couples_list[2])
plt.loglog(mplot_range,mplot_range*0+10**Mmicr_lim,label='Microscope M',color='gray',alpha=0.5)
plt.loglog(mplot_range,mplot_range*0+10**Mmicr_lim2,label='Microscope Me',color='gray',alpha=0.9)
plt.axvline(x=2*np.pi*3.15e7**(-1)*hbar/3,c='grey',label='m=1/(3yr)',linestyle='dashed')
plt.title('DM constraints for the model '+r"$d=\frac{\sqrt{2\rho}}{m*M}\cos(m*t+\phi)$")
plt.xlabel(r"m [eV]")
plt.ylabel(r"M [Mpl]")
plt.legend(loc='upper right')
plt.grid()
plt.show()