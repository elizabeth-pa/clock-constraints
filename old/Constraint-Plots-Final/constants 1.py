import numpy as np

Ts=86400*365
rhoDE = (2.4e-3)**4 # energy density of Dark Energy in units of eV^4
rhoDM=(2.6e-6)**4 # energy density of Dark Matter in units of eV^4
hbar = 4.135667696e-15 /(2 * np.pi)
Mpl = 1.22e28/np.sqrt(8 * np.pi) # Planck mass in eV

##### Microscope #####

Msun=1.12e66
AUev=1.21e17

Rmic=AUev/(1.496e11/7e6)
Mearth=Msun*(5.972/1988000)
eta_microscope = np.sqrt(2.3**2 + 1.5**2) * 1e-15

gToeV = 1e9 / (1.8e-24)
Mearth = 6e27 * gToeV

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


###### Clocks Noise default ######

h0a = 3e-33
hm1a = 7.2e-37
Ka=0
h0b = 8e-28
hm1b = 2.8e-32
Kb=1