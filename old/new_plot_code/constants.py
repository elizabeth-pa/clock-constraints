import numpy as np

hbar = 4.135667696e-15 /(2 * np.pi) # eV s
c = 3e10    # cm / s
Mpl = 1.22e28 / np.sqrt(8 * np.pi) # Planck mass in eV
gToeV = 1e9 / (1.8e-24)

rhoDE = (2.4e-3)**4 # energy density of Dark Energy in units of eV^4
rhoDM = 2.6e-6 # energy density of Dark Matter in units of eV^4

# BE: I recommend converting Ts to natural units
# here, so as to minimize unit conversions in other parts of the code
Ts = 3600*24*365
Msun = 2e33 * gToeV
AUev = 7.59e17
Mearth = Msun * (5.972/1988000)
eps = 0.0167    # Earth's orbit eccentricity

##### Microscope #####

Rmic=AUev/(1.496e11/7e6)
eta_microscope = np.sqrt(2.3**2 + 1.5**2) * 1e-15

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
