# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 17:39:02 2017

@author: gawe
"""



import scipy as _scipy
import numpy as _np
import matplotlib.pyplot as _plt
from pybaseutils import speed_of_light

import cmath


"""
Reference 1:  Electrical properties of metal loaded radomes - Robinson (1960)


"""

freq = 1e9*_np.linspace(100.0, 250.0, 250)
#freq = 1e9*_np.linspace(10.0, 200.0, 200-10)

# ====================================================================== #
th = 45    # [deg], angle of incidence to dichroic plate (measured in free space, between ray and normal to plate surface) 
#l3 = 2.4e-3 # [m], plate thickness

## My prototype
##   fco = 125.5 GHz
#D = 1.4e-3 #[m], diameter of guide holes
#S = 1.5e-3 #[m] spacing of guide holes
#l3 = 2.4e-3
#l3 = 3e-3
#S = 2e-3

##   fco = 121.1 GHz
D = 1.45e-3 #[m], diameter of guide holes
S = 1.5e-3 #[m] spacing of guide holes
l3 = 3e-3

##   fco = 110 GHz
#D = 1.6e-3 #[m], diameter of guide holes
#S = 1.8e-3 #[m] spacing of guide holes
# l3 = 2.7e-3 # [m], plate thickness

##   fco = 113 GHz
#D = 1.55e-3 #[m], diameter of guide holes
#S = 1.7e-3 #[m] spacing of guide holes

## Radome A from reference
#D = 0.3125 #[in], diameter of guide holes
#S = 0.360  #[in] spacing of guide holes
#th = 0.0    # [deg], angle of incidence to dichroic plate (measured in free space, between ray and normal to plate surface) 
#l3 = 0.250 # [in], plate thickness
#
## Convert between metric and imperial
#D *= 25.4e-3
#S *= 25.4e-3
#l3 *= 25.4e-3

#l3 *= 1.5
#l3 *= 10

# ====================================================================== #

##### Material 1 - Free space
cc, mu0, eps0 = speed_of_light()

# =============== #

##### Material 2 - Material in contact with the metal plate
eps2 = 1.0006 # relative permeability of material in waveguide

# =============== #

##### Material 3 - Material filling the metal plate (guide holes)

Ltot = l3   # [m], total physical thickness of plate

#mur = 0.999991 # Copper, relative permeability of cavity walls
#rho = 1.724e-8 # Copper, ohm-meter, resistivity of walls of cavity 
#
#mur = 1.00002 # Aluminum, relative permeability of cavity walls
#rho = 2.65e-8 # Aluminum, ohm-meter, resistivity of walls of cavity

mur = 1.05 # Brass, relative permeability of cavity walls
rho = 6.39e-8 # Brass, ohm-meter, resistivity of walls of cavity

# ============== Air filled guide ============== #
#       My prototype
eps3 = 1.0006  # relative permeability of material in the guide, air
loss_tangent = 0.0  # loss tangent of material in guide

# ============== Polystyrene filled guide ============== #
#eps3 = 2.4  # relative permeability of material in the guide, air
#loss_tangent = 0.0  # loss tangent of material in guide

## ============== Polyamide filled guide ============== #
#eps3 = 4.3  # relative permeability of material in the guide
#loss_tangent = 0.004  # loss tangent of material in guide
#
## ============== Mica filled guide ============== #
#eps3 = 5.7  # relative permeability of material in the guide
#loss_tangent = 0.000  # loss tangent of material in guide
#
## ============== Teflon (PTFE) filled guide ============== #
#eps3 = 2.1  # relative permeability of material in the guide
#loss_tangent = 0.001  # loss tangent of material in guide
#
## ============== Sapphire filled guide ============== #
#eps3 = 10.0  # relative permeability of material in the guide
#loss_tangent = 0.000  # loss tangent of material in guide
#
## ============== Fused Quartz filled guide ============== #
#eps3 = 3.78  # relative permeability of material in the guide
#loss_tangent = 0.000  # loss tangent of material in guide

# ============== Alumina ceramic ============== #
#       Randome A
#eps3 = 8.66  # relative permeability of material in the guide
#loss_tangent = 0.0018  # loss tangent of material in guide

# ====================================================================== #

# ============== Guide parameters ============== #

#A = 1   # square guide array 
A = 0.5 * _np.sqrt(3.0)  # hexagonal guide array

fco = 1e-9*1.841*cc/(_np.pi*D)  # [GHz], designed lower cut-off frequency
fcd = 1e-9*cc/(S*A)             # [GHz], diffraction limited upper cut-off frequency
  
wlco = cc/(1e9*fco)/_np.sqrt(eps3)
  
wavelength = cc/freq
wl_3 = wavelength/_np.sqrt(eps3)
guide_wl = wl_3/(1.0-(wl_3/(1.706*D))**2.0)   # guide wavelength for the TE11 mode

# ====== #

# Hole spacing must be small enough to ensure first grating lobe response lies
# at an angle of 90 degrees to the direction of propagation.  
#    For an air-dielectric:
maxS = 1.0*wavelength/A/(_np.sqrt(eps2) + _np.sin(th*_np.pi/180.0))
# maxS = 1.1547*wavelength/(1+_np.sin(th))   # reference 2

# ====== #

# Electrical length of each cavity
phi3 = 2*_np.pi*l3/guide_wl

# This is a correction that is not always applicable, maybe don't use it
# If you replace the physical length of the cavities with the electrical length:
tauOgl = 0.0022 + 0.0055*(D/wl_3)
phi3 -= 2.0*_np.pi * 2.0*tauOgl # second order correction ... not always applicable

# ====== # 

# Attenuation constant due to dielectric
alphd = _np.pi*(guide_wl/wavelength)*loss_tangent/wavelength # np/m

# Attenuation constant due to dissipation in conducting cavity walls
rhoe = 1.724e-8 # ohm-meter = resistivity of copper
alphc = 1.5e-4 * _np.sqrt(mur*rho/rhoe)*_np.sqrt(eps3/wavelength) * (guide_wl/(D*wavelength))
alphc *= 0.420 + (wavelength/(1.706*D))**2.0

# Attenuation constant 
alph3 = alphc + alphd

# Propagation constant (multiplied by plate thickness)
gl = alph3*l3 + 1j*phi3

#_plt.figure()
#_ax1 = _plt.subplot(2,1,1)
#_ax2 = _plt.subplot(2,1,2)
#_ax1.plot(1e-9*freq, 1.0-_np.exp(-alphc*l3), 'b-')
#_ax1.plot(1e-9*freq, 1.0-_np.exp(-alphd*l3), 'g-')
#_ax1.plot(1e-9*freq, 1.0-_np.exp(-alph3*l3), 'r-')[ii]
#_ax2.plot(1e-9*freq, 1.0-_np.exp(-1.0*_np.abs(gl)), 'r-')

# ======================================= #
# Reference 1:
#
## Admittance of free space for waves polarized perpendicular and parallel to plane of incidence
#Y1_perp = _np.cos(th*_np.pi/180.0)/377.0
#Y1_parr = 1.0/(377.0*_np.cos(th*_np.pi/180.0))
#
## ==== # 
#
## Shunt susceptance of zero-thickness perforated metal plate:
#Bs = (S/D)**2.0 * (wavelength/D) * (1.0-(1.706*D/wavelength)**2.0) # siemens
#Bs *= -1.0 * (A/377.0) * (3.0/(2.0*_np.pi))  # = 1.096e-3, hexagonal array
#
## ==== # 
#
## Characteristic admittance of the metal plate containing dielectric filled cavities
#C3 = 1.522   # constant of proportionality attributed to Marcuvitz referenced in reference 1
#
#J1prime = _scipy.special.jvp(v=1, z=_np.pi*D/(4.0*S), n=1) # 1st Derivative of bessel function of 1st kind, of order 1
#Y3 = (1.0 - (0.426*D/S)**2.0 )/(2.0*J1prime)
#Y3 = Y3**2.0
#Y3 *= (S/D)**2.0 * (wavelength/guide_wl) 
#Y3 *= A*C3/377.0  # siemens
#
## Circuit parameter propagation
#T2_perp = _np.zeros_like(freq)
#T2_parr = _np.zeros_like(freq)
#R2_perp = _np.zeros_like(freq)
#R2_parr = _np.zeros_like(freq)
#ph_perp = _np.zeros_like(freq)
#ph_parr = _np.zeros_like(freq)
#for ii in range(len(freq)):
#    ABCD1 = _np.zeros( (2,2), dtype=complex)
#    ABCD2 = _np.zeros( (2,2), dtype=complex)
#
#    ABCD1[0,0] = 1.0
#    ABCD1[0,1] = 0.0    
#    ABCD1[1,0] = 1j*Bs[ii]
#    ABCD1[1,1] = 1.0
#[ii]
#    ABCD2[0,0] = _np.cosh(gl[ii])
#    ABCD2[0,1] = _np.sinh(gl[ii])/Y3[ii]
#    ABCD2[1,0] = Y3[ii]*_np.sinh(gl[ii])
#    ABCD2[1,1] = _np.cosh(gl[ii])
#    
#    ABCD = _np.dot(ABCD1, _np.dot(ABCD2, ABCD1))
#
#    perp = ABCD[0,0]+ABCD[0,1]*Y1_perp+ABCD[1,0]/Y1_perp+ABCD[1,1]
#    parr = ABCD[0,0]+ABCD[0,1]*Y1_parr+ABCD[1,0]/Y1_parr+ABCD[1,1]
#    
#    # Power transmission coefficient
#    T2_perp[ii] = 4.0/_np.abs(perp.copy())**2.0
#    T2_parr[ii] = 4.0/_np.abs(parr.copy())**2.0
#
#    # Power reflection coefficient
##    R2_perp[ii] = ((ABCD[0,0]+ABCD[0,1]*Y1_perp-ABCD[1,0]/Y1_perp-ABCD[1,1])/perp)**2.0
##    R2_parr[ii] = ((ABCD[0,0]+ABCD[0,1]*Y1_parr-ABCD[1,0]/Y1_parr-ABCD[1,1])/parr)**2.0
#    
#    # Insertion delay - Phase delay caused by guide (degrees)
#    ph_perp[ii] = _np.arctan(_np.imag(perp) /_np.real(perp)) - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength[ii]  # degrees
#    ph_parr[ii] = _np.arctan(_np.imag(parr) /_np.real(parr)) - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength[ii]
## end for
#
#R2_perp = 1.0-T2_perp
#R2_parr = 1.0-T2_parr
#
# ======================================= #
## Reference 2:
## Above cutoff, power transmission for normal incidence
#Y1_perp = 2.652e-3 # siemens (mho = inverse Ohm), free space admittance
#Y1_parr = 2.652e-3 # siemens (mho = inverse Ohm), free space admittance
#Bs = 1.096e-3*(S/D)**2.0 * (wavelength/D) * (1.0-(1.706*D/wavelength)**2.0) # siemens
#
#
#J1prime = _scipy.special.jvp(v=1, z=_np.pi*D/(4.0*S), n=1)
#Y2 = 3.496e-3*(S/D)**2.0 * ( (1.0 - (0.426*D/S)**2.0 )/(2.0*J1prime) )*(wavelength/guide_wl) # siemens
#
### Above cut-off the power transmission for NORMAL incidence:
#beta = gl.copy()
#
#Cb = _np.cos(beta*Ltot)
#Sb = _np.sin(beta*Ltot)
#
#T2_perp = (Y1_perp/Y2)*Sb + 2.0*Bs*Cb/Y1_perp + Y2*Sb/Y1_perp - (Bs**2.0)*Sb/(Y1_perp*Y2)
#T2_perp = T2_perp**2.0 + 4.0*(Cb-Bs*Sb/Y2)**2.0
#T2_perp *= 4.0
#
#T2_parr = (Y1_parr/Y2)*Sb + 2.0*Bs*Cb/Y1_parr + Y2*Sb/Y1_parr - (Bs**2.0)*Sb/(Y1_parr*Y2)
#T2_parr = T2_parr**2.0 + 4.0*(Cb-Bs*Sb/Y2)**2.0
#T2_parr *= 4.0
#
#ph_perp = _np.zeros_like(T2_perp)
#ph_parr = _np.zeros_like(T2_parr)
#
#R2_perp = 1.0-T2_perp
#R2_parr = 1.0-T2_parr
#    
# ======================================= #
# Reference 3:  Chen
# Circular openings with Equilateral triangular lattice

#       For 0.5*D>0.28*S   and S<0.57 * wavelength
J1prime = _scipy.special.jvp(v=1, z=4.0*_np.pi*0.5*D/(_np.sqrt(3)*S), n=1)
A = 12.0 * _np.sqrt(_np.asarray(4.0/3.0 * (wavelength/S)**2.0 - 1.0, dtype=complex)) \
    * (J1prime/(1.0-(4*_np.pi*0.5*D/(1.841*_np.sqrt(3.0)*S))**2.0))**2.0
A -= 12.0/_np.sqrt(_np.asarray(4.0/3.0 * (wavelength/S)**2.0 - 1.0, dtype=complex)) \
    * (J1prime/(4.0*_np.pi*0.5*D/(_np.sqrt(3.0)*S)))**2.0    

B = 0.33*(S/(0.5*D))**2.0 * _np.sqrt(_np.asarray((0.293*wavelength/(0.5*D))**2.0 - 1.0, dtype=complex) )

#beta = (0.293*wavelength/(0.5*D))**2.0 - 1.0
#beta[beta>=0], beta[beta<0] = _np.sqrt( beta[beta>=0] ), 1j*_np.sqrt( -1*beta[beta<0] )
#beta *= (2.0*_np.pi/wavelength)

beta = (2.0*_np.pi/wavelength)*_np.sqrt(_np.asarray((0.293*wavelength/(0.5*D))**2.0 - 1.0, dtype=complex))

def coth(val):
    return 1.0/cmath.tanh(val)

R2 = _np.zeros( (len(freq),), dtype=complex)
T2 = _np.zeros_like(R2)
for ii in range(len(freq)):    
    R2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*l3))) + 1.0/(1.0-1j*(A[ii]+B[ii]*coth(beta[ii]*l3))) - 1.0
    T2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*l3))) - 1.0/(1.0-1j*(A[ii]+B[ii]*coth(beta[ii]*l3)))
    print(_np.abs(R2[ii]), _np.abs(1-T2[ii]))

# For oblique incidence, there is a correction here:
porosity = _np.pi*D**2.0 / (2*_np.sqrt(3)*S**2.0)
T2_perp = T2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
T2_parr = T2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))

#T2_perp = -20.0*_np.log10(T2_perp)
#T2_parr = -20.0*_np.log10(T2_parr)
#T2_perp = 10**(T2_perp/20.0)
#T2_parr = 10**(T2_parr/20.0)

ph_perp = _np.zeros_like(T2_perp)
ph_parr = _np.zeros_like(T2_parr)

R2_perp = _np.zeros_like(T2_perp)
R2_parr = _np.zeros_like(T2_parr)
R2_perp[_np.abs(T2_perp)<1] = 1.0-_np.abs(T2_perp[_np.abs(T2_perp)<1])
R2_parr[_np.abs(T2_parr)<1] = 1.0-_np.abs(T2_parr[_np.abs(T2_parr)<1])

#R2_perp = (1.0-_np.sqrt(T2_perp)**2.0)
#R2_parr = (1.0-_np.sqrt(T2_parr)**2.0)

# ======================================= #

_plt.figure()
_plt.plot(1e-9*freq, 20*_np.log10(_np.abs(T2_perp)), '--')
_plt.plot(1e-9*freq, 20*_np.log10(_np.abs(T2_parr)), '-')

# ======================================= #
    
_plt.figure()

_ax1 = _plt.subplot(3,2,1)
_ax1.plot(1e-9*freq, _np.abs(T2_perp))
_ax1.plot(1e-9*freq, _np.abs(T2), 'r--')

_ax1.set_ylabel('Pow. Trans. Coeff.')
_ax1.set_title('Perpendicular incidence')
_ax1.axvline(x=fco, linestyle='--')
_ax1.axvline(x=fcd, linestyle='--')

_ax2 = _plt.subplot(3,2,3, sharex=_ax1)
_ax2.plot(1e-9*freq, _np.abs(R2_perp))
_ax2.plot(1e-9*freq, _np.abs(R2), 'r--')
_ax2.set_ylabel('Pow. Refl Coeff')
_ax2.axvline(x=fco, linestyle='--')
_ax2.axvline(x=fcd, linestyle='--')

_ax3 = _plt.subplot(3,2,5, sharex=_ax1)
_ax3.plot(1e-9*freq, ph_perp)
_ax3.set_ylabel('Phase Delay')
_ax3.set_xlabel('Freq [GHz]')
_ax3.axvline(x=fco, linestyle='--')
_ax3.axvline(x=fcd, linestyle='--')

_ax4 = _plt.subplot(3,2,2, sharex=_ax1)
_ax4.plot(1e-9*freq, _np.abs(T2_parr))
_ax4.plot(1e-9*freq, _np.abs(T2), 'r--')
_ax4.set_title('Parrallel incidence')
_ax4.axvline(x=fco, linestyle='--')
_ax4.axvline(x=fcd, linestyle='--')

_ax5 = _plt.subplot(3,2,4, sharex=_ax1)
_ax5.plot(1e-9*freq, _np.abs(R2_parr))
_ax5.plot(1e-9*freq, _np.abs(R2), 'r--')
_ax5.axvline(x=fco, linestyle='--')
_ax5.axvline(x=fcd, linestyle='--')

_ax6 = _plt.subplot(3,2,6, sharex=_ax1)
_ax6.plot(1e-9*freq, ph_parr)
_ax6.set_xlabel('Freq [GHz]')
_ax6.axvline(x=fco, linestyle='--')
_ax6.axvline(x=fcd, linestyle='--')

# ======================================= #

#Qo = # quality of resonant guide cavity 
#Qd = 1.0/loss_tangent  # Q due to dissipation in the guide dielectric 
## Qc = 1.0/Qo - 1.0/Qd  # quality of resonant cavity due to dissipation loss in the metal wall
#
#Ql = # loaded Q, after taking into account all losses
#
#Qe = # Q t´due to energey coupled out of the guides into space
#
#
#T = 1.0-Ql/Qo















