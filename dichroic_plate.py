# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 17:39:02 2017

@author: gawe
"""



#import scipy as _scipy
from scipy.special import jvp
import numpy as _np
import os as _os
import matplotlib.pyplot as _plt
from pybaseutils import speed_of_light
import cmath
#
#from matplotlib.patches import Ellipse

"""
Reference 1:  Electrical properties of metal loaded radomes - Robinson (1960)


"""

#wd = _os.path.abspath(_os.path.curdir)
#wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
wd = _os.path.join('G://','Workshop','ECE','QMF','OP2','Dichroic Plate')

freq = 1e9*_np.linspace(100.0, 250.0, 250)
#freq = 1e9*_np.linspace(10.0, 200.0, 200-10)

# ====================================================================== #
th = 45    # [deg], angle of incidence to dichroic plate (measured in free space, between ray and normal to plate surface)
#th = 90    # [deg], angle of incidence to dichroic plate (measured in free space, between ray and normal to plate surface)
#l3 = 2.4e-3 # [m], plate thickness

#l3 = 1e-3
#l3 = 3e-3
l3 = 5e-3
#l3 = 15e-3
#l3 = 20e-3
#l3 = 3.0e-3
#l3 = _np.round(1e3*0.5*l3)/1e3
thickness = l3

# My prototype for OP2 CECE protection
#D = 1.30e-3 #[m], diameter of guide holes, 135.1 GHz
D = 1.27e-3 #[m], diameter of guide holes, 135.1 GHz
#D = 1.25e-3 #[m], diameter of guide holes, 140.5 GHz
#S = 1.4e-3 #[m] spacing of guide holes
S = 1.45e-3 #[m] spacing of guide holes
#l3 = 3.0e-3

##   fco = 146.4 GHz
#D = 1.20e-3 #[m], diameter of guide holes
#S = 1.4e-3 #[m] spacing of guide holes
##l3 = 3.0e-3

##   fco = 125.5 GHz
#D = 1.4e-3 #[m], diameter of guide holes
#S = 1.5e-3 #[m] spacing of guide holes
#S = 1.6e-3 #[m] spacing of guide holes
#S = 1.8e-3 #[m] spacing of guide holes
#S = 2.0e-3 #[m] spacing of guide holes
#l3 = 3.0e-3


###   fco = 121.1 GHz
#D = 1.45e-3 #[m], diameter of guide holes
#S = 1.5e-3 #[m] spacing of guide holes
#l3 = 3e-3

#####   fco = 117.1 GHz
#D = 1.5e-3 #[m], diameter of guide holes
##S = 1.6e-3 #[m] spacing of guide holes
#S = 1.7e-3 #[m] spacing of guide holes  - second prototype
##l3 = 5e-3

## My prototype for OP1.2a Reflectometry protection + ECE signal diplexing
#   fco = 110 GHz
#D = 1.6e-3 #[m], diameter of guide holes
#S = 1.8e-3 #[m] spacing of guide holes, 0.1 mm wall spacing too small. 0.2 acceptable by shop
#S = 2.0e-3 #[m] spacing of guide holes, 0.1 mm wall spacing too small. 0.2 acceptable by shop
#l3 = 3.0e-3 # [m], plate thickness

##   fco = 113 GHz
#D = 1.55e-3 #[m], diameter of guide holes
#S = 1.7e-3 #[m] spacing of guide holes
#l3 = 3.0e-3 # [m], plate thickness

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
matname = "Air"
eps3 = 1.0006  # relative permeability of material in the guide, air
loss_tangent = 0.0  # loss tangent of material in guide

## ============== Polystyrene filled guide ============== #
#matname = "Polystyrene"
#eps3 = 2.4  # relative permeability of material in the guide, air
#loss_tangent = 0.0  # loss tangent of material in guide

## ============== Polyamide filled guide ============== #
#matname = "Polyamide"
#eps3 = 4.3  # relative permeability of material in the guide
#loss_tangent = 0.004  # loss tangent of material in guide
#
## ============== Mica filled guide ============== #
#matname = "Mica"
#eps3 = 5.7  # relative permeability of material in the guide
#loss_tangent = 0.000  # loss tangent of material in guide
#
## ============== Teflon (PTFE) filled guide ============== #
#matname = "PTFE"
#eps3 = 2.1  # relative permeability of material in the guide
#loss_tangent = 0.001  # loss tangent of material in guide
#
## ============== Sapphire filled guide ============== #
#matname = "Saphire"
#eps3 = 10.0  # relative permeability of material in the guide
#loss_tangent = 0.000  # loss tangent of material in guide
#
## ============== Fused Quartz filled guide ============== #
#matname = "Quartz"
#eps3 = 3.78  # relative permeability of material in the guide
#loss_tangent = 0.000  # loss tangent of material in guide

# ============== Alumina ceramic ============== #
#       Randome A
#matname = "Alumina"
#eps3 = 8.66  # relative permeability of material in the guide
#loss_tangent = 0.0018  # loss tangent of material in guide

# ============== Macor ceramic ============== #
#matname = "Macor"
#eps3 = 5.67  # relative permeability of material in the guide
#loss_tangent = 0.0071  # loss tangent of material in guide

# ====================================================================== #

# ============== Guide parameters ============== #

#A = 1   # square guide array
A = 0.5 * _np.sqrt(3.0)  # hexagonal guide array

fco = 1e-9*1.841*cc/(_np.pi*D)  # [GHz], designed lower cut-off frequency
fcd = 1e-9*cc/(S*A)             # [GHz], diffraction limited upper cut-off frequency

wlco = cc/(1e9*fco)/_np.sqrt(eps3)
wl140 = cc/(140e9)/_np.sqrt(eps3)

wavelength = cc/freq
wl_3 = wavelength/_np.sqrt(eps3)
guide_wl = _np.ones((len(freq),), dtype=complex)
guide_wl *= wl_3/(1.0-(wl_3/(1.706*D))**2.0)   # guide wavelength for the TE11 mode in cylindrical waveguide
guide_140 = wl140/(1.0-(wl140/(1.706*D))**2.0)   # guide wavelength for the TE11 mode

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
#tauOgl = 0.0022 + 0.0055*(D/wl_3)
#phi3 -= 2.0*_np.pi * 2.0*tauOgl # second order correction ... not always applicable

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

## ======================================= #
## Reference 1:
#
## Admittance of free space for waves polarized perpendicular and parallel to plane of incidence
#Y1_perp = _np.cos(th*_np.pi/180.0)/377.0
#Y1_parr = 1.0/(377.0*_np.cos(th*_np.pi/180.0))
#
## ==== #
#
## Shunt susceptance of zero-thickness perforated metal plate:
#Bs = _np.ones( (len(freq),), dtype=complex)
#Bs *= (S/D)**2.0 * (wavelength/D) * (1.0-(1.706*D/wavelength)**2.0) # siemens
#Bs *= -1.0 * (A/377.0) * (3.0/(2.0*_np.pi))  # = 1.096e-3, hexagonal array
#
## ==== #
#
## Characteristic admittance of the metal plate containing dielectric filled cavities
#C3 = 1.522   # constant of proportionality attributed to Marcuvitz referenced in reference 1
#
##J1prime = jvp(v=1, z=_np.pi*D/(4.0*S), n=1) # 1st Derivative of bessel function of 1st kind, of order 1
#J1prime = jvp(v=1, z=4.0*_np.pi*0.5*D/(_np.sqrt(3)*S), n=1)
#Y3 = _np.ones((len(freq),), dtype=complex)
#Y3 *= (1.0 - (0.426*D/S)**2.0 )/(2.0*J1prime)
#Y3 = Y3**2.0
#Y3 *= (S/D)**2.0 * (wavelength/guide_wl)
#Y3 *= A*C3/377.0  # siemens
#
## Circuit parameter propagation
#T2_perp = _np.zeros( (len(freq),), dtype=complex)
#T2_parr = _np.zeros( (len(freq),), dtype=complex)
#R2_perp = _np.zeros( (len(freq),), dtype=complex)
#R2_parr = _np.zeros( (len(freq),), dtype=complex)
#ph_perp = _np.zeros( (len(freq),), dtype=complex)
#ph_parr = _np.zeros( (len(freq),), dtype=complex)
#for ii in range(len(freq)):
#    ABCD1 = _np.zeros( (2,2), dtype=complex)
#    ABCD2 = _np.zeros( (2,2), dtype=complex)
#
#    ABCD1[0,0] = 1.0
#    ABCD1[0,1] = 0.0
#    ABCD1[1,0] = 1j*Bs[ii]
#    ABCD1[1,1] = 1.0
#
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
#    R2_perp[ii] = ((ABCD[0,0]+ABCD[0,1]*Y1_perp-ABCD[1,0]/Y1_perp-ABCD[1,1])/perp)**2.0
#    R2_parr[ii] = ((ABCD[0,0]+ABCD[0,1]*Y1_parr-ABCD[1,0]/Y1_parr-ABCD[1,1])/parr)**2.0
#
#    # Insertion delay - Phase delay caused by guide (degrees)
#    ph_perp[ii] = _np.arctan(_np.imag(perp) /_np.real(perp)) - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength[ii]  # degrees
#    ph_parr[ii] = _np.arctan(_np.imag(parr) /_np.real(parr)) - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength[ii]
## end for
#
##R2_perp = 1.0-T2_perp
##R2_parr = 1.0-T2_parr

# ======================================= #
## Reference 2:
## Above cutoff, power transmission for normal incidence
#Y1_perp = 2.652e-3 # siemens (mho = inverse Ohm), free space admittance
#Y1_parr = 2.652e-3 # siemens (mho = inverse Ohm), free space admittance
#Bs = 1.096e-3*(S/D)**2.0 * (wavelength/D) * (1.0-(1.706*D/wavelength)**2.0) # siemens
#
#J1prime = jvp(v=1, z=_np.pi*D/(4.0*S), n=1)
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
J1prime = jvp(v=1, z=4.0*_np.pi*0.5*D/(_np.sqrt(3)*S), n=1)
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
#ph = _np.zeros( (len(freq),), dtype=float)
for ii in range(len(freq)):

    AA = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*l3)))
    BB = 1.0/  (1.0 - 1j*(A[ii]+B[ii]*      coth(beta[ii]*l3)))
    # Reflection
    R2[ii] = AA.copy() + BB.copy() - 1.0

    # Transmission
    T2[ii] = AA.copy() - BB.copy()

    # Insertion delay - Phase delay caused by guide (degrees)
#    ph[ii] = _np.arctan(_np.imag(T2[ii]) /_np.real(T2[ii])) - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength[ii]  # degrees
    print(_np.abs(R2[ii]), _np.abs(T2[ii])) #, ph[ii])
#    print(_np.abs(R2[ii]), _np.abs((1-_np.sqrt(T2[ii]))**2.0))

# For oblique incidence, there is a correction here:
porosity = _np.pi*D**2.0 / (2*_np.sqrt(3)*S**2.0)
T2_perp = T2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
T2_parr = T2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))

R2_perp = R2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
R2_parr = R2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))

ph_perp = _np.arctan(_np.imag(T2_perp) /_np.real(T2_perp))*180.0/_np.pi - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength  # degrees
ph_parr = _np.arctan(_np.imag(T2_parr) /_np.real(T2_parr))*180.0/_np.pi - 360.0*Ltot*_np.cos(th*_np.pi/180.0)/wavelength  # degrees

T2_perp = _np.abs(T2_perp)
T2_parr = _np.abs(T2_parr)

R2_perp = _np.abs(R2_perp)
R2_parr = _np.abs(R2_parr)

#T2_perp *= -1
#T2_parr *= -1

#ph_perp = _np.zeros_like(T2_perp)
#ph_parr = _np.zeros_like(T2_parr)

#R2_perp = _np.zeros_like(T2_perp)
#R2_parr = _np.zeros_like(T2_parr)
#R2_perp[_np.abs(T2_perp)<1] = 1.0-_np.abs(T2_perp[_np.abs(T2_perp)<1])
#R2_parr[_np.abs(T2_parr)<1] = 1.0-_np.abs(T2_parr[_np.abs(T2_parr)<1])

R2_perp = 1.0-T2_perp
R2_parr = 1.0-T2_parr

#R2_perp = (1.0-_np.sqrt(T2_perp)**2.0)
#R2_parr = (1.0-_np.sqrt(T2_parr)**2.0)

T2_perp_log = 20*_np.log10(T2_perp)
T2_parr_log = 20*_np.log10(T2_parr)
por_log = 10*_np.log10(porosity)

T2_perp_140 = _np.interp(140,1e-9*freq,T2_perp_log)
T2_parr_140 = _np.interp(140,1e-9*freq,T2_parr_log)

# ======================================= #

# sketch
length = 8.3e-2 # cm
width = 6e-2 # cm
offset = thickness

amaj = 0.5*3.9e-2
bmin = 0.5*2.8e-2
box_start = [1.4e-2, 2.0e-2]

Nvert = 2.0*bmin / (S*_np.sin(60.0*_np.pi/180.0))
Nhoriz = 2.0*amaj / S

Nvert = int(_np.round(Nvert))
Nhoriz = int(_np.round(Nhoriz))

print(Nvert, Nhoriz)

# =========================================================================== #

def hexagon_generator(edge_length, offset):
  """Generator for coordinates in a hexagon."""
  npts = 6
  x = _np.zeros((npts,), dtype=float)
  y = _np.zeros((npts,), dtype=float)
  angle = _np.linspace(30, 390, 7)
  for ii in range(npts):
    x[ii] = offset[0] + edge_length*_np.cos(angle[ii]*_np.pi/180.0)
    y[ii] = offset[1] + edge_length*_np.sin(angle[ii]*_np.pi/180.0)
  return x, y

def closest_approach(shape1, shape2):
    minDist = 100
    n1 = len(shape1[:,0])
    n2 = len(shape2[:,0])
    for ii in range(n1):
        for jj in range(n2):
            dist = _np.sqrt( (shape1[ii,0]-shape2[jj,0])**2.0 + (shape1[ii,1]-shape2[jj,1])**2.0 )
            minDist = min(minDist, dist)
        # end for
    # end for
    return minDist


# =========================================================================== #

angle = _np.linspace(0, 2*_np.pi, 180)
elli = _np.vstack((amaj*_np.cos(angle), bmin*_np.sin(angle))).T
Dscrew = 0.004

hfig = _plt.figure()

# 4 plate walls - plan view
_plt.plot((-1e3*length/2, -1e3*length/2), (-1e3*width/2, 1e3*width/2), 'k-')
_plt.plot(( 1e3*length/2,  1e3*length/2), (-1e3*width/2, 1e3*width/2), 'k-')
_plt.plot((-1e3*length/2,  1e3*length/2), (-1e3*width/2,-1e3*width/2), 'k-')
_plt.plot((-1e3*length/2,  1e3*length/2), ( 1e3*width/2, 1e3*width/2), 'k-')

# 4 x 4mm bolt holes
_plt.plot(1e3*length/2-1e3*1.3e-2 +1e3*0.5*Dscrew*_np.cos(angle), 1e3*width/2-1e3*0.85e-2+1e3*0.5*Dscrew*_np.sin(angle), 'k-')
_plt.plot(1e3*length/2-1e3*1.3e-2 +1e3*0.5*Dscrew*_np.cos(angle), -1e3*width/2+1e3*0.85e-2+1e3*0.5*Dscrew*_np.sin(angle), 'k-')
_plt.plot(-1e3*length/2+1e3*1.3e-2 +1e3*0.5*Dscrew*_np.cos(angle), 1e3*width/2-1e3*0.85e-2+1e3*0.5*Dscrew*_np.sin(angle), 'k-')
_plt.plot(-1e3*length/2+1e3*1.3e-2 +1e3*0.5*Dscrew*_np.cos(angle), -1e3*width/2+1e3*0.85e-2+1e3*0.5*Dscrew*_np.sin(angle), 'k-')
_plt.plot(1e3*elli[:,0], 1e3*elli[:,1], 'k--')  # ellipse
#_plt.axvline(x= 1e3*length/2-1e3*1.3e-2, color='k', linestyle='--')
#_plt.axvline(x=-1e3*length/2+1e3*1.3e-2, color='k', linestyle='--')
#_plt.axhline(y= 1e3*width/2-1e3*0.85e-2, color='k', linestyle='--')
#_plt.axhline(y=-1e3*width/2+1e3*0.85e-2, color='k', linestyle='--')

# 4 plate walls - lower projection (side view)
_plt.plot((-1e3*length/2, -1e3*length/2), (-1e3*thickness/2-1e3*offset-1e3*0.5*width, 1e3*thickness/2-1e3*offset-1e3*0.5*width), 'k-')
_plt.plot(( 1e3*length/2,  1e3*length/2), (-1e3*thickness/2-1e3*offset-1e3*0.5*width, 1e3*thickness/2-1e3*offset-1e3*0.5*width), 'k-')
_plt.plot((-1e3*length/2,  1e3*length/2), (-1e3*thickness/2-1e3*offset-1e3*0.5*width,-1e3*thickness/2-1e3*offset-1e3*0.5*width), 'k-')
_plt.plot((-1e3*length/2,  1e3*length/2), ( 1e3*thickness/2-1e3*offset-1e3*0.5*width, 1e3*thickness/2-1e3*offset-1e3*0.5*width), 'k-')

# 4 plate walls - right projection (side view)
_plt.plot((-1e3*thickness/2+1e3*(offset+length/2), -1e3*thickness/2+1e3*(offset+length/2)), (-1e3*width/2, 1e3*width/2), 'k-')
_plt.plot(( 1e3*thickness/2+1e3*(offset+length/2),  1e3*thickness/2+1e3*(offset+length/2)), (-1e3*width/2, 1e3*width/2), 'k-')
_plt.plot((-1e3*thickness/2+1e3*(offset+length/2),  1e3*thickness/2+1e3*(offset+length/2)), (-1e3*width/2,-1e3*width/2), 'k-')
_plt.plot((-1e3*thickness/2+1e3*(offset+length/2),  1e3*thickness/2+1e3*(offset+length/2)), ( 1e3*width/2, 1e3*width/2), 'k-')

xrow = S*_np.cos(60.*_np.pi/180.)
ycol = S*_np.sin(60.*_np.pi/180.)

#   1.6 x 1.80 mm
# odd - odd : 272
# even - even : 272
# odd - even : 281
# even -odd : 281
if Nvert%2>0: # odd
    # hole spacing symmetric about x=0 line, no point at x=0
    voffset = S*_np.sin(60.*_np.pi/180.)
else:
    # point at x=0
    voffset = 0.0
# endif

if Nhoriz%2==0: # even
    hoffset = 0.5*S
else: #odd
    hoffset = 0.0
# endif

ncircles = 0
centers = list()
for ii in range(Nvert):
    for jj in range(Nhoriz):
        xcen = S*(jj-Nhoriz/2)+hoffset
        ycen = S*_np.sin(60.*_np.pi/180.)*(ii-Nvert/2) + voffset
        if ii%2>0:
            xcen += xrow
        # end if
        circ = _np.vstack((xcen+0.5*D*_np.cos(angle), ycen+0.5*D*_np.sin(angle))).T
        ybound = (_np.abs(ycen)+0.5*D)<_np.abs(bmin/amaj)*_np.sqrt(_np.abs(amaj**2.0-xcen**2.0))
        xbound = _np.abs(xcen)<_np.abs(amaj-0.5*D)
        if ybound and xbound and closest_approach(_np.atleast_2d([xcen,ycen]), elli) >= 0.5*D:
            xhex, yhex = hexagon_generator(S*_np.tan(30.0*_np.pi/180.0), (xcen,ycen))

#            _plt.plot(1e3*xhex, 1e3*yhex, 'k-')
            _plt.plot(1e3*circ[:,0], 1e3*circ[:,1], 'k-')

            centers.append([xcen, ycen])
            ncircles += 1
    # end for
# end for
centers = _np.asarray(centers)

#ax.set_xlim((-1e3*(0.5*length+0.1), 1e3*(0.5*length+0.1)))
#ax.set_ylim((-1e3*(0.5*width+0.1), 1e3*(0.5*width+0.1)))
#_plt.axis((-1e3*(0.5*length+0.1), 1e3*(0.5*length+0.1),-1e3*(0.5*width+0.1), 1e3*(0.5*width+0.1)))
#_plt.axis('equal')
_plt.xlim((-1e3*(0.5*length+offset), 1e3*(0.5*length+2*offset+thickness)))
_plt.ylim((-1e3*(0.5*width+2*offset+thickness), 1e3*(0.5*width+offset)))
print(ncircles)

#_plt.axis('off')
_plt.title('%0.1f mm Dichroic Plate:  %3.1f GHz < f < %3.1f GHz \n S=%0.2f mm, D=%0.2f mm, N=%i holes'%(1e3*thickness, fco, fcd, 1e3*S, 1e3*D, ncircles))
#hfig.savefig(_os.path.join(wd,'DichroicPlate_drawing_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.png'%(matname, fco,1e3*D,1e3*S,1e3*l3)), dpi=200, transparent=True)

# ======================================= #

delimiter = '\n'
hdr = "Dichroic plate characteristics: Filled with %s"%(matname,) + delimiter
hdr += "Hexagonal hole pattern (%i holes): diameter=%2.2f mm, spacing=%2.2f mm, thickness=%2.1f mm"%(ncircles, 1e3*D, 1e3*S, 1e3*l3) + delimiter
hdr += "filter cut-offs: %3.1f<f<%3.1f GHz"%(fco, fcd) + delimiter
hdr += "Power transmission (perpendicular): %3.1f dB@%3.0f GHz"%(T2_perp_140, 140) + delimiter
hdr += "Power transmission (parallel): %3.1f dB@%3.0f GHz"%(T2_parr_140, 140) + delimiter
hdr += "Porosity limit (%0.2f): %3.1f dB"%(porosity, por_log) + delimiter

print(hdr)

filnam = _os.path.join(wd,'DichroicPlate_holes_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fco,1e3*D,1e3*S,1e3*l3))
#_np.savetxt(filnam, 1e3*centers, fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n%6s 6%s'%('x[mm]', 'y[mm]') )

filnam = _os.path.join(wd,'DichroicPlate_Transmission_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fco,1e3*D,1e3*S,1e3*l3))
#_np.savetxt(filnam, (freq,T2_parr,T2_perp), fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n %8s %8s %8s'%('freq[GHz]','T2[parr]', 'T2[perp]'))
#_np.savetxt(filnam, _np.asarray((freq,T2_parr,T2_perp), dtype=float).T, fmt='%7.3e %6.3f %6.3f', delimiter=' ', newline='\n', header=hdr + '\n %8s %8s %8s'%('freq[GHz]','T2[parr]', 'T2[perp]'))

# ======================================= #

#hfig = _plt.figure('thickness_scan')
hfig = _plt.figure()

_plt.plot(1e-9*freq, T2_perp_log, '-')
_plt.plot(1e-9*freq, T2_parr_log, '--')
xlims = _plt.xlim()
xlims = (xlims[0],210)
ylims = _plt.ylim()
#ylims = (ylims[0], 0.0)
ylims = (-30, 0.0)
_plt.xlim(xlims)
_plt.ylim(ylims)

_plt.xlabel('frequency [GHz]')
_plt.ylabel(r'|T$^2$| [dB]')
_plt.title(r'Power Transmission Coefficient: f$_{c,o}$<%3.1f, f$_{c,d}$<%3.1f GHz'%(fco,fcd) )
_plt.axvline(x=fco, linestyle='--', color='k')
_plt.axvline(x=fcd, linestyle='--', color='k')
_plt.axhline(y=por_log, linestyle='--', color='k')

_plt.text(x=fco+5, y=-15, s='Hexagonal hole pattern: \n diameter=%2.2f mm, \n spacing=%2.2f mm, \n thickness=%2.1f mm'%(1e3*D, 1e3*S, 1e3*l3))
#_plt.text(x=fco+5, y=ylims[1]-15, s='Hexagonal hole pattern: \n diameter=%2.2f mm, \n spacing=%2.2f mm'%(1e3*D, 1e3*S))
#_plt.text(x=xlims[0]+5.0, y=ylims[1]-20, s=' thickness=%2.1f mm'%(1e3*l3,))

# ==== #

#hfig.savefig(_os.path.join(wd,'DichroicPlate_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.png'%(matname, fco,1e3*D,1e3*S,1e3*l3)), dpi=200, transparent=True)

# ======================================= #

hfig = _plt.figure(figsize=(8,3.5))

_ax1 = _plt.subplot(1,2,1)
_ax1.set_position([ 0.125,  0.15, 0.35,  0.75])
_ax1.plot(1e-9*freq, _np.abs(T2_perp))
_ax1.plot(1e-9*freq, _np.abs(T2), 'r--')

_ax1.set_ylabel('Pow. Trans. Coeff.')
_ax1.set_title('Perpendicular Polarization')
_ax1.set_ylim((0,1))
_ax1.axvline(x=fco, linestyle='--')
_ax1.axvline(x=fcd, linestyle='--')
_ax1.set_xlabel('Freq [GHz]')

_ax4 = _plt.subplot(1,2,2, sharex=_ax1)
_ax4.set_position([ 0.55,  0.15, 0.35,  0.75])

_ax4.plot(1e-9*freq, _np.abs(T2_parr))
_ax4.plot(1e-9*freq, _np.abs(T2), 'r--')
_ax4.set_title('Parrallel Polarization')
_ax4.set_ylim((0,1))
_ax4.axvline(x=fco, linestyle='--')
_ax4.axvline(x=fcd, linestyle='--')
_ax4.set_xlabel('Freq [GHz]')

#_ax2 = _plt.subplot(2,2,3, sharex=_ax1)
#_ax2.plot(1e-9*freq, _np.abs(R2_perp))
#_ax2.plot(1e-9*freq, _np.abs(R2), 'r--')
#_ax2.set_ylabel('Pow. Refl Coeff')
#_ax2.set_ylim((0,1))
#_ax2.axvline(x=fco, linestyle='--')
#_ax2.axvline(x=fcd, linestyle='--')
#_ax2.set_xlabel('Freq [GHz]')
#
#_ax5 = _plt.subplot(2,2,4, sharex=_ax1)
#_ax5.plot(1e-9*freq, _np.abs(R2_parr))
#_ax5.plot(1e-9*freq, _np.abs(R2), 'r--')
#_ax5.set_ylim((0,1))
#_ax5.axvline(x=fco, linestyle='--')
#_ax5.axvline(x=fcd, linestyle='--')
#_ax5.set_xlabel('Freq [GHz]')

#_ax3 = _plt.subplot(2,2,5, sharex=_ax1)
#_ax3.plot(1e-9*freq, ph_perp.real)
#_ax3.set_ylabel('Phase Delay')
#_ax3.set_xlabel('Freq [GHz]')
#_ax3.set_ylim((0,1))
#_ax3.axvline(x=fco, linestyle='--')
#_ax3.axvline(x=fcd, linestyle='--')

#
#_ax6 = _plt.subplot(3,2,6, sharex=_ax1)
#_ax6.plot(1e-9*freq, ph_parr.real)
#_ax6.set_xlabel('Freq [GHz]')
#_ax6.set_ylim((0,1))
#_ax6.axvline(x=fco, linestyle='--')
#_ax6.axvline(x=fcd, linestyle='--')

#hfig.savefig(_os.path.join(wd,'DichroicPlate_AngleEffect_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.png'%(matname, fco,1e3*D,1e3*S,1e3*l3)), dpi=200, transparent=True)

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















