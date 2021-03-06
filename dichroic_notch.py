# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 17:39:02 2017

@author: gawe
"""



import scipy as _scipy
import numpy as _np
import os as _os
import matplotlib.pyplot as _plt
from pybaseutils import speed_of_light

import cmath


"""
Reference 1:  Electrical properties of metal loaded radomes - Robinson (1960)


"""

freq = 1e9*_np.linspace(100.0, 250.0, 250)
#freq = 1e9*_np.linspace(10.0, 200.0, 200-10)

##### Material 1 - Free space
cc, mu0, eps0 = speed_of_light()

wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Notch')

# ====================================================================== #
th = 45    # [deg], angle of incidence to dichroic plate (measured in free space, between ray and normal to plate surface) 

## My prototype
##   fco = 125.5 GHz
#D = 1.4e-3 #[m], diameter of guide holes
#S = 1.5e-3 #[m] spacing of guide holes

#fco = _np.array([139e9, 141e9], dtype=float)
fco = _np.array([135e9, 141☺e9, 145e9], dtype=float)
#fco = _np.array([125e9, 139.5e9], dtype=float)
#fco = _np.array([150e9, 130e9], dtype=float)
D = (1.841*cc)/(fco*_np.pi)  # [m], diameter of guide holes
#S = _np.round(10e3*D)/10e3
#S = _np.array([1.5e-3, 1.5e-3], dtype=float)
#S = D+0.4e-3
S = D+1e-9
thickness = 5e-3

# ====================================================================== #

##### Material 1 - Free space
cc, mu0, eps0 = speed_of_light()

# ====================================================================== #

wavelength = cc/freq

# ======================================= #
# Reference 3:  Chen
# Circular openings with Equilateral triangular lattice

def coth(val):
    return 1.0/cmath.tanh(val)

def dichroic_plate(radius, spacing, thickness):    
    #       For radius>0.28*spacing   and spacing<0.57 * wavelength
    A = 0.5 * _np.sqrt(3.0)  # hexagonal guide array  
    fc1 = 1e-9*1.841*cc/(_np.pi*2.0*radius)  # [GHz], designed lower cut-off frequency
    fc2 = 1e-9*cc/(spacing*A)             # [GHz], diffraction limited upper cut-off frequency  
    # wlco = cc/(1e9*fc1)/_np.sqrt(eps3)

    J1prime = _scipy.special.jvp(v=1, z=4.0*_np.pi*radius/(_np.sqrt(3)*spacing), n=1)
    A = 12.0 * _np.sqrt(_np.asarray(4.0/3.0 * (wavelength/spacing)**2.0 - 1.0, dtype=complex)) \
        * (J1prime/(1.0-(4*_np.pi*radius/(1.841*_np.sqrt(3.0)*spacing))**2.0))**2.0
    A -= 12.0/_np.sqrt(_np.asarray(4.0/3.0 * (wavelength/spacing)**2.0 - 1.0, dtype=complex)) \
        * (J1prime/(4.0*_np.pi*radius/(_np.sqrt(3.0)*spacing)))**2.0    
    
    B = 0.33*(spacing/radius)**2.0 * _np.sqrt(_np.asarray((0.293*wavelength/radius)**2.0 - 1.0, dtype=complex) )
    
    beta = (2.0*_np.pi/wavelength)*_np.sqrt(_np.asarray((0.293*wavelength/radius)**2.0 - 1.0, dtype=complex))

    R2 = _np.zeros( (len(freq),), dtype=complex)
    T2 = _np.zeros_like(R2)
    for ii in range(len(freq)):    
        R2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness))) + 1.0/(1.0-1j*(A[ii]+B[ii]*coth(beta[ii]*thickness))) - 1.0
        T2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness))) - 1.0/(1.0-1j*(A[ii]+B[ii]*coth(beta[ii]*thickness)))
        # print(_np.abs(R2[ii]), _np.abs(1-T2[ii]))
    # end for

    # Porosity
    por = _np.pi*(2.0*radius)**2.0 / (2.0*_np.sqrt(3)*spacing**2.0)

    T2 = _np.abs(T2)
    R2 = _np.abs(R2)
        
    print("Dichroic plate characteristics: ")
    print("Hexagonal hole pattern: diameter=%2.2f mm, spacing=%2.2f mm, thickness=%2.2f mm"%(1e3*2.0*radius, 1e3*spacing, 1e3*thickness))
    print("filter cut-offs: %3.1f<f<%3.1f GHz"%(fc1, fc2))
#    return T2perp, T2parr, por, fc1, fc2
    return T2, R2, por, fc1, fc2

def correct_angle(T2, radius, spacing, por, th):
    # For oblique incidence, there is a correction here:

    T2perp = T2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-por))
    T2parr = T2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-por))
    return T2perp, T2parr


T2perp, R2perp = [_np.zeros((len(freq),3), dtype=float) for ii in range(2)]
T2parr, R2parr = [_np.zeros((len(freq),3), dtype=float) for ii in range(2)]
porosity, fco, fcd = [_np.zeros((3,), dtype=float) for ii in range(3)]

#T2perp, R2perp = [_np.zeros((len(freq),2), dtype=float) for ii in range(2)]
#T2parr, R2parr = [_np.zeros((len(freq),2), dtype=float) for ii in range(2)]
#porosity, fco, fcd = [_np.zeros((2,), dtype=float) for ii in range(3)]

T2, R2, porosity[0], fco[0], fcd[0] = dichroic_plate(0.5*D[0].copy(), S[0].copy(), thickness)
T2perp[:,0], T2parr[:,0] = correct_angle(T2, 0.5*D[0].copy(), S[0].copy(), porosity[0], th)
R2perp[:,0], R2parr[:,0] = correct_angle(R2, 0.5*D[0].copy(), S[0].copy(), porosity[0], th)

T2, R2, porosity[1], fco[1], fcd[1] = dichroic_plate(0.5*D[1].copy(), S[1].copy(), thickness)
T2perp[:,1], T2parr[:,1] = correct_angle(T2, 0.5*D[1].copy(), S[1].copy(), porosity[1], th)
R2perp[:,1], R2parr[:,1] = correct_angle(R2, 0.5*D[1].copy(), S[1].copy(), porosity[1], th)

T2, R2, porosity[2], fco[2], fcd[2] = dichroic_plate(0.5*D[2].copy(), S[2].copy(), thickness)
T2perp[:,2], T2parr[:,2] = correct_angle(T2, 0.5*D[2].copy(), S[2].copy(), porosity[2], th)
R2perp[:,2], R2parr[:,2] = correct_angle(R2, 0.5*D[2].copy(), S[2].copy(), porosity[2], th)

#    T2perp = 20*_np.log10(T2perp)
#    T2parr = 20*_np.log10(T2parr)    
    
# Notch with 1 big dichroic filter, 1 med. dichroic filter, 1 small dichroic filter, and 1 solid plate
T2perp = R2perp[:,0]*R2perp[:,1] + T2perp[:,0]*T2perp[:,2]
T2parr = R2parr[:,0]*R2parr[:,1] + T2parr[:,0]*T2parr[:,2]
    
#T2_perp_log = -T2_perp_log[:,0]+T2_perp_log[:,1]
#T2_parr_log = -T2_parr_log[:,0]+T2_parr_log[:,1]

#T2_perp_log = T2_perp_log[:,0]*(1.0-T2_perp_log[:,1])
#T2_parr_log = T2_parr_log[:,0]*(1.0-T2_parr_log[:,1])

#T2_perp_log = T2_perp_log[:,0]*T2_perp_log[:,1]
#T2_parr_log = T2_parr_log[:,0]*T2_parr_log[:,1]

#T2_perp_log = _np.abs(1-T2_perp_log[:,0])*_np.abs(1-T2_perp_log[:,1])
#T2_parr_log = _np.abs(1-T2_parr_log[:,0])*_np.abs(1-T2_parr_log[:,1])

#T2_perp_log = _np.abs(1-T2_perp_log)
#T2_parr_log = _np.abs(1-T2_parr_log)

T2_perp_log = 20*_np.log10(T2perp)
T2_parr_log = 20*_np.log10(T2parr)

#T2_perp_log = 20*_np.log10(_np.abs(T2_perp_log[:,0])) + 20*_np.log10(_np.abs(T2_perp_log[:,1]))
#T2_parr_log = 20*_np.log10(_np.abs(T2_parr_log[:,0])) + 20*_np.log10(_np.abs(T2_parr_log[:,1]))
    
# ======================================= #

hfig = _plt.figure()
_plt.plot(1e-9*freq, T2_perp_log, '-')
_plt.plot(1e-9*freq, T2_parr_log, '--')
_plt.xlim((105,180))

_plt.xlabel('frequency [GHz]')
_plt.ylabel(r'|T$^2$| [dB]')
#_plt.title(r'Power Transmission Coefficient: f$_{c,o}$<%3.1f, f$_{c,d}$<%3.1f GHz'%(fco,fcd) )
_plt.axvline(x=fco[0], linestyle='--', color='k')
_plt.axvline(x=fco[1], linestyle='--', color='k')


# ==== # 
import os as _os
#wd = _os.path.abspath(_os.path.curdir)
wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
hfig.savefig(_os.path.join(wd,'DichroicNotch_%3.1fGHz_%3.1fGHz.png'%(fco[0],fco[1])), dpi=200, transparent=True)

# ======================================= #
    

