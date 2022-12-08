# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 17:28:13 2021

@author: gawe
"""

import numpy as _np


mu0 = 4.0*_np.pi*1e-7        # [H/m], approximate permeability of free space
cc = 299792458               # [m/s], defined speed of light, weird flex: 3e8 m/s is fine
eps0 = 1.0/(mu0*(cc**2.0))   # [F/m], permittivity of free space


"""
Note:
    In "rectangular" waveguide generally "a" is measured parallel to the
    traditional H-field (B-field in free space) direction for a fundamental
    TE10 mode (transverse electric wave travelling in the z-direction).
        TEmn where m, n are mode numbes
    In this sense "b" is measured parallel to the traditional D-field
    (E-field in free space) direction for a fundamental TE10 mode.

    if a>b (rectangular waveguide)
        TE10 is the dominant mode
    if a=b (square waveguide)
        TE10 and TE01 are both dominant (degenerate)
    if a<b (rectangular waveguide with long dimension parallel to the E-field
            direction)
        TE01 is the dominant mode (Transverse magnetic wave
"""

# ----

def beta2(freq, mu=None, eps=None):
    """
    freespace wavenumber squared
    """
    if mu is None: mu = mu0 # end if
    if eps is None: eps = eps0 # end if
    return (2*_np.pi*freq)**2.0*(eps*mu)

def betax(m, a):
    """
    wavenumber in the x-direction
    """
    return m*_np.pi/a

def betay(n, b):
    """
    wavenumber in the y-direction
    """
    return n*_np.pi/b

def cutoff_wavenumber(a, b, m=1, n=0):
    return _np.sqrt( betax(m, a)**2.0 + betay(n, b)**2.0 )

def cutoff_frequency(a, b, m=1, n=0, mu=None, eps=None):
    if mu is None: mu = mu0 # end if
    if eps is None: eps = eps0 # end if
    return cutoff_wavenumber(a, b, m, n)/(2.0*_np.pi*_np.sqrt(eps*mu))

# ----

def betaz(freq, a, b, m=1, n=0, mu=None, eps=None):
    """
    wavenumber in the z-direction

    Note: betaz is the guide value
    """
    bet2 = beta2(freq, mu=mu, eps=eps)
#    betax2 = betax(m, a)**2.0
#    betay2 = betay(n, b)**2.0
#    return _np.sqrt(bet2-betax2-betay2)
    return _np.sqrt(bet2-cutoff_wavenumber(a, b, m, n)**2.0)





# ----

def lambdax_TE(m, a):
    return 2*a/m

def lambday_TE(n, b):
    return 2*b/n

def lambda_freespace(freq, mu=None, eps=None):
    return 2.0*_np.pi/_np.sqrt(beta2(freq, mu=mu, eps=eps))

def lambdaz_TE(freq, a, b, m=1, n=0, mu=None, eps=None):
    """
    Guide wavelength: wavelength in the z-direction (along the waveguide)
    """
    #ilambdax2 = 1.0/lambdax_TE(m,a)**2.0
    #ilambday2 = 1.0/lambday_TE(n,b)**2.0
    #ilambda2 = 1.0/lambda_freespace(freq, mu=mu, eps=eps)**2.0
    #return 1.0/_np.sqrt( ilambda2 - ilambdax2 - ilambday2 )
    return 2.0*_np.pi/betaz_TE(freq, a, b, m=m, n=n, mu=mu, eps=eps)

# ----

def TE10_cutoff(a, mu=None, eps=None):
    if mu is None: mu = mu0 # end if
    if eps is None: eps = eps0 # end if
    return 1.0/(2.0*a*_np.sqrt(mu*eps))  # (Fc)10

def cutoff_TE(a, b, m=1, n=0, mu=None, eps=None):
    fc10 = TE10_cutoff(a, mu=mu, eps=eps)
    return fc10*_np.sqrt(m**2.0 + (n*a/b)**2.0)

# ----

def TM11_cutoff(a, b, mu=None, eps=None):
    if mu is None: mu = mu0 # end if
    if eps is None: eps = eps0 # end if
    fc10 = TE10_cutoff(a, mu=mu, eps=eps)
    return fc10*_np.sqrt(1 + (a/b)**2.0)   # (Fc)11

def cutoff_TM(a, b, m=1, n=1, mu=None, eps=None):
    fc10 = TE10_cutoff(a, mu=mu, eps=eps)
    return fc10*_np.sqrt(m**2.0 + (n*a/b)**2.0)

# ----

def wave_impedance_TEmn(freq, a, b, m=1, n=0, mu=None, eps=None):
    """
    Z_w in the +z-direction
        f>fc: real and greater than the intrinsic impedance of the medium
                in the guide

        f=fc: infinite wave impedance

        f<fc: reactively inductive    (imaginary)


        ---> rectangular waveguide is an inductive storage element for
             TEmn waves travelling in the +z direction for f<fc
    """
    if mu is None: mu = mu0 # end if
    return 2.0*_np.pi*freq*mu/betaz_TE(freq, a, b, m=m, n=n, mu=mu, eps=eps)









# ====================================================================== #
