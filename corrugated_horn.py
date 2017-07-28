# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 17:42:08 2017

@author: gawe
"""


"""
 // ex from 
 A Novel Quasi-Optical Frequency Multiplier Design for Millimeter and Submillimeter Wavelengths
    by J. W. Archer
  DOI: 10.1109/TMTT.1984.1132693     
 Published in: IEEE Transactions on Microwave Theory and Techniques ( Volume: 32, Issue: 4, Apr 1984 ) 



Theoretical properties of the Corrugated horn as a function of operating frequency
 a - aperture radius
 d - slot depth
 delt - slot width / pitch ratio
 th_o - horn semi-angle
 gamma - mode content factor (ratio of longitudinal fields of TE11 and TM11 components)
 Delta - aperture phase deviation, edge to axis

corrugation slot period should be close to a half-wavelength
 - better cross-polarization performance observed with thinner wall between corrugations
 
 H-plane 1/e full-width is essentially frequency independent at 29.3 degrees
"""

import numpy as _np
try:
    from pybaseutils.utils import speed_of_light
except:
    from ..utils import speed_of_light  
# end try      

a = 0.295    # aperture radius
Nslots = 38  # number of corrugations
Dslots = 0.012 # depth of corrugations 
th_o = 22    # horn semi-angle
 
freq = _np.linspace(50, 180, 1)  # GHz 
cc, _, _ = speed_of_light()
wavelength = cc / freq  # free space wavelength


tmp1 = 2.0*_np.pi*a*th_o/wavelength
tmp2 = 2.0*_np.pi*d / wavelength

Delta = lambda wl, a, th_o: a/wl * _np.tan(th_o*_np.pi/360.0)
