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

# ============== #

DEBUG = True
# DEBUG = False

# ============== #

wd = _os.path.abspath(_os.path.curdir)
#wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
# wd = _os.path.join('G://','Workshop','ECE','QMF','OP2','Dichroic Plate')

OD = 28e-3    # m, outer dimension of waveguide
ID = 27.8e-3  # m, inner dimension of waveguide


# ====================================================================== #

# Number of blades in shutter
Nslipblades = int(2) # fixed number of offset slip-blades (always 2)
Nblades = int(6)  # total number of blades = 2 * number of slots and pins per half-plate

# Blade inner diameter
bladeID = ID  # m, the inner diameter of each blade should match the waveguide diameter

# ==== half-plate design section

# alignment ridge depth
ridge_depth = 3e-3 # m

# Slot design
Nslots = int(Nblades // 2)
Npins  = int(Nblades // 2)

slot_depth = 2e-3 # m

#
# Half-plate thickness (total)
delta = 2e-3  # m, extra thickness between slots and wg input


Ltot = delta + slot_depth + ridge_depth +    # [m], total physical thickness of half-plate


# sketch
length = 8.3e-2 # cm
width = 6e-2 # cm
offset = thickness

amaj = 0.5*ID
bmin = 0.5*ID
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

if not DEBUG:
    delimiter = '\n'
    hdr = "Dichroic plate characteristics: Filled with %s"%(matname,) + delimiter
    hdr += "Hexagonal hole pattern (%i holes): diameter=%2.2f mm, spacing=%2.2f mm, thickness=%2.1f mm"%(ncircles, 1e3*D, 1e3*S, 1e3*l3) + delimiter
    hdr += "filter cut-offs: %3.1f<f<%3.1f GHz"%(fco, fcd) + delimiter
    hdr += "Power transmission (perpendicular): %3.1f dB@%3.0f GHz"%(T2_perp_140, 140) + delimiter
    hdr += "Power transmission (parallel): %3.1f dB@%3.0f GHz"%(T2_parr_140, 140) + delimiter
    hdr += "Porosity limit (%0.2f): %3.1f dB"%(porosity, por_log) + delimiter

    print(hdr)

    filnam = _os.path.join(wd,'DichroicPlate_holes_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fco,1e3*D,1e3*S,1e3*l3))
    _np.savetxt(filnam, 1e3*centers, fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n%6s 6%s'%('x[mm]', 'y[mm]') )

    filnam = _os.path.join(wd,'DichroicPlate_Transmission_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fco,1e3*D,1e3*S,1e3*l3))
    _np.savetxt(filnam, (freq,T2_parr,T2_perp), fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n %8s %8s %8s'%('freq[GHz]','T2[parr]', 'T2[perp]'))
    #_np.savetxt(filnam, _np.asarray((freq,T2_parr,T2_perp), dtype=float).T, fmt='%7.3e %6.3f %6.3f', delimiter=' ', newline='\n', header=hdr + '\n %8s %8s %8s'%('freq[GHz]','T2[parr]', 'T2[perp]'))
# end if
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

if not DEBUG:
    hfig.savefig(_os.path.join(wd,'DichroicPlate_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.png'%(matname, fco,1e3*D,1e3*S,1e3*l3)), dpi=200, transparent=True)
# end if

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

if not DEBUG:
    hfig.savefig(_os.path.join(wd,'DichroicPlate_AngleEffect_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.png'%(matname, fco,1e3*D,1e3*S,1e3*l3)), dpi=200, transparent=True)
# end if

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















