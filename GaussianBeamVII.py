# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 19:00:57 2023

@author: gawe

Based on GaussianBeamVI
    'GaussianBeamVI - copyright @ wwang2@caltech.edu'

A nice little utility for visualizing the effect of two lens

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import text

# Define plots and range
fig, ax = plt.subplots()
fig.canvas.manager.set_window_title('GaussianBeamVII')
plt.subplots_adjust(bottom=0.50)

# beam waist function
def beam(w,lamda,start,stop):
    """ Trace the beam waist in the forward direction

    """
    x = np.arange(start, stop, res)
    return w*(1+ (lamda*(x-start)/(np.pi*w**2))**2)**0.5

# Plot the back trace of the beam
def backbeam(m,w,lamda,start,stop):
    """ Trace the beam waist from the end to the beginning

    """
    zr = np.pi*(w*m)**2/lamda
    x = np.arange(start, stop, res)
    return m*w*(1+ ((stop-x)/zr)**2)**0.5

# Compute the waist position
def imageposition(w,lamda,f,s):
    """ Calculate the imagine position of a thin lens

    """
    zr = np.pi*w**2/lamda
    return 1/(1/f - 1/(s+zr**2/(s-f)))

# Compute the magnification
def mag(w,lamda,f,s):
    """ Magnification of output image from a thin lens

    inputs
        w - [m] - input beam waist
        lambda - [m] - wavelength
        f - [m]  - focal length of optical element
        s - [m] - position of optical element

    outputs
        mag - Magnification. Essentially wout/win
    """
    zr = np.pi*w**2/lamda
    return 1/(((1-s/f)**2+(zr/f)**2)**0.5)

# Reset
def reset(event):
    spos1.reset()
    sf1.reset()
    sw0.reset()


# Define initial setuplength and plot resolution
length = 800.0
res = 0.1

cc = 299792458  # [m/s], speed of light in vacuum
# lambda_max = 300 # [mm]   1 mm => 300 GHz, 300 mm => 1 GHz, or micro-meters: 1000 GHz

um_units = True
if um_units:
    # Working units in micro-meters
    #
    # initial beam and lens position
    lamda0 = 782e-6   # wavelength in free-space
    pos1ini = 51      # initial position of lens 1
    f1ini = 50.0      # initial focal length of lens 1
    wini = 0.5        # initial waist input to lens 1

    # 2nd lens initial parameter
    f2ini = 50.0
    pos2ini = 151.0


else:
    # Working units in micro-meters
    #
    # initial beam and lens position
    freq0 = 140.0e9     # [Hz], frequency of interest
    lamda0 = cc/(freq0)   # wavelength in free-space
    pos1ini = 51      # initial position of lens 1
    f1ini = 50.0      # initial focal length of lens 1
    wini = 0.5        # initial waist input to lens 1

    # 2nd lens initial parameter
    f2ini = 50.0
    pos2ini = 151.0

m1ini = mag(wini, lamda0, f1ini, pos1ini)



# Compute Initial magnification and beamwaist position
im1posini = imageposition(wini,lamda0,f1ini,pos1ini)
im2posini = imageposition(wini*m1ini,lamda0,f2ini,pos2ini-im1posini-pos1ini)
w1ini = m1ini*wini
m2ini = mag(w1ini,lamda0,f2ini,pos2ini-im1posini-pos1ini)
w2ini = m2ini*w1ini

# Initialize Plots with initial parameters
beam1down ,= plt.plot(np.arange(0,pos1ini,res),-beam(wini,lamda0,0.0,pos1ini),color = "Blue")
beam1up ,= plt.plot(np.arange(0,pos1ini,res),beam(wini,lamda0,0.0,pos1ini),color = "Blue")

beam2down ,= plt.plot(np.arange(pos1ini+im1posini,pos2ini,res),
                      -beam(w1ini,lamda0,pos1ini+im1posini,pos2ini),color = "Blue")
beam2up ,= plt.plot(np.arange(pos1ini+im1posini,pos2ini,res),
                    beam(w1ini,lamda0,pos1ini+im1posini,pos2ini),color = "Blue")

bbeam1down ,= plt.plot(np.arange(pos1ini,pos1ini+im1posini,res),
                       -backbeam(m1ini,wini,lamda0,pos1ini,pos1ini+im1posini),color = "Blue")
bbeam1up ,= plt.plot(np.arange(pos1ini,pos1ini+im1posini,res),
                     backbeam(m1ini,wini,lamda0,pos1ini,pos1ini+im1posini),color = "Blue")

beam3down ,= plt.plot(np.arange(pos2ini+im2posini,length,res),
                      -beam(w2ini,lamda0,pos2ini+im2posini,length),color = "Blue")
beam3up ,= plt.plot(np.arange(pos2ini+im2posini,length,res),
                      beam(w2ini,lamda0,pos2ini+im2posini,length),color = "Blue")

bbeam2down ,= plt.plot(np.arange(pos2ini,pos2ini+im2posini,res),
                       -backbeam(m2ini,w1ini,lamda0,pos2ini,pos2ini+im2posini),color = "Blue")

bbeam2up ,= plt.plot(np.arange(pos2ini,pos2ini+im2posini,res),
                       backbeam(m2ini,w1ini,lamda0,pos2ini,pos2ini+im2posini),color = "Blue")

im1 ,= plt.plot([pos1ini+im1posini,pos1ini+im1posini], [-w1ini,w1ini])
im2 ,= plt.plot([pos2ini+im2posini,pos2ini+im2posini], [-w2ini,w2ini])

lens1 ,= plt.plot([pos1ini,pos1ini],[-2,2])
lens2 ,= plt.plot([pos2ini,pos2ini],[-2,2])

plt.axis([0, length, -2, 2])

ax.xaxis.set_ticks(np.arange(0, length, 50))
ax.yaxis.set_ticks(np.arange(-2.5, 2.5, 0.25))

ax.tick_params(labeltop=True, labelright=True)

axcolor = 'lightgoldenrodyellow'

# Define wavelength Slider
axlamda  = plt.axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)
slamda = Slider(axlamda, 'wavelent', 200, 1200, valinit=782)

# Define lens position slider
axpos1  = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
axpos2 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

spos1 = Slider(axpos1, 'position1', 0.0, length, valinit=pos1ini)
spos2 = Slider(axpos2, 'position2', 0.0, length, valinit=pos2ini)

# Define initial beam wasit slider
axw0 = plt.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
sw0 = Slider(axw0, 'beam waist', 0.0, 2.0, valinit=wini)

# Define lens1 focus slider
axf1 = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
sf1 = Slider(axf1, 'lens 1 focus', 0.0, 300, valinit=f1ini)  # 300 mm -> 1 GHz, 2 mm -> 140 GHz, 1 mm --> 300 GHz

# Define lens2 foucs slider
axf2 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sf2 = Slider(axf2, 'lens 2 focus', 0.0, 300, valinit=f2ini)

# Update plots in response to sliders
def update(val):

    if um_units:
        # Convert wavelength to micrometers
        lamda = slamda.val* 0.000001
    else:
        # Convert wavelength to mm
        lamda = slamda.val * 1e-3

    # print(lamda)
    pos1 = spos1.val
    pos2 = spos2.val
    w0 = sw0.val
    f1 = sf1.val
    f2 = sf2.val

    m1 = mag(w0,lamda,f1,pos1)
    im1pos = imageposition(w0,lamda,f1,pos1)
    w1 = m1*w0

    m2 = mag(w1,lamda,f2,pos2-im1pos-pos1)
    im2pos = imageposition(w1,lamda,f2,pos2-im1pos-pos1)
    w2 = m2*w1

    # print(w2)

    # Plot the lenses
    lens1.set_data([pos1, pos1], [-2, 2])
    lens2.set_data([pos2, pos2], [-2, 2])

    # Plot the ray-traces of the beam
    beam1up.set_data(np.arange(0.0, pos1,res),np.array(beam(w0,lamda,0.0,pos1)))
    bbeam1up.set_data(np.arange(pos1,pos1+im1pos,res),np.array(backbeam(m1,w0,lamda,pos1,pos1+im1pos)))

    beam1down.set_data(np.arange(0.0, pos1,res),-np.array(beam(w0,lamda,0.0,pos1)))
    bbeam1down.set_data(np.arange(pos1,pos1+im1pos,res),-np.array(backbeam(m1,w0,lamda,pos1,pos1+im1pos)))

    beam2up.set_data(np.arange(pos1+im1pos,pos2,res),np.array(beam(w1,lamda,pos1+im1pos,pos2)))
    beam2down.set_data(np.arange(pos1+im1pos,pos2,res),-np.array(beam(w1,lamda,pos1+im1pos,pos2)))

    beam3up.set_data(np.arange(pos2+im2pos,length,res),np.array(beam(w2,lamda,pos2+im2pos,length)))
    beam3down.set_data(np.arange(pos2+im2pos,length,res),-np.array(beam(w2,lamda,pos2+im2pos,length)))

    bbeam2up.set_data(np.arange(pos2,pos2+im2pos,res),
                       backbeam(m2,w1,lamda,pos2,pos2+im2pos))
    bbeam2down.set_data(np.arange(pos2,pos2+im2pos,res),
                       -backbeam(m2,w1,lamda,pos2,pos2+im2pos))

    im1.set_data([im1pos+pos1,im1pos+pos1],[-w1,w1])
    im2.set_data([pos2+im2pos,pos2+im2pos], [-w2,w2])

    fig.canvas.draw_idle()

# Call the update function once slider changes
sw0.on_changed(update)
spos1.on_changed(update)
sf1.on_changed(update)
spos2.on_changed(update)
sf2.on_changed(update)
slamda.on_changed(update)

# Define Reset Button
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
button.on_clicked(reset)

plt.show()