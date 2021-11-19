# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:01:28 2021

@author: gawe

Functions dealing with rectangular patch antenna.

"""
import math
import numpy as np
from math import cos, sin, sqrt, pi, log10, atan2, acos, radians

from scipy import integrate
import scipy.integrate
import json


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# from pybaseutils.utils import sph2cart as sph2cart1
# from pybaseutils.utils import cart2sph as cpatchart2sph1

# constants
light_velocity = 299792458
impedance = 50

# ======================================== #

# import plotly
# from plotly.offline import iplot
# import plotly.graph_objs as go
# plotly.offline.init_notebook_mode(connected=True)


def S_i(a):
    temp = scipy.integrate.quad(lambda x:sin(x)/x,0,a)
    return temp[0]


def J0(s):
    temp = scipy.integrate.quad(lambda x:cos(s*sin(x)),0,pi)
    temp = (1/pi)*temp[0]
    return temp


def get_k(f):
    lamda_0 = light_velocity/f
    k0 = (2*pi)/lamda_0
    return k0


def getG1 (W, f):
    k0 = get_k (f)
    X = k0 * W
    I1 = -2 + cos(X) + X*S_i(X) + sin(X)/X
    G1 = I1 / ( 120 * pi**2 )
    return G1


def getG12 (W, k0, L):
    temp = scipy.integrate.quad(lambda x: (((sin(k0*W*cos(x)/2)/cos(x))**2)*J0(k0*L*sin(x))*sin(x)**3), 0, pi)
    G12 = (1/(120*pi**2))*temp[0]
    return G12


def getGs(f, W, L):
    G1 = getG1(W, f)
    k0 = get_k(f)
    G12 = getG12(W, k0, L)
    return G1, G12


def input_impedance (f, W, L):
    k0 = get_k (f)
    G1, G12 = getGs(f, W, L)
    Rin = 1/(2*(G1+G12))
    print("Input Impedance:", Rin, "ohms")
    return Rin


def inset_feed_position(Rin, L):
    # R = 50.0
    R = impedance
    y0 = (L/pi)*(math.acos(sqrt(R/Rin)))
    return y0


def get_directivity(G1, G12, W, f, I1, I2):
    lamda_0 = light_velocity/f
    g_12 = G12/G1
    D_AF = 2/(1+g_12)
    D0 = ((2*pi*W)/lamda_0)**2*(1/I1)
    D2 = D0 * D_AF
    DIR_1 = 10*log10(D2)
    D_2 = ((2*pi*W)/lamda_0) ** 2 * (pi/I2)
    DIR_2 = 10 * log10(D_2)
    return DIR_1, DIR_2

# ======================================== #


def PatchFunction(thetaInDeg, phiInDeg, Freq, W, L, h, Er):
    """
    Taken from Design_patchr
    Calculates total E-field pattern for patch as a function of theta and phi
    Patch is assumed to be resonating in the (TMx 010) mode.
    E-field is parallel to x-axis

    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate

    Refrence C.A. Balanis 2nd Edition Page 745
    """
    lamba = light_velocity / Freq

    theta_in = math.radians(thetaInDeg)
    phi_in = math.radians(phiInDeg)

    ko = 2 * math.pi / lamba

    xff, yff, zff = sph2cart1(999, theta_in, phi_in)                            # Rotate coords 90 deg about x-axis to match array_utils coord system with coord system used in the model.
    xffd = zff
    yffd = xff
    zffd = yff
    r, thp, php = cart2sph1(xffd, yffd, zffd)
    phi = php
    theta = thp

    if theta == 0:
        theta = 1e-9                                                              # Trap potential division by zero warning

    if phi == 0:
        phi = 1e-9

    Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * (1 + 12 * (h / W)) ** -0.5        # Calculate effictive dielectric constant for microstrip line of width W on dielectric material of constant Er

    F1 = (Ereff + 0.3) * (W / h + 0.264)                                        # Calculate increase length dL of patch length L due to fringing fields at each end, giving total effective length Leff = L + 2*dL
    F2 = (Ereff - 0.258) * (W / h + 0.8)
    dL = h * 0.412 * (F1 / F2)

    Leff = L + 2 * dL

    Weff = W                                                                    # Calculate effective width Weff for patch, uses standard Er value.
    heff = h * sqrt(Er)

    # Patch pattern function of theta and phi, note the theta and phi for the function are defined differently to theta_in and phi_in
    Numtr2 = sin(ko * heff * cos(phi) / 2)
    Demtr2 = (ko * heff * cos(phi)) / 2
    Fphi = (Numtr2 / Demtr2) * cos((ko * Leff / 2) * sin(phi))

    Numtr1 = sin((ko * heff / 2) * sin(theta))
    Demtr1 = ((ko * heff / 2) * sin(theta))
    Numtr1a = sin((ko * Weff / 2) * cos(theta))
    Demtr1a = ((ko * Weff / 2) * cos(theta))
    Ftheta = ((Numtr1 * Numtr1a) / (Demtr1 * Demtr1a)) * sin(theta)

    # Due to groundplane, function is only valid for theta values :   0 < theta < 90   for all phi
    # Modify pattern for theta values close to 90 to give smooth roll-off, standard model truncates H-plane at theta=90.
    # PatEdgeSF has value=1 except at theta close to 90 where it drops (proportional to 1/x^2) to 0

    rolloff_factor = 0.5                                                       # 1=sharp, 0=softer
    theta_in_deg = theta_in * 180 / math.pi                                          # theta_in in Deg
    F1 = 1 / (((rolloff_factor * (abs(theta_in_deg) - 90)) ** 2) + 0.001)       # intermediate calc
    PatEdgeSF = 1 / (F1 + 1)                                                    # Pattern scaling factor

    UNF = 1.0006                                                                # Unity normalisation factor for element pattern

    if theta_in <= math.pi / 2:
        Etot = Ftheta * Fphi * PatEdgeSF * UNF                                   # Total pattern by pattern multiplication
    else:
        Etot = 0

    return Etot

def patch_function(theta_in_deg, phi_in_deg, freq, w, l, h, er):
    """
    Taken from Design_patchr
    Calculates total E-field pattern for patch as a function of theta and phi
    Patch is assumed to be resonating in the (TMx 010) mode.
    E-field is parallel to x-axis
    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate
    Refrence C.A. Balanis 2nd Edition Page 745
    """
    lambda_ = light_velocity / freq
    theta_in = math.radians(theta_in_deg)
    phi_in = math.radians(phi_in_deg)

    ko = 2 * math.pi / lambda_

    xff, yff, zff = sph2cart1(999, theta_in, phi_in)                            # Rotate coords 90 deg about x-axis to match array_utils coord system with coord system used in the model.
    xffd = zff
    yffd = xff
    zffd = yff
    r, thp, php = cart2sph1(xffd, yffd, zffd)
    phi = php
    theta = thp

    if theta == 0:
        # Trap potential division by zero warning
        theta = 1e-9

    if phi == 0:
        phi = 1e-9

    # Calculate effective dielectric constant for micro_strip line of width W on dielectric material of constant Er
    e_ref = ((er + 1) / 2) + ((er - 1) / 2) * (1 + 12 * (h / w)) ** -0.5

    # Calculate increase length dL of patch length L due to fringing fields at each end,
    # giving total effective length Leff = L + 2*dL

    f1 = (e_ref + 0.3) * (w / h + 0.264)
    f2 = (e_ref - 0.258) * (w / h + 0.8)
    d_l = h * 0.412 * (f1 / f2)

    l_eff = l + 2 * d_l

    # Calculate effective width Weff for patch, uses standard Er value.
    w_eff = w
    h_eff = h * sqrt(er)

    # Patch pattern function of theta and phi,
    # Note the theta and phi for the function are defined differently to theta_in and phi_in
    num_tr_2 = sin(ko * h_eff * cos(phi) / 2)
    dem_tr_2 = (ko * h_eff * cos(phi)) / 2
    f_phi = (num_tr_2 / dem_tr_2) * cos((ko * l_eff / 2) * sin(phi))

    num_tr_1 = sin((ko * h_eff / 2) * sin(theta))
    dem_tr_1 = ((ko * h_eff / 2) * sin(theta))
    num_tr_1a = sin((ko * w_eff / 2) * cos(theta))
    dem_tr_1a = ((ko * w_eff / 2) * cos(theta))
    f_theta = ((num_tr_1 * num_tr_1a) / (dem_tr_1 * dem_tr_1a)) * sin(theta)

    # Due to groundplane, function is only valid for theta values :   0 < theta < 90   for all phi
    # Modify pattern for theta values close to 90 to give smooth roll-off, standard model truncates H-plane at theta=90.
    # PatEdgeSF has value=1 except at theta close to 90 where it drops (proportional to 1/x^2) to 0

    # 1=sharp, 0=softer
    roll_off_factor = 0.5
    # theta_in in Deg
    theta_in_deg = theta_in * 180 / math.pi
    # intermediate calc
    f1 = 1 / (((roll_off_factor * (abs(theta_in_deg) - 90)) ** 2) + 0.001)
    # Pattern scaling factor
    pat_edge_sf = 1 / (f1 + 1)
    # Unity normalisation factor for element pattern
    UNF = 1.0006

    # Total pattern by pattern multiplication
    if theta_in <= math.pi / 2:
        e_tot = f_theta * f_phi * pat_edge_sf * UNF
    else:
        e_tot = 0

    return e_tot


def GetPatchFields(PhiStart, PhiStop, ThetaStart, ThetaStop, Freq, W, L, h, Er):
    """"
    Calculates the E-field for range of thetaStart-thetaStop and phiStart-phiStop
    Returning a numpy array of form - fields[phiDeg][thetaDeg] = eField

    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate
    """
    fields = np.ones((PhiStop, ThetaStop))                                      # Create initial array to hold e-fields for each position

    for phiDeg in range(PhiStart, PhiStop):
            for thetaDeg in range(ThetaStart, ThetaStop):                       # Iterate over all Phi/Theta combinations
                eField = PatchFunction(thetaDeg, phiDeg, Freq, W, L, h, Er)     # Calculate the field for current Phi, Theta
                fields[phiDeg][thetaDeg] = eField                               # Update array with e-field

    return fields



def get_patch_fields(phi_start, phi_stop, theta_start, theta_stop, freq, w, l, h, er):
    """"
    Calculates the E-field for range of thetaStart-thetaStop and phiStart-phiStop
    Returning a numpy array of form - fields[phiDeg][thetaDeg] = eField
    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate
    """
    # Create initial array to hold e-fields for each position
    fields = np.ones((phi_stop, theta_stop))
    # Iterate over all Phi/Theta combinations
    for phiDeg in range(phi_start, phi_stop):
        for thetaDeg in range(theta_start, theta_stop):
            # Calculate the field for current Phi, Theta
            eField = patch_function(thetaDeg, phiDeg, freq, w, l, h, er)
            # Update array with e-field
            fields[phiDeg][thetaDeg] = eField

    return fields


def PatchEHPlanePlot(Freq, W, L, h, Er, isLog=True):
    """
    Plot 2D plots showing E-field for E-plane (phi = 0°) and the H-plane (phi = 90°).
    """

    fields = GetPatchFields(0, 360, 0, 90, Freq, W, L, h, Er)                                                   # Calculate the field at each phi, theta

    Xtheta = np.linspace(0, 90, 90)                                                                             # Theta range array used for plotting

    if isLog:                                                                                                   # Can plot the log scale or normal
        plt.plot(Xtheta, 20 * np.log10(abs(fields[90, :])), label="H-plane (Phi=90°)")                          # Log = 20 * log10(E-field)
        plt.plot(Xtheta, 20 * np.log10(abs(fields[0, :])), label="E-plane (Phi=0°)")
        plt.ylabel('E-Field (dB)')
    else:
        plt.plot(Xtheta, fields[90, :], label="H-plane (Phi=90°)")
        plt.plot(Xtheta, fields[0, :], label="E-plane (Phi=0°)")
        plt.ylabel('E-Field')

    plt.xlabel('Theta (degs)')                                                                                  # Plot formatting
    plt.title("Patch: \nW=" + str(W) + " \nL=" + str(L) +  "\nEr=" + str(Er) + " h=" + str(h) + " \n@" + str(1e-9*Freq) + "GHz")
    plt.ylim(-40)
    plt.xlim((0, 90))

    start, end = plt.xlim()
    plt.xticks(np.arange(start, end, 5))
    plt.grid(b=True, which='major')
    plt.legend()
    plt.show()                                                                                                  # Show plot

    return fields                                                                                               # Return the calculated fields


def patch_eh_plane_plot(freq, w, l, h, er, is_log=True):
    """
    Plot 2D plots showing E-field for E-plane (phi = 0) and the H-plane (phi = 90).
    """

    fields = get_patch_fields(0, 360, 0, 90, freq, w, l, h, er)

    Xtheta = np.linspace(0, 90, 90)

    if is_log:
        # Log = 20 * log10(E-field)# Can plot the log scale or normal
        plt.plot(Xtheta, 20 * np.log10(abs(fields[90, :])), label="H-plane (Phi=90)")
        plt.plot(Xtheta, 20 * np.log10(abs(fields[0, :])), label="E-plane (Phi=0)")
        plt.ylabel('E-Field (dB)')
    else:
        plt.plot(Xtheta, fields[90, :], label="H-plane (Phi=90)")
        plt.plot(Xtheta, fields[0, :], label="E-plane (Phi=0)")
        plt.ylabel('E-Field')

    plt.xlabel('Theta (degs)')
    plt.title("EH Plane - Theta ")
    plt.ylim(-40)
    plt.xlim((0, 90))

    start, end = plt.xlim()
    plt.xticks(np.arange(start, end, 5))
    plt.grid(b=True, which='major')
    plt.legend()
    plt.show()
    return fields



def SurfacePlot(Fields, Freq, W, L, h, Er):
    """Plots 3D surface plot over given theta/phi range in Fields by calculating cartesian coordinate equivalent of spherical form."""

    print("Processing SurfacePlot...")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    phiSize = Fields.shape[0]                                                                                   # Finds the phi & theta range
    thetaSize = Fields.shape[1]

    X = np.ones((phiSize, thetaSize))                                                                           # Prepare arrays to hold the cartesian coordinate data.
    Y = np.ones((phiSize, thetaSize))
    Z = np.ones((phiSize, thetaSize))

    for phi in range(phiSize):                                                                                  # Iterate over all phi/theta range
        for theta in range(thetaSize):
            e = Fields[phi][theta]

            xe, ye, ze = sph2cart1(e, math.radians(theta), math.radians(phi))                                   # Calculate cartesian coordinates

            X[phi, theta] = xe                                                                                  # Store cartesian coordinates
            Y[phi, theta] = ye
            Z[phi, theta] = ze

    ax.plot_surface(X, Y, Z, color='b')                                                                         # Plot surface
    plt.ylabel('Y')
    plt.xlabel('X')                                                                                             # Plot formatting
    plt.title("Patch: \nW=" + str(W) + " \nL=" + str(L) +  "\nEr=" + str(Er) + " h=" + str(h) + " \n@" + str(1e-9*Freq) + "GHz")
    plt.show()


# def surface_plot(fields, is_note_book=False):
def surface_plot_go(fields, is_note_book=False):
    """Plots 3D surface plot over given theta/phi range in Fields by calculating cartesian
    coordinate equivalent of spherical form."""

    print("Processing SurfacePlot...")
    # Finds the phi & theta range
    phiSize = fields.shape[0]
    thetaSize = fields.shape[1]
    # Prepare arrays to hold the cartesian coordinate data.
    X = np.ones((phiSize, thetaSize))
    Y = np.ones((phiSize, thetaSize))
    Z = np.ones((phiSize, thetaSize))

    for phi in range(phiSize):
        for theta in range(thetaSize):
            e = fields[phi][theta]

            xe, ye, ze = sph2cart1(e, math.radians(theta), math.radians(phi))

            X[phi, theta] = xe
            Y[phi, theta] = ye
            Z[phi, theta] = ze
    surface = go.Surface(x=X, y=Y, z=Z)
    data = [surface]

    layout = go.Layout(
        title='Surface Plot of EH Plane',
        scene=dict(
            xaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            yaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            zaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            )
        )
    )

    fig = go.Figure(data=data, layout=layout)
    if is_note_book:
        iplot(fig)
    else:
        plotly.offline.plot(fig)


def DesignPatch(Er, h, Freq):
    """
    Returns the patch_config parameters for standard lambda/2 rectangular microstrip patch. Patch length L and width W are calculated and returned together with supplied parameters Er and h.

    Returned values are in the same format as the global patchr_config variable, so can be assigned directly. The patchr_config variable is of the following form [Er,W,L,h].
    Usage: patchr_config=design_patchr(Er,h,Freq)
    Er.....Relative dielectric constant
    h......Substrate thickness (m)
    Freq...Frequency (Hz)

    e.g. patchr_config=design_patchr(3.43,0.7e-3,2e9)
    """
    Eo = 8.854185e-12

    lambd = light_velocity / Freq
    lambdag = lambd / sqrt(Er)

    W = (light_velocity / (2 * Freq)) * sqrt(2 / (Er + 1))

    Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * (1 + 12 * (h / W)) ** -0.5                              # Calculate effictive dielectric constant for microstrip line of width W on dielectric material of constant Er

    F1 = (Ereff + 0.3) * (W / h + 0.264)                                                 # Calculate increase length dL of patch length L due to fringing fields at each end, giving actual length L = Lambda/2 - 2*dL
    F2 = (Ereff - 0.258) * (W / h + 0.8)
    dL = h * 0.412 * (F1 / F2)

    lambdag = lambd / sqrt(Ereff)
    L = (lambdag / 2) - 2 * dL

    print('Rectangular Microstrip Patch Design')
    print("Frequency (GHz): " + str(1e-9*Freq))
    print("Dielec Const, Er : " + str(Er))
    print("Patch Width,  W: " + str(W) + "m")
    print("Patch Length,  L: " + str(L) + "m")
    print("Patch Height,  h: " + str(h) + "m")

    return W, L, h, Er


def design_patch(er, h, freq):

    lambda_ = light_velocity / freq
    w = (light_velocity / (2 * freq)) * sqrt(2 / (er + 1))
    temp = 1 + 12*(h/w)

    e_ref = ((er + 1) / 2) + ((er - 1) / 2) * temp ** -0.5

    f1 = (e_ref + 0.3) * (w / h + 0.264)
    f2 = (e_ref - 0.258) * (w / h + 0.8)
    d_l = h * 0.412 * (f1 / f2)

    lambda_g = lambda_ / sqrt(e_ref)
    L = (lambda_g / 2) - 2 * d_l

    print('Rectangular Microstrip Patch Design')
    print("Frequency: (GHz) " + str(1e-9*freq))
    print("Dielec Const, Er : " + str(er))
    print("Patch Width,  W: " + str(w) + "m")
    print("Patch Length,  L: " + str(L) + "m")
    print("Patch Height,  h: " + str(h) + "m")
    return w, L



def exampleRectPatch():
    """Some example patches with various thickness & Er."""
    # print("Patch.py")

    freq = 14e9
    Er = 3.66           # RO4350B
    # Er = 2.1           #

    # h = 1.0e-3
    h = 0.101e-3
    W, L, h, Er = DesignPatch(Er, h, freq)
    fields = PatchEHPlanePlot(freq, W, L, h, Er)
    SurfacePlot(fields, freq, W, L, h, Er)

    h = 1.524e-3
    # h = 2e-3
    W, L, h, Er = DesignPatch(Er, h, freq)     # RO4350B
    fields = PatchEHPlanePlot(freq, W, L, h, Er)
    SurfacePlot(fields, freq, W, L, h, Er)

    # fields = PatchEHPlanePlot(freq, 10.7e-3, 10.47e-3, 3e-3, 2.5)                # Random
    # SurfacePlot(fields, freq, W, L, h, Er)

# ======================================== #


def ArrayFactor(ElementArray, Freq):
    """
    Summation of field contributions from each element in array, at frequency freq at theta 0°-95°, phi 0°-360°.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    Returns arrayFactor[theta, phi, elementSum]
    """

    arrayFactor = np.ones((360, 95))

    Lambda = light_velocity / Freq

    for theta in range(95):
        for phi in range(360):                                                                                                      # For all theta/phi positions
            elementSum = 1e-9 + 0j

            for element in ElementArray:                                                                                            # Summation of each elements contribution at theta/phi position.
                relativePhase = CalculateRelativePhase(element, Lambda, math.radians(theta), math.radians(phi))                     # Find relative phase for current element
                elementSum += element[3] * math.e ** ((relativePhase + element[4]) * 1j)                                            # Element contribution = Amp * e^j(Phase + Phase Weight)

            arrayFactor[phi][theta] = elementSum.real

    return arrayFactor


def CalculateRelativePhase(Element, Lambda, theta, phi):
    """
    Incident wave treated as plane wave. Phase at element is referred to phase of plane wave at origin.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    theta & phi in radians
    See Eqn 3.1 @ https://theses.lib.vt.edu/theses/available/etd-04262000-15330030/unrestricted/ch3.pdf
    """
    phaseConstant = (2 * math.pi / Lambda)

    xVector = Element[0] * math.sin(theta) * math.cos(phi)
    yVector = Element[1] * math.sin(theta) * math.sin(phi)
    zVector = Element[2] * math.cos(theta)

    phaseOfIncidentWaveAtElement = phaseConstant * (xVector + yVector + zVector)

    return phaseOfIncidentWaveAtElement

# ======================================== #


def FieldSumPatch(ElementArray, Freq, W, L, h, Er):
    """
    Summation of field contributions from each patch element in array, at frequency freq for theta 0°-95°, phi 0°-360°.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    Returns arrayFactor[theta, phi, elementSum]
    """

    arrayFactor = np.ones((360, 95))

    Lambda = light_velocity / Freq

    for theta in range(95):
        for phi in range(360):                                                                                                      # For all theta/phi positions
            elementSum = 1e-9 + 0j

            xff, yff, zff = sph2cart1(999, math.radians(theta), math.radians(phi))                                                  # Find point in far field

            for element in ElementArray:                                                                                            # For each element in array, find local theta/phi, calculate field contribution and add to summation for point
                xlocal = xff - element[0]
                ylocal = yff - element[1]                                                                                           # Calculate local position in cartesian
                zlocal = zff - element[2]

                r, thetaLocal, phiLocal = cart2sph1(xlocal, ylocal, zlocal)                                                         # Convert local position to spherical

                patchFunction = PatchFunction(math.degrees(thetaLocal), math.degrees(phiLocal), Freq, W, L, h, Er)            # Patch element pattern for local theta, phi

                if patchFunction != 0:                                                                                              # Sum each elements contribution
                    relativePhase = CalculateRelativePhase(element, Lambda, math.radians(theta), math.radians(phi))                 # Find relative phase for current element
                    elementSum += element[3] * patchFunction * math.e ** ((relativePhase + element[4]) * 1j)                        # Element contribution = Amp * e^j(Phase + Phase Weight)

            arrayFactor[phi][theta] = elementSum.real

    return arrayFactor


def FieldSumHorn(ElementArray, Freq):
    """
    Summation of field contributions from each horn element in array, at frequency freq for theta 0°-95°, phi 0°-360°.
    Horn pattern estimate using cos q(theta) function.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    Returns arrayFactor[theta, phi, elementSum]
    """

    arrayFactor = np.ones((360, 95))

    Lambda = light_velocity / Freq

    for theta in range(95):
        for phi in range(360):                                                                                                      # For all theta/phi positions
            elementSum = 1e-9 + 0j

            xff, yff, zff = sph2cart1(999, math.radians(theta), math.radians(phi))                                                  # Find point in far field

            for element in ElementArray:                                                                                            # For each element in array, find local theta/phi, calculate field contribution and add to summation for point
                if theta > 90:
                    hornFunction = 0                                                                                                # Assume no radiation behind horn for simplification
                else:
                    xlocal = xff - element[0]                                                                                       # Calculate local position in cartesian
                    ylocal = yff - element[1]
                    zlocal = zff - element[2]

                    r, thetaLocal, phiLocal = cart2sph1(xlocal, ylocal, zlocal)                                                     # Convert local position to spherical

                    # TODO: This is a random horn antenna pattern: cos^28(theta). You can use a real one or input one later
                    hornFunction = math.cos(thetaLocal) ** 28                                                                       # Horn q function, q = 28

                if hornFunction != 0:                                                                                               # Sum each elements contribution
                    relativePhase = CalculateRelativePhase(element, Lambda, math.radians(theta), math.radians(phi))                 # Find relative phase for current element
                    elementSum += element[3] * hornFunction * math.e ** ((relativePhase + element[4]) * 1j)                         # Element contribution = Amp * e^j(Phase + Phase Weight)

            arrayFactor[phi][theta] = elementSum.real

    return arrayFactor


# ======================================== #


"""
Returns the efficiency of a rectangular microstrip patch as a percentage. Based on ArrayCalc calc_patchr_eff.m.
References:
Microstrip Antennas, I.J Bahl and P.Bhartia, Published Atrech House, Page 60
Advances in Microstrip and Printed Antennas", Lee and Chen (Ch5)

Some useful numbers :

                CONDUCTORS                                      DIELECTRICS

        Material             Conductivity S/m           Material         Er     Tand

        Perfect              9.90E+99 (lossless)        FR4_Epoxy        4.4    0.02
        Silver               6.29E+07                   Arlon 25FR       3.43   0.0035
        Copper               5.80E+07                   Arlon AD300      3.00   0.003
        Pure Alumin.         3.77E+07                   Arlon AR1000    10.00   0.0035
        Al. 6063-T832        3.08E+07                   Rogers RO3003    3.00   0.0013
        Al. 6061-T6          2.49E+07                   Rogers RO3006    6.15   0.0025
        Brass                1.56E+07                   Rogers RO3010   10.20   0.0035
        Phospor bronze       9.09E+06                   Rogers RO4350    3.48   0.004
        Stainless Steel 302  1.39E+06                   Glass             5.5   0.000
                                                        Plexiglass        3.4   0.001
                                                        Polyamide         4.3   0.004
                                                        Polyester         3.2   0.003
                                                        Polyethylene      2.25  0.001
"""

def CalculatePatchEff(Er, W, L, h, tand, sigma, Freq, VSWR):
    """
    Er: Relative dielectric constant
    W: Patch width (m)
    L: Patch length (m)
    h: dielectric thickness (m)
    tand: Loss tangent of dielectric (units)
    sigma: Conductivity of patch (Siemens/m)
    Freq: Frequency (Hz)
    VSWR: VSWR for bandwidth estimate (ratio). http://www.antenna-theory.com/definitions/vswr.php describes how well the antenna is impedance matched to the line it is connected to.
    """

    if Er <= 1.000001:
        Er = 1.000001

    if tand <= 0.000001:
        tand = 0.000001

    Eo = 8.854185e-12                                                                               # Free space dielectric constant
    Ee = Eo * Er                                                                                    # Effective dielectric constant

    lamba = light_velocity / Freq

    """
    % Calculation for space and surface wave efficiency factor, gives roughly the same results.
    % Reference :  "Advances in Microstrip and Printed Antennas" Lee and Chen(Ch 5)
    % Efficiency due to surface wave component, dominant for larger h/lambda values
    """
    Mur = 1
    n1 = sqrt(Er * Mur)
    ko = 2 * pi * Freq * (sqrt((8.854e-12) * (pi * 4e-7)))
    Lo = lamba

    Psur = (1 / Lo ** 2) * ((ko * h) ** 3) * (60 * pi ** 3 * Mur ** 3 * (1 - 1 / n1 ** 2) ** 3)    # Power radiated as surface wave

    c1 = 1 - 1 / n1 ** 2 + (2 / 5) / n1 ** 4
    Pr = (1 / Lo ** 2) * (ko * h) ** 2 * (80 * pi ** 2 * Mur ** 2 * c1)                                   # Total power radiated

    Effsw = Pr / (Pr + Psur)                                                                        # Efficiency factor for surface wave losses

    """
    % Efficiency due to ohmic and dielectric losses, dominant for smaller h/lambda values
    % ***********************************************************************************
    % Reference : "Microstrip Antennas" Bahl and Bartia
    """
    if W < lamba:
        Rr = 90 * lamba ** 2 / W ** 2                                                               # Radiation resistance for W<lambda
    else:
        Rr = 120 * lamba / W                                                                        # Radiation resistance for W>lambda

    Qr = (light_velocity * sqrt(Ee)) / (2 * (Freq / 1e6) * h)                                                    # Quality factor, modified by me, not sure freq was in Ghz, more like MHz !?
    Rc = (1 / sigma) * 0.5 * sqrt(Freq) * (L / W) * Qr ** 2                                         # Equivalent resistance for conductor loss (ohms)
    Rd = (30 * tand / Er) * ((h * lamba) / (L * W)) * Qr ** 2                                       # Equivalent resistance for dielectric loss (ohms)

    Rtot = Rr + Rd + Rc                                                                             # Total resistance (ohms)
    Effcd = Rr / Rtot                                                                               # Efficiency factor for combined dielectric and ohmic losses

    Eff1 = Effsw * Effcd
    Eff = Eff1 * 100                                                                                # Total efficiency including ohmic, dielectric and surface wave losses (percent)

    Qt = Qr * Eff1 / (pi)                                                                           # Ref Balanis p762  ( Qtotal = Qradiated*Efficiency ) Not the pi factor, I added that, seems necassary get sensible results using Qr from above !?

    BW = (VSWR - 1) / (Qt * sqrt(VSWR))

    BWf = BW * Freq / 1e6                                                                           # Bandwidth as a frequency span in MHz

    print("Rectangular patch overall efficency " + str(Eff) + "%")
    print("Surface wave efficiency factor " + str(Effsw))
    print("Ohmic and dielectric efficiency factor " + str(Effcd))
    print("BW=" + str(BWf) + "MHz for VSWR=" + str(VSWR) + " at Fo=" + str(Freq / 1e6) + " MHz")

    return Eff

def examplePatchEfficiency():
    """
    Calculates efficiencies for two patches @14GHz, one with FR4 and one RO4350.
    """
    freq = 14e9                     # Same parameters for both patches
    h = 1.524e-3
    VSWR = 2.0
    sigma = 5.8e7                   # Copper

    # FR4_Epoxy
    Er = 4.4
    tand = 0.02

    print("\n\nCalculating for FR4 patch.")
    W, L, h, Er = DesignPatch(Er, h, freq)
    eff = CalculatePatchEff(Er, W, L, h, tand, sigma, freq, VSWR)
    CalcDirectivity(eff, PatchFunction, freq, W, L, h, Er)

    # Rogers RO4350
    print("\n\nCalculating for RO4350 patch.")
    Er = 3.48
    tand = 0.004
    W, L, h, Er = DesignPatch(Er, h, freq)
    eff = CalculatePatchEff(Er, W, L, h, tand, sigma, freq, VSWR)
    CalcDirectivity(eff, PatchFunction, freq, W, L, h, Er)
# end def


# ======================================== #

"""
Function to calculate peak directivity.
Also includes some examples that are used to check result.
"""


def SqrtSinPattern(Theta, Phi, *args):
    """
    See Fig1 @ http://www.antenna-theory.com/basics/directivity.php
    Expect Directivity to be 1.05dB.
    """
    return sqrt(sin(radians(Theta)))


def SinPowerPattern(Theta, Phi, *args):
    """
    See Fig1 @ http://www.antenna-theory.com/basics/directivity.php
    Expect Directivity to be 2.707dB.
    """
    return sin(radians(Theta)) ** 5


def IsotropicPattern(Theta, Phi, *args):
    """
    Isotropic directional pattern. i.e. radiation is same in all directions.
    Expect directivity to be 0dB.
    """
    return 1


def xfrange(start, stop, step):
    """
    Creates range of float values.
    """
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1


def CalcDirectivity(Efficiency, RadPatternFunction, *args):
    """
    Based on calc_directivity.m from ArrayCalc.
    Calculates peak directivity in dBi value using numerical integration.
    If the array efficiency is set to below 100% then the returned value is referred to as Gain (dB).

    Usage: ThetaMax, PhiMax = CalcDirectivity(RadPatternFunction, Efficiency)

    RadPatternFunction - antennas radiation pattern function. F(Theta, Phi)
    Efficiency - Efficiency of antenna in %. Default 100%.

    Returned values:
    ThetaMax - Theta value for direction of maximum directivity (Deg)
    PhiMax - Phi value for direction of maximum directivity (Deg)

    Integration is of the form :
    %
    %       360   180
    %     Int{  Int{  (E(theta,phi)*conj(E(theta,phi))*sin(theta) d(theta) d(phi)
    %        0     0
    %
    %         z
    %         |-theta   (theta 0-180 measured from z-axis)
    %         |/
    %         |_____ y
    %        /\
    %       /-phi       (phi 0-360 measured from x-axis)
    %      x
    %
    """
    print("Calculating Directivity for " + RadPatternFunction.__name__)

    deltheta = 2                                                                # Step value of theta (Deg)
    delphi = 2                                                                  # Step value for phi (Deg)

    dth = radians(deltheta)
    dph = radians(delphi)

    Psum = 0
    Pmax = 0
    Thmax = 0
    Phmax = 0

    for phi in xfrange(0, 360, delphi):                                                                     # Phi Integration Loop 0-360 degrees
        for theta in xfrange(0, 180, deltheta):                                                             # Theta Integration Loop 0-180 degrees
            eField = RadPatternFunction(theta, phi, *args)                                       # Total E-field at point
            Pthph = eField * np.conjugate(eField)                                                                             # Convert to power

            if Pthph > Pmax:
                Pmax = Pthph                                                                                # Store peak value
                Thmax = theta                                                                               # Store theta value for the maximum
                Phmax = phi                                                                                 # Store phi value for the maximum

            # print(str(theta) + "," + str(phi) + ": " + str(Pthph))
            Psum = Psum + Pthph * sin(radians(theta)) * dth * dph                                           # Summation

    Pmax = Pmax * (Efficiency / 100)                                                                        # Apply antenna efficiency

    directivity_lin = Pmax / (Psum / (4 * pi))                                                              # Directivity (linear ratio)
    directivity_dBi = 10 * log10(directivity_lin)                                                           # Directivity (dB wrt isotropic)

    if Efficiency < 100:                                                                                    # Gain case
        dBdiff = 10 * log10(abs(100 / Efficiency))                                                          # Difference between gain and directivity
        print("Directivity = " + str(directivity_dBi + dBdiff) + "dBi")                                     # Display what directivity would be for ref.
        print("Efficiency = " + str(Efficiency) + "%")
        print("Gain = " + str(directivity_dBi) + "dB")
    else:                                                                                                   # Directivity case
        print("Directivity = " + str(directivity_dBi) + "dBi")

    print("At Theta = " + str(Thmax) + ", Phi = " + str(Phmax))

    return Thmax, Phmax


def exampleDirectivity():
    CalcDirectivity(100, SqrtSinPattern)
    print("\n\n")
    CalcDirectivity(90, SinPowerPattern)
    print("\n\n")
    CalcDirectivity(100, IsotropicPattern)

    print("\n\n")

    freq = 14e9
    Er = 3.66                                                           # RO4350B

    h = 0.101e-3
    W, L, h, Er = DesignPatch(Er, h, freq)
    CalcDirectivity(100, PatchFunction, freq, W, L, h, Er)
    fields = PatchEHPlanePlot(freq, W, L, h, Er)
    SurfacePlot(fields, freq, W, L, h, Er)

    W = 10.7e-3
    L = 10.47e-3
    h = 3e-3
    Er = 2.5

    print("\n\n")
    CalcDirectivity(100, PatchFunction, freq, W, L, h, Er)
    fields = PatchEHPlanePlot(freq, W, L, h, Er)
    SurfacePlot(fields, freq, W, L, h, Er)
# end def


# ======================================== #


def PatchFunction(thetaInDeg, phiInDeg, Freq, W, L, h, Er):
    """
    Taken from Design_patchr
    Calculates total E-field pattern for patch as a function of theta and phi
    Patch is assumed to be resonating in the (TMx 010) mode.
    E-field is parallel to x-axis

    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate

    Refrence C.A. Balanis 2nd Edition Page 745
    """
    lamba = light_velocity / Freq

    theta_in = math.radians(thetaInDeg)
    phi_in = math.radians(phiInDeg)

    ko = 2 * math.pi / lamba

    xff, yff, zff = sph2cart1(999, theta_in, phi_in)                            # Rotate coords 90 deg about x-axis to match array_utils coord system with coord system used in the model.
    xffd = zff
    yffd = xff
    zffd = yff
    r, thp, php = cart2sph1(xffd, yffd, zffd)
    phi = php
    theta = thp

    if theta == 0:
        theta = 1e-9                                                              # Trap potential division by zero warning

    if phi == 0:
        phi = 1e-9

    Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * (1 + 12 * (h / W)) ** -0.5        # Calculate effictive dielectric constant for microstrip line of width W on dielectric material of constant Er

    F1 = (Ereff + 0.3) * (W / h + 0.264)                                        # Calculate increase length dL of patch length L due to fringing fields at each end, giving total effective length Leff = L + 2*dL
    F2 = (Ereff - 0.258) * (W / h + 0.8)
    dL = h * 0.412 * (F1 / F2)

    Leff = L + 2 * dL

    Weff = W                                                                    # Calculate effective width Weff for patch, uses standard Er value.
    heff = h * sqrt(Er)

    # Patch pattern function of theta and phi, note the theta and phi for the function are defined differently to theta_in and phi_in
    Numtr2 = sin(ko * heff * cos(phi) / 2)
    Demtr2 = (ko * heff * cos(phi)) / 2
    Fphi = (Numtr2 / Demtr2) * cos((ko * Leff / 2) * sin(phi))

    Numtr1 = sin((ko * heff / 2) * sin(theta))
    Demtr1 = ((ko * heff / 2) * sin(theta))
    Numtr1a = sin((ko * Weff / 2) * cos(theta))
    Demtr1a = ((ko * Weff / 2) * cos(theta))
    Ftheta = ((Numtr1 * Numtr1a) / (Demtr1 * Demtr1a)) * sin(theta)

    # Due to groundplane, function is only valid for theta values :   0 < theta < 90   for all phi
    # Modify pattern for theta values close to 90 to give smooth roll-off, standard model truncates H-plane at theta=90.
    # PatEdgeSF has value=1 except at theta close to 90 where it drops (proportional to 1/x^2) to 0

    rolloff_factor = 0.5                                                       # 1=sharp, 0=softer
    theta_in_deg = theta_in * 180 / math.pi                                          # theta_in in Deg
    F1 = 1 / (((rolloff_factor * (abs(theta_in_deg) - 90)) ** 2) + 0.001)       # intermediate calc
    PatEdgeSF = 1 / (F1 + 1)                                                    # Pattern scaling factor

    UNF = 1.0006                                                                # Unity normalisation factor for element pattern

    if theta_in <= math.pi / 2:
        Etot = Ftheta * Fphi * PatEdgeSF * UNF                                   # Total pattern by pattern multiplication
    else:
        Etot = 0

    return Etot

def sph2cart1(r, th, phi):
  x = r * cos(phi) * sin(th)
  y = r * sin(phi) * sin(th)
  z = r * cos(th)

  return x, y, z

def cart2sph1(x, y, z):
  r = sqrt(x**2 + y**2 + z**2) + 1e-15
  th = acos(z / r)
  phi = atan2(y, x)

  return r, th, phi


# ========================================================================= #
# ========================================================================= #

class Result:
    def __init__(self):
        self.frequency = None
        self.patch_width = None
        self.patch_length = None
        self.feeder_width = None
        self.feeder_length = None
        self.inset_gap_width = None
        self.inset_length = None
        self.ground_length = None
        self.ground_width = None
        self.input_edge_impedance = None


def design_string(resonant_frequency, dielectric_constant, thickness):
    return json.dumps(design_result(resonant_frequency, dielectric_constant, thickness).__dict__, indent=4)


def design_result(resonant_frequency, dielectric_constant, thickness):
    return design(resonant_frequency, dielectric_constant, thickness).get_result()


def design(resonant_frequency, dielectric_constant, thickness):
    """calculates length and width of patch antenna from dielectric constant, thickness and resonant frequency"""
    return PatchDesigner(resonant_frequency, dielectric_constant, thickness)


class PatchDesigner:
    """All parameter calculations"""
    freq = None
    er = None
    h = None
    patch_length = None
    patch_lengthl_eff = None
    patch_width = None
    feeder_length = None
    feeder_width = None
    inset_gap = None
    e_eff = None
    delta_l = None
    wavelength = None
    electrical_length = None
    ground_length = None
    ground_width = None
    inset_length = None
    input_impedance = None

    def __init__(self, freq, er, h):
        """
        Designs the patch parameters
        Parameters:
            freq (float): Resonant frequency in Hz.
            er (float): Dielectric constant of the cavity material.
            h (float): Thickness of the cavity in m.
        """
        if not 10 ** 6 <= freq <= 100 * 10 ** 9:
            raise ValueError("Frequency value should be in between 1MHz to 100 GHz")

        if not 0 < er <= 10**5:
            raise ValueError("Dielectric constant value should be in greater than 0 and smaller or equals 100,000")

        if not 0 < h <= 1:
            raise ValueError("Thickness value should be in greater than 0 and smaller or equals 1 meter")

        self.freq = freq
        self.er = er
        self.h = h
        self.set_wavelength()
        self.set_length_width_e_eff()
        self.set_feeder_width_length()

    def set_wavelength(self):
        self.wavelength = light_velocity / self.freq

    def set_length_width_e_eff(self):
        self.patch_width = (light_velocity / (2 * self.freq)) * sqrt(2 / (self.er + 1))
        temp = 1 + 12*(self.h / self.patch_width)
        self.e_eff = ((self.er + 1) / 2) + ((self.er - 1) / 2) * temp ** -0.5
        f1 = (self.e_eff + 0.3) * (self.patch_width / self.h + 0.264)
        f2 = (self.e_eff - 0.258) * (self.patch_width / self.h + 0.8)
        self.delta_l = self.h * 0.412 * (f1 / f2)
        self.patch_lengthl_eff = (self.wavelength / sqrt(self.e_eff)) / 2
        self.patch_length = self.patch_lengthl_eff - 2 * self.delta_l

    def set_feeder_width_length(self):
        self.feeder_length = (light_velocity / (4 * self.freq)) * (sqrt(1 / self.e_eff))
        self.feeder_width = self.patch_width / 5
        self.inset_gap = self.patch_width / 5
        self.set_input_impedance()
        self.inset_length = (self.patch_length / pi) * (math.acos(sqrt(impedance / self.input_impedance)))
        self.ground_length = self.patch_length + self.feeder_length + self.get_fringing_l()
        self.ground_width = self.patch_width + self.feeder_width + self.get_fringing_l()

    def get_result(self):
        result = Result()
        result.frequency = self.freq
        result.patch_width = self.patch_width
        result.patch_length = self.patch_length
        result.feeder_width = self.feeder_width
        result.feeder_length = self.feeder_length
        result.inset_gap_width = self.inset_gap
        result.inset_length = self.inset_length
        result.ground_length = self.ground_length
        result.ground_width = self.ground_width
        result.edge_impedance = self.input_impedance
        return result

    def get_fringing_l(self):
        return 6 * self.h

    def get_k(self):
        k0 = (2*pi)/self.wavelength
        return k0

    def S_i(self, a):
        temp = integrate.quad(lambda x: sin(x)/x, 0, a)
        return temp[0]

    def getG1 (self):
        k0 = self.get_k()
        X = k0 * self.patch_width
        I1 = -2 + cos(X) + X * self.S_i(X) + sin(X)/X
        G1 = I1 / (120 * pi**2)
        return G1

    def J0(self, s):
        temp = integrate.quad(lambda x: cos(s*sin(x)), 0, pi)
        return (1/pi) * temp[0]

    def getG12 (self):
        k0 = self.get_k()
        temp = integrate.quad(lambda x: (((sin(k0 * self.patch_width * cos(x) / 2) / cos(x)) ** 2) * self.J0(k0 * self.patch_length * sin(x)) * sin(x) ** 3), 0, pi)
        G12 = (1/(120*pi**2))*temp[0]
        return G12

    def set_input_impedance(self):
        G1, G12 = self.getG1(), self.getG12()
        self.input_impedance = 1 / (2 * (G1 + G12))


def m_to_inch(val):
    return 39.3701 * val


def get_gerber_str(d, feed_type):
    fl = m_to_inch(d.feeder_length)
    fw = m_to_inch(d.feeder_width)
    pl = m_to_inch(d.patch_length)
    pw = m_to_inch(d.patch_width)
    fringing_l = m_to_inch(d.get_fringing_l())
    gerber = get_inset_feed_gerber(fl, fw, pl, pw, fringing_l, d) if feed_type == 'inset' else \
        get_normal_feed_gerber(fl, fw, pl, pw, fringing_l)
    return gerber


def get_normal_feed_gerber(fl, fw, pl, pw, fringing_l):
    init_x = "{:.4f}".format((fl/2) + fringing_l).replace('.', '')
    init_y = "{:.4f}".format(fringing_l).replace('.', '')
    patch_x = "{:.4f}".format(fl + fringing_l + (pl/2)).replace('.', '')
    gerber_format = f"""
G04 ===== Begin FILE IDENTIFICATION =====*
G04 File Format:  Gerber RS274X*
G04 ===== End FILE IDENTIFICATION =====*
%FSLAX24Y24*%
%MOIN*%
%SFA1.0000B1.0000*%
%OFA0.0B0.0*%
%ADD14R,{fl}X{fw}*%
%ADD15R,{pl}X{pw}*%
%LNcond*%
%IPPOS*%
%LPD*%
G75*
D14*
X{init_x}Y{init_y}D03*
D15*
X{patch_x}*
M02*
    """
    return gerber_format


def get_inset_feed_gerber(fl, fw, pl, pw, fringing_l, d):
    inset_l = m_to_inch(d.inset_length)
    inset_g = m_to_inch(d.inset_gap)
    pl_s = pl - inset_l
    init_x = "{:.4f}".format((fl/2) + fringing_l).replace('.', '')
    init_y = "{:.4f}".format(fringing_l).replace('.', '')
    patch_x = "{:.4f}".format(fl + fringing_l + inset_l + (pl_s/2)).replace('.', '')
    inset_x = "{:.4f}".format(fl + fringing_l + (inset_l/2)).replace('.', '')
    inset_top_y = "{:.4f}".format(fw/2 + inset_g + (inset_g/2) + fringing_l).replace('.', '')
    inset_y = "{:.4f}".format(fringing_l).replace('.', '')
    inset_down_y = "{:.4f}".format(fringing_l - (fw/2 + inset_g + (inset_g/2))).replace('.', '')
    gerber_format = f"""
G04 ===== Begin FILE IDENTIFICATION =====*
G04 File Format:  Gerber RS274X*
G04 ===== End FILE IDENTIFICATION =====*
%FSLAX24Y24*%
%MOIN*%
%SFA1.0000B1.0000*%
%OFA0.0B0.0*%
%ADD14R,{fl}X{fw}*%
%ADD15R,{pl_s}X{pw}*%
%ADD16R,{inset_l}X{inset_g}*%
%LNcond*%
%IPPOS*%
%LPD*%
G75*
D14*
X{init_x}Y{init_y}D03*
D15*
X{patch_x}*
D16*
X{inset_x}Y{inset_top_y}*
D16*
X{inset_x}Y{inset_y}*
D16*
X{inset_x}Y{inset_down_y}*
M02*
    """
    return gerber_format


def write_gerber(resonant_frequency, dielectric_constant, thickness, file_name, feed_type):
    """Calculate design values in inch"""
    d = PatchDesigner(resonant_frequency, dielectric_constant, thickness)
    write_gerber_design(d, file_name, feed_type)


def write_gerber_design(design_: PatchDesigner, file_name, feed_type="normal"):
    content = get_gerber_str(design_, feed_type)
    with (open(file_name, 'w')) as f:
        f.write(content)


# ======================================= #


def plotPatch3D():
    #
    #
    #
    # Under Construction
    #
    #
    #
    import plotly
    from plotly.offline import init_notebook_mode
    import plotly.graph_objs as go
    plotly.offline.init_notebook_mode(connected=True)

    intensity = [0, 0.142857142857143, 0.285714285714286, 0.428571428571429, 0.571428571428571, 0.714285714285714, 0.85714257142857, 1]

    i = [7, 0, 0, 0, 4, 4, 2, 6, 4, 0, 3, 7]
    j = [3, 4, 1, 2, 5, 6, 5, 5, 0, 1, 2, 2]
    k = [0, 7, 2, 3, 6, 7, 1, 2, 5, 5, 7, 6]


    cavitycolor = [[0, 'rgb(0, 100, 0)']*6]

    coppercolor = [[0, 'rgb(139, 69, 19)'],
                   [1, 'rgb(139, 69, 19)'],
                   [2, 'rgb(139, 69, 19)'],
                   [3, 'rgb(139, 69, 19)'],
                   [4, 'rgb(139, 69, 19)'],
                   [5, 'rgb(139, 69, 19)'],
                   ]

    ct = 0.05 # copper_thickness
    pl = float(5) # patch length
    pw = float(5) # patch width
    fl = float(2) # feeder length
    fw = float(2) # feeder width

    # height
    ch = float(1)

    tl = pl + fl # total length

    tfp = (pw / 2) + (fw / 2)  # top feeder point
    bfp = (pw / 2) - (fw / 2)  # bottom feeder point

    data = [
        go.Mesh3d(
            x = [0, 0, tl, tl, 0, 0, tl, tl],
            y = [0, pw, pw, 0, 0, pw, pw, 0],
            z = [0, 0, 0, 0 ] + ([ct]*4),
            colorbar = go.ColorBar(
                title='ground'
            ),
            facecolor = coppercolor,
            intensity = intensity,
            i = i,
            j = j,
            k = k,
            name = 'ground',
            showscale = True
        ),
        go.Mesh3d(
            x = [0, 0, pl, pl, 0, 0, pl, pl],
            y = [0, pl, pl, 0, 0, pl, pl, 0],
            z = ([ch]*4) + ([ch + ct]*4),
            colorbar = go.ColorBar(
                title='patch_top'
            ),
            facecolor = coppercolor,
            intensity = intensity,
            i = i,
            j = j,
            k = k,
            name = 'patch_top',
            showscale = True
        ),
        go.Mesh3d(
            x = [pl, pl, tl, tl, pl, pl, tl, tl],
            y = [tfp, bfp]*4,
            z = ([ch]*4) + ([ch+ct]*4),
            colorbar = go.ColorBar(
                title='feeder_top'
            ),
            facecolor= coppercolor,
            i = i,
            j = j,
            k = k,
            name = 'feeder_top',
            showscale = True
        ),
        go.Mesh3d(
            x = [0, 0, tl, tl, 0, 0, tl, tl],
            y = [0, pw, pw, 0, 0, pw, pw, 0],
            z = ([0 + ct] * 4) + ([ch] * 4),
            colorbar = go.ColorBar(
                title='cavity'
            ),
            facecolor = cavitycolor,
            intensity = intensity,
            i = i,
            j = j,
            k = k,
            name='cavity',
            showscale=True
        ),
        go.Mesh3d()
    ]

    layout = go.Layout(
        xaxis=go.XAxis(
            title='x'
        ),
        yaxis=go.YAxis(
            title='y'
        )

    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig)
# end def

def exampleDesignPatch():
    freq = 2.4e9
    Er = 4.4
    h = 1.6 * 10 ** -3
    # v = 3 * 10 ** 8

    W, L = design_patch(Er, h, freq)

    Rin = input_impedance(freq, W, L)
    print('Inset Feed Position : ', inset_feed_position(Rin, L))

    G1, G12 = getGs(freq, W, L)
    print('G1 : ', G1)
    print('G12 : ', G12)

    I1 = 1.863
    I2 = 3.59801

    d1, d2 = get_directivity(G1, G12, W, freq, I1, I2)

    print('Directivity : ', d1, ' dB')
    print('Directivity : ', d2, ' dB')

    fields = patch_eh_plane_plot(freq, W, L, h, Er)
    # surface_plot(fields)
    SurfacePlot(fields, freq, W, L, h, Er)
# end def

def exampleSimpleDesign():
    import pytest

    # resonant frequency in Hz
    freq = 2.4 * 10 ** 9

    # dielectric constant
    er = 4.4

    # thickness of the cavity in meter
    h = 1.6 * 10 ** -3

    # Quick print result
    print(design_string(freq, er, h))


    # Using result object
    simpleresult = design_result(freq, er, h)
    print(simpleresult.patch_width)

    # Quick write gerber file
    #
    # normal feed
    write_gerber(freq, er, h, 'patch_normal_design.gbr', 'normal')

    # inset feed
    write_gerber(freq, er, h, 'patch_inset_design.gbr', 'inset')

    # Custom change design parameters
    # Using design object
    simpledesign = design(freq, er, h)

    # Changing feeder length and feeder width
    simpledesign.feeder_length *= 1.25
    simpledesign.feeder_width *= 1.10


    # To write as gerber file for both types

    # normal feed
    write_gerber_design(simpledesign, 'patch_normal_design.gbr', 'normal')

    # inset feed
    write_gerber_design(simpledesign, 'patch_inset_design.gbr', 'inset')
# end def

# ======================================= #

def test_frequency_limit():
    import pytest

    with pytest.raises(ValueError) as execinfo:
        design(0, 0, 0)

    return execinfo.value.args[0] == 'Frequency value should be in between 1MHz to 100 GHz'


def test_dielectric_limit():
    import pytest

    with pytest.raises(ValueError) as execinfo:
        design(10 ** 9, 0, 0)

    return execinfo.value.args[0] == 'Dielectric constant value should be in greater than 0 and smaller or equals 100,000'


def test_thickness_limit():
    import pytest

    with pytest.raises(ValueError) as execinfo:
        design(10 ** 9, 1, 0)

    return execinfo.value.args[0] == 'Thickness value should be in greater than 0 and smaller or equals 1 meter'

# ======================================= #

if __name__ == "__main__":
    exampleRectPatch()

    examplePatchEfficiency()

    exampleDirectivity()

    exampleDesignPatch()

    exampleSimpleDesign()

    test_frequency_limit()
    test_dielectric_limit()
    test_thickness_limit()
# end if



# ======================================== #