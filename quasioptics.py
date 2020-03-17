# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:44:33 2017

@author: gawe
"""

# ========================================================================== #
# ========================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
#import scipy as _scipy
import os as _os
import matplotlib.pyplot as _plt
import cmath

# Pybaseutils
from pybaseutils.Struct import Struct
from pybaseutils.utils import speed_of_light #, factorial

try:
    from numpy import matmul  # not available in numpy 1.9
except:
    from numpy import dot as matmul  # kluge for numpy 1.9
# end try


# ========================================================================== #


# Base Material - Free space
cc, mu0, eps0 = speed_of_light()

class qoptics_params(Struct):

    # ====================================================================== #
    # ====================================================================== #
    # auto-updating properties

    # qz is the main control property. Everything is set by qz:
    @property
    def qz(self):
        return self._qz
    @qz.setter
    def qz(self, value):
        """
        q = z + i*zr = z + i * pi*N*wo^2/lambda0
            or
        1/qz = q/Rz - i*lambda0/(pi*Nref*wz^2)
        """
        # value = _np.asarray(value)
        self._qz = value
        self._Rz = _np.real(1.0/self._qz)
        self._zr = _np.atleast_1d(_np.imag(self._qz))

        tmp = _np.imag(-1.0/self._qz)
        tmp = self.lambda0/self.Nrefr/(_np.pi*tmp)
        self._wz = _np.sqrt(tmp)
        self._wo = self._wz/_np.sqrt( 1.0+(self.zz/self._zr)**2.0 )
        self.fwhm = self._wz/_np.sqrt(2*_np.log(2))
    @qz.deleter
    def qz(self):
        del self._qz

    # =============================== #

    @property
    def wo(self):
        return self._wo
    @wo.setter
    def wo(self, value):
        # value = _np.asarray(value)
        self._wo = value
        # self.RayleighRange()
    @wo.deleter
    def wo(self):
        del self._wo

    @property
    def zr(self):
        return self._zr
    @zr.setter
    def zr(self, value):
        self._zr = value
#        self._wo = _np.sqrt(self.zr*self.lambda0/(_np.pi*self.Nrefr))
    @zr.deleter
    def zr(self):
        del self._zr

    @property
    def wz(self):
        return self._wz
    @wz.setter
    def wz(self, value):
        # value = _np.asarray(value)
        self._wz = value
    @wz.deleter
    def wz(self):
        del self._wz

    @property
    def Rz(self):
        return self._Rz
    @Rz.setter
    def Rz(self, value):
        self._Rz = value
    @Rz.deleter
    def Rz(self):
        del self._Rz

    # ====================================================================== #
    # ======================== Grid properties ============================= #


    # Grid properties:
    @property
    def zz(self):
        return self._zz
    @zz.setter
    def zz(self, value):
        # value = _np.asarray(value)
        self._zz = value
        self.nz = len(self._zz)
    @zz.deleter
    def zz(self):
        del self._zz

    @property
    def rr(self):
        return self._rr
    @rr.setter
    def rr(self, value):
        # value = _np.asarray(value)
        self._rr = value
        self.nr = len(self._rr)
    @rr.deleter
    def rr(self):
        del self._rr

    @property
    def Nrefr(self):
        return self._Nrefr
    @Nrefr.setter
    def Nrefr(self, value):
        # value = _np.asarray(value)
        self._Nrefr = value
#        self.WaveVector()  # update the wavevector in media
    @Nrefr.deleter
    def Nrefr(self):
        del self._Nrefr

    @property
    def kvec(self):
        return self._kvec
    @kvec.setter
    def kvec(self, value):
        # value = _np.asarray(value)
        self._kvec = value
#        self.RefractiveIndex() # update the refractive index if the wavevector is set
    @kvec.deleter
    def kvec(self):
        del self._kvec

    # ====================================================================== #
    # ===================== Quasi-optical formula ========================== #


    def WaveVector(self):
        if not hasattr(self,'k0'):
            self.k0 = 2.0*_np.pi/self.lambda0
        # endif
        self._kvec = self.Nrefr*self.k0

    def RefractiveIndex(self):
        self._Nrefr = _np.sqrt(cc*self.kvec/self.omega)

    def ComplexBeamParameter(self):
        self.qz = self.zz + 1j*self.zr

    def RayleighRange(self):
        self._zr = _np.pi*self.Nrefr*self.wo**2.0 / self.lambda0

    def BeamWaist(self):
        self._wo = _np.sqrt(2.0*self.zr/self.kvec)

    def RadiusCurvaturePhaseFront(self):
        self.Rz = self.zz*( 1.0 + (self.zr/self.zz)**2.0)

    def SpotSize(self):
        self._wz = self.wo * _np.sqrt( 1.0 + (self.zz/self.zr)**2.0)

    def GouyPhaseShift(self):
        self._gphz = _np.arctan(self.zz/self.zr)

    def PlaneWavePhaseShift(self):
        self._plphz = self.kvec*self.zz

    def WavefrontBendingPhaseShift(self):
        nr = len(self.rr)
        nz = len(self.zz)
        kvec = _np.copy(self.kvec)
        rr = _np.copy(self.rr)
        Rz = _np.copy(self.Rz)
        self._wbphz = _np.zeros((nz, nr), dtype=_np.float64)
        for ii in range(nr):
            # wavefront bending
            self._wbphz[:,ii] = kvec*rr[ii]**2.0/(2.0*Rz)
        # end for

    def power(self, rr=None, P0=1.0):
        """
        Stores the beam power as a function of beam radius and distance

            for circle of radius r, the fraction of power transmitted
                r=w(z)       -> is 0.865
                r=1.07*w(z)  -> is 0.90
                r=1.224*w(z) -> is 0.95
                r=1.52*w(z)  -> is 0.99
        """
        return P0*(1.0-_np.exp(-2.0*rr**2.0/self.wz**2.0))


    def intensity(self, P0=1.0):
        """
        Stores the beam intensity as a function of beam radius and distance
        """
        nz = len(self.zz)
        nr = len(self.rr)

        self.I0 = 2.0*P0/(_np.pi*self.wo**2.0)
#        self.Iaxis = 2.0*P0/(_np.pi*self.wz**2.0)

        self.Irz = _np.zeros((nz, nr), dtype=_np.float64)
        for ii in range(nr):
            self.Irz[:,ii] = self.I0*(self.wo/self.wz)**2.0 * _np.exp(-2.0*self.rr[ii]**2.0/self.wz**2.0)
        return self.Irz

    def GaussianPhase(self):
        self.GouyPhaseShift()
        self.PlaneWavePhaseShift()

        if hasattr(self, 'rr'):
            self.WavefrontBendingPhaseShift()

            nr = len(self.rr)
            nz = len(self.zz)
            self.phase = _np.zeros((nz, nr), dtype=_np.float64)

            for ii in range(nr):
                # add in the wavefront bending
                self.phase[:, ii] = self._plphz - self._gphz + self._wbphz[:, ii]
            # end for
        else:
            self.phase = self._plphz-self._gphz
        # end if

#    def ComplexAmplitude(self):
        # Uv(x,y,z) = Av*exp(-1j*Phase(x,y,z))

    # ======= #

    def RadiusCurvatureAtDistance(self, d1, wo=None):
        if wo is None:
            wo = self.wo
        self.WaveVector()
        return d1*( 1.0+(0.5*self._kvec*wo**2.0/d1)**2.0 )

    # ======= #

#    @staticmethod
#    def coth(val):
#        return 1.0/cmath.tanh(val)

    # ====================================================================== #
    # ======================= Common Materials ============================= #

    @staticmethod
    def _mats(material, metal=False):
        # ============== Vacuum filled guide ============== #
        if material.upper().find("VACUUM")>-1 or material.upper().find("NONE")>-1 or material is None:
            epsr = 1.0
            loss_tangent = 0.0

        # ============== Air filled guide ============== #
        elif material.upper().find("AIR")>-1:
            epsr = 1.0006  # relative permeability of material in the guide, air
            loss_tangent = 0.0  # loss tangent of material in guide

        # ============== Polystyrene filled guide ============== #
        elif material.lower().find("polystyrene")>-1 or material.lower().find('ps')>-1:
            epsr = 2.4  # relative permeability of material in the guide, air
            loss_tangent = 0.0  # loss tangent of material in guide

        ## ============== Polyamide filled guide ============== #
        elif material.lower().find("polyamide")>-1 or material.lower().find("kapton")>-1:  # is this actually for 'polyimide'?
            epsr = 4.3  # relative permeability of material in the guide
            loss_tangent = 0.004  # loss tangent of material in guide

        # ============== Mica filled guide ============== #
        elif material.lower().find("mica")>-1 or material.lower().find("glimmer")>-1:
            epsr = 5.7  # relative permeability of material in the guide
            loss_tangent = 0.000  # loss tangent of material in guide

        # ============== Teflon (PTFE) filled guide ============== #
        elif material.lower().find("ptfe")>-1 or material.lower().find("teflon")>-1:
            epsr = 2.1  # relative permeability of material in the guide
            loss_tangent = 0.001  # loss tangent of material in guide

        # ============== Sapphire filled guide ============== #
        elif material.lower().find("saphire")>-1:
            epsr = 10.0  # relative permeability of material in the guide
            loss_tangent = 0.000  # loss tangent of material in guide

        # ============== Fused Quartz filled guide ============== #
        elif material.lower().find("quartz")>-1:
            epsr = 3.78  # relative permeability of material in the guide
            loss_tangent = 0.000  # loss tangent of material in guide

        # ============== Alumina ceramic ============== #
        elif material.lower().find("alumina")>-1:
            epsr = 8.66  # relative permeability of material in the guide
            loss_tangent = 0.0018  # loss tangent of material in guide

        # ============== Macor ceramic ============== #
        elif material.lower().find("macor")>-1:
            epsr = 5.67  # relative permeability of material in the guide
            loss_tangent = 0.0071  # loss tangent of material in guide
        # endif

        if metal:
            if material.lower().find("copper")>-1 or material.lower().find('cu')>-1:
                mur = 0.999991 # Copper, relative permeability of cavity walls
                rho = 1.724e-8 # Copper, ohm-meter, resistivity of walls of cavity
            elif material.lower().find("aluminum")>-1 or material.lower().find('al')>-1:
                mur = 1.00002 # Aluminum, relative permeability of cavity walls
                rho = 2.65e-8 # Aluminum, ohm-meter, resistivity of walls of cavity
            elif material.lower().find("brass")>-1 or material.lower().find('bra')>-1:
                mur = 1.05 # Brass, relative permeability of cavity walls
                rho = 6.39e-8 # Brass, ohm-meter, resistivity of walls of cavity
            # end if
            return mur, rho
        # end if
        return epsr, loss_tangent
    # end def materials

    # ====================================================================== #
    # ====================================================================== #

# end class qoptics_base

class qoptics_abcd(qoptics_params):

    # ============= Built-in geometric optics ray matrices ================= #

    def propagate_ray(self, ABCD, r1, th1):
        """
        Propagate a ray at an angle through an ABCD matrix

        from a ray perspective (r is distance from beam-axis, theta is angle
        of approach)
        """
        arr = matmul(ABCD, _np.array([r1, th1], dtype=float) )
        return arr[0], arr[1]   # r2, th2

    def propagate_beamparameter(self, ABCD, qz):
        """
        use the ray tranfser method to propagate a Gaussian beam parameter

        Note:
            q(z) = z + i*zr
                zr = imag(q(z))
        and
            1/q(z) = 1/R(z) - i * lambda0 / (pi*N*w(z)^2)
                R(z) = 1.0/real( 1/q(z) )
                w(z) can be determined from the imaginary part
                (this is done automatically)

        """
        #self.qz = matmul(ABCD, _np.asarray([self.qz, 1.0], dtype=complex) )
        return (ABCD[0,0]*qz + ABCD[0,1])/(ABCD[1,0]*qz + ABCD[1,1])


    def ABCD(self, A, B, C, D):
        """
        Return the ray matrix from geometric optics given its parts
        """
        return _np.array( [[A, B],[C, D]], dtype=complex)

    def uniform_media(self, d):
        """
        Propagation through uniform media: ABCD matrix for a Gaussian beam
        Input
            d - [m], distance to propagate
        Output
            ABCD - [m], ABCD matrix for beam propagation through uniform media a distance d
        """
        return self.ABCD(1.0, d, 0.0, 1.0)

    def freespace(self, d):
        """
        Free space propagation ABCD matrix for a Gaussian beam
        Input
            d - [m], distance to propagate through free space
        Output
            ABCD - [m], ABCD matrix for beam propagation through free space a distance d
        """
        return self.uniform_media(d)

    def refraction_planarinterface(self, n1, n2):
        """
        Refraction at a flat interface between two media
            initial media index of refraction - n1
            final media index of refraction - n2
        """
        return self.ABCD(1.0, 0.0, 0.0, n1/n2)

    def refraction_curvedinterface(self, R, n1, n2):
        """
        Refraction at a curved interface between two media (spherical boundary)
            R - radius of curvature,
                R>0 for convex  (center of curvature after interface)
                R<0 for concave  (center of curvature in front of interface)
            initial media index of refraction - n1
            final media index of refraction - n2
        """
        return self.ABCD(1.0, 0.0, (n1-n2)/(R*n2), n1/n2)

    def thinlens(self, flen):
        """
        ABCD matrix for a Gaussian beam propagation through a thin lens
        Input
            flen - [m], focal length of the thin lens
        Output
            ABCD - [m], ABCD matrix for beam propagation through free space a distance d
        """
        return self.ABCD(1.0, 0.0, -1.0/flen, 1.0)

    def planarmirror(self):
        """
        Ray transfer matrice for reflection from a flat mirror

        Only valid for mirrors perpendicular to the ray.
            (identity matrix)
        """
        return self.ABCD(1.0, 0.0, 0.0, 1.0)

    def sphericalmirror(self, R, th, plane='horizontal'):
        """
        Ray transfer matrice for reflection from a curved (spherical) mirror
            (valid in the paraxial approximation)
            R is the radius of curvature of the spherical mirror

            theta in radians, the mirror angle of incidence in the horizontal plane.

            concave:  R>0
            convex:   R<0
        """
        if plane.lower == 'horizontal':
            # Effective radius of curvature in the tangential plane (horizontal direction)
            Re = R*_np.cos(th)
        else:
            # Effective radius of curvature in the sagittal plane (vertical direction)
            Re = R / _np.cos(th)
        # end if

        return self.ABCD(1.0, 0.0, -2.0/Re, 1.0)

    def thicklens(self, n1, n2, R1, R2, t):
        """
        ABCD matrix for a Gaussian beam propagation through a thick lens
        Input
            n1 = refractive index outside of the lens
            n2 = refractive index of the lens itself (inside the lens)
            R1 = Radius of curvature of First surface (interface)
            R2 = Radius of curvature of Second surface (interface)
            t = center thickness of lens
        Output
            ABCD - [m], ABCD matrix for beam propagation through free space a distance d

        Note that the implemented formulation is more general, but there is an
        analytic check that works:

          flat surface on left, curved surface, curvature (R), on right
                ____________
               |           )
               |            )
               |     n2      )     n1
               |            )
               |___________)
                    t

                    Mlens = |   1.0          t*n1/n2     |
                            |  -1/f      1.0-n1*t/(n2*f) |
                        f = (n2-1)*(1/R2-1/R1) + (n-1)**2.0 * t/(n2*R1*R2)
                          = thicklens_focallength(R1, R2, n2)
        """
        ABCD = matmul(self.refraction_curvedinterface(R2, n2, n1),
                                      self.uniform_media(t))
        return matmul( ABCD, self.refraction_curvedinterface(R1, n1, n2) )

    def refraction_thickplate(self, n1, n2, t):
        """
        ABCD matrix for a Gaussian propagating through a thick plate
            n1 = refractive index outside of the plate
            n2 = refractive index of the plate itself (inside the lens)
            t = center thickness of plate
        Output
            ABCD - [m], ABCD matrix for beam propagation through free space a distance d
        """
        ABCD = matmul(self.refraction_planarinterface(n2, n1), self.uniform_media(t))
        return matmul( ABCD, self.refraction_planarinterface(n1, n2))

    def rightangleprism(self, th, psi, d, n):
        """
        k=cos(psi)/cos(th) is the beam expansion factor
                th is angle of incidence
                psi is the angle of refraction
                d is the prism path length
                n is the refractive index in the prism material
            applies for orthogonal beam exit
        """
        k = _np.cos(psi)/_np.cos(th)
        return self.ABCD(k, d/(n*k), 0.0, 1.0/k)

    def multiple_prism(self, M, B):
        """
        M is the total beam magnification: M=k1*k2*k3...*kr
        B is the total optical proapgation distance of the multiple prism extender
        """
        return self.ABCD(M, B, 0.0, 1.0/M)

    def thinlens_freespace(self, flen, d):
        """
        propagate into a thin lens with focal length 'flen'
        then propagate through free space a distance d
        """
        # RTM for lens followed by free-space
        # A = 1-d/f
        # B = d
        # C = -1/f
        # D = 1
        return matmul(self.freespace(d), self.thinlens(flen))

    def freespace_thinlens(self, flen, d):
        """
        propagate through free space a distance d
        then propagate into a thin lens with focal length 'flen'
        """
        # RTM for free space followed by a lens
        # A = 1
        # B = d
        # C = -1/f
        # D = 1-d/f
        return matmul(self.thinlens(flen), self.freespace(d))

    def freespace_thinlens_freespace(self, d1, flen, d2):
        """
        propagate a distance 'd1' into a thin lens with focal length 'flen'
        then propagate through free space a distance d2
        """
        return matmul(self.freespace(d2), self.freespace_thinlens(flen, d1))
        # return matmul(self.thinlens_freespace(flen, d2), self.freespace(d1))

    # ============= Basic stuff ================= #

    def thinlens_focallength(self, r1, r2, n=1.0):
        """
        inputs
            r1 - radius of curvature of first surface
            r2 - radius of curvature of second surface
            n  - index of refraction in the lens  (default: 1)
        outputs
            flen - focal length of the thin lens
        """
        inverse_flen = 1.0/r2 + 1.0/r1  # flen = R1*R2/(R2+R1) = 1/(1/R1+1/R2)
#        inverse_flen *= n-1.0
        return 1.0/inverse_flen

    def thicklens_focallength(self, r1, r2, n, t):
        """
        inputs
            r1 - radius of curvature of first surface
            r2 - radius of curvature of second surface
            n  - index of refraction in the lens  (default: 1)
            t  - thickness of the lens
        outputs
            flen - focal length of the thick lens
        """
        inverse_flen = (n-1.0)*(1.0/r2 - 1.0/r1) + (n-1.0)**2.0*t/(n*r1*r2)
        return 1.0/inverse_flen

    def elliptic_mirror(self, r1, r2):
        """

        """
        foc = self.thinlens_focallength(r1, r2)
        return self.thinlens(foc)

    def rectangular_horn(self, freq, distance, R_A, A_h, B_h, awg):
        """
        inputs
            freq
            distance - [m], requested position
            R_A      - [m], Dist. from Horn aperture to WG apert.
            A_h      - [m], E-plane aperture dimension
            B_h      - [m], H-plane aperture dimension
            awg      - [m], waveguide aperture in the E-field direction
        outputs
            wo       - [m,m], output beamwaist at the requested position (E-plane / H-plane)
            ro       - [m,m], output radius of curvature at the requested position (E-plane / H-plane)

        """
        wavelength = cc/freq
        R_a = R_A*A_h/(A_h-awg)      #The length to the throat of the horn
        L_a = _np.sqrt( R_a^2 + 0.25*A_h^2.0 ) #The horn edge length in the E-field direction
        L_b = _np.sqrt( R_a^2 + 0.25*B_h^2.0 ) # the horn edge length in the H-field direction

#        alpha = _np.arctan(A_h/(2*R_a))     #The Horn flare angles
#        beta  = _np.arctan(B_h/(2*R_a))

        # Optimum coupling to the horn occurs for a beam waist of diameter:
        w_h = 0.35*(A_h/2)   # Beam radius at the horn aperture for best coupling
        w_e = 0.50*(B_h/2)   # pg 163, and 165 or Table 7.1 pg 169

        # The beam waist within the horn antenna is then (pg 167, eqn 7.3a):
        woh = w_h/_np.sqrt(1+( (_np.pi*w_h**2)/(wavelength*L_a))**2)
        woe = w_e/_np.sqrt(1+( (_np.pi*w_e**2)/(wavelength*L_b))**2)
        wo_horn = [woe, woh]

        # The Rayleigh range, or confocal parameter (pg 22, eqn 2.41)
        z_rayH = _np.pi*woh**2/wavelength;    qo_H = 1j*z_rayH #beam parameter at waist
        z_rayE = _np.pi*woe**2/wavelength;    qo_E = 1j*z_rayE
#        qo_horn = [qo_E,qo_H]

        # The waist of the beam occurs offset from the horn aperture by:
        zh = L_a/( 1 + ( L_a/z_rayH)**2 ) # loc. of beam radius w/in apert.
        ze = L_b/( 1 + ( L_b/z_rayE)**2 ) # (pg 167 eqn 7.30b)
        zz_horn = [ze, zh]

        qoh = qoptics_abcd();   qoe = qoptics_abcd()
        qoh.qz = qoh.propagate_beamparameter(qoh.freespace(distance+zh), qo_H)
        qoe.qz = qoe.propagate_beamparameter(qoe.freespace(distance+ze), qo_E)

        return [qoe.wz, qoh.wz], [qoe.Rz, qoh.Rz], wo_horn, zz_horn


    def corrugated_horn(self, freq, distance, DH, LH):
        """
        inputs
            freq
            distance - [m], requested position
            DH       - [m], diameter of the corrugated horn (aperture)
            LH       - [m], length of the corrugated horn
        outputs
            wo       - [m], output beamwaist at the requested position
            ro       - [m], output radius of curvature at the requested position
        """
        wi = 0.32*DH                            # beam radius at the horn aperture
        ri = _np.sqrt(LH**2.0 + 0.25*DH**2.0)   # radius of curvature at the horn aperture

        wavelength = cc/freq
        wo = wi*_np.sqrt((1.0+distance/ri)**2.0 (wavelength*distance/(_np.pi*wi**2.0))**2.0)
        ro = distance/( 1.0 - (1.0+distance/ri)/( wo/wi )**2.0 )
        return wo, ro

    def two_thin_lens(self, f1, d, f2):
        return self.ABCD(-f2/f1, f1+f2, 0.0, -f1/f2)

    @staticmethod
    def elliptic_mirror_design(freq, win, wout, d1, d2, grazing_angle=45, **kwargs):
        """

        - Radii of curvature match the beams at input and output

        The equation of an ellipse in polar form from one of the focii to any
        point on the ellipse is
            r(theta) = amaj*(1-ecc**2.0)/( 1.0 +- ecc*cos(theta) )

        Note - ellipsoidal mirrors can be used to prevent/correct astigmastism
            to do that you have to know the beam radii in each plane

            example//

            # E-plane first
            f_e, (a_e, b_e), elength = elliptic_mirror_design(freq, wein, wout, d1, d2, grazing_angle, trunclevel)

            # Now calculate H-plane parameters to have identical output beam radii at d2
            f_h, (a_h, b_h), hlength = elliptic_mirror_design(freq, whin, wout, d1, d2, grazing_angle, trunclevel)


        """
        plotit = kwargs.setdefault('plotit', True)
        trunclevel = kwargs.setdefault('trunclevel', 5.0)
        plotit = kwargs.setdefault('plotit', True)

        inc_angle = _np.pi-grazing_angle

        wavelength = cc/freq
        R1 = d1*( 1.0+(_np.pi*win**2.0 /(wavelength*d1))**2.0 )
        R2 = d2*( 1.0+(_np.pi*wout**2.0/(wavelength*d2))**2.0 )

        # required dimension of the mirror to meet the truncation level
        qo = qoptics(freq=freq, wo=win, zz=_np.atleast_1d(d1))
        qo.RayleighRange()
        qo.ComplexBeamParameter()
        wproj = trunclevel*qo.wz/_np.sin(0.5*_np.pi-inc_angle)

        foc = qo.thinlens_focallength(R1, R2)
        M = R2/R1    # optical magnification

        #  amaj = the arithmetic mean of the distance between a point on the ellipse and both focii
        amaj = 0.5*(R1 + R2)

        # use eccentricity to calculate the ellipse parameters
        ecc = _np.sqrt(R1**2.0+R2**2.0-2.0*R1*R2*_np.cos(_np.pi-2.0*grazing_angle))/(2.0*amaj)
        bmin1 = amaj*_np.sqrt(1.0-ecc**2.0)

#        # this is analytically identical to the combination of the above:
#        bmin2 = _np.sqrt(0.5*R1*R2*(1.0+_np.cos(_np.pi-2.0*grazing_angle)))
#
#        print("bmin1 %6.4f, bmin2 %6.4f"%(bmin1,bmin2))

        print('foc %6.4f, R1 %6.4f, R2 %6.4f, d1 %6.4f, d2 %6.4f, amaj %6.4f, bmin %6.4f, ecc %6.4f'%(foc, R1, R2, d1, d2, amaj, bmin1, ecc))
        if plotit:
            qo.elliptic_mirror_design_plot(R1, R2, amaj, bmin1, ecc, mirlength=2*wproj,**kwargs)
        # end if
        return foc, M, (amaj, bmin1), wproj

    @staticmethod
    def rotated_ellipse(amaj, ecc, theta, phi):
#        return amaj*(1.0-ecc**2.0)/(1.0+ecc*_np.cos(theta - phi))
        return amaj*(1.0-ecc**2.0)/(1.0-ecc*_np.cos(theta - phi))

    @staticmethod
    def mirangle_on_ellipse(amaj, ecc, R1, phi):
#        return _np.arccos((amaj*(1.0-ecc**2.0)/R1-1)/ecc) + phi
        return _np.arccos((1.0-amaj*(1.0-ecc**2.0)/R1)/ecc) + phi

    @staticmethod
    def elliptic_mirror_design_plot(R1, R2, amaj, bmin, ecc, **kwargs):

        _ax1 = kwargs.setdefault('_ax1', None)
        phi = kwargs.setdefault('phi', 0.0)
        mirlength = kwargs.setdefault('mirlength', 0)

        import matplotlib.pyplot as _plt
        if _ax1 is None:
            hfig1, _ax1 = _plt.subplots(1, 1, sharex=True, squeeze=True)
        # end if

        qo = qoptics_abcd()

        # ellipse properties:
        angle = _np.linspace(0, 2*_np.pi, 180)
        c = ecc*amaj
        theta = qo.mirangle_on_ellipse(amaj, ecc, R1, phi)

        # Mirror position on ellipse
        rmir = qo.rotated_ellipse(amaj, ecc, theta, phi)
#        rmir = qo.rotated_ellipse(bmin, ecc, theta, phi)

        pmir = _np.vstack((rmir*_np.cos(theta), rmir*_np.sin(theta))).T
        pmir = _np.atleast_2d(pmir)

        # ================ #

        # rotated ellipse with zero at origin focus:
        relli = qo.rotated_ellipse(amaj, ecc, angle, phi)
#        relli = qo.rotated_ellipse(bmin, ecc, angle, phi)
        elli = _np.vstack((relli*_np.cos(angle), relli*_np.sin(angle))).T

        # ================ #

        # origin and target focus:
        _ax1.plot([0, 2*c*_np.cos(phi)], [0, 2*c*_np.sin(phi)], 'ko')
        _ax1.plot([0, 2*c*_np.cos(phi)], [0, 2*c*_np.sin(phi)], 'k--')

        # rotated ellipse:
        _ax1.plot(elli[:,0], elli[:,1], 'k--')  # ellipse

        # plot mirror position on ellipse
        _ax1.plot(pmir[:, 0], pmir[:, 1], 'ko')

        _ax1.plot([0, R1*_np.cos(theta)], [0, R1*_np.sin(theta)], 'r-')
        _ax1.plot([rmir*_np.cos(theta), 2*c*_np.cos(phi)], [rmir*_np.sin(theta), 2*c*_np.sin(phi)], 'r-')

        _ax1.set_aspect('equal', 'box')
#        _min = 1.1*_np.min(_np.min(elli))
#        _max = 1.1*_np.max(_np.max(elli))
#        _ax1.axhline(xmin=_min, xmax=_max, y=0.0, linestyle='-', color='k')
#        _ax1.axvline(ymin=_min, ymax=_max, x=0.0, linestyle='-', color='k')
#        _ax1.set_xlim((_min, _max))
#        _ax1.set_ylim((_min, _max))
        # end if


    # ====================================================================== #
    # ====================================================================== #
# end class


class qoptics(qoptics_abcd):
    def __init__(self, **kwargs):
        self.init(**kwargs)
    # end def __init__

    def init(self, **kwargs):
        # Calculate the free-space wavelength
        self.verbose = kwargs.get('verbose', True)
        self.freq = kwargs.get('freq', None)
        self.omega = 2.0*_np.pi*self.freq
        self.lambda0 = cc/self.freq

        # ====== #

        self.Nrefr = kwargs.get('Nrefr', None) # index of refraction in material
        self.kvec = kwargs.get('kvec', None)   # vector of wavenumbers in material

        self.setMaterial()

        # ====== #

        # beam waist
        if 'wo' in kwargs:
            self.wo = kwargs.get('wo')
            self.RayleighRange()
        # end if


        #  Grid parameters (assume on beam-axis at beam-waist if no input)
        if 'rr' in kwargs:
            self.rr = kwargs.get('rr')
        # end if

        if 'zz' in kwargs:
            self.zz = kwargs.get('zz')
            if 'wo' in kwargs:
                self.ComplexBeamParameter()
            # end if
        # end if
    # end def

    def setMaterial(self, Nrefr=None, kvec=None):
        # Material specification:
        # You  must input either a refractive index, or a wave-vector to get a material
        if Nrefr is not None:  self.Nrefr = Nrefr # endif
        if kvec is not None:   self.kvec = kvec # endif

        if self.Nrefr is None and self.kvec is None:
            if self.verbose:
                print('assuming free-space propagation')
            # endif
            self.Nrefr = 1.0
        #end if
    # end def setMaterial

    def add_element(self, ABCD, zinsert, lbl=None):
        if not hasattr(self, '_ABCD'):
            self._ABCD = []
            self._zins = []
            self._lbls = []
            self.nelements = 0
        # end if
        self.nelements += 1
        self._ABCD.append(ABCD)
        self._zins.append(zinsert)
        self._lbls.append(lbl)
    # end def

#    def add_element_params(self, **kwargs):
#
#        if 'incangle' in kwargs and not hasattr(self, '_incangle'):
#            self._incangle = []
#        # end if
#
#        if 'incangle' in kwargs:
#            self._incangle.append(incangle)
#        # end if
#    # end def

    def reset_elements(self):
        delattr(self, '_ABCD')
        delattr(self, '_zins')
        delattr(self, '_lbl')
        self.nelements = 0

        if hasattr(self, '_incangle'): delattr(self, '_incangle') # end if
    # end def

    def BeamPropagation(self):
        # ==== #
        if not hasattr(self, '_ABCD'):
            self.nelements = 0
            self._zins = []
        # end if

        zz = self.zz.copy()
        nz = len(zz)
        # append the end of the line to make the looping easier
        self._zins.append(zz[-1]-zz[0])

        # ==== #
        ABCDz = []
        qz = _np.zeros((nz,), dtype=_np.complex128)
        if hasattr(self, 'rays'):
            nrays = _np.atleast_2d(self.rays).shape[0]
            rayz = _np.zeros((nrays, 2, nz), dtype=_np.complex128)
        # end if

        for iz, z in enumerate(zz):
            # first build the ABCD matrix for this position

            dz = (z - zz[0])   # [m], distance along beam center

            if dz<self._zins[0]:
                # if the position is before the first element
                ABCD = self.freespace(dz)   # free-space up to this position
            else:
                # if after the first element
                ABCD = self.freespace(self._zins[0]) # free-space up to first element

                # loop over each element and test current position
                for nn in range(self.nelements):

                    # if the current position is greater than this element
                    if dz>self._zins[nn]:
                        # add in the element
                        ABCD = matmul(self._ABCD[nn], ABCD) # compound the matrices

                        # check current position against next element position
                        if dz<self._zins[nn+1]:
                            # if between this element and the next one
                            ABCD = matmul(self.freespace(dz-self._zins[nn]), ABCD) # add free-space between them
                        else:
                            # if position is greater than next element
                            ABCD = matmul(self.freespace(self._zins[nn+1]-self._zins[nn]), ABCD) # add free-space between them
                        # end if
                    # end if
                # end for
            # end if
            ABCDz.append(ABCD)
            qz[iz] = self.propagate_beamparameter(ABCD, self.qz[0])
            if hasattr(self, 'rays'):
                rays = _np.atleast_2d(self.rays)
                for ir, ray in enumerate(rays):
                    rayz[ir, 0, iz], rayz[ir, 1, iz] = self.propagate_ray(ABCD, ray[0], ray[1])
                # end for
            # end if
        # end for
        self.qz = qz
        self.ABCDz = ABCDz
        if len(self._zins)>0:
            self._zins = self._zins[:-1]
        else:
            delattr(self, '_zins')
        # end if

        if hasattr(self, 'rays'):
            self.rayz = _np.real(rayz)
        # end if
    # end def

    def plot(self, _ax1=None, _ax2=None):
        import matplotlib.pyplot as _plt
        if _ax1 is None and _ax2 is None:
            hfig1, (_ax1, _ax2) = _plt.subplots(2,1, sharex=True, squeeze=True)
        else:
            if _ax1 is None:
                hfig1, _ax1 = _plt.subplots(1, 1, sharex=True, squeeze=True)
            if _ax2 is None:
                hfig2, _ax2 = _plt.subplots(1, 1, sharex=True, squeeze=True)
        # end if

        _ax1.plot(self.zz, self.wz, 'k-')
        # _ax1.set_xlabel('Distance along beam-axis: z [m]')
        _ax1.set_ylabel('Beam Radius: w(z) [m]')
        _ax1.set_title('Beam propagation:')
        # _ax1.axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
        # _ax1.axvline(x=zantenna+0.40, linewidth=1.0, color='k', linestyle='--')
#        xlims = _ax1.get_xlim()
        ylims = _ax1.get_ylim()
        _ax1.set_ylim((0.0, 1.3*_np.max(self.wz)))
        ylims = _ax1.get_ylim()

        # xmax = zz[-1]
        # self.dvrg = (180.0/_np.pi)*self.lambda0/(_np.pi*self.Nrefr*self.wo)
        # ymax = _np.tan(self.dvrg*_np.pi/180.0)*xmax

        _ax2.plot(self.zz, self.Rz, 'k-')
        _ax2.set_xlabel('distance along beam-axis')
        _ax2.set_ylabel('Phase Curvature: R(z) [m]')

        if hasattr(self, '_ABCD'):
            for ii in range(self.nelements):
                _ax1.axvline(x=self._zins[ii]+self.zz[0], linewidth=1.0, color='k', linestyle='--')
                _ax2.axvline(x=self._zins[ii]+self.zz[0], linewidth=1.0, color='k', linestyle='--')

                if self._lbls[ii] is not None:
                    _ax1.text(x=(self._zins[ii]+self.zz[0])+0.01, y=0.8*ylims[1], s=self._lbls[ii])
                # end if
            # end for
        # end if
        return _ax1, _ax2
    # end def plot

    def plotrays(self, _ax1=None):
        import matplotlib.pyplot as _plt
#        if _ax1 is None and _ax2 is None:
#            hfig1, (_ax1, _ax2) = _plt.subplots(2,1, sharex=True, squeeze=True)
#        else:
        if 1:
            if _ax1 is None:
                hfig1, _ax1 = _plt.subplots(1, 1, sharex=True, squeeze=True)
#            if _ax2 is None:
#                hfig2, _ax2 = _plt.subplots(1, 1, sharex=True, squeeze=True)
        # end if

        for ray in self.rayz:
            _ax1.plot(self.zz, ray[0, :], 'k-')
        # _ax1.set_xlabel('Distance along beam-axis: z [m]')

        _ax1.set_ylabel('Ray Radius: r(z) [m]')
        _ax1.set_title('Ray propagation:')

#        _ax2.set_ylabel('Ray angle: th(z) [rad]')
#        _ax2.set_title('Ray propagation:')

        ylims = _ax1.get_ylim()

        if hasattr(self, '_ABCD'):
            for ii in range(self.nelements):
                _ax1.axvline(x=self._zins[ii]+self.zz[0], linewidth=1.0, color='k', linestyle='--')

                if self._lbls[ii] is not None:
                    _ax1.text(x=(self._zins[ii]+self.zz[0])+0.01, y=0.8*ylims[1], s=self._lbls[ii])
                # end if
            # end for
        # end if
        return _ax1
    # end def plot

    # ====================================================================== #


    def Efield_TEM00(self):
        """"
        returns the wave electric field from the fundamental Gaussian beam (TEM00)
        """
        rr = self.rr.copy()
        zz = self.zz.copy()
        self.nr = len(rr)
        self.nz = len(zz)
#        rv, zv = _np.meshgrid(rr, zz)
#
#        rr = _np.atleast_2d(rr)
#        zz = _np.atleast_2d(zz)
#        if _np.size(rr, axis=1) == 1:   rr = rr.T  # end if
#        if _np.size(zz, axis=0) == 1:   zz = zz.T  # end if
#        nr = self.nr
#        nz = self.nz

        # Update the beam radius, curvature of phase fronts and Guoy phase shift
        self.RaylieghRange()
        self.SpotSize()
        self.RadiusCurvaturePhaseFront()
        self.GaussianPhase()

        # === #

        # although the matrix way might be smarter, it is harder to debug
        self.Erz = _np.zeros((self.nz, self.nr), dtype=_np.complex128)
        for ii in range(self.nr):
            self.Erz[:,ii] = (self.wo/self.wz)*_np.exp(-rr[ii]**2.0/self.wz**2.0)*_np.exp(-1j*self.phase[:,ii])
        # end for
        return self.Erz

#    def elliptic_mirror_design(self):
        """
        a = semi-major radius of ellipse
        b = semi-minor radius of ellipse
        ellipsoid:  (x**2.0 + y**2.0)/b**2.0 + z**2.0 / a**2.0 = 1.0

        eccentricity:        e = sqrt(1.0-b**2.0 / a**2.0)

        Distance between any point on ellipse and focii:  P is point on ellipse
            (1) r1 + r2 = a (r1 is distance to focus 1 from P, r2 is distance
                             to focus 2 from P)

        Focal distance:
            (2) c = 2*e*a  Linear distance between the two focii, F1 / F2

        Waist of input / output beam (w1 / w2) are not at F1 / F2, but at
            distances d1 / d2 from a point on the ellipse, P

        Ideal design has radius of curvature of phase fronts equal to the
            distance from a point on ellipsoidal surface to each focii
            (3) r1 = R(d1) / r2 = R(d2)

        Treat ellipsoidal mirror as a thin lens:
            (4) 1/f = 1/r1 + 1/r2

        Ellipse design requires choosing w01/w02, then setting radius of
            curvature of the phase fronts on the elllipsoidal surface (3)

            That determines the focal lengths of the equivalent thin lens
            (ellisoidal mirror), as well as the semi-major axis of the
            ellipse. (1) + (4)

            The eccentricity (and semi-minor axis of the ellipse) is set by the
            required distance between the two Focii. (2)
        """

    # ====================================================================== #
    # ===================== Higher order modes ============================= #

    """
    A thin wrapper around numpy's hermite polynomial module
    """
    def HermDecomp(self, nmodes=3, mmodes=3):
        """
        Decompose the gaussian beam into the Gauss-Hermite basis set in the
        cartesian coordinate system
        """
#        self.nmodes = nmodes
#        self.mmodes = mmodes
#        self.nx = len(self.xx)
#        self.ny = len(self.yy)
#        self.nz = len(self.zz)
#        psi = _np.zeros((self.nmodes, self.mmodes, self.nz, self.nx, self.ny), dtype=_np.complex128)
#        for ii, zz in enumerate(self.zz):
#            for nn in range(self.nmodes):
#                for mm in range(self.mmodes):
#                    psi[nn,mm,ii,:,:] = _np.sqrt(2**(0.5-nn)/(_np.sqrt(_np.pi)*self.wz[ii]*factorial(nn)))
#                    psi[nn,mm,ii,:,:] *=H
        raise NotImplementedError

    herm = _np.polynomial.hermite
    def HermPoly(coef, domain=None, window=None):
        raise NotImplementedError
        # pass

#        _np.polynomial.hermite.Hermite
#        These are the methods available in that class:
#            __call__(arg)
#            basis(deg[, domain, window]) 	Series basis polynomial of degree deg.
#            cast(series[, domain, window]) 	Convert series to series of this class.
#            convert([domain, kind, window]) 	Convert series to a different kind and/or domain and/or window.
#            copy() 	Return a copy.
#            cutdeg(deg) 	Truncate series to the given degree.
#            degree() 	The degree of the series.
#            deriv([m]) 	Differentiate.
#            fit(x, y, deg[, domain, rcond, full, w, window]) 	Least squares fit to data.
#            fromroots(roots[, domain, window]) 	Return series instance that has the specified roots.
#            has_samecoef(other) 	Check if coefficients match.
#            has_samedomain(other) 	Check if domains match.
#            has_sametype(other) 	Check if types match.
#            has_samewindow(other) 	Check if windows match.
#            identity([domain, window]) 	Identity function.
#            integ([m, k, lbnd]) 	Integrate.
#            linspace([n, domain]) 	Return x, y values at equally spaced points in domain.
#            mapparms() 	Return the mapping parameters.
#            roots() 	Return the roots of the series polynomial.
#            trim([tol]) 	Remove trailing coefficients
#            truncate(size) 	Truncate series to length size.
#
    lagu = _np.polynomial.laguerre
    def LauguerrePoly(coef, domain=None, window=None):
        raise NotImplementedError

class qoptics2d(qoptics):

    def add_element(self, coords, style='reflector', **kwargs):
        """
        add each element by spatial coordinate
        inputs:
            coords
            style  - [s], 'reflector', 'splitter'

        """
        pass
    # end def

    def determine_paths(self):
        """
        calculate the unique paths
        """
        self.npaths = 2*self.ncombiner+self.nsplitter+1

#        for coord in coords:

        # check type of each element along the ray
        #   - track phase along each path independently (phase purposes)
        #   - if reflector, continue on this path
        #   - if splitter: copy path up to this point, then finish current path
        #       (for loop after copying)
        #   - check if any paths overlap
        #        where paths overlap combine beams
        #
        # combining beams:
        #      for two beams, there is an approximation
        #           use a weighted sum of their beam parametrs
        #           - 1/qe =  1.0/(|a1|+|a2|) (|a1|/q1 + |a2|/q2)
        #              M. A. Arain and G. Mueller, On the Interference of
        #              two Gaussian beams and their ABCD Matrix Representation
        #              optics express 19183 vol 17, No. 21, 2009
        #      in general, you have to decompose into a decomposition of higher
        #      order modes and combine there
        #
        # may be necessary:
        # -  define coordinates of each point on each element
        #  1) physical optics: calculate currents generated on each element by incident Gaussian beam
        #  2) ray trace to next optical element
        #  3) calculate scattered fields
        #  4) represent as Gaussian beam / higher order modes
        #      (interaction integral of of each mode with reflector current)
        # -  repeat 1-4 for each element
        pass
# end class


class dichroic(qoptics_abcd):
    """
    Dichroic plate as a quasi-optical high-pass filter / mirror

        A dichroic plate is essentially a plate with a dense hole pattern.
          - Each hole is basically a fundamental circular waveguide
              (transmiting higher-frequencies through the plate)

          - The transmission efficiency of the plate is basically the ratio of
            holes/solid that the incident wave samples

          - The transmitted waves generate an interference pattern due to the
            dense spacing of identical sources. Each having an Airy-function
            like antenna-pattern.

          - The transmitted wave-mode is dependent on the incident wave and
            the hole pattern itself.
           (An advanced design could be used as a quasi-optical mode-converter)

          - wavelengths on the same order as the hole spacing are reflected,
            but are also diffracted if the plate is not purely perpendicular

          - The dense pattern of holes appears to be a smooth surface to
            wavelengths that are much greater than the hole spacing, and these
            frequencies are reflected.

    """

    # ====================================================================== #
    # ======================= Dichroic filters ============================= #

    def example_fco135p1GHz(self, wd=None):
        """
        Example of a 135.1 GHz cut-off highpass filter
            Material 1   - Free space
            Material 2   - Material in contact with the metal plate
            Material 3   - material filling the metal plate (guide holes)
        """
        if wd is None:
            wd = _os.path.join('G://','Workshop','ECE','QMF','OP2','Dichroic Plate')
        # end if
        th = 45    # [deg], angle of incidence to dichroic plate (measured in free space, between ray and normal to plate surface)
        l3 = 5e-3  # [m], plate thickness
        D = 1.27e-3 # [m], diameter of guide holes, 135.1 GHz
        S = 1.45e-3 #[m] spacing of guide holes

        # frequency range to plot response
        freq = 1e9*_np.linspace(100.0, 250.0, 250)

        # Material 1 - free space
        eps1, loss_tangent1 = self._mats('Vacuum', metal=False)

        # Material 2 - Material in contact with the metal plate
        # relative permeability of material in waveguide
        eps2, loss_tangent = self._mats('Air', metal=False)

        # Material 3 - Material filling the metal plate (guide holes)
        # relative permeability / resistivity of cavity walls
        Ltot = l3   # [m], total physical thickness of plate
        mur, rho = self._mats('Brass', metal=True)

#        outputs = self.plate(freq=freq, diameter=D, spacing=S, thickness=Ltot, th=th, shape='hex')
        self.response(freq=freq, diameter=D, spacing=S, thickness=Ltot, th=th, plotit=True, wd=None)

    @staticmethod
    def max_hole_spacing(freq, epsr_plate, A, th):
        """
        Hole spacing must be small enough to ensure first grating lobe response
        lies at an angle of 90 degrees to the direction of propagation.
        """
        wavelength = cc/freq

        #    For an air-dielectric:
        return 1.0*wavelength/A/(_np.sqrt(epsr_plate) + _np.sin(th*_np.pi/180.0))
#        return 1.1547*wavelength/(1+_np.sin(th))   # reference 2

    @staticmethod
    def guide_cutoff_lower(D):
        """
        Designed cylindrical waveguide cutoff frequency (lower)

        inputs
            D - [m], diameter of waveguide
        output
            fco - [GHz], designed lower cut-off frequency
        """
        return 1e-9*1.841*cc/(_np.pi*D)

    @staticmethod
    def guide_cutoff_upper(S, A=1):
        """
        Designed cutoff frequency due to diffraction from the hole pattern (upper)

        inputs
            S - Hole spacing
            A - Effective Guide array (A=1, default, for square guide array,
                                   A=0.5*sqrt(3) for hexagonal guide array)
        output
            fco - [GHz], diffraction limited upper cut-off frequency
        """
        return 1e-9*cc/(S*A)

    @staticmethod
    def guide_wavelength(freq, epsr, D):
        """
        Guide wavelength in cylindrical waveguide
        """
        wavelength = cc/freq
        wl = wavelength/_np.sqrt(epsr)

        return wl/(1.0-(wl/(1.706*D))**2.0)

    @staticmethod
    def guide_electrical_length(freq, epsr, D, thickness):
        """
        Electrical length of each cavity
        """
        guide_wl = dichroic.guide_wavelength(freq, epsr, D)

        return 2.0*_np.pi*thickness/guide_wl  # [np/m]

    @staticmethod
    def attenuation_dielectric(freq, epsr, loss_tangent, D):
        """
        Attenuation due to the dielectric in the waveguide
        """
        wavelength = cc/freq
        guide_wl = dichroic.guide_wavelength(freq, epsr, D)

        return _np.pi*(guide_wl/wavelength)*loss_tangent/wavelength  # [np/m]

    @staticmethod
    def attenuation_wgwall(freq, epsd, rho, mur, D):
        """
        Attenuation constant due to dissipation in conducting cavity walls
        """
        wavelength = cc/freq
        guide_wl = dichroic.guide_wavelength(freq, epsd, D)

        rhoe = 1.724e-8 # ohm-meter = resistivity of copper
        alphc = 1.5e-4 * _np.sqrt(mur*rho/rhoe)*_np.sqrt(epsd/wavelength) * (guide_wl/(D*wavelength))
        alphc *= 0.420 + (wavelength/(1.706*D))**2.0
        return alphc  # [np/m]

    @staticmethod
    def square_guide_array():
        return 1.0 # square guide array

    @staticmethod
    def hexagonal_guide_array():
        return 0.5 * _np.sqrt(3.0)  # hexagonal guide array

    def plate(self, freq, diameter, spacing, thickness, th, shape='hex'):
        """
        A dichroic plate function
        inputs:
            diameter  - [m], diameter of guide holes (sets cut-off frequency)
            spacing   - [m], spacing of guide holes (sets porosity)
            thickness - [m], plate thickness (sets tunneling efficiency and should be multiple of guide wavelength)
            th        - [deg], angle of incidence to the dichroic plate
            (measured in free space, between ray and normal to plate surface)
            shape     - [str],  "square" or "hexagonal" sets array effeiciency
        outputs -

        """
        from scipy.special import jvp
        wavelength = cc/freq
        radius = 0.5*diameter

        #       For radius>0.28*spacing   and spacing<0.57 * wavelength
        if shape.lower().find('hex')>-1:
            A = dichroic.hexagonal_guide_array()
        elif shape.lower().find('sq')>-1:
            A = dichroic.square_guide_array()
        # end if

        # Cut-off frequencies
        fc1 = self.guide_cutoff_lower(diameter)       # [GHz], designed lower cut-off frequency
        fc2 = self.guide_cutoff_upper(spacing, A=A)   # [GHz], diffraction limited upper cut-off frequency

        # wlco = cc/(1e9*fc1)/_np.sqrt(eps3)

        J1prime = jvp(v=1, z=4.0*_np.pi*radius/(_np.sqrt(3)*spacing), n=1)
        A = 12.0 * _np.sqrt(_np.asarray(4.0/3.0 * (wavelength/spacing)**2.0 - 1.0, dtype=complex)) \
            * (J1prime/(1.0-(4*_np.pi*radius/(1.841*_np.sqrt(3.0)*spacing))**2.0))**2.0
        A -= 12.0/_np.sqrt(_np.asarray(4.0/3.0 * (wavelength/spacing)**2.0 - 1.0, dtype=complex)) \
            * (J1prime/(4.0*_np.pi*radius/(_np.sqrt(3.0)*spacing)))**2.0

        B = 0.33*(spacing/radius)**2.0 * _np.sqrt(_np.asarray((0.293*wavelength/radius)**2.0 - 1.0, dtype=complex) )

        beta = (2.0*_np.pi/wavelength)*_np.sqrt(_np.asarray((0.293*wavelength/radius)**2.0 - 1.0, dtype=complex))

        R2 = _np.zeros( (len(freq),), dtype=complex)
        T2 = _np.zeros_like(R2)
        for ii in range(len(freq)):

            AA = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness)))
            BB = 1.0/  (1.0 - 1j*(A[ii]+B[ii]* 1.0/cmath.tanh(beta[ii]*thickness)))

            # Reflection
            R2[ii] = AA.copy() + BB.copy() - 1.0

            # Transmission
            T2[ii] = AA.copy() - BB.copy()

            # R2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness))) + 1.0/(1.0-1j*(A[ii]+B[ii]*1.0/cmath.tanh(beta[ii]*thickness))) - 1.0
            # T2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness))) - 1.0/(1.0-1j*(A[ii]+B[ii]*1.0/cmath.tanh(beta[ii]*thickness)))
            # if self.verbose:    print(_np.abs(R2[ii]), _np.abs(1-T2[ii]))

            # Insertion delay - Phase delay caused by guide (degrees)

        # end for

        # For oblique incidence, there is a correction here:
        porosity = _np.pi*(2.0*radius)**2.0 / (2.0*_np.sqrt(3)*spacing**2.0)
        T2perp = T2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
        T2parr = T2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))

        # Does it work the same for reflection?
        # R2_perp = R2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
        # R2_parr = R2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))

#        phperp = _np.arctan(_np.imag(T2perp) /_np.real(T2perp))*180.0/_np.pi - 360.0*thickness*_np.cos(th*_np.pi/180.0)/wavelength  # degrees
#        phparr = _np.arctan(_np.imag(T2parr) /_np.real(T2parr))*180.0/_np.pi - 360.0*thickness*_np.cos(th*_np.pi/180.0)/wavelength  # degrees

        T2perp = _np.abs(T2perp)
        T2parr = _np.abs(T2parr)

#        R2perp = _np.abs(R2perp)
#        R2parr = _np.abs(R2parr)

#        R2perp = 1.0-T2perp
#        R2parr = 1.0-T2parr

        if self.verbose:
            print("Dichroic plate characteristics: ")
            print("Hexagonal hole pattern: diameter=%2.2f mm, spacing=%2.2f mm, thickness=%2.2f mm"%(1e3*2.0*radius, 1e3*spacing, 1e3*thickness))
            print("filter cut-offs: %3.1f<f<%3.1f GHz"%(fc1, fc2))
        # end if

        return T2perp, T2parr, porosity, fc1, fc2

    def response(self, freq, diameter, spacing, thickness, th, plotit=True, wd=None):
         """

         """

         if freq is None and hasattr(self, 'freq'):
             freq=self.freq
         elif freq is None:
             freq = 1e9*_np.linspace(100.0, 250.0, 250)
         # endif

         T2_perp_140, T2_parr_140, _, _, _ = self.plate(freq, diameter=diameter, spacing=spacing, thickness=thickness, th=th, shape='hex')
         T2perp, T2parr, porosity, fc1, fc2 = self.plate(freq, diameter=diameter, spacing=spacing, thickness=thickness, th=th, shape='hex')

         T2_perp_log = 20*_np.log10(T2perp)
         T2_parr_log = 20*_np.log10(T2parr)
         por_log = 10*_np.log10(porosity)

         # ======================================= #

         delimiter = '\n'
         hdr = "Dichroic plate characteristics: Filled with %s"%(matname,) + delimiter
         hdr += "Hexagonal hole pattern (%i holes): diameter=%2.2f mm, spacing=%2.2f mm, thickness=%2.1f mm"%(ncircles, 1e3*diameter, 1e3*spacing, 1e3*thickness) + delimiter
         hdr += "filter cut-offs: %3.1f<f<%3.1f GHz"%(fc1, fc2) + delimiter
         hdr += "Power transmission (perpendicular): %3.1f dB@%3.0f GHz"%(T2_perp_140, 140) + delimiter
         hdr += "Power transmission (parallel): %3.1f dB@%3.0f GHz"%(T2_parr_140, 140) + delimiter
         hdr += "Porosity limit (%0.2f): %3.1f dB"%(porosity, por_log) + delimiter

         print(hdr)

         filnam = _os.path.join(wd,'DichroicPlate_holes_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fc1,1e3*diameter,1e3*spacing,1e3*thickness))
         _np.savetxt(filnam, 1e3*centers, fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n%6s 6%s'%('x[mm]', 'y[mm]') )

         filnam = _os.path.join(wd,'DichroicPlate_Transmission_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fc1,1e3*diameter,1e3*spacing,1e3*thickness))
         _np.savetxt(filnam, (freq,T2parr,T2perp), fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n %8s %8s %8s'%('freq[GHz]','T2[parr]', 'T2[perp]'))




         if plotit:
             # ======================================= #

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
             _plt.title(r'Power Transmission Coefficient: f$_{c,o}$<%3.1f, f$_{c,d}$<%3.1f GHz'%(fc1,fc2) )
             _plt.axvline(x=fc1, linestyle='--', color='k')
             _plt.axvline(x=fc2, linestyle='--', color='k')
             _plt.axhline(y=por_log, linestyle='--', color='k')

             _plt.text(x=fc1+5, y=-15, s='Hexagonal hole pattern: \n diameter=%2.2f mm, \n spacing=%2.2f mm, \n thickness=%2.1f mm'%(1e3*diameter, 1e3*spacing, 1e3*thickness))



         if wd is not None and plotit and _os.path.exists(wd):
             # wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
             hfig.savefig(_os.path.join(wd,'DichroicNotch_%3.1fGHz_%3.1fGHz.png'%(fc1[0],fc2[1])), dpi=200, transparent=True)


#
#    def notch_response(self, freq=None, diameter=1.6e-3, spacing=2e-3, thickness=5e-3, th=45, plotit=True, wd=None):
#        if freq is None and hasattr(self, 'freq'): freq=self.freq # endif
#
#        T2perp, T2parr, por, fc1, fc2 = self.dichroic_plate(freq, diameter, spacing, thickness, th)
#
#        T2_perp_log, T2_parr_log = [_np.zeros((len(freq),2), dtype=float) for ii in range(2)]
#        porosity, fco, fcd = [_np.zeros((2,), dtype=float) for ii in range(3)]
#        T2_perp_log[:,0], T2_parr_log[:,0], porosity[0], fco[0], fcd[1] = \
#            self.dichroic_plate(0.5*diameter.copy(), spacing.copy(), thickness.copy())
#        T2_perp_log[:,1], T2_parr_log[:,1], porosity[1], fco[0], fcd[1] = \
#            self.dichroic_plate(0.5*diameter.copy(), spacing.copy(), thickness.copy())
#
#        T2_perp_log = T2_perp_log[:,0]*(1-T2_perp_log[:,1])
#        T2_parr_log = T2_parr_log[:,0]*(1-T2_parr_log[:,1])
#
#        T2_perp_log = 20*_np.log10(T2_perp_log)
#        T2_parr_log = 20*_np.log10(T2_parr_log)
#
#        if plotit:
#            # ======================================= #
#
#            hfig = _plt.figure()
#            _plt.plot(1e-9*freq, T2_perp_log, '-')
#            _plt.plot(1e-9*freq, T2_parr_log, '--')
#            _plt.xlim((105,180))
#
#            _plt.xlabel('frequency [GHz]')
#            _plt.ylabel(r'|T$^2$| [dB]')
#            #_plt.title(r'Power Transmission Coefficient: f$_{c,o}$<%3.1f, f$_{c,d}$<%3.1f GHz'%(fco,fcd) )
#            _plt.axvline(x=fco[0], linestyle='--', color='k')
#            _plt.axvline(x=fco[1], linestyle='--', color='k')
#
#        if wd is not None and plotit and _os.path.exists(wd):
#            # wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
#            hfig.savefig(_os.path.join(wd,'DichroicNotch_%3.1fGHz_%3.1fGHz.png'%(fco[0],fco[1])), dpi=200, transparent=True)

    @staticmethod
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

    @staticmethod
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

# end class special





# ========================================================================== #
# ========================================================================== #

def Basic_antenna_plot(freq, w0x, w0y, z0x, z0y):
    # axis along which to determine data: antenna position
    zantenna = -z0x if _np.abs(-z0x)>_np.abs(-z0y) else -z0y
    zz = zantenna + _np.linspace(0, 1.00, num=250)

    tst_x = qoptics(freq=freq, wo=w0x, zz=zz)
    tst_x.RayleighRange()
    tst_x.SpotSize()  # calculates beam radius along assigned vector "zz" asssuming free space

    tst_y = qoptics(freq=freq, wo=w0y, zz=zz)
    tst_y.RayleighRange()
    tst_y.SpotSize()  # calculates beam radius along assigned vector "zz" asssuming free space

#    dvrg_x = (180.0/_np.pi)*tst_x.lambda0/(_np.pi*tst_x.Nrefr*tst_x.wo)
#    dvrg_y = (180.0/_np.pi)*tst_y.lambda0/(_np.pi*tst_y.Nrefr*tst_y.wo)

    _plt.figure()
    _plt.plot(zz, tst_x.wz, 'k--')
    _plt.plot(zz, tst_y.wz, 'k-')
    _plt.xlabel('distance along beam-axis')
    _plt.ylabel('beam radius (x-- / y-)')
    _plt.axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')

    xlims = _plt.xlim()
    ylims = _plt.ylim()
#    xmax = zz[-1]
    # ymax = tst.wz[-1]
#    ymax = _np.nanmax((_np.tan(dvrg_x*_np.pi/180.0),_np.tan(dvrg_y*_np.pi/180.0)))*xmax
#    _plt.plot([0, xmax], [0.0, ymax], 'k--')
    _plt.xlim(xlims)
    _plt.ylim([0, ylims[1]])
# end def

def Basic_antenna_thinlens_plot(freq, w0x, w0y, z0x, z0y, d1=0.1, flen=0.5):

    # axis along which to determine data: antenna position
    zantenna = -z0x if _np.abs(-z0x)>_np.abs(-z0y) else -z0y
    zz = zantenna + _np.linspace(0, 1.00, num=250)

    tst_x = qoptics(freq=freq, wo=w0x, zz=zz)
    tst_y = qoptics(freq=freq, wo=w0y, zz=zz)

    tst_x.add_element(tst_x.thinlens(flen), d1, 'Mirror 1')
    tst_y.add_element(tst_y.thinlens(flen), d1)

    tst_x.BeamPropagation()
    tst_y.BeamPropagation()

    _ax1, _ax2 = tst_x.plot()
    tst_y.plot(_ax1=_ax1, _ax2=_ax2)
# end def

def qme_telescope_design():
    """
    """
    freq = 140e9

    # axis along which to determine data: antenna position
    zantenna = 0.0
    zz = zantenna + _np.linspace(0, 3.00, num=500)

    tst = qoptics(freq=freq, wo=3.2e-3, zz=zz)

    tst.add_element(tst.thinlens(171.37e-3), 175e-3, '1 Ellipsoidal\nMirror')
    tst.add_element(tst.thinlens(1078.6e-3), 1250e-3, '2 Plane\nMirror')
    tst.add_element(tst.thinlens(1078.6e-3), 1079e-3, '3 Ellipsoidal\nMirror')
    tst.add_element(tst.thinlens(1078.6e-3), 1250e-3, '4 Plane\nMirror')

    tst.BeamPropagation()
    _ax1, _ax2 = tst.plot()
# end def


def qme_telescope_op11():
    """
    """
    freqs, w0, z0, hom, purity = __qme_op11_antenna__()

    freq = 140e9
    wavelength = cc/freq
#    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
#    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3

    # axis along which to determine data: antenna position
#    zantenna = -z0x if _np.abs(-z0x)>_np.abs(-z0y) else -z0y
#    zantenna = 0.0
    zantenna = -z0y
    zz = zantenna + _np.linspace(0, 3.00, num=500)

    tst = qoptics(freq=freq, wo=w0y, zz=zz)

    # d23 = f1 + f2 = 1250e-3î
    d1 = 175e-3;       f1 = 171.4e-3;        # R1=176.34mm, R2=6080.8mm; a=3.1286m, b=0.73222 m
    d2 = d1 + 180e-3; # planar mirror
    d3 = d2 + 1071e-3; f3 = 1078.6e-3  # ellipsoidal mirror,
    d4 = d3 + 189.3e-3 # planar mirror

    # type = reflector, combiner, splitter
    #  number of unique paths through system ~ 1 + splitter + 2*combiner
    element1 = {'incangle': 45.00, 'shape':'ellipsoid', 'type':'reflector', 'dim':127.3e-3, "amaj":3.1286, "bmin":0.73222}
    element2 = {'incangle':-45.00, 'shape':'planar',    'type':'reflector', 'dim':127.3e-3}
    element3 = {'incangle':-26.45, 'shape':'ellipsoid', 'type':'reflector', 'dim':100.0e-3} #, "amaj":, "bmin":}
    element4 = {'incangle': 29.74, 'shape':'planar',    'type':'reflector', 'dim':100.0e-3}

    tst.add_element(tst.thinlens(f1),   d1, '1 Ellipsoidal\nMirror')
    tst.add_element(tst.planarmirror(), d2, '2 Plane\nMirror')
    tst.add_element(tst.thinlens(f3),   d3, '3 Ellipsoidal\nMirror')
    tst.add_element(tst.planarmirror(), d4, '4 Plane\nMirror')

#    tst.add_element_params(**element1)
#    tst.add_element_params(**element2)
#    tst.add_element_params(**element3)
#    tst.add_element_params(**element4)

#    d1 = 175e-3; f1 = 171.37e-3; th1=45  # [deg], R1=176.34mm, R2=6080.8mm; a=3.128.6m, b=0.73222 m      # analysis:ignore
#    d2 = d1+189e-3;  th2 = 45            # [deg], planar mirror                                   # analysis:ignore
#    d3 = d2 + 1250e-3; f3=1078.6e-3; th3 = 37.1  # 90-52.9 degrees from horizonal, ellipsoidal mirror,  # analysis:ignore
#    d4 = d3 + 208e-3;  th4 = 13.16      # [deg]                       # analysis:ignore

    tst.rays = []
    angrange = wavelength/(_np.pi*w0y)*_np.linspace(0, 1, num=5, endpoint=True)
    for ii in range(5):
#        tst.rays.append([ii*0.001, 0.0])
        tst.rays.append([0.0, angrange[ii]])
    # end for
    tst.BeamPropagation()
    _ax1, _ax2 = tst.plot()
    tst.plotrays()
#    tst.layout()
# end def

def qme_telescope_op12():
    """
    """
    freqs, w0, z0, hom, purity = __qme_op12_antenna__()

    freq = 140e9
    wavelength = cc/freq
#    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
#    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3

    # axis along which to determine data: antenna position
#    zantenna = -z0x if _np.abs(-z0x)>_np.abs(-z0y) else -z0y
    zantenna = -z0y
    zz = zantenna + _np.linspace(0, 3.00, num=500)

    tst = qoptics(freq=freq, wo=w0y, zz=zz)

    # d23 = f1 + f2 = 1250e-3î
    # 175e-3 - 53e-3 = 122e-3 to aperture of antenna
    d1 = 122e-3;  f1 = 171.4e-3;  th1 = 45.0    # [deg], R1=176.34mm, R2=6080.8mm; a=3.128.6m, b=0.73222 m      # analysis:ignore
    d2 = d1 + 180e-3;   f2 = 1e8;       th2 = 45.0    # [deg], planar mirror                                   # analysis:ignore
    d3 = d2 + 1071e-3;  f3 = 1078.6e-3; th3 = 90.0-63.55  # ellipsoidal mirror,  # analysis:ignore
    d4 = d3 + 189.3e-3; f4 = 1e8;       th4 = 13.16   # [deg]                       # analysis:ignore

    tst.add_element(tst.thinlens(f1), d1, '1 Ellipsoidal\nMirror')
    tst.add_element(tst.thinlens(f2), d2, '2 Plane\nMirror')
    tst.add_element(tst.thinlens(f3), d3, '3 Ellipsoidal\nMirror')
    tst.add_element(tst.thinlens(f4), d4, '4 Plane\nMirror')

    tst.rays = []
    angrange = wavelength/(_np.pi*w0y)*_np.linspace(0, 1, num=5, endpoint=True)
    for ii in range(5):
#        tst.rays.append([ii*0.001, 0.0])
        tst.rays.append([0.0, angrange[ii]])
    # end for
    tst.BeamPropagation()
    _ax1, _ax2 = tst.plot()
    tst.plotrays()
# end def

def __Plaum_20200212__():
    freqs, wx, wy, zx, zy, hom, purity = [[] for _ in range(7)]
    freqs.append(129e9)
    wy.append(16.43e-3)
    wx.append(17.59e-3)
    zy.append(325e-3)
    zx.append(258e-3)
    hom.append(18.5e-2)   #  % as in 18.5 percent higher order modes (non-Gaussian)
    purity.append(None)

    freqs.append(138e9)
    wy.append(18.45e-3)
    wx.append(18.23e-3)
    zy.append(308e-3)
    zx.append(403e-3)
    hom.append(19.61e-2)   #  % as in 18.5 percent higher order modes (non-Gaussian)
    purity.append(None)

    freqs = _np.asarray(freqs)
    w0 = _np.vstack((wx, wy))
    z0 = _np.vstack((zx, zy))
    hom = _np.asarray(hom)
    purity = _np.asarray(purity)
    return freqs, w0, z0, hom, purity

def __qme_op11_antenna__():
    freqs, wx, wy, zx, zy, hom, purity = [[] for _ in range(7)]
    freqs.append(129e9)
    wy.append(3.40e-3)
    wx.append(3.04e-3)
    zy.append(-4.06e-3)
    zx.append(-2.90e-3)
    hom.append(0.79e-2)   #  % as in
    purity.append(None)

    freqs.append(138e9)
    wy.append(3.32e-3)
    wx.append(3.06e-3)
    zy.append(-1.22e-3)
    zx.append(-1.44e-3)
    hom.append(0.44e-2)   #  % as in
    purity.append(None)

    freqs = _np.asarray(freqs)
    w0 = _np.vstack((wx, wy))
    z0 = _np.vstack((zx, zy))
    hom = _np.asarray(hom)
#    purity = _np.asarray(purity)
    return freqs, w0, z0, hom, purity

def __qme_op12_antenna__():
    """
    antenna parameters from 1-QME11-Q0004.0
    OP12 antennas from the QME diagnostic
    """
    freqs, wx, wy, zx, zy, hom, purity = [[] for _ in range(7)]

    offset = 105e-3
    freqs.append(50e9)
    wy.append(7.41e-3)
    wx.append(6.25e-3)
    zy.append(offset-0.21051)
    zx.append(offset-0.23370)
    hom.append(0.87e-2)
    purity.append(0.9913)

    freqs.append(62e9)
    wy.append(6.41e-3)
    wx.append(5.47e-3)
    zy.append(offset-0.14214)
    zx.append(offset-0.15406)
    hom.append(1.10e-2)
    purity.append(0.9890)

    freqs.append(90e9)
    wy.append(4.94e-3)
    wx.append(3.96e-3)
    zy.append(offset-0.15048)
    zx.append(offset-0.15848)
    hom.append(1.43e-2)
    purity.append(0.9857)

    freqs.append(113e9)
    wy.append(4.23e-3)
    wx.append(3.40e-3)
    zy.append(offset-0.15114)
    zx.append(offset-0.15719)
    hom.append(2.56e-2)
    purity.append(0.9744)

    freqs.append(140e9)
    wy.append(4.31e-3)
    wx.append(3.44e-3)
    zy.append(offset-0.14735)
    zx.append(offset-0.16142)
    hom.append(1.55e-2)
    purity.append(0.9845)

    freqs.append(170e9)
    wy.append(3.23e-3)
    wx.append(2.85e-3)
    zy.append(offset-0.15475)
    zx.append(offset-0.15046)
    hom.append(3.18e-2)
    purity.append(0.9682)

    freqs.append(193e9)
    wy.append(2.48e-3)
    wx.append(2.27e-3)
    zy.append(offset-0.15053)
    zx.append(offset-0.15944)
    hom.append(2.28e-2)
    purity.append(0.9771)

    freqs = _np.asarray(freqs)
    w0 = _np.vstack((wx, wy))
    z0 = _np.vstack((zx, zy))
    hom = _np.asarray(hom)
    purity = _np.asarray(purity)
    return freqs, w0, z0, hom, purity

def Plaum_20200212():
    # start by calculating the free space propagation from a horn a gaussian
    # beam of specified waist at a specfied position
    freqs, w0, z0, hom, purity = __Plaum_20200212__()

    freq = 138e9
    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3

    Basic_antenna_plot(freq, w0x, w0y, z0x, z0y)
# end def

def qme_op11_antenna():
    """
    propagate a quasi-optical beam through free-space
    """
    freqs, w0, z0, hom, purity = __qme_op11_antenna__()

    freq = 138e9
    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3

    Basic_antenna_plot(freq, w0x, w0y, z0x, z0y)
# end def

def qme_op12_antenna():
    """
    propagate a quasi-optical beam through free-space
    """
    freqs, w0, z0, hom, purity = __qme_op12_antenna__()

    freq = 140e9
    # freq = 132e9
#    freq = 113e9
    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3

    Basic_antenna_plot(freq, w0x, w0y, z0x, z0y)
# end def

def abcd_prop_test():
    freqs, w0, z0, hom, purity = __qme_op12_antenna__()

    freq = 140e9
    wavelength = cc/freq
#    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
#    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3

    # axis along which to determine data: antenna position
#    zantenna = -z0x if _np.abs(-z0x)>_np.abs(-z0y) else -z0y
    zantenna = -z0y
    zz = zantenna + _np.linspace(0, 1.00, num=250)

    tst_y = qoptics(freq=freq, wo=w0y, zz=zz)
#    tst_y.rays = [0.02, 0*_np.pi/180.0]
    tst_y.rays = []
    angrange = wavelength/(_np.pi*w0y)*_np.linspace(0, 1, num=5, endpoint=True)
    for ii in range(5):
#        tst_y.rays.append([0.025*ii, 0*_np.pi/180.0])
        tst_y.rays.append([0.0, angrange[ii]])
    # end for

#    tst_y.BeamPropagation()
#    _ax1, _ax2 = tst_y.plot()
#    tst_y.plotrays()

    # curved mirror at 45 degrees
    th = 45*_np.pi/180.0
    Rmirror  = 0.10   # for an elliptic mirror flen = R1*R2/(R2+R1), if R2>>R1, then flen = R1
    tst_y.add_element(tst_y.sphericalmirror(Rmirror, th, 'vertical'), 0.15, 'Mirror 1')

#    tst_y.BeamPropagation()
#    _ax1, _ax2 = tst_y.plot()

    # thin lens
    tst_y.add_element(tst_y.thinlens(0.1), 0.35, 'Mirror 2')

#    tst_y.BeamPropagation()
#    _ax1, _ax2 = tst_y.plot()

    # planarinterface
    eps_r, loss_tangent = tst_y._mats(material='mica')
    tst_y.add_element(tst_y.refraction_thickplate(1.0, _np.sqrt(eps_r), 2.0e-3), 0.55, 'Mica Window')

#    tst_y.BeamPropagation()
#    _ax1, _ax2 = tst_y.plot()

    # thick lens
    eps_r, loss_tangent = tst_y._mats(material='teflon')
    tst_y.add_element(tst_y.thicklens(1.0, _np.sqrt(eps_r), 0.07, 0.7, 0.05), 0.85, 'Teflon lens')


    tst_y.BeamPropagation()
    _ax1, _ax2 = tst_y.plot()
    tst_y.plotrays()
    return tst_y
# end def
# end def

#def qme_op11_antenna_mirror():
#    """
#    propagate a quasi-optical beam through free-space, to a mirror (thin lens), and onwards
#    """
#    from matplotlib.patches import Arc, Ellipse
#    mirdim = 0.20  # 10 cm mirror dimension
#    truncation_level = 3.0  # the mirror should be 5 x larger than the beam waist at its position (no. standard deviations)
#    inc_angle = 45   # deg, angle of incidence on ellipsoidal mirror
#    mirror2plasma = 0.40  # distance from mirror center to plasma focal point
#
#    # freq = 129e9
#    # w0x = 3.04e-3
#    # w0y = 3.40e-3
#    # z0x = -2.90e-3    # waist position from aperture
#    # z0y = -4.06e-3
#
#    freq = 138e9
#    w0x = 3.06e-3
#    w0y = 3.32e-3
#    z0x = -1.44e-3    # waist position from aperture
#    z0y = -1.22e-3
#
#    # axis along which to determine data
#    # zantenna = -z0x   # aperture of antenna
#    zantenna = -z0y   # aperture of antenna
#    zz = zantenna + _np.linspace(0, 0.50, num=250)
#
#    # start by calculating the free space propagation from a horn a gaussian
#    # beam of specified waist at a specfied position
#    # tst = qoptics(freq=freq, wo=w0x, zz=zz)
#    tst = qoptics(freq=freq, wo=w0y, zz=zz)
#    # tst.wo = w0x   # assignment also calculated Rayleigh range
#    # tst.zz = zz
#    tst.RayleighRange()         # calculate the confocal parameter / depth of focus
#    tst.ComplexBeamParameter()
#    tst.SpotSize()  # calculates beam radius along assigned vector "zz" asssuming free space
#
#    # ======================== #
#
#    # equations of an ellipse
#    #  R1 - input radii of curvature of wavefront. R1 is the distance between focii F1 (antenna) and the mirror center
#    #  R2 - output radii of curvature of wavefront. R2 is the distance between focii F2 (in plasma) and the mirror center
#    #  R1+R2 = 2bmaj;   where x^2/a^2 + y^2/b^2 = 1 and a<b
#    #  flens = R1*R2/(R1+R2)=R1*R2/(2bmaj)   is the focal length
#    #     R2 = 2*bmaj*flens / R1
#    bmaj = 0.5*mirdim
#
#    # use the calculated beam radius to determine the necessary mirror position
#    ilens = _np.argwhere(tst.wz>=bmaj/truncation_level).squeeze()[0]
#    # flens = 10e-2
#    # flens = 11.8e-2
#    # flens = _np.copy(tst.zz[ilens])
#    zlens = _np.copy(tst.zz[ilens])
#    Rparr = _np.copy(tst.Rz[ilens])     # Radius of curvature in parallel to the plane of incidence
#    Rperp = bmaj**2.0/Rparr             # Radius of curvature in perpendicular to the plane of incidence
#    amin = Rperp
#
#    # F1 = _np.sqrt(
#    # ================ #
#    angle = _np.linspace(0, 2*_np.pi, 180)
#    elli = _np.vstack((bmaj*_np.cos(angle), amin*_np.sin(angle))).T
#
#    _plt.figure()
#    _plt.plot(elli[:,0], elli[:,1], 'k--')  # ellipse
#    _plt.plot(bmaj*_np.linspace(0.0, 1.0, num=30), _np.zeros((30,)), 'k-')
#    _plt.plot(_np.zeros((30,)), amin*_np.linspace(0.0, 1.0, num=30), 'k-')
#    _plt.text(x=0.5*bmaj, y=0.1*amin, s=r'$b_{maj}$=%4.2f'%(bmaj,))
#    _plt.text(x=0.1*bmaj, y=0.5*amin, s=r'$a_{min}$=%4.2f'%(amin,))
#    # _plt.plot(
#    # ================ #
#
#    fparr = 0.5*Rparr*_np.cos(inc_angle*_np.pi/180.0)
#    fperp = 0.5*Rperp/_np.cos(inc_angle*_np.pi/180.0)
#    # flens = fparr
#    flens = Rparr*mirror2plasma/(Rparr + mirror2plasma)
#
#    # qx = _np.zeros((len(zz),), dtype=complex)
#    qx = _np.copy(tst.qz)
#    for ii, zx in enumerate(zz):
#
#        if zx<=zlens:
#            # ABCD = tst.freespace(zx)
#            pass
#        elif zx>zlens:
#            ABCD = tst.freespace_thinlens_freespace(zlens, flens, zx-zlens)
#
#            qx[ii] = tst.propagate_beamparameter(ABCD, tst.qz[0])
#        # end if
#
#        # qx[ii] = tst.propagate_beamparameter(ABCD, tst.qz[0])
#    # end for
#    tst.qz = _np.copy(qx)
#
#    # ================== #
#
#    dvrg = (180.0/_np.pi)*tst.lambda0/(_np.pi*tst.Nrefr*tst.wo)
#
#    _plt.figure()
#    _plt.plot(zz, tst.wz, 'k-')
#    _plt.xlabel('distance along beam-axis')
#    _plt.ylabel('beam radius')
#    _plt.axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
#
#    xlims = _plt.xlim()
#    ylims = _plt.ylim()
#    _plt.xlim(xlims)
#    _plt.ylim([0, ylims[1]])
#
#    # ================== #
#
#    # fwhm = _np.sqrt(2*_np.log(2))*tst.wz        # FWHM of the quasi-optical beam
#    # kmax_antenna = 2.0*_np.pi/fwhm
#    kmax_antenna = 2.0*_np.pi/(2.0*tst.wz)      # attenuation begins here for rectangular antenna pattern
#    # kmax_antenna = 2.0*_np.sqrt(2)/(2.0*tst.wz)   # 3 dB point in wavenumber sensitivity
#    rhos = _np.sqrt(1.67e-27*1.6e-19*500)/(1.6e-19*2.4)
#
#    _, _ax = _plt.subplots(2,1, sharex='all', num='kres')
#    _ax[0].plot(zz, 1e-2*kmax_antenna, 'k-')
#    _ax[0].set_ylabel('Max. K resolved [cm-1]')
#    _ax[0].axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
#
#    _ax[1].plot(zz, kmax_antenna*rhos, 'k-')
#    _ax[1].set_ylabel(r'($k_\bot\rho_s$)_max: $\rho_s$=1e-3')
#    _ax[1].set_xlabel('distance along beam-axis')
#    _ax[1].axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
#
#    xlims = _ax[0].get_xlim()
#    ylims = _ax[0].get_ylim()
#    _ax[0].set_xlim(xlims)
#    _ax[0].set_ylim([0, ylims[1]])
#    ylims = _ax[1].get_ylim()
#    _ax[1].set_ylim([0, ylims[1]])
#
#    print((amin, bmaj, flens, fperp, fparr))
## end def
#
#def qme_op12_antenna_mirror():
#    """
#    propagate a quasi-optical beam through free-space, to a mirror (thin lens), and onwards
#    """
#    from matplotlib.patches import Arc, Ellipse
#    mirdim = 0.20  # 10 cm mirror dimension
#    truncation_level = 3.0  # the mirror should be 5 x larger than the beam waist at its position (no. standard deviations)
#    inc_angle = 45   # deg, angle of incidence on ellipsoidal mirror
#    mirror2plasma = 0.40  # distance from mirror center to plasma focal point
#
#    freqs, w0, z0, hom, purity = __qme_op12_antenna__()
#
#    freq = 140e9
#    w0x = _np.interp(freq, freqs, w0[0,:]) # 4.31e-3
#    w0y = _np.interp(freq, freqs, w0[1,:]) # 3.44e-3
#    z0x = _np.interp(freq, freqs, z0[0,:]) # 147.35e-3    # waist position from aperture
#    z0y = _np.interp(freq, freqs, z0[1,:]) # 161.42e-3
#
#    # axis along which to determine data
#    # zantenna = -z0x   # aperture of antenna
#    zantenna = -z0y   # aperture of antenna
#    zz = zantenna + _np.linspace(0, 0.50, num=250)
#
#    # start by calculating the free space propagation from a horn a gaussian
#    # beam of specified waist at a specfied position
#    # tst = qoptics(freq=freq, wo=w0x, zz=zz)
#    tst = qoptics(freq=freq, wo=w0y, zz=zz)
#    # tst.wo = w0x   # assignment also calculated Rayleigh range
#    # tst.zz = zz
#    tst.RayleighRange()         # calculate the confocal parameter / depth of focus
#    tst.ComplexBeamParameter()
#    tst.SpotSize()  # calculates beam radius along assigned vector "zz" asssuming free space
#
#    # ======================== #
#
#    # equations of an ellipse
#    #  R1 - input radii of curvature of wavefront. R1 is the distance between focii F1 (antenna) and the mirror center
#    #  R2 - output radii of curvature of wavefront. R2 is the distance between focii F2 (in plasma) and the mirror center
#    #  R1+R2 = 2bmaj;   where x^2/a^2 + y^2/b^2 = 1 and a<b
#    #  flens = R1*R2/(R1+R2)=R1*R2/(2bmaj)   is the focal length
#    #     R2 = 2*bmaj*flens / R1
#    bmaj = 0.5*mirdim
#
#    # use the calculated beam radius to determine the necessary mirror position
#    ilens = _np.argwhere(tst.wz>=bmaj/truncation_level).squeeze()[0]
#    # flens = 10e-2
#    # flens = 11.8e-2
#    # flens = _np.copy(tst.zz[ilens])
#    zlens = _np.copy(tst.zz[ilens])
#    Rparr = _np.copy(tst.Rz[ilens])     # Radius of curvature in parallel to the plane of incidence
#    Rperp = bmaj**2.0/Rparr             # Radius of curvature in perpendicular to the plane of incidence
#    amin = Rperp
#
#    # F1 = _np.sqrt(
#    # ================ #
#    angle = _np.linspace(0, 2*_np.pi, 180)
#    elli = _np.vstack((bmaj*_np.cos(angle), amin*_np.sin(angle))).T
#
#    _plt.figure()
#    _plt.plot(elli[:,0], elli[:,1], 'k--')  # ellipse
#    _plt.plot(bmaj*_np.linspace(0.0, 1.0, num=30), _np.zeros((30,)), 'k-')
#    _plt.plot(_np.zeros((30,)), amin*_np.linspace(0.0, 1.0, num=30), 'k-')
#    _plt.text(x=0.5*bmaj, y=0.1*amin, s=r'$b_{maj}$=%4.2f'%(bmaj,))
#    _plt.text(x=0.1*bmaj, y=0.5*amin, s=r'$a_{min}$=%4.2f'%(amin,))
#    # _plt.plot(
#    # ================ #
#
#    fparr = 0.5*Rparr*_np.cos(inc_angle*_np.pi/180.0)
#    fperp = 0.5*Rperp/_np.cos(inc_angle*_np.pi/180.0)
#    # flens = fparr
#    flens = Rparr*mirror2plasma/(Rparr + mirror2plasma)
#
#    # qx = _np.zeros((len(zz),), dtype=complex)
#    qx = _np.copy(tst.qz)
#    for ii, zx in enumerate(zz):
#
#        if zx<=zlens:
#            # ABCD = tst.freespace(zx)
#            pass
#        elif zx>zlens:
#            ABCD = tst.freespace_thinlens_freespace(zlens, flens, zx-zlens)
#
#            qx[ii] = tst.propagate_beamparameter(ABCD, tst.qz[0])
#        # end if
#
#        # qx[ii] = tst.propagate_beamparameter(ABCD, tst.qz[0])
#    # end for
#    tst.qz = _np.copy(qx)
#
#    # ================== #
#
#    dvrg = (180.0/_np.pi)*tst.lambda0/(_np.pi*tst.Nrefr*tst.wo)
#
#    _plt.figure()
#    _plt.plot(zz, tst.wz, 'k-')
#    _plt.xlabel('distance along beam-axis')
#    _plt.ylabel('beam radius')
#    _plt.axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
#
#    xlims = _plt.xlim()
#    ylims = _plt.ylim()
#    _plt.xlim(xlims)
#    _plt.ylim([0, ylims[1]])
#
#    # ================== #
#
#    # fwhm = _np.sqrt(2*_np.log(2))*tst.wz        # FWHM of the quasi-optical beam
#    # kmax_antenna = 2.0*_np.pi/fwhm
#    kmax_antenna = 2.0*_np.pi/(2.0*tst.wz)      # attenuation begins here for rectangular antenna pattern
#    # kmax_antenna = 2.0*_np.sqrt(2)/(2.0*tst.wz)   # 3 dB point in wavenumber sensitivity
#    rhos = _np.sqrt(1.67e-27*1.6e-19*500)/(1.6e-19*2.4)
#
#    _, _ax = _plt.subplots(2,1, sharex='all', num='kres')
#    _ax[0].plot(zz, 1e-2*kmax_antenna, 'k-')
#    _ax[0].set_ylabel('Max. K resolved [cm-1]')
#    _ax[0].axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
#
#    _ax[1].plot(zz, kmax_antenna*rhos, 'k-')
#    _ax[1].set_ylabel(r'($k_\bot\rho_s$)_max: $\rho_s$=1e-3')
#    _ax[1].set_xlabel('distance along beam-axis')
#    _ax[1].axvline(x=zantenna, linewidth=1.0, color='k', linestyle='--')
#
#    xlims = _ax[0].get_xlim()
#    ylims = _ax[0].get_ylim()
#    _ax[0].set_xlim(xlims)
#    _ax[0].set_ylim([0, ylims[1]])
#    ylims = _ax[1].get_ylim()
#    _ax[1].set_ylim([0, ylims[1]])
#
#    print((amin, bmaj, flens, fperp, fparr))
## end def

# ========================================================================== #
# ========================================================================== #


if __name__=="__main__":
    qme_op11_antenna()
#    qme_op12_antenna()
    Plaum_20200212()

#    qo = abcd_prop_test()

#    qme_telescope_op11()
#    qme_telescope_op12()
#
#    qo = qoptics_abcd()
#    foc, M, (amaj, bmin), wproj = qo.elliptic_mirror_design(140.0e9, 30e-3, 20e-3, 200e-3, 300e-3, grazing_angle=15*_np.pi/180.0, trunclevel=5.0, phi=0)
#    foc, M, (amaj, bmin), wproj = qo.elliptic_mirror_design(140.0e9, 3.2e-3, 20e-3, 200e-3, 300e-3, grazing_angle=45*_np.pi/180.0, trunclevel=5.0, phi=-_np.pi/6.0)
#    foc, M, (amaj, bmin), wproj = qo.elliptic_mirror_design(140.0e9, 3.2e-3, 20e-3, 200e-3, 300e-3, grazing_angle=60*_np.pi/180.0, trunclevel=5.0, phi= _np.pi/3.0)
#    foc, M, (amaj, bmin), wproj = qo.elliptic_mirror_design(140.0e9, 3.2e-3, 20e-3, 100e-3, 500e-3, grazing_angle=30*_np.pi/180.0, trunclevel=5.0, phi=-_np.pi/3.0)
#    foc, M, (amaj, bmin), wproj = qo.elliptic_mirror_design(142.75e9, 3.2e-3, 34.83e-3, 175e-3, 601.06e-3, grazing_angle=45*_np.pi/180.0, trunclevel=5.0, phi=0.0)
# end if



# ========================================================================== #
# ========================================================================== #
