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
import scipy as _scipy
import os as _os
import matplotlib.pyplot as _plt
import cmath

# Pybaseutils
try:
    from ..Struct import Struct
    from ..utils import speed_of_light
except:
    from pybaseutils.Struct import Struct
    from pybaseutils.utils import speed_of_light
# end try

# ========================================================================== #


##### Material 1 - Free space
cc, mu0, eps0 = speed_of_light()    



class qoptics(Struct):

    def __init__(self, **kwargs):       
        self.cc, self.mu0, self.eps0 = speed_of_light()

        self.init(kwargs)
    # end def __init__
        
    def init(self, **kwargs):
        # Calculate the free-space wavelength
        self.freq = kwargs.get('freq', None)
        self.omega = 2.0*_np.pi*self.freq
        self.lambda0 = self.cc/self.freq        

        #  Grid parameters (assume on beam-axis at beam-waist if no input)        
        self.zz = kwargs.get('zz', 0.0)
        self.nz = len(self.zz)
        
        self.rr = kwargs.get('rr', 0.0)
        self.nr = len(self.rr)

        # ====== #
        
        self.Nrefr = kwargs.get('Nrefr', None) # index of refraction in material
        self.kvec = kwargs.get('kvec', None)   # vector of wavenumbers in material
        self.wo = kwargs.get('wo', None)
        self.zr = kwargs.get('zr', None)

                    
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
            

#    def BeamPropagation(self, ABCD):
#
#        qz = _np.zeros((nz,), dtype=complex)        
#        for ii in range(nz):
#            qz[ii] = propagate_beamparameter(self, ABCD, qz[ii-1])
            
    # ====================================================================== #            
    # ===================== Quasi-optical formula ========================== #

    def WaveVector(self):
        if not hasattr(self,'k0'):
            self.k0 = 2.0*_np.pi/self.lambda0
        # endif
        self._kvec = self.Nrefr*self.k0

    def RefractiveIndex(self):
        self._Nrefr = _np.sqrt(self.cc*self.kvec/self.omega)
        
    def BeamWaist(self):
        self._wo = _np.sqrt(2.0*self.zr/self.kvec)

    def RayleighRange(self):
        self._zr = _np.pi*self.Nrefr*self.wo**2.0 / self.lambda0

    def SpotSize(self):
        self.wz = self.wo * _np.sqrt( 1.0 + (self.zz/self.zr)**2.0)

    def RadiusCurvaturePhaseFront(self):
        self.Rz = self.zz*( 1.0 + (self.zr/self.zz)**2.0)
        
    def GouyPhaseShift(self):
        self.gphz = _np.arctan(self.zz/self.zr)

    def ComplexBeamParameter(self):
        self._qz = self.zz + 1j*self.zr

    # ======= #
                
    def Efield_TEM00(self):
        """"
        returns the wave electric field from the fundamental Gaussian beam (TEM00)        
        """
        rr = self.rr.copy()
        rr = _np.atleast_2d(rr)
        if _np.size(rr, axis=1) == 1:   rr = rr.T  # end if
        nr = self.nr
        
        zz = self.zz.copy()
        zz = _np.atleast_2d(zz)        
        if _np.size(zz, axis=0) == 1:   zz = zz.T  # end if
        nz = self.nz
            
        # Update the beam radius, curvature of phase fronts and Guoy phase shift
        self.RaylieghRange()
        self.SpotSize()
        self.RadiusCurvaturePhaseFront()
        self.GouyPhaseShift()        
        
        kvec = self.kvec.copy()
        wz = self.wz.copy()
        Rz = self.Rz.copy()
        gphz = self.gphz.copy()

        # === #
        rr = _np.dot( _np.ones( (nz,1), dtype=float), rr)
        zz = _np.dot( zz, _np.ones( (1, nr), dtype=float))
        
        kvec = _np.dot(kvec.reshape(nz,1), _np.ones((1,nr), dtype=float))
        wz = _np.dot(wz.reshape((nz,1), _np.ones((1, nr), dtype=float)))
        Rz = _np.dot(Rz.reshape((nz,1), _np.ones((1, nr), dtype=float)))
        gphz = _np.dot(gphz.reshape((nz,1), _np.ones((1, nr), dtype=float)))
                
        Exyz = (self.wo/wz) * _np.exp(-rr**2.0 / wz**2.0) 
        Exyz *= _np.exp(-1j*(kvec*rr**2.0)/(2.0*Rz))
        Exyz *= _np.exp(-1j*(kvec*zz-gphz))
        
    # ======= #
                
    def power(self, rr):
        """
        Stores the beam intensity as a function of beam radius and distance        
        """
        if hasattr(self, 'I0'):
            I0 = self.I0
        else:
            I0 = 1.0
        # endif
        rr = _np.atleast_2d(rr)
        
        return I0*(self.wo/self.wz)**2.0 * _np.exp(-2.0*rr**2.0/self.wz**2.0)            


    # ============= Basic stuff ================= #

    def elliptic_mirror(self, r1, r2):
        foc = 1.0/ (1.0/r1 + 1.0/r2)
        return self.thinlens(foc)

    def two_thin_lens(self, f1, d, f2):
        return self.ABCD(-f2/f1, f1+f2, 0.0, -f1/f2)
            
#    def elliptic_mirror_design(self):
        """        
        a = semi-major radius of ellipse
        b = semi-minor radius of ellipse
        ellipsoid:  (x**2.0 + y**2.0)/b**2.0 + z**2.0 / a**2.0 = 1.0

        eccentricity:        e = sqrt(1.0-b**2.0 / a**2.0)

        Distance between any point on ellipse and focii:  P is point on ellipse  
            (1) r1 + r2 = a (r1 is distance to focus 1 from P, r2 is distance to focus 2 from P)

        Focal distance:      
            (2) c = 2*e*a  Linear distance between the two focii, F1 / F2                    
            
        Waist of input / output beam (w1 / w2) are not at F1 / F2, but at distances
            d1 / d2 from a point on the ellipse, P            
           
        Ideal design has radius of curvature of phase fronts equal to the distance 
            from a point on ellipsoidal surface to each focii
            (3) r1 = R(d1) / r2 = R(d2)
            
        Treat ellipsoidal mirror as a thin lens:   
            (4) 1/f = 1/r1 + 1/r2
        
        Ellipse design requires choosing w01/w02, then setting radius of 
            curvature of the phase fronts on the elllipsoidal surface (3)
            
            That determines the focal lengths of the equivalent thin lens (ellisoidal mirror), 
            as well as the semi-major axis of the ellipse. (1) + (4)
            
            The eccentricity (and semi-minor axis of the ellipse) is set by the 
            required distance between the two Focii. (2)
        """
        
    # ============= Built-in geometric optics ray matrices ================= #

    def propagate_ray(self, ABCD, r1, th1):
        # from a ray perspective (r is distance from beam-axis, theta is angle of approach)
        arr = _np.outer(ABCD, _np.array([r1, th1], dtype=float) )
        return arr[0], arr[1]   # r2, th2

    def propagate_beamparameter(self, ABCD, qz):        
        #self.qz = _np.outer(ABCD, _np.asarray([self.qz, 1.0], dtype=complex) )
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
        
    def refraction_flatinterface(self, n1, n2):
        """
        Refraction at a flat interface between two media
            initial media index of refraction - n1
            final media index of refraction - n2
        """
        return self.ABCD(1.0, 0.0, 0.0, n1/n2)
        
    def refraction_curvedinterface(self, R, n1, n2):
        """
        Refraction at a curved interface between two media
            R - radiuso f curvature, R>0 for convex (center of curvature  after interface)
            initial media index of refraction - n1
            final media index of refraction - n2
        """
        return self.ABCD(1.0, 0.0, (n1-n2)/(R*n2), n1/n2)
        
    def reflection_flatmirror(self):
        """
        Ray transfer matrice for reflection from a flat mirror
        
        Only valid for mirrors perpendicular to the ray.
        """
        return self.ABCD(1.0, 0.0, 0.0, 1.0)

    def curvedmirror(self, R, th, plane='horizontal'):
        """
        Ray transfer matrice for reflection from a curved mirror
            R is the radius of curvature (R>0) for concave, valid in the paraxial approximation
            theta in radians, the mirror angle of incidence in the horizontal plane.
        
        """

        if plane.lower == 'horizontal':
            # Effective radius of curvature in the tangential plane (horizontal direction)
            Re = R*_np.cos(th)    
        else:
            # Effective radius of curvature in the sagittal plane (vertical direction)
            Re = R / _np.cos(th)
        # end if

        return self.ABCD(1.0, 0.0, -2.0/Re, 1.0)        

    def thinlens(self, flen):
        """
        ABCD matrix for a Gaussian beam propagation through a thin lens
        Input
            flen - [m], focal length of the thin lens 
        Output 
            ABCD - [m], ABCD matrix for beam propagation through free space a distance d
        """
        return self.ABCD(1.0, 0.0, -1.0/flen, 1.0)
        
    def thicklens(self, n1, n2, R1, R2, t):
        """
        ABCD matrix for a Gaussian beam propagation through a thin lens
        Input
            n1 = refractive index outside of the lens
            n2 = refractive index of the lens itself (inside the lens)
            R1 = Radius of curvature of First surface (interface)
            R2 = Radius of curvature of Second surface (interface)
            t = center thickness of lens
        Output 
            ABCD - [m], ABCD matrix for beam propagation through free space a distance d
        """       
        ABCD = self.outer(self.refraction_curvedinterface(self, R2, n2, n1), 
                                      self.ABCD(1.0, t, 0.0, 1.0))
        return self.outer( ABCD, self.refraction_curvedinterface(self, R1, n1, n2) )
        
    def thinlens_freespace(self, flen, d):
        """        
        propagate into a thin lens with focal length 'flen' 
        then propagate through free space a distance d
        """
        return _np.outer(self.freespace(d), self.thinlens(flen))
        
    def freespace_thinlens_freespace(self, d1, flen, d2):
        """
        propagate a distance 'd1' into a thin lens with focal length 'flen' 
        then propagate through free space a distance d2
        """
        return _np.outer(self.thinlens_freespace(flen, d2), self.ABCD_freespace(d1))

    def rightangleprism(self, th, psi, d, n):
        """
        k is the beam expansian factor, 
            th is the angle of incidence, 
            psi is the angle of refraction, 
        d is the prism path length
        n is the refractive index of prism material
        ... 
        only applies for orthogonal beam exit
        """
        k = _np.cos(psi)/_np.cos(th)
        return self.ABCD(k, d/(n*k), 0.0, 1.0/k)

    # Reference 3:  Chen
    # Circular openings with Equilateral triangular lattice

    # ====================================================================== #
    # ======================= Dichroic filters ============================= #
    
    def dichroic_plate(self, freq, diameter, spacing, thickness, th):           

        wavelength = cc/freq
        radius = 0.5*diameter
        
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

            AA = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness)))
            BB = 1.0/  (1.0 - 1j*(A[ii]+B[ii]*      coth(beta[ii]*thickness))) 
            # Reflection
            R2[ii] = AA.copy() + BB.copy() - 1.0
        
            # Transmission
            T2[ii] = AA.copy() - BB.copy()
            
            # R2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness))) + 1.0/(1.0-1j*(A[ii]+B[ii]*coth(beta[ii]*thickness))) - 1.0
            # T2[ii] = 1.0 / (1.0 - 1j*(A[ii]+B[ii]*cmath.tanh(beta[ii]*thickness))) - 1.0/(1.0-1j*(A[ii]+B[ii]*coth(beta[ii]*thickness)))
            # if self.verbose:    print(_np.abs(R2[ii]), _np.abs(1-T2[ii]))
        
        # For oblique incidence, there is a correction here:
        porosity = _np.pi*(2.0*radius)**2.0 / (2.0*_np.sqrt(3)*spacing**2.0)
        T2perp = T2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
        T2parr = T2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))

        # Does it work the same for reflection?
        # R2_perp = R2*_np.cos(th*_np.pi/180.0)**(2.0*(1.0-porosity))
        # R2_parr = R2*_np.cos(th*_np.pi/180.0)**(1.5*(1.0-porosity))
                
        phperp = _np.arctan(_np.imag(T2perp) /_np.real(T2perp))*180.0/_np.pi - 360.0*thickness*_np.cos(th*_np.pi/180.0)/wavelength  # degrees
        phparr = _np.arctan(_np.imag(T2parr) /_np.real(T2parr))*180.0/_np.pi - 360.0*thickness*_np.cos(th*_np.pi/180.0)/wavelength  # degrees        

        T2perp = _np.abs(T2perp)
        T2parr = _np.abs(T2parr)
        
        R2perp = 1.0-T2perp
        R2parr = 1.0-T2parr

        if self.verbose:
            print("Dichroic plate characteristics: ")
            print("Hexagonal hole pattern: diameter=%2.2f mm, spacing=%2.2f mm, thickness=%2.2f mm"%(1e3*2.0*radius, 1e3*spacing, 1e3*thickness))
            print("filter cut-offs: %3.1f<f<%3.1f GHz"%(fc1, fc2))

        return T2perp, T2parr, porosity, fc1, fc2
        
    def dichroic_notch_response(self, freq=None, diameter=1.6e-3, spacing=2e-3, thickness=5e-3, th=45, plotit=True, wd=None):
        if freq is None and hasattr(self, 'freq'): freq=self.freq # endif
        
        T2perp, T2parr, por, fc1, fc2 = self.dichroic_plate(freq, diameter, spacing, thickness, th)
        
        T2_perp_log, T2_parr_log = [_np.zeros((len(freq),2), dtype=float) for ii in range(2)]
        porosity, fco, fcd = [_np.zeros((2,), dtype=float) for ii in range(3)]
        T2_perp_log[:,0], T2_parr_log[:,0], porosity[0], fco[0], fcd[1] = \
            self.dichroic_plate(0.5*diameter.copy(), spacing.copy(), thickness.copy())
        T2_perp_log[:,1], T2_parr_log[:,1], porosity[1], fco[0], fcd[1] = \
            self.dichroic_plate(0.5*diameter.copy(), spacing.copy(), thickness.copy())        
        
        T2_perp_log = T2_perp_log[:,0]*(1-T2_perp_log[:,1])
        T2_parr_log = T2_parr_log[:,0]*(1-T2_parr_log[:,1])
        
        T2_perp_log = 20*_np.log10(T2_perp_log)
        T2_parr_log = 20*_np.log10(T2_parr_log)
        
        if plotit:
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

        if wd is not None and plotit and _os.path.exists(wd):
            # wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
            hfig.savefig(_os.path.join(wd,'DichroicNotch_%3.1fGHz_%3.1fGHz.png'%(fco[0],fco[1])), dpi=200, transparent=True)            

    def dichroic_plate_response(self, freq=None, diameter=1.6e-3, spacing=2e-3, thickness=5e-3, th=45, plotit=True, wd=None):
        if freq is None and hasattr(self, 'freq'): 
            freq=self.freq 
        elif freq is None:
            freq = 1e9*_np.linspace(100.0, 250.0, 250)            
        # endif
        
        T2perp, T2parr, por, fc1, fc2 = self.dichroic_plate(freq, diameter, spacing, thickness, th)
                
        T2_perp_log = 20*_np.log10(T2perp)
        T2_parr_log = 20*_np.log10(T2parr)


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
        _np.savetxt(filnam, 1e3*centers, fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n%6s 6%s'%('x[mm]', 'y[mm]') )
        
        filnam = _os.path.join(wd,'DichroicPlate_Transmission_%s_%3.1fGHz_d%0.2f_s%0.2f_t%0.1f.txt'%(matname, fco,1e3*D,1e3*S,1e3*l3))
        _np.savetxt(filnam, (freq,T2_parr,T2_perp), fmt='%6.3f', delimiter=' ', newline='\n', header=hdr + '\n %8s %8s %8s'%('freq[GHz]','T2[parr]', 'T2[perp]'))



        
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
            _plt.title(r'Power Transmission Coefficient: f$_{c,o}$<%3.1f, f$_{c,d}$<%3.1f GHz'%(fco,fcd) )
            _plt.axvline(x=fco, linestyle='--', color='k')
            _plt.axvline(x=fcd, linestyle='--', color='k')
            _plt.axhline(y=por_log, linestyle='--', color='k')
            
            _plt.text(x=fco+5, y=-15, s='Hexagonal hole pattern: \n diameter=%2.2f mm, \n spacing=%2.2f mm, \n thickness=%2.1f mm'%(1e3*D, 1e3*S, 1e3*l3))



        if wd is not None and plotit and _os.path.exists(wd):
            # wd = _os.path.join('G://','Workshop','QMB','Documentation','Design','Dichroic Plate')
            hfig.savefig(_os.path.join(wd,'DichroicNotch_%3.1fGHz_%3.1fGHz.png'%(fco[0],fco[1])), dpi=200, transparent=True)            
            
    # ====================================================================== #
    # ===================== Higher order modes ============================= #

    """
    A thin wrapper around numpy's hermite polynomial module
    """
    herm = _np.polynomial.hermite
    def HermPoly(coef, domain=None, window=None):
        pass


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
    # ====================================================================== #            
    # ====================================================================== #            
                
    @property
    def wo(self):
        return self._wo        
    @wo.setter
    def wo(self, value):
        self._wo = _np.asarray(value)
        self.RayleighRange()
    @wo.deleter
    def wo(self):
        del self._wo        

    @property
    def zr(self):
        return self._zr
    @zr.setter
    def zr(self, value):
        self._zr = _np.asarray(value)
        self.BeamWaist()
    @zr.deleter
    def zr(self):
        del self._zr        

    @property
    def qz(self):
        return self._qz
    @qz.setter
    def qz(self, value):
        self._qz = _np.asarray(value)
        self.Rz = _np.real(1.0/self._qz)
        
        tmp = _np.imag(-1.0/self._qz)
        tmp = self.lambda0/self.Nrefr/(_np.pi*tmp)
        self.wz = _np.sqrt(tmp) 
    @qz.deleter
    def qz(self):
        del self._qz        
        
    @property
    def zz(self):
        return self._zz
    @zz.setter
    def zz(self, value):
        self.nz = len(value)
        self._zz = _np.asarray(value)
        if len(self.Nrefr) == 1:  self.Nrefr *= _np.ones( (self.nz,), dtype=complex)
    @zz.deleter
    def zz(self):
        del self._zz                

    @property
    def rr(self):
        return self._rr
    @rr.setter
    def rr(self, value):
        self.nr = len(value)
        self._rr = _np.asarray(value)
    @rr.deleter
    def rr(self):
        del self._rr        
        
    @property
    def Nrefr(self):
        return self._Nrefr        
    @Nrefr.setter
    def Nrefr(self, value):
        self._Nrefr = _np.asarray(value)
        self.WaveVector()
    @Nrefr.deleter
    def Nrefr(self):
        del self._Nrefr        

    @property
    def kvec(self):
        return self._kvec        
    @kvec.setter
    def kvec(self, value):
        self._kvec = _np.asarray(value)
        self.RefractiveIndex()
    @kvec.deleter
    def kvec(self):
        del self._Nrefr        
        
    # ====================================================================== #            
    # ======================= Common Materials ============================= #            


    def _mats(self, material):           
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
        
        return epsr, loss_tangent
    # ====================================================================== #            
    # ====================================================================== #            
        
# end class quasioptics
    
def coth(val):
    return 1.0/cmath.tanh(val)    