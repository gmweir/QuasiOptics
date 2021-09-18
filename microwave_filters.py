# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 11:55:55 2021

@author: gawe
"""

import numpy as _np
import matplotlib.pyplot as _plt

from FFT import hilbert

"""

General filter terminology:
    Low Pass Filter - allows frequency below a given frequency to pass
                      (to be transmitted or received) while rejecting
                      frequencies above the given frequency.

    High pass filter - allows frequencies above a given frequency to pass
                       through the filter, while rejecting frequencies
                       above the given frequency (opposite of low pass).

    Band pass filter - Passes signal between a lower and an upper frequency
                       while rejecting signal at frequencies outside of
                       that range

    Band stop (or Band Reject or Notch) - Prevent all signal between an
                       upper and lower frequency from passing while
                       allowing all others to pass (opposite of band pass)

    Attenuation - measured in dB, degree by which a signal sees a loss in
                  amplitude after passing through a filter

    Bandwidth - width of the passband of a bandpass filter and is expressed
                  as the frequency difference between lower and upper 3 dB
                  points

    Cutoff - Usually the point at which the response of the filter has
             fallen by 3 dB from passband level

    Group delay - Group delay is a measure of how different components of a
                  modulated signal (which is a sum of sine waves at various
                  frequencies) would propagate through the filter. Measured
                  in units of time (seconds) and is a derivative of the
                  filter's phase with respect to frequency

    Passband - the portion of the frequency spectrum that the filter allows
               to be transmitted

    Stopband -  the band (frequency) where the filter has reached the
                required out-of-band rejection

    Shape factor - The ratio of a filters stop band to pass band. The
                   higher the shape factor, typically the closer the filter
                   is to theoretical performance.

    Selectivity - Measurement of the capability of the filter to pass or
                  reject specific frequencies relative to the center
                  frequency of the filter. Selectivity is typically stated
                  as the loss through a filter that occurs at some
                  specified distance from the center frequency. A filter
                  with high selectivity exhibits high slope in the
                  transition from pass to stop - selectivity is crucial in
                  environments where adjacent channels are close together
                  and high selectivity enables designers to make good use
                  of available bandwidth.

    Temperature stability - Concerns how the temperatures performance
                            varies with temperature. An approach is to
                            define in ppm/degree C the shift of the
                            filters, cutoff, passband, etc. in frequency
                            as temperature varies.

    Insertion loss - ratio of a signal level in a test configuration
                     without a filter present (|V1|) to that when the
                     filter is present (|V2|). When discussed this is
                     typically referencing the loss in the passband.

    Ripple - Ripple is a measure of the variation of insertion loss within
             the passband and is measured in dB

    S11 - Scattering parameter that represents the reflection coefficient
          (Gamma) at the input. Related to Return Loss.
                  RL [dB] = -20*log10|S11|

    S21 - Scattering parameter that represents the loss of signal while
          passing from input to output ports. When the measurement ports
          are at the same impedance, it is a measure of the insertion loss.
                  IL [dB] = -20*log10|S21|

    VSWR - A measure of the filters match to a given impedance (ex// to a
           50 Ohm system), Voltage Standing Wave Ratio is calculated from
           S11 (Gamma).
                       VSWR = (1+|Gamma|)/(1-|Gamma|)

    Return loss - Relative amount of power reflected by an input signal:
                  Measure of the amount of signal that is returned or
                  reflected by the filter. measured in dB, it is the
                  negativ eof hte magnitude of the reflection coefficient
                  expressed as power. Return loss is expressed as a
                  positive number.

                  RL [dB] = -20*log10|Gamma|
                      For example:
                        for -3 dB of return loss
                         3 dB = (- 20 dB ) * log10 | Gamma |
                         | Gamma | = 10^(-3/20) = 0.7079 --> 70% reflection

                        for 10% reflection
                         -20*log10(0.10) = 20 dB --> "20 dB return loss"
                                                  or "-20 dB reflected"

                Note:  Return loss is mathematically positive, but
                colloquially negative... "19 dB of loss"
                                      or "Return loss of 19 dB"

    Passband return loss - return loss in the filters passband

    Percent bandwidth - Common relative figure of merit that compares
                        bandwidth with carrier frequency.

                        normalized measure of how much frequency variation
                        a component / system can handle.

                        Commonly calculated as (3 dBW)/(Center Frequency)

                        BW% = BW/Fc
                            BW is the absolute bandwidth
                            Fc is the center frequency
                        OR

                        BW%=2(F_H/F_L - 1)/(F_H/F_L+1)
                            F_H - Upper frequency of the passband
                            F_L - Lower frequency of the passband

                        1 filter with 1 GHz passband centered at 10 GHz
                        has 10% bandwidth

    Q Factor - The quality factor (Q) of a resonator is expressed as the
               ratio of stored versus lost energy per oscillation cycle.
               Overall losses through a resonator increase as Q factor
               drops and will increase more rapidly with frequency for
               lower values of resonator Q. As a result, the edges of the
               passband become more rounded and badnwidth narrows as the
               Q-factor decreases.

    Rejection - Attenuation of signals outside the passband. Typically
                measured in dB or dBc if referenced from insertion loss of
                passband

General filter technologies:

    Crystal filters - make use of a quartz crystal as the resonant element.
                      The high Q of a quartz resonator makes for a very
                      steep band-pass. These filters are usually
                      implemented at IF frequencies  in the range of 10 MHz
                      and Q factors fall in the range of 10e3 to 100e3

    SAQ/BAQ (Surface Acoustic Wave and Bulk Acoustic Wave) -
                  Acoustic filters cover a range of frequencies up to 6 GHz
                  and offer a good performance/cost tradeoff, making them
                  the dominant off chip filter approach in mobile devices.

    Ceramic Filters - Cover a range of ~100 MHz to ~8 GHz. Offer similar
                      performance to discrete lumped element inductor-
                      capacitor (LC) designs but can be implemented in
                      small form factor surface mount packages. Performance
                      and package thickness can be a limiting factor when
                      comparing ceramic filters with SAW/BAW

    Lumped Element - Discrete LC approaches provide a low-cost approach to
                     implement a filter, but the attainable Q factors are
                     limited in such devices. Discrete lumped element
                     filters are usually used around 30 MHz to 300 MHz
                     range but can in principle be built for applications
                     up to 40 GHz. At mmWave frequencies though discrete
                     lumped element filters are very hard to implement
                     because of the dimensional limitations imposed by the
                     frequency, since the filter elements must be much
                     smaller than the wavelength of the transmission lines.
                     Discrete LC designs are perforamnce and repeatability
                     limited by the tolerances of the discrete components.

    Cavity Filters - Cavity filters are a common approach in the 40 MHz to
                     960 MHz frequency range and can offer high selectivity
                     under high power. They can achieve good performance,
                     but are physically large, and usualyl only seen in
                     infrastructure applications, such as for additional
                     filtering at a cell site.

    Planar filters - Such as Microstrip filters - are manufactured using a
                     thin-film process, and depending on the filter
                     topology can offer high-Q and a reasonable approach to
                     achieving performanc ein a small footprint when
                     comapred with discrete lumped element designs. In a
                     thin film Lumped Element approach, the filter's
                     transmission lines are printed in various
                     configurations, depending on the frequired perofrmance
                     and filter lements are realized through discrete
                     resistive, capactivie, and inductive elements. Planar
                     distributed Element Filters rely on carefully
                     distributed transmission lines to create resonant
                     structures and can be designed to tighter tolerances
                     than a lumped element filter. Distributed Element
                     designs are more practical than Lumped Element designs
                     at increased frequencies.

    Waveguide Filters - Characterized by high power handling capability,
                        leading to their wide adoption in radar
                        applications, high selectivity and rejection and
                        low loss givne that the waveguide itself is a low
                        loss medium.




---> Rule of thumb and fundamental issue:

 Passband loss is inversely proportional to filter bandwidth.
  ->Very narrow bands require very high Q-factors to achieve low passband loss.
      Air filled waveguide resonators enable Q factors from 5 to 20e3
      Dielectric resonators as the resonant elements can push Q up to 50e3


"""

def percBW(freqc):
    """
    Percent bandwidth (BW) - [dBW]
    """
    return

def ReturnLoss(S11):
    """
    mathematically it is the negative of the magnitude of the reflection
    coefficient expressed as power (dB)
        Lr = -20*log10|Gamma|
        |Gamma| = S11 = (VSWR+1)/(VSWR-1)
            Voltage Standing Wave Ratio -- VSWR
    in a perfect transmission line, the return loss is infinite
    """
    return -20*_np.log10(_np.abs(S11))


def VSWR(S11):
    return (1.0+_np.abs(S11))/((1.0-_np.abs(S11)))


def InsertionLoss(S21):
    """
    Relative amount of power absorbed/dissipated while passing through the filter
    """
    return -20.0*_np.log10(_np.abs(S21))



# ==========================================================================
"""
Impedance inverters

    An alternative realization of ladders networks
    --> lossless, reciprocal, frequency-independent, two-port network
           ______
    o-----|     |------
          |     |      |
    Zin   |  K  |      Zload
          |     |      |
    o-----|_____|-------

        Zin(p) = (A*Zload + B)/(C*Zload+D)
               = jK/(j*Zload/K) = K^2/Zload
    ---> K=1 ---> Zin = 1/Zload



Maximally flat prototype (Butterworth approximation)
    |S12(1j*omega)|^2 = 1/( 1+omega^(2N) )    (3.3, page 49)


"""

def impedance_inverter(Kimp):
    """
    K is the characteristic impedance or admittance of the inverter
    """
    return _np.matrix([[0, 1j*Kimp],[1j/Kimp, 0]])

# ==========================================================================

"""
ideal lowpass filter
    magnitude of gain is unity in passband, zero in the stopband
       |H(jomega)| = 1 |omega|<|omegac|
       |H(jomega)| = 0 |omega|>|omegac|
    phase response is linear in the passband:
        psi(omega)=k*omega
    group delay is the derivative
        Tg = -d[psi(omega)]/[domega] = k
    --> linear phase --> constant group delay
        --> zero phase distortion for finite bandwidth signals

    impulse response
       h(t) = (1/pi) * (sin[(t-k)*omegac])/(t-k)
       for omegac = 1
          h(t)_(omegac=1) = (1/pi) = sinc(t-k)
            --> zeros at t-k = m*pi,  m=+-1, +-2, etc.
               t= m*pi+k, peak of sync function at +k
            --> noncausal... physically unrealizable
             (infinite selectivity-->infinite group delay-->infinite elements)

      we make filters by truncating the impulse response and tailoring the
      frequency spectrum so the selectivity isn't too horrendous


Lumped element filter with transfer function
    S12(p) = N(p) / D(p)
        p - complex frequency variable (nu in most signals and systems books)

--> Minimum phase network if there are no poles or zeros in the right half p plane
        N(p) /= 0,   D(p) /= 0,   Re[p]>0
            Defining N(p) and D(p) == Hurwitz polynomials
       when energy transfer between input and output can only take one path,
       a physical system is a minimum phase system.

   Transfer function of minimum phase network:
       H(jomega) = exp[-alpha(omega) - j*psi(omega)]
         alpha(omega) - magnitude
         psi(omega)   - phase
       --> magnitude and phase are a Hilbert transform pair
         psi(omega) = (omega/pi) * int_-inf^inf[ alpha(y)/(y^2-omega^2) dy ]
         alpha(omega) = alpha(0) + (omega^2/pi) * int_-inf^inf[ psi(y)/(y*(y^2-omega^2)) dy ]


Approximation for the amplitude:
    |H(jomega)| = 1.0 / An(omega^2);
        An(omega^2) is a polynomial degree N in omega^2
    --> "All pole" transfer function (this class of ladder network)
    Equiripple characteristics provides optimum selectivity for a given degree


In a linear, time-invariant system, the relationship between voltage and
current can be expressed in terms of their Laplace transforms:
    V(p) = Z(p)*I(p)

    Z(p) is the input impedance of the network. For a finite lumped network
    (not infinite in extent), Z is a rational function of p and may be
    expressed as the ratio of two polynomials:
            Z(p) = N(p)/D(p)
"""



# ==========================================================================


class ChebyshevLPF(object):
    """
    low-pass prototype network based on a Chebyshev filter:
        two-port lumped element network
        with angular cutoff = 1rad/s
        and operating in a 1-Ohm system

    The simplest 1 degree low-pass filter is a series resistor + parallel
    capacitor. An N-degree low-pass filter just cascades these lumped elements.
    In this class the series resistor is modeled as a parallel admittance.

    Initialization Inputs
        N - [-] - Order of the prototype Chebyshev filter
        La - [dB] - stopband insertion loss
        Lr - [dB] - passband return loss
        S  - [-]  - Shape factor (ratio of stopband to passband frequencies)
    """

    def __init__(self, La = 3, Lr = 20, S = 1.2, N= None):
        if N is None:
            N = ChebyshevLPF.degree_Chebyshev(La, Lr, S)
        # end if
        self.N = N    # order of filter
        self.La = La  # rejection [dB]
        self.Lr = Lr  # return loss [dB]
        self.S = S    # Shape factor

        self.Krn = _np.zeros((N,), dtype=_np.float64) # Admittance of each stage of the filter
        self.Crn = _np.zeros_like(self.Krn) # Capacitance of each stage of filter

        for nn in range(N):
            rr = nn+1
            #rn = rr+1

            # This sets the admittance and capacitance of each stage to be different.
            # they should be rearranged a bit
            self.Krn[nn] = self.Kcheb(self.N, self.Lr, rr)
            self.Crn[nn] = self.Ccheb(self.N, self.Lr, rr)
        # end for
    # end def


    @staticmethod
    def degree_Chebyshev(La, Lr, S):
        """
        Formula for calculating the minimum degree of a Chebyshev low-pass filter

        Inputs
            La - [dB] - stopband insertion loss
            Lr - [dB] - passband return loss
            S  - [-]  - Shape factor (ratio of stopband to passband frequencies)

        Outputs
            Nmin - minimum degree of the Chebyshev filter
                N >= (La + Lr + 6)/(20*log10[S+sqrt(S^2-1)])
        """
        return (La + Lr + 6.0)/(20.0*_np.log10(S+_np.sqrt(S*S-1.0))) # formula (1)

    @staticmethod
    def __eta_prototype(N, Lr):
        """
        intermediate variable
        """
        epsl = (10.0**(Lr/10.0) - 1.0)**(-0.5)   # formula (2)
        return _np.sinh(_np.arcsinh(1.0/epsl)/N)  # formula (3)

    @staticmethod
    def Kcheb(N, Lr, r):
        """
        Element value K_(r,r+1) for a Chebyshev low-pass prototype filter
            rth admittance in the prototype network
        """
        eta = ChebyshevLPF.__eta_prototype(N, Lr)
        return _np.sqrt( eta*eta + (_np.sin(r*_np.pi/N)**2.0) )/eta   # formula (4)

    @staticmethod
    def Ccheb(N, Lr, r):
        """
        Element value C_Lr for a Chebyshev low-pass prototype filter
            rth capacitor in the prototype network
        """
        eta = ChebyshevLPF.__eta_prototype(N, Lr)
        return (2.0/eta)*_np.sin((2.0*r-1.0)*_np.pi/(2.0*N))   # formula (5)
# end class

# ==========================================================================

"""
The design of this tuneable cavity filter loosely is based on the paper:
    Design of Low-Loss Coaxial Cavity Bandpass Filter with Post-Manufacturing
    Tuning Capabilities
    2012 IEEE Symposium on Business, Engineering and Industrial Applications
    by Z. Zakaria, A. Sabah, and W. Y. Sam

-----

The majority of the math / physics comes from:
    Theory and Design of Microwave Filters, published 2001
    (IEE electromagnetic waves series; no.48), by Ian. Hunter.
    (ISBN: 0852967772 (or 978-0-85296-777-5)
"""

def Cheb2CombLine_Admittance(alpha, theta, Ccheb):
    """
    inputs
        alpha - [-] - Bandwidth scaling factor
        theta - [radians] - electrical length of the resonators at the
                            center frequency omega0 of the filter
        Ccheb - [Farad?] - Capacitance of the rth capacitor in the prototype network
    outputs
        Yr    - [Siemans] - Characteristic admittance of the short-circuited
                            stub with capacitance Ccheb in a Combline filter
    """
    return alpha*Ccheb*_np.tan(theta)                       # formula (6)


def Cheb2CombLine_Beta(omega0, theta):
    """
    inputs
        omega0 - [rad/s] - resonant cyclic frequency of filter
        theta - [radians] - electrical length of the resonators at the
                            center frequency omega0 of the filter
    outputs
        beta - [s/rad] - scaling parameter between rth admittance / capacitance
    """
    return 1.0/(omega0*_np.tan(theta))   # formula (8)


def Cheb2CombLine_Capacitance(omega0, alpha, theta, Ccheb):
    """
    inputs
        omega0 - [rad/s] - resonant cyclic frequency of filter
        alpha - [-] - Bandwidth scaling factor
        theta - [radians] - electrical length of the resonators at the
                            center frequency omega0 of the filter
        Ccheb - [Farad?] - Capacitance of the rth capacitor in the prototype network
    outputs
        Cr    - [Farrad] - Equivalent rth capacitance in Combline filter
                        Cr = beta*Yr
    """
    beta = Cheb2CombLine_Beta(omega0, theta)
    Yr = Chem2CombLine_Admittance(alpha, theta, Ccheb)
    return beta*Yr    # formula (7)

# ==========================================================================
# ==========================================================================



def prototype_LPF():
    La = 30 # [dB], stopband insertion loss
    Lr = 3  # [dB], passband return loss


#def stopbad