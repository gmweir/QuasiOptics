# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 11:55:55 2021

@author: gawe
"""

import numpy as _np
import matplotlib.pyplot as _plt

# from FFT import hilbert

try:
    from . import rectangular_waveguide as _rw
except:
    from QO import rectangular_waveguide as _rw
# end try


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

"""

def transfer_losslessTEM(beta, l, Zo, Yo):
    """
    Transfer function of a lossless transmission line to the TEM mode
    (assuming sinusoidal excitation)

    beta - [rad/m] - propagation constant of the line (guide wavenumber)
    l - [m] - length of the transmission line
    Zo - [Ohms] - characteristic impedance of the ilne
    Yo - [Mohrs] - characteristic admittance of the line

    """
    return _np.matrix([[_np.cos(beta*l), 1j*Zo*_np.sin(beta*l)], [1j*Yo*_np.sin(beta*l), _np.cos(beta*l)]])


def richards_transform():
    """
    s -> j*omega -> j*alpha*tan(theta)
    """
    pass


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


Reactance functions and circuit synthesis (lossless networks):

In a linear, time-invariant system, the relationship between voltage and
current can be expressed in terms of their Laplace transforms:
    V(p) = Z(p)*I(p)

    Z(p) is the input impedance of the network. For a finite lumped network
    (not infinite in extent), Z is a rational function of p and may be
    expressed as the ratio of two polynomials:
            Z(p) = N(p)/D(p)

    For passive networks, since all physical networks give rise to real responses for real input, and conservation of energy applies, Z(p) is real for p real. Re{Z(p)}>0 for Re{p}>0. Then the coefficients of N(p) and D(p) are real, and Z(p) has no poles or zeros in the right half-plane (N,P are Hurwitz polynomials).


=======

    The reflection coefficient in terms of input impedance:
        Gamma(p) = +- (Z(p)-1)/(Z(p)+1); Gamma(p) is real for p real
                                         0 <= |Gamma(p)| <= 1 for Re{p} > 0

    Note that the input impedance can be split into its real and imaginary parts:
        Z = R + j*X
        Gamma = +- (R+j*X - 1) / (R+j*X + 1)
        |Gamma|^2 = 1 - 4R/( X^2 + (1+R)^2 )

A lossless network is comprised entirely of reactive components (no resistance)
    Z(p) |_(p=j*omega) = Z(j*omega) = R(omega) + j*X(omega)
        R(omega) = Re{Z(j*omega)} = 0 (lossless network)
        X(omega) = Im{Z(j*omega)}

 Z(p) = (m1 + n1)/(m2+n2)
    m1, n1 are the even, odd parts of N(p); m2, n2 even/odd parts of D(p)

 Split the impedance into even / odd functions: Z(p) = Ev{Z(p)}+Odd{Z(p)}
  even polynomials -> even powers -> p=j*omega -> Ev{Z(j*omega)} = Real
  odd polynomials -> odd powers -> p = j*omega -> Odd{Z(j*omega)} = Imaginary

In a lossless network:
    R-> 0:   Ev{Z(p)} = 0 = ( Z(p) + Z(-p) ) / 2
        or
        (m1 + n1)/(m2+n2) + (m1 - n1)/(m2 - n2) = 0
    --> m1/n1 = n2/m2

    Z(p) = m1/n2  or  n1/m2
            If the numerator is even, the denominator is odd and visa-versa

    ---> These types of impedances are known as "reactance functions"

Generally:
    Z(p) = Ainf*p + A0/p + sum_1^m ( 2*Ai*p / (p^2+omega_i^2) )

if p = j*omega
    Z(j*omega) = j*X(omega)
       = Ainf*omega - A0/omega + sum_1^m ( 2*Ai*omega / (omega_i^2-omega^2) )

And the derivative in frequency space:
    dX(omega)/domega > 0
    = Ainf + A0/omega^2 + sum_1^m ( 2*Ai*( omega_i^2+omega^2) / (omega_i^2-omega^2)^2 )

---> The impedance of a parallel tuned circuit is identical to the reactance
     function for a lossless network ... this is synthesis

           ____inductor (L)____
   Z(p)   |                    |
    o-----|                    |------o
          |___             ____|
              Capacitor (C)

--> Foster synthesis: circuit is derived by partial fraction decomposition
    (expansion) of the impedance circuit, then cascading parallel tuned
    circuits following a series inductor + capacitor as below:

                     ___L____         ___L____            ___L____
     Z(p)           |        |       |        |          |        |
      o---L---| |-o-|        |---o---|        |--- ...---|        |---o
                    |___ ____|       |___ ____|          |___ ____|
                        C                C                   C

--> Ladder synthesis: circuit is derived by alternating reactance / admittances
    (Lossless two-port network: series inductor, parallel capacitor) and
    subtracting the value from the reactance function until complete

          _____      _____      _____         _____
    o----|    |-----|    |-----|    |--...---|    |---o
         | Z1 |     | Z2 |     | Z3 |        | ZN |
    o----|____|-----|____|-----|____|--...---|____|---o
   ex//

           _______         _____       _____          _____
     Z(p) |       |       |    |      |    |         |    |
      o---| Z1(p) |---o---| Z2 |--o---| Z3 |--o...o--| Zn |---o----o
          |_______|   |   |____|  |   |____|         |____|   |
                     _|_         _|_                         _|_
                    |  |        |  |                        |  |
                    |Y1|        |Y2|                        |Yn|
                    |__|        |__|                        |__|
                     |           |                           |
                     0           0                           0

--> Darlington synthesis: any positive rael function can be synthesized as
    the input impedance of a lossless passive reciprocal two-port network
    which is terminated in al oad resistor (Rload>=0).
               _______
         o----|      |---o
              | Z(p) |   Rload>=0
         o----|______|---o

    ==========================================================
    Example Synthesis using Ladder networks:
        For an impedance in partial fraction form of (and combined after)
            Z(p) = 2*p + 1/p + 2p/(p^2+1)
                 = ( 2*p^4 + 2*p^2 + p^2 + 1 + 2*p^2 ) / (p*(p^2+1))
                 = ( 2*p^4 + 5*p^2 + 1 ) / (p^3+p))

    Evaluate the residue as p->infinity (Ainf above)
        Z(p->infinity) -> 2*p    (from inspection of the first line, it is obvious)
            Z1(p) = 2*p

    Z(p) = 2*p + 1/p + 2p/(p^2+1) = ( 2*p^4 + 5*p^2 + 1 ) / (p*(p^2+1))
         simple pole of order 1 at 0,
         and simple pole at +-j on imaginary axis of order 1

    First order pole: no derivative. Residue at infinity evaluated by
        Res( Z(p), infinity) = lim (p->infinity) 1/p * Z(p)
           = lim (p->infinity) ( 2*p^4 + 5*p^2 + 1 ) / (p^2*(p^2+1))

           complete the square on the numerator
           = lim (p->infinity) ( (2*p^2+1)*(p^2+1) + 2*p^2) / (p^2*(p^2+1))

            factor the removable poles from each term
           = lim (p->infinity) ( (2*p^2+1)/p^2 + 2/(p^2+1) )

            second term limits to 0 as p tends to infinity
            first term need L'hospitals rules to evaluate
           = lim (p->infinity)  (4*p)/(2p) = 2

            therefore
        1/p * Z1(p) = Res( Z(p), infinity) = lim (p->infinity) 1/p * Z(p) = 2

    Subtract off the reactance:
         In a lossless network, this series reactance represents an inductor
         of value L=2

        Z1(p) = Z(p) - 2p
              = ( 2*p^4 + 5*p^2 + 1 ) / (p*(p^2+1)) - 2p
        Use long-division on the first term to make drop the order of the numerator
                                      ____________(2p+0)
                         (p*(p^2+1)) |2*p^4 + 5*p^2 + 1
                                      2*p^4 + 2*p^2
                                      _________________
                                              3*p^2 + 1
                Z(p) = 2p + (3*p^2+1)/(p^3+p)

        Z1(p) = Z(p) - 2p = (3*p^2+1)/(p^3+p)

    Now invert Z1 to form the first admittance Y1(p)
        Y1(p) = p*(p^2+1) / (3*p^2+1)

    Determine the residue as p tends to infinity
       1/p*Y1(p)|_(p=infinity) = lim(p->infinity) (p^2+1) / (3*p^2+1)
           L'Hospitals
           = lim(p->infinity) (2*p) / (6*p) = 1/3

    Therefore, the first admittance is a shunt capacitor of value 1/3
        Subtract it off from the admittance
     1/Z2 = 1/Z1 - Y1(p) = p*(p^2+1) / (3*p^2+1) - p/3
          = (3*p*(p^2+1) - 3*p^3-p)/(3*(3*p^2+1))
          = 2*p/(9*p^2+3) = (2/3) * (p/(3p^2+1))

    Invert to form an impedance, and take the residue as p->infinity
      1/p*Z2|_(p=infinity) = lim (p-> infinity) (1/p) * (3/2p) * (3*p^2+1)
          = lim (p-> infinity) (3/2p^2) * (3*p^2+1)
          = lim (p-> infinity) (9/2 + 3/2p^2)
          = 9/2

    So subtract off a series inductance of 9/2 and invert to form the 2nd admittance:
       1/Y2(p) = Z2(p)- 9*p/2
             = (3/2p) * (3*p^2+1) - 9*p/2
             = (9*p/2 + 3/2p - 9*p/2) = 3/2p

       Y2(p) = 2p/3
    Take the residue as p tends to infinity
        1/p * Y2(p) |_(p=infinity) = lim (p-> infinity) (1/p) * (2*p/3)
              = 2/3
    Subtract off a shunt capacitor of value 2/3

        1/Z3 = Y2(p) - 2*p/3 = 0  ... so we are done

    For an input impedance of
            Z(p) = 2*p + 1/p + 2p/(p^2+1)
                 = ( 2*p^4 + 2*p^2 + p^2 + 1 + 2*p^2 ) / (p*(p^2+1))
                 = ( 2*p^4 + 5*p^2 + 1 ) / (p^3+p))
                      poles at p=0 and +-j
                      zeros at p^2 = -5/4 +- 1/4 * sqrt(25 - 8)
                                   = -5/4 +- sqrt(17/16) = -0.219.. or -2.28..
                               p = +- 0.468...*j and  p = +- 1.510...*j

        The 2-port lossless network, ladder circuit synthesis consists of
            - series inductor:  L=2
            - shunt capacitor:  C=1/3
            - series inductor:  L=9/2
            - shunt capacitor:  C=2/3


    ==========================================================

    A lossless tw-oport network may be decomposed intoa  cascade of first-,
    second- and fourth-degree networks depending on the locations of the zeros
    of the even part of Z(p). These are the transmission zeros of the network
        transmission zeros: Ev{Z(p)}=0
        Zero on jomega-axis corresponds to zero transmission at real omega
            --> measured zero in frequency response

    Lossless network driven by 1-Ohm generator and terminated in 1-Ohm load:
        Zin(p) = (m1+n1) / (m2+n2)
        Ev{Zin} = 0.5*(Z(p) + Z(-p)) = (m1*m2 - n1*n2) / (m2^2 - n2^2)

        Input power:    Pin = |Iin(jomega)|^2 * Re{Zin(jomega)}
                = 0.5*Iin(jomega)*Iin(-jomega)*[ Zin(jomega)+Zin(-jomega) ]
                = Vg*conj(Vg)*[ Zin(p)+Zin(-p) ]/[(1+Zin(p))*(1+Zin(-p))]

   Power transmission zeros at [Zin(p)+Zin(-p)]/[(1+Zin(p))*(1+Zin(-p))] = 0
   ---> zeros of the even part of the input impedance

   ---> Additional zeros at poles of Zin(p) and Zin(-p) (at d.c., infinity or
        on the imaginary axis) can be removed as elements of a reactance
       function by Foster synthesis

   ---> Remaining transmission zeros are not poles of Zin(p) and Zin(-p) and
       can be removed by second-order or fourth-order networks
      Second-order - Brune section: Finite real-frequency transmission zeros
                                    extracted in complex conj. pairs

      Second-order - Darlington C section: Remove transmission zeros on the real-axis

      Fourth-order - Darlington D section: Remove complex transmission zeros


    Cascade Synthesis:  Synthesize an imput impedance as a cascade of
                        Brune, C, and D sections terminated in a resitor
   (assuming any zeros that are also poles of Z(p) and Z(-p) have been removed)
       Zeros occur where m1*m2-n1*n2 = 0.
           Assume it is a perfect square
           m1*m2-n1*n2 = {
                          product_(i=1^q)[1+p^2/omegai^2] *
                          product_(i=1^r)[1-p^2/sigmai^2] *
                          product_(i=1^s)[p^4+2(omegai^2-sigmai^2)*p^2+(omegai^2+sigmai^2)^2]
                          }^2
              ----> three types of transmission zeros
              (i) imaginary axis pair
              (ii) real axis pair
              (iii) complex quadruplet
              (asymmetric zeros come later)

"""

"""
    ==========================================================

Quick refresher:                            Evaluating residues
    A mathematically rigorous way to evaluate the residue

    In a Laurent series:   f(z) = sum_(n=-infty^infty) a_n*(z-zo)^n
        of f(z) about a point z0, a_-1 is the residue of f(z)

        if f(z) is analytic at z0, the residue is 0

        a_-1 = 1/(m-1)! * d^(m-1) / dz^(m-1) [(z-z0)^m * f(z)]_z=z0
            the residue of f(z) as z approaches z0, is related to the (m-1)th
            weighted derivative of f near that point for an order m pole


    simple pole of the function, f, given by c:
        Res(f, c) = lim(z->c) of (z-c)*f(z)

        if f(z) = g(z) / h(z) (holomorphic functions near c)
                Note that holomorphic is a synonym for analytic
        and h(c) = 0 and h'(c) /= 0, then use L'Hospitals rule
        Res(f, c) = lim(z->c) of (z-c)*f(z)
                  = lim(z->c) of (g(z) + z*g'(z) - c*g'(z)) / h'(z)
                  = lim(z->c) of g(z) / h'(z)

    higher order poles of function, f, given by c (order of pole: n)
        Res(f, c) = 1/(n-1)! lim(z->c) of d^(n-1)/dz^(n-1) ( (z-c)^n f(z) )
                                        (n-1)th derivative

    Special case of z-> infinity:

    Res( f(z), infinity) = - Res(1/z^2*f(1/z), 0)

    (i) if lim(z->infinity) f(z)=0
        then
        Res( f(z), infinity) = - lim(z->infinity) z*f(z)

    (ii) if lim(z->infinity) f(z)= c /= 0
        then
        Res( f(z), infinity) = lim( |z|->infinity) z^2*f'(z)

    Special case where part or all of f(z) can be expanded into a
    Taylor or Laurent series: it is easier


    ==========================================================
"""


"""
    ==========================================================
    Scaling the 1-Ohm prototypes to arbitrarty impedance / cutoff frequency

    Most microwave filters operate at 50-Ohms. Historically, 50-Ohms was
    chosen as a compromise between the losses and power handling capacity
    of coax cable.

    To convert from 1-Ohm to Zo-Ohms, scale the impedance of each element
        Inductors:    Z = L*p -->  Zo*L*p = (Zo*L)*p;     L-> Zo*L
        Capacitors:   Z = 1/(C*p) --> Zo/Cp = 1/(C/Zo)p;  C-> C/Zo
        Impedance invertors: K--> Zo*K  (characteristic impedance)

        Ladder coupled filters have series inductors + shunt-capacitors
        Admittance inverter coupled filters have shunt capacitors and impedance inverters

    ==========================================================
    Converting the lowpass prototype to an abritrary cutoff frequency, omega_c

    Lowpass prototype networks have a cut-off of omega = 1 rad/s
        lowpass transmission characteristic:   |S12(jomega)|^2 = 1/(1+Fn(omega)^2)
            omega --> omega/omega_c
            |S12(jomega)|^2 = 1/(1+Fn(omega/omega_c)^2)

            Fn(omega=1) = Fn(omega/omega_c=1)

        Inductors: Z = L*p; Z(jomega)=jomegaL --> jomegaL/omega_c; L-> L/omega_c
        Capacitors: Z=1/C*p; Z(jomega)=-j/omegaC --> -j/(omega/omega_c)C; C-> C/omega_c
        Impedance inverters are frequency independent

    ==========================================================
    Converting the lowpass prototype to a highpass filter with arbitrary omega_c

    |S12(jomega)|^2 = 1/(1+Fn(omega)^2)
        transform omega--> -omega_c/omega to map dc to infinite frequency (and vice versa)
    |S12(jomega)|^2 = 1/(1+Fn(-omega_c/omega)^2)
        inductors: Z(jomega)=jomegaL--> -jomega_c*L/omega
                                   = -j/(omega*(1/omeca_cL))
                                   = -j/(omega*Cprime);
                            Cprime = (1/omega_C)/L

        capacitors: Z(jomega)=-j/omegaC--> -jomega/(omega_c*C) = jomega*Lprime
                            Lprime = 1/(omega_c*C)

        Impedance inverters are unaffected.

        Shifts transmission zeros from omega=infinity to omega=0

    Note:  first convert designed load impedance to match system (50 Ohms),
           then convert to high pass filter

    ==========================================================
    Converting the lowpass prototype to a bandpass filter with arbitrary omega_1,2
        map omega = +-1 in the lowpass to omega_1 and omega_2
        transmission zeros in lowpass --> both omgea=0 and omega=infinity
        midband of the lowpass prototype (omega=0)
                --> center of passband in bandpass: omega_o

    Transformation:
        omega --> alpha*(omega/omega_o - omega_o/omega)
            omega=-1:
                -1 = alpha*(omega_1/omega_o - omega_o/omega_1)
            omega=+1:
                +1 = alpha*(omega_2/omega_o - omega_o/omega_2)

        Solving:
            omega_o = _np.sqrt(omega_1*omega_2)   --- geometric mean
            alpha = omega_o/(omega_2-omega_1)     --- bandwidth scaling factor

    Inductor: Z=jomegaL--> j*alpha*L*(omega/omega_o-omega_o/omega)
                   =j*(alpha*L/omega_o)*omega - j/(omega*(1/(alpha*L*omega_o))
             --> A series connected LC circuit
           o---L---o ==>  o---Lprime--| Cprime |----o
                   Lprime = alpha*L/omega_o
                   Cprime = 1/(alpha*L*omega_o)

    Capacitor: Y=jomegaC--> j*alpha*C*(omega/omega_o-omega_o/omega)
                  = j*(alpha*C/omega_o)*omega - j/(omega*(1/(alpha*C*omega_o)))
            --> admittance of a parallel connected LC circuit


            o---|C|---o ==>   Lprime = 1 / (alpha*C*omega_o)
                              Cprime = alpha*C/omega_o

                           ____Lprime___
                   Z(p)   |            |
                    o-----|            |------o
                          |____Cprime__|


    Impedance inverter:  invariant under the frequency transformation

    ==========================================================
    Realizing an impedance inverter in a circuit:

        Use a pi-network of reactances to simulate an impedance inverter

                o--------- Y=jB ------------o
                      |               |
                    Y=-jB           Y=-jB
                      |               |
                o---------------------------o

       Transfer matrix
               |  1    0  ||  1  -j/B ||  1    0  |     | 0   -j/B |
       [T] =   |          ||          ||          |  =  |          |
               | -jB   1  ||  0   1   || -jB   1  |     | -jB   0  |

               | 0   j/K  |
           =   |          |     where K=-B
               | jK   0   |
      --->  a pi-network of reactance elements equates eactly to an inverter
            of characteristic admittance K=-B.
      ---> Approximate the ideal reactive element jB with a series capacitor
           of admittance Y=jB=jomegaC
                       Then K = -omega*C
                 o - no longer frequency independent, but if the filter is
                     sufficiently narrowband, then it is okay.
                 o - negative sign on capacitance doesn't matter
                     ... just flips the phase response of the filter.
                 o - realize inverters with capacitive pi-sections, and you
                     can use the positive (filter) shunt capacitances to absorb
                     negative capacitances required by the pi-network

    Then the lowpass prototype ladder with impedance inverters and shunt capacitors
           ______       ______       ______         __________
  o---o---|     |--o---|     |--o---|     |--...---|          |--o-----o
      C1  | K12 |  C2  | K23 |  C3  | K34 |        | K(n-1,n) |  CN    Load (1-Ohm)
  o---o---|_____|--o---|_____|--o---|_____|--...---|__________|--o-----o

    becomes a bandpass prototype with impedance inverters and parallel shunt inductors/capacitors
              ______           ______                        __________
  o---o--o---|     |---o--o---|     |---o--o---...---o--o---|          |---o--o--------o
     L1',C1' | K12 |  L2',C2' | K23 |  L3',C3'   L',C'(n-1)'| K(n-1,n) |  Ln',Cn'    Load (1-Ohm)
  o---o--o---|_____|---o--o---|_____|---o--o---...---o--o---|__________|---o--o--------o

        where L1' = 1/(alpha*C1*omega_o),    C1'=(alpha*C1)/omega_o , etc.

  And has a pi-network of capacitors inserted to replace the inverters

              C12          C23                        C(n-1,nn)
  o---o--o---| |---o--o---| |---o--o---...---o--o-------| |---o--o--------o
     L11,C11     L22,C22       L33,C33     L,C(n-1,n-1)      Lnn,Cnn    Load (1-Ohm)
  o---o--o---------o--o---------o--o---...---o--o-------------o--o--------o

        The rth shunt inductor is L_rr = 1/(alpha*C_r*omegao)
        The rth shunt capacitor is C_rr = alpha C_r/omega_o - C_(r-1,r) - C_(r,r+1)
                                            bandpass xform     pi-net1    pi-net2

                     and finally, C_(r,r+1) = K_(r,r+1)/omega_o

    -----> The same could be achieved by inductively coupling the resonators
           (make an admittance inverter from inductors)
           or by alternating inductors/capacitors
    -----> Note that for very small bandwidths, the admittance of the filter
           may need to be scaled by 1/alpha to make the components realizable
           ...
           admittance of the rth shunt resonator prior to forming capacitive inverters:
           --> Y_r(jomega) = j*[alpha*Cr/omega_o)*omega - 1/(omega/(alpha*Cr*omega_o))]
                alpha = omega_o/(delta_omega)
                delta_omega = omega_2 - omega_1
           for delta_omega very small compared to omega_o, alpha is very large

           the inductance of the rth shunt inductor is
               Lrr = 1/(alpha*Cr*omega_o)
               Lrr is very small if alhpa is very large (bad for manufacturing)

         ... the system impedance then needs to be transformed to match the load / generator

     Insert the impedance transformer between the filter and its terminations:
         Y(jomega) = jomegaC_a + 1/(1-j/(omegaC_b))
                   = jomegaC_a + (1.0+j/omegaC_b)/(1+1/(omega^2C_b^2))

    Real part must have denominator equal to 1/alpha at omega_o:
            Re{Y(jomega)} = 1.0/(1.0+1.0/(omega^2*C_b^2))
                force equal to 1/alpha at omega=omega_o (our scaling factor)
            1.0 + 1.0/(omega^2*C_b^2) = alpha
            -->  Cb = 1.0/(omega_o * _np.sqrt(alpha-1))
    Imaginary part must
            Im{Y(jomega)} = omega*Ca + 1/(omega*Cb) / (1+1/(omega^2*Cb^2))
                must be equal to zero at omega=omega_o
            omega*Ca = - 1/(omega*Cb) / (1+1/(omega^2*Cb^2))
                Ca = -sqrt(alpha-1) / (omega_o * alpha)
   ---> Cb is the first and last series capacitor coupling into/out of network
                   Cb = 1.0/(omega_o * _np.sqrt(alpha-1))
   ---> -Ca is absorbed into the capacitance of the first and last resonators
                   Ca = -sqrt(alpha-1) / (omega_o * alpha)

    C01         C12          C23                        C(n-1,nn)       C(nn,nn+1)
  o-| |--o--o---| |---o--o---| |---o--o---...---o--o-------| |---o--o---| |--o
       L11,C11     L22,C22       L33,C33     L,C(n-1,n-1)      Lnn,Cnn    Load (1-Ohm)
  o------o--o---------o--o---------o--o---...---o--o-------------o--o--------o

    C01 = C(n, n+1) = 1.0/(omega_o*sqrt(alpha-1))

    C(r,r+1) = K(r,r+1)/(alpha*omega_o)
    C11 = C1/omega_o - sqrt(alpha-1)/(omega_o*alpha) - C12
    Cnn = Cn/omega_o - sqrt(alpha-1)/(omega_o*alpha) - C(n+1,n)

    Crr = Cr/omega_o - C(r-1,r) - C(r,r+1)    (r=2 ... N-1)
    Lrr = 1/(Cr*omega_o)


        The rth shunt inductor is L_rr = 1/(alpha*C_r*omegao)
        The rth shunt capacitor is C_rr = alpha C_r/omega_o - C_(r-1,r) - C_(r,r+1)
                                            bandpass xform     pi-net1    pi-net2

                     and finally, C_(r,r+1) = K_(r,r+1)/omega_o


"""

# ==========================================================================


class prototypeFilter(object):
    """

    Initialization Inputs
        La    - [dB] - stopband insertion loss
        Lr    - [dB] - passband return loss
        BW    - [MHz] - passband bandwidth (3 dB)
        delta - [MHz] - stopband bandwidth (point at which stopband insertion loss is achieved)

    """
    def __init__(self, La = 3, Lr = 20, BW = 100, delta = None):
        if delta is None:
            delta = 2*BW
        # end if

        pass
    # end def __init__


# end class


# ==========================================================================


class ButterworthLPF(object):
    """
    low-pass prototype network based on a Butterworth filter:
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
            N = ButterworthLPF.degree(La, Lr, S)
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
            self.Krn[nn] = self.Kinverter(self.N, self.Lr, rr)
            self.Crn[nn] = self.Cshunt(self.N, self.Lr, rr)
        # end for
    # end def


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
            N = ChebyshevLPF.degree(La, Lr, S)
        # end if
        self.N = N    # order of filter
        self.La = La  # rejection [dB]
        self.Lr = Lr  # return loss [dB]
        self.S = S    # Shape factor

        self.Krn = _np.zeros((N,), dtype=_np.float64) # Admittance of each stage of the filter
        self.Crn = _np.zeros_like(self.Krn) # Capacitance of each stage of filter
        self._eta = __eta_prototype(N, Lr)

        for nn in range(N):
            rr = nn+1
            #rn = rr+1

            # This sets the admittance and capacitance of each stage to be different.
            # they should be rearranged a bit
            self.Krn[nn] = self.Kinverter(self.N, self.Lr, rr)
            self.Crn[nn] = self.Cshunt(self.N, self.Lr, rr)
        # end for
    # end def


    @staticmethod
    def degree(La, Lr, S, returnint=True):
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
        Nmin = (La + Lr + 6.0)/(20.0*_np.log10(S+_np.sqrt(S*S-1.0))) # eqn 3.71
        if returnint:
            return int(_np.ceil(Nmin))
        else:
            return Nmin
        # end Nmin
    # end def

    @staticmethod
    def __eta_prototype(N, Lr):
        """
        intermediate variable
        """
        epsl = (10.0**(Lr/10.0) - 1.0)**(-0.5)   # formula (2)
        return _np.sinh(_np.arcsinh(1.0/epsl)/N)  # formula (3)

    @staticmethod
    def Kinverter(N, Lr, r):
        """
        Element value K_(r,r+1) for a Chebyshev low-pass prototype filter
            rth admittance in the prototype network
        """
        eta = ChebyshevLPF.__eta_prototype(N, Lr)
        return _np.sqrt( eta*eta + (_np.sin(r*_np.pi/N)**2.0) )/eta   # formula (4)

    @staticmethod
    def Cshunt(N, Lr, r):
        """
        Element value C_Lr for a Chebyshev low-pass prototype filter
            rth shunt-capacitor in the prototype network
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

def Cheb2CombLine_Admittance(alpha, theta, Cshunt):
    """
    inputs
        alpha - [-] - Bandwidth scaling factor
        theta - [radians] - electrical length of the resonators at the
                            center frequency omega0 of the filter
        Cshunt - [Farad?] - Capacitance of the rth capacitor in the prototype network
    outputs
        Yr    - [Siemans] - Characteristic admittance of the short-circuited
                            stub with capacitance Cshunt in a Combline filter
    """
    return alpha*Cshunt*_np.tan(theta)                       # formula (6)


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


def Cheb2CombLine_Capacitance(omega0, alpha, theta, Cshunt):
    """
    inputs
        omega0 - [rad/s] - resonant cyclic frequency of filter
        alpha - [-] - Bandwidth scaling factor
        theta - [radians] - electrical length of the resonators at the
                            center frequency omega0 of the filter
        Cshunt - [Farad?] - Capacitance of the rth capacitor in the prototype network
    outputs
        Cr    - [Farrad] - Equivalent rth capacitance in Combline filter
                        Cr = beta*Yr
    """
    beta = Cheb2CombLine_Beta(omega0, theta)
    Yr = Cheb2CombLine_Admittance(alpha, theta, Cshunt)
    return beta*Yr    # formula (7)

# ==========================================================================
# ==========================================================================



def prototype_LPF(La=30, Lr=3, S=4, N=None):
    # La = 30 # [dB], stopband insertion loss
    # Lr = 3  # [dB], passband return loss

    # 1-Ohm low pass filter based ona  Chebyshev filter characteristic
    chebLPF = ChebyshevLPF(La=La, Lr=Lr, S=S, N=N)



#def stopbad



def example451():
    """
    example from page 114 of "Theory nad Design of Microwave Filters - I. Hunter"

    Design a capacitively coupled Chebyshev bandpass filter to meet the
    following specifications:
        Center frequency:           1 GHz
        Passband bandwidth:         50 MHz
        Passband return loss:       >= 20 dB
        Stopband insertion loss:    >40 dB at fo +-100 MHz
        System impedance:           50 Ohms

    - First design the lowpass filter prototype
    - Next convert to a bandpass filter
    - Then check the value of the inductors / capacitors
    - Scale impedances by 1/alpha to make them reasonable
    - Use capacitive impedance transformer to match system impedance of 50 Ohms


    (i) Chebyshev prototype
     - inverter coupled and
     |S12(jomega)|^2 = 1/(1+eps^2*T_N^2(omega));
         T_N is the Chebyshev polynomial, determined by condition
         for equiripple response (Cheb) s. eqn. 3.54-3.59

        C_N determined so that T_N(omega) is an nth degree polynomial in omega
            T_N(omega) = cos(C_N*theta)
            is zero when C_N*theta = (2*r-1)*pi/2 (r=1,2,...)
                --> theta = (2*r-1)*pi/(2*C_N)  therefore C_N = N
                T_N(omega) = cos(N*arccos(omega))

     Generating function for T_N(omega):
         T_(N+1)(omega) = 2*omega*T_N(omega) - T_(N-1)(omega)
     Initial conditions:
         T_0(omega) = 1 and T_1(omega) = omega

    4th order Chebyshev prototype:
         T_0 = 1
         T_1 = omega
         T_2 = 2*omega*T_1 - T_0 = 2*omega^2-1                          eqn 3.62
         T_3 = 2*omega*T_2 - T_1 = 2*omega*(2*omega^2-1) - omega
                                 = 4*omega^3 - 3*omega                  eqn 3.63
         T_4 * 2*omega*T_3 - T_2 = 2*omega*(4*omega^3 - 3*omega) - 2*omega^2 + 1
                                 = 8*omega^4 - 6*omega^2 - 2*omega^2 + 1
                                 = 8*omega^4 - 8*omega^2 + 1


      |S12(jomega)|^2 = 1/(1+eps^2*T_N^2(omega))
             = 1/(1+eps^2*cos^2(N*arccos(omega)))

            ---> poles at T_N^2(omega) = -1/eps^2
                 cos^2(N*arccos(omega)) = -1/eps^2      eqn 3.74

            solved using
                 eta = sinh( arcsinh(1/eps)/N )     eqn 3.75
            or   1/eps = sinh(N*arcsinh(eta))       eqn 3.76


     |S11(jomega)|^2 = 1 - |S12(jomega)|^2
                     = eps^2*T_N^2(omega)/(1+eps^2*T_N^2(omega))

    """

    def genTN(N, omega):
        TN = [1.0, omega]
        if N>1:
            for nn in range(2,N):
                TN += [2*omega*TN[-1] - TN[-2]]
            # end for
#            return TN
#        else:
#            return TN[:(N+1)]
        # end if
        return TN
    # end def

    def S12_squared(eps, omega, N):
        return 1.0/(1.0+eps*eps*(_np.cos(N*_np.arccos(omega)))**2.0)

    def S11_squared(eps, omega, N):
        return 1.0-S12_squared(eps,omega, N)

    def minN(La, Lr, S):
        """ Return the minimum order required for a Chebyshev filter  """
        return (La+Lr+6.0)/(20.0*_np.log10(S+_np.sqrt(S*S-1.0)))   # eqn 3.71

    def selectivity(delta, BW):
        """ Return the selectivity of the filter: ratio of stopband to passband """
        return delta/BW

    def ripple(Lr):
        """ passband ripple determined from return loss"""
        # eqn. 3.47; passband return loss at ripple level
        return 1.0/_np.sqrt(_np.power(10.0, Lr/10.0) - 1.0)

    def chebyshev_eta(eps, N):
        """ eta parameter for a Chebyshev filter """
        # used in coordinate transformation for filter design
        return _np.sinh( _np.arcsinh(1.0/eps)/N )

    # --- General Nth degree Chebyshev prototype network with series
    # inductors / impedance transformers

    def Kcheb_rr1(r, eta, N):
        """ Characteristic impedance of the inverters pg 63"""
        return _np.sqrt(eta*eta+_np.sin(r*_np.pi/N)**2.0)/eta

    def Lcheb_r(r, eta, N):
        """ inductance of series inductors pg 63"""
        return (2.0/eta)*_np.sin((2.0*r-1.0)*_np.pi/(2.0*N))

    # --- General Nth degree Chebyshev prototype network with shunt
    # capacitors / impedance transformers
    # Note that the formula is identical to that for inductance. This is admittance
    def Ccheb_r(r, eta, N):
        """ capacitors of shunt capacitors """
        #return (2.0/eta)*_np.sin((2.0*r-1.0)*_np.pi/(2.0*N))
        return Lcheb_r(r, eta, N)

    def create_cheb_prototype(eta, N):
        """
        A loop that generates the prototype values of capacitance and
        characteristic impedance of the inverters
        """
        _Klpf = _np.zeros((N-1,), dtype=_np.float64)
        _Clpf = _np.zeros((N,), dtype=_np.float64)
        for ii in range(N):
            if ii+1<N:
                _Klpf[ii] = Kcheb_rr1(ii+1, eta, N)
            # end if
            _Clpf[ii] = Ccheb_r(ii+1, eta, N)
        # end for
        return _Klpf, _Clpf

    def alpha(fc, BW):
        """ bandwidth scaling factor """
        return fc/BW


    # ========================== #
    La = 40 # [dB], stopband insertion loss
    Lr = 20 # [dB], passband return loss

    fc = 1e9       # [Hz], center frequency of filter
    BW = 50e6      # [Hz], passband bandwidth
    delta = 200e6  # [Hz], stopband bandwidth (points at which the insertion loss meets requirements)

    # First evaluate the degree of the lowpass prototype
    # N>= (La+Lr+6)/(20*log10(S+sqrt(S^2-1)))
    #  S = 200/50 = 4
    #  N>=3.68 --> N>= 4
    S = selectivity(delta, BW)     # [-], selectivity of the filter: ratio of stopband to passband
    minN = minN(La, Lr, S)
    N = int(_np.ceil(minN))  # 4

    # Ripple level determined from passband return loss: epsrl=0.1005...
    #    Lr = 10*log10(1+ 1/eps^2)   # eqn. 3.47; passband return loss at ripple level
    epsrl = ripple(Lr)

    # eta parameter of filter
    eta = chebyshev_eta(epsrl, N)  # 0.8201...

    # Scattering parameters:
    # S12_squared = 1.0/(1.0+epsrl*epsrl)        # 0.99; transmission
    # S11_squared = epsrl*epsrl/(1+epsrl*epsrl)  # 0.01; reflection

    # ============================ #
    #
    # plot the response of the general low pass filter
    # ff = _np.linspace(-1, 1.0, num=100, endpoint=True)
    ff = _np.linspace(0, 2.5, num=100, endpoint=True)
    S12 = S12_squared(epsrl, ff, N)
    S11 = S11_squared(epsrl, ff, N)

    _plt.figure()
    _ax1 = _plt.subplot(2,1,1)
    _ax1.plot(ff, 10*_np.log10(S12), '-')
    _ax2 = _plt.subplot(2,1,2, sharex=_ax1)
    _ax2.plot(ff, 10*_np.log10(S11), '-')
#    _ax1.set_xlim((0,2))
    # _ax1.set_xlim((-1,1))
    # _ax1.set_xlim((0, 1))
    _ax1.set_ylabel('|S12| [dB]')
    _ax2.set_ylabel('|S11| [dB]')
    _ax1.set_title('Prototype LPF')

    # ============================ #

    # Create the prototype Chebyshev lowpass filter using impedance inverters
    # and paralell caps
    _Klpf, _Clpf = create_cheb_prototype(eta, N)

    # ------------------ transform to bandpass filter
    # Use the shunt capacitor admittance example
    # (same value as series inductor impedance)
    # scale the lowpass prototype to a bandpass filter
    alpha = alpha(fc, BW)    # bandwidth scaling factor
    omega_o = 2.0*_np.pi*fc  # [rad/s], cyclic frequency at center of bandpass filter


    # Admittance, Y=jomegaC
    #   --> admittance of a parallel connected LC circuit
    #
    #        o---|C|---o ==>   Lprime = 1 / (alpha*C*omega_o)
    #                          Cprime = alpha*C/omega_o
    #
    #                      ____Lprime___
    #              Z(p)   |            |
    #               o-----|            |------o
    #                     |____Cprime__|
    #
    # These are the values of the shunt inductors and shunt capacitors used
    # for the bandpass Chebyshev filter: in parallel to the impedance inverters
    #  --> this would be the intermediary result, but we need to add impedance
    #      transformers at the input andoutput to make the values independent
    #      of bandwidth
    # Lprime = 1/ (alpha*_Clpf*omega_o)
    # Cprime = alpha*_Clpf/omega_o

    def scaled_impedance_inverter_caps(_Klpf, alpha, omega_o):
        """
        return the scaled capacitances necessary to realize an impedance inverter
        """
        return _Klpf/(alpha*omega_o)

    def scaled_bpf_caps(_Clpf, omega_o):
        """
        return the scaled shunt capacitances necessary for this band pass filter
        """
        return _Clpf/omega_o

    def scaled_bpf_inductors(_Clpf, omega_o):
        """
        return the scaled shunt inductors for this band pass filter
        """
        return 1.0/(_Clpf*omega_o)   # Inductors get put in parallel

    def impedance_transformer_caps(alpha, omega_o):
        """
        return the capacitances required make a capacitively coupled
        impedance transformer on the input and output terminations of the filter
        """
        _Ca = 1.0 * _np.sqrt(alpha-1)/(omega_o*alpha)  # negative in application to others
        _Cb = 1.0/(omega_o*_np.sqrt(alpha-1))
        return _Ca, _Cb


    # ===================== #

    def cap_coupled_bpf(N, omega_o, _Clpf, _Ckk, _Cyy, _Ca, _Cb, verbose=1):
        """
        the capcitors (parallel / series by columns) necessary to generate a
        capacitively coupled chebyshev band pass filter
        ... combine these with shunt inductors generated by scaled_bpf_inductors
        """
        Lbpf = scaled_bpf_inductors(_Clpf, omega_o)  # Inductors get put in parallel, but don't change with cap coupling
        Cbpf = _np.zeros((N+1,2), dtype=_np.float64) # rows obvious, columsn parallel and series

        # First and last series capacitors are from the impedance transformers
        Cbpf[0,1] += _Cb
        Cbpf[N,1] += _Cb

        # First and last shunt capacitors absorb the neg. shunt capacitance from the imp. transformers
        Cbpf[1,0] -= _Ca
        Cbpf[N,0] -= _Ca

        # Series capacitors in the middle come from the impedance inverter realization
        Cbpf[1:N  , 1] += _Ckk

        # Shunt capacitors in the middle come from the original filter
        Cbpf[1:N+1, 0] += _Cyy

        # First and last shunt capacitor only are adjacent to one impedance
        # inverter realization, so they absorb neg. capacitance from those individually
        # ... the others in the middle are adjacent to two impedance inverters
        Cbpf[1:N  , 0] -= _Ckk[:]
        Cbpf[2:N+1, 0] -= _Ckk[:]

        if verbose:
            print(Lbpf)
            print(Cbpf)
        # end if
        return Lbpf, Cbpf

    def scale_impedances_to_system(Lbpf, Cbpf, Zo=50.0):
        """
        scale the components to match the required system impedance (input and output)
        """
        # Zo = 50.0  # [Ohms]

        # Inductors:    Z = L*p -->  Zo*L*p = (Zo*L)*p;     L-> Zo*L
        Lbpf *= Zo
        # Capacitors:   Z = 1/(C*p) --> Zo/Cp = 1/(C/Zo)p;  C-> C/Zo
        Cbpf /= Zo

        return Lbpf, Cbpf

    # ===================== #

    # Realize each impedance inverter by a pi-network connected set of capacitors
    # 4th order network with capacitive coupling impedance transformers on
    # input / output
    _Ca, _Cb = impedance_transformer_caps(alpha, omega_o)
    _Ckk = scaled_impedance_inverter_caps(_Klpf, alpha, omega_o)
    _Cyy = scaled_bpf_caps(_Clpf, omega_o)

    Lbpf, Cbpf = cap_coupled_bpf(N, omega_o, _Clpf, _Ckk, _Cyy, _Ca, _Cb, verbose=1)

    # Now scale the impedances to match 50 ohms at the input and output
    Lbpf, Cbpf = scale_impedances_to_system(Lbpf, Cbpf, Zo=50.0)

    return Lbpf, Cbpf
# end def

def cot(x):
    """
    cotangent
    """
    return 1.0/_np.tan(x)


def physical_filter():
    C1 = C4 = 15.9894e-12 # pF
    C2 = C3 = 38.1670e-12 # pF

    Y1 = Y4 = 9.1327 # mhor
    Y2 = Y3 = 22.0488 # mhor

    K12 = K34 = 0.026
    K23 = 0.031

    fo = 2.5e9
    l = 12e-3
    b = 24e-3
    d = 8e-3

    # d = 2*d

    for Kij in [K12, K23, K34]:
        Sij = Scomb_ij(fo, l, b, d, Kij)
        print(Sij)
        # 0.07346126821028004
        # 0.07212307820975238
        # 0.07346126821028004

    # end for

    for Cij in [C1, C2]:
        # Mg = Mgap(d, Cij)
        # print(Mg)
        # -0.002130268362364526
        # -0.002130268588630083

        Mg = Mgap(d, 1e12*Cij)
        print(Mg)
    # end for
# end def




def loading_capacitance(theta0, Zs, Zr, fc):
    """
    Loading capacitance: Cl

    Zs - system impedance (50 Ohms)
    Zr - Resonator impedance (equation 11)
    theta0 - electrical length --> must be less than 90 degrees (less than pi/4)
    fc - center frequency
    """
    return Zs*_np.tan(theta0)/(Zs*Zr*2.0*_np.pi*fc)

def inverse_loadCap(Cl, Zs, Zr, fc):
    """
    return the electrical length based on the loading capacitance
    """
    return _np.arctan( Cl*(Zs*Zr*2.0*_np.pi*fc)/Zs )

def ewall(b,d):
    """
    return the distance to the wall from the resonator edge that matches the cavity diameter

    equation 12
    """
    return 0.5*(b+d)

def Mgap(d, Cij):
    """
    returns the resonator gap between the lid and resonator to provide the necessary capacitance:
        given the capacitance and resonator diameter

    equation 13
    """
    return 0.695*d*d/(100*Cij-2.61*d)


def Scomb_ij(fo, l, b, d, Kij):
    """
    Return the distance between each resonaor (i to j) based on
    Kij - the admittance invertor value,
    electrical length (frequency and resonator length, l)
    b - cavity diameter
    d - resonator diameter

    equation 14
    """
    ftheta = f_theta(theta(fo, l))
    Sij = (b/1.37) * ((0.91*b/d) + 0.048 - _np.log10( 4*ftheta*Kij/_np.pi ))
    return Sij

def theta(fo, l):
    """
    equation 16
    electrical length: ... should probably use guide wavelength here
    """
    wavelength = 3e8 / fo  # free space wavelength
    return 2.0*_np.pi*l/wavelength

def f_theta(theta):
    """ equation 15 """
    return 0.5*(1.0+2*theta/_np.sin(2*theta))



if __name__ == '__main__':
    # Lbpf, Cbpf = example451()

    physical_filter()
# end if