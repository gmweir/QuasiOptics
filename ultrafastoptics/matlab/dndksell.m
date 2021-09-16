function [n, dnk, ddnk] = dndksell(ks, S)
%DNSELL  Calculate index and derivatives from Sellmeir coefficients
% [n, dn, ddn] = dnsell(y, S) is d^(k) n(k) for k = 0...2 at wavelength y,
% given the Sellmier coefficients given in matrix S, where the first column
% is the A coefficients and the second are the B coefficients.

[n, dny, ddny] = dnsell(2*pi./ks, S);
dnk = -2*pi*dny./ks.^2;
ddnk = 4*pi*dny./ks.^3 + 4*pi^2*ddny./ks.^4;
