function [R, gd, Rgrad, gdgrad] = stackgdgradez(ys, ds, s1, s2, ssub, theta, pol)
% STACKGDDEZ  
%  [R, gd, gdd] = stackgddez(ys, ds, s1, s2, ssub, theta, pol)
%
% All arguments are the same as that for stackgdd, with the exception of
% s1, s2, and ssub, which are sellmeier coefficients for the first, second
% and substrate materials, respectively. They are of the form:
% 
% s = [A1 A2 A3...; B1 B2 B3...], where An are the numerator coefficient and
% Bn are the denomator coefficients, in units of microns squared. Each term
% in the Sellmeier series is given by: A*lambda^2./(lambda^2 - B).
%
% This function is provided for convenience. In optimization, one
% should precompute all of the indices and use stackgdd directly.

ks = 2*pi./ys;

[n1, n1d] = dndksell(ks, s1);
[n2, n2d] = dndksell(ks, s2);
[nsub, nsubd] = dndksell(ks, ssub);

[R, gd, Rgrad, gdgrad] = ...
  stackgdgrad(ks, ds, 1, [n1 n2 nsub], [n1d n2d nsubd], theta, pol);

%[R, gd, Rgrad, gdgrad] = ...
%  stackgdgradmex(ks.', ds.', 1, [n1 n2 nsub].', [n1d n2d nsubd].', theta, pol);
