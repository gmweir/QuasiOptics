function [R, gd, gdd] = stackgddez(ys, ds, s1, s2, ssub, theta, pol, approx)
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

if nargin < 8
  approx = 'none';
end

n = length(ds);
np = floor(n/2);
ks = 2*pi./ys;

[n1, n1d, n1dd] = dndksell(ks, s1);
[n2, n2d, n2dd] = dndksell(ks, s2);
[nsub, nsubd, nsubdd] = dndksell(ks, ssub);

if mod(n,2)  % odd
  ns = [repmat([n1, n2], 1, np), n1, nsub];
  dns = [repmat([n1d, n2d], 1, np), n1d, nsubd];
  ddns = [repmat([n1dd, n2dd], 1, np), n1dd, nsubdd];
else
  ns = [repmat([n1, n2], 1, np), nsub];
  dns = [repmat([n1d, n2d], 1, np), nsubd];
  ddns = [repmat([n1dd, n2dd], 1, np), nsubdd];
end

[R, gd, gdd] = ...
  stackgdd(ks, ds, 1, ns, dns, ddns, theta, pol, approx);
