function [T, R] = stacktrez(ys, ds, s1, s2, ssub, theta, pol)
% STACKTREZ  
%  [T, R] = stacktrez(ys, ds, s1, s2, ssub, theta, pol)
%
% All arguments are the same as that for stacktr, with the exception of
% s1, s2, and ssub, which are sellmeier coefficients for the first, second
% and substrate materials, respectively. They are of the form:
% 
% s = [A1 A2 A3...; B1 B2 B3...], where An are the numerator coefficient and
% Bn are the denomator coefficients, in units of microns squared. Each term
% in the Sellmeier series is given by: A*lambda^2./(lambda^2 - B).
%
% This function is provided for convenience. In optimization, one
% should precompute all of the indices and use stackgdd directly.

n = length(ds);
np = floor(n/2);
ks = 2*pi./ys;

n1 = nsell(ys, s1);
n2 = nsell(ys, s2);
nsub = nsell(ys, ssub);

if mod(n,2)  % odd
  ns = [repmat([n1, n2], 1, np), n1, nsub];
else
  ns = [repmat([n1, n2], 1, np), nsub];
end

[T, R] = stacktr(ks, ds, 1, ns, theta, pol);
