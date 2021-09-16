% STACKGDFASTMEX  Calculate complex reflection, group delay and GDD of
% arbitrary dielectric structure using approximate analytic methods,
% including only first order material dispersion. Assumes two indices.
%
% NB: The variable conventions used in this file are transposed from that
% used in the MATLAB function stackgdfast. (This is because the MATLAB
% version is parallelized over wavenumber, whereas the MEX version does
% each wavenumber independently. If you want to use the MEX function with
% the same arguments as for the MATLAB version, use stackgdfastmexwrap.)
%
% [r, gd, gdd] = stackgdfastmex(ks, ds, n0, ns, dns, theta, pol)
%
% Input:
% ks: wavenumber row vector (1/um)
% ds: thickness column vector (um)
% ns(3,length(ks)): zeroth index matrix, from first layer to substrate
% dns(3,length(ks)): first order index
% n0: initial index, assumed to be dispersionless
% theta: incidence angle (radians)
% pol: {'TE','TM'}
%
% Output:
% r(ks): complex reflection coefficient
% gd(ks): group delay (fs)
