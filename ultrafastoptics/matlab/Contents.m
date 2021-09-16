% Ultrafast Optics Toolbox
% Version 1.0.0 07-March-2006
%
% Dispersion computation.
%   stackgd - Compute group delay and reflectivity using exact analytic methods.
%   stackgdfast - Compute group delay using approximate analytic methods.
%   stackgdfastmex - MEX (C) version of stackgdfast (uses different interface).
%   stackgdfastmexwrap - Call MEX version of stackgdfast using same interface.
%   stackgdd - Compute GDD, GD, and R using exact analytic methods.
%   stackgddfast - Compute dispersion using approximate analytic methods.
%   stackgddfastmex - MEX (C) version of stackgddfast (uses different interface).
%   stackgddfastmexwrap - Call MEX version of stackgddfast using same interface.
%   stackgddez - Compute dispersion using simple (but inefficient) interface.
% 
% Demonstrations.
%   dispcompdemo - Demonstrate dispersion computation functions.
%
% Utility.
%   nsell - indices from Sellmeier coefficients.
%   dnsell - derivatives of index with respect to wavelength.
%   dndksell - derivatives of index with respect to wavenumber.
%   makemex - Compile mex files.
%   
% Copyright 2004-2006  Jonathan R. Birge and MIT.
