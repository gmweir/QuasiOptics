function gdd = gddsell(y, L, S)
%GDDSELL  Calculate GDD of material from Sellmeir coefficients
% gdd = dnsell(y, L, S)
% The base units are um and fs.

% Initialization.
c = 0.299792458;  % um/fs

% Compute index derivatives wrt wavelength.
[n, dn, ddn] = dnsell(y, S);

% Compute GDD.
gdd = L*y.^3.*ddn/(2*c^2*pi);