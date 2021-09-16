function gd = gdsell(y, L, S)
%GDDSELL  Calculate GD of material from Sellmeir coefficients
% gd = gdsell(y, L, S)
% The base units are um and fs.

% Initialization.
c = 0.299792458;  % um/fs

% Compute index derivatives wrt wavelength.
[n, dn, ddn] = dnsell(y, S); %#ok<NASGU>

% Compute GD.
gd = L*(n - y.*dn)/c;
