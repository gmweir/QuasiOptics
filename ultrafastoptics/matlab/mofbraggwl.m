function ms = mofbraggwl(ybs, Ds, y0, kappa)
% MOFKBCALC  Calculate m values for each Bragg wavenumber.
%  Given a vector of Bragg wavelengths and a vector of dispersion
%  Taylor coefficents (D2...) about k0, returns a vector of normalized stack
%  locations with m(max(kbs)) = 0. Kappa is dimensionless coupling
%  coefficient. The code is in no way optimized, nor is any optimization
%  really needed, however. This implements formula (3.79) in Matuschek.
%
%  ybs: nm, Ds: fs, y0: nm, kappa: 1.

% Initialization.
c = 0.2997924580;  % um/fs
n = length(ybs);
vn = length(Ds);
kbs = 2*pi*1000./ybs;  % Bragg wavenumber in 1/um
kbmax = max(kbs);
k0 = 2*pi*1000/y0;
ms = zeros(1, n);

% Calculate m over the expanded Taylor series.
% (Implicitly parallel over kb vector.)
for v = 0:(vn-1),
  for u = 0:v,
    ms = ms + ...
      1/(factorial(v)*(u+2))*nchoosek(v,u)*Ds(v+1)*c^(v+2)*(-k0)^(v-u) * ...
      (1-kappa/pi)^(u+1)*(kbmax^(u+2) - kbs.^(u+2));
  end
end
ms = ms/2/pi;