function ms = mofkb(kbs, Ds, k0, kappa)
% MOFKBCALC  Calculate m values for each Bragg wavenumber.
%  Given a vector of Bragg wavenumbers and a vector of dispersion
%  Taylor coefficents (D2...) about k0, returns a vector of normalized stack
%  locations with m(max(kbs)) = 0. Kappa is dimensionless coupling
%  coefficient. The code is in no way optimized, nor is any optimization
%  really needed, however. This implements formula (3.79) in Matuschek.
%
%  kbs: 1/um, Ds: fs, k0: 1/um, kappa: 1.

% Initialization.
n = length(kbs);
vn = length(Ds);
kbmax = max(kbs);
ms = zeros(1, n);
c = 0.2997924580;  % um/fs

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