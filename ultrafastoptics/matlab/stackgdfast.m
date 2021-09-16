function [r, gd] = stackgdfast(ks, ds, n0, ns, dns, theta, pol)
% STACKGDFAST  Calculate complex reflection, group delay and GDD of
% arbitrary dielectric structure using approximate analytic methods,
% including only first order material dispersion. Assumes two indices.
%
% [r, gd, gdd] = stackgdd(ks, ds, n0, ns, dns, ddns, theta, pol)
%
% Input:
% ks: wavenumber column vector (1/um)
% ds: thickness row vector (um)
% ns(length(ks),3): zeroth index matrix, from first layer to substrate
% dns(length(ks),3): first order index
% n0: initial index, assumed to be dispersionless
% theta: incidence angle (radians)
% pol: {'TE','TM'}
%
% Output:
% r(ks): complex reflection coefficient
% gd(ks): group delay (fs)

% Implementation:
%
% T matrices are represented by a matrix where each rows holds the top row
% elements of the T matrix for a given wavenumber.


%%% Initialization and input handling.
c = .2997924580;  % um/fs
if (nargin < 7), pol = 'TE'; end
if (nargin < 6), theta = 0; end
nk = length(ks);  % number of wavelengths
ds = [ds, 0];  % create dummy layers for substrate
n = length(ds);  % number of layers
if (dns == 0)
  dns = zeros(nk, 3);
end

%%% Initialize loop variables.
T = [ones(nk,1), zeros(nk,1)];  % initial matrix
dT = zeros(nk,2);  % initial 1st derivative matrix

%%% Precalculate material variables.
% index variables: [n1, n2, nsub]
% p variables: [p01, p12, p21, p2sub]

% Effective indices.
neffs = sqrt(ns.^2 - (n0*sin(theta))^2);
neff0 = n0*cos(theta);
dneffs = ns.*dns./neffs;  % not really needed for nsub
% Calculate pTE (needed regardless of polarization).
pTEs(:,1) = neff0./neffs(:,1); % n1 from air
pTEs(:,2) = neffs(:,1)./neffs(:,2); % n2 from n1
pTEs(:,3) = 1./pTEs(:,2);  % n1 from n2
pTEs(:,4) = neffs(:,2)./neffs(:,3);  % nsub from n2
% Calculate p = pTM or p = pTE as needed.
if (pol == 'TM')
  % Calculate p = pTM from p0 and pTE.
  % p0 is just p (as above) for normal incidence.
  % p0 = ns(:,L - 1)./nL;
  % p = pTE./p0.^2;
  p0s(:,1) = n0./ns(:,1);
  p0s(:,2) = ns(:,1)./ns(:,2);
  p0s(:,3) = 1./p0s(:,2);
  p0s(:,4) = ns(:,2)./ns(:,3);
  ps = pTEs./p0s.^2;
else  % pol == 'TE'
  ps = pTEs;
end
% Common subexpressions.
dterms = (neffs + [ks,ks,ks].*dneffs);
kneffs = [ks,ks,ks].*neffs;
pps = 1 + ps;
pms = 1 - ps;

%%% Step through structure, calculating all wavelengths in parallel.
for L = 1:n,
  % Initialize local variables.
  d = ds(L);
  id = 1i*d;  % layer thickness times i
    
  % Select appropriate material parameter indices
  if (L == n)  % last layer
    nx = 3; px = 4;
  elseif (L == 1)  % first layer
    nx = 1; px = 1;
  elseif (mod(L,2))  % odd layer, n2->n1
    nx = 1; px = 3;
  else  % even layer, n1->n2
    nx = 2; px = 2;
  end
    
  % Calculate new matrices for full structure up to L.
  % Only three matrix multiplys are required due to the approximation.
  % [98% of the total execution time is spent here.]
  ephi = exp(-id*kneffs(:,nx))/2;  % phasor/2 of layer
  D = -id*dterms(:,nx);
  TL1 = ephi.*pps(:,px); TL2 = ephi.*pms(:,px);  % current layer elements
  T = [TL1.*T(:,1) + TL2.*conj(T(:,2)), TL1.*T(:,2) + TL2.*conj(T(:,1))];  % new T
  dT = [D.*T(:,1) + TL1.*dT(:,1) + TL2.*conj(dT(:,2)), ...
	D.*T(:,2) + TL1.*dT(:,2) + TL2.*conj(dT(:,1))];
end

% Calculate complex reflectance and group delay from total T matrices.
% [This isn't as numerically stable as it could be.]
r = -T(:,2)./T(:,1);  % complex reflection coefficient
dr = (T(:,2).*dT(:,1) - dT(:,2).*T(:,1))./T(:,1).^2;
dphi = (imag(dr).*real(r) - real(dr).*imag(r))./abs(r).^2;
gd = dphi/c;


function c = tmatrixprod(b, a)
% Compute the top row of the products of two T matrices given their top
% rows, as in c = b*a. Parallelized along columns.
c = [a(:,1).*b(:,1) + conj(a(:,2)).*b(:,2), ...
      a(:,2).*b(:,1) + conj(a(:,1)).*b(:,2)];
