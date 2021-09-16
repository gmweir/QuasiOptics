function [r2s, gds, r2grads, gdgrads] = stackgdgrad(ks, ds, n0, ns, dns, theta, pol, approx)
% STACKGDGRAD  Calculate complex reflection, group delay and gradients
% thereof using exact analytic methods. Assumes two indices.
%
% [r, gd, gdd] = stackgdgrad(k, ds, n0, ns, dns, ddns, theta, pol)
%
% Input:
% ds: thickness vector (um)
% ns(nk,3): n1, n2, ns
% dns(nk,3): first order index
% n0: initial index, assumed to be dispersionless
% ks: wavenumbers (1/um)
% theta: incidence angle (radians)
% pol: {'TE','TM'}
%
% Output:
% r2(nk): reflectivity
% gd(nk): group delay (fs)
% r2grad(nk,n)
% gdgrad(nk,n)

% Implementation:
% T matrix rows are the top elements of transfer matrix with
% layer changing along columns.
%
% [Better comments are needed.]
% Currently only works for even number of layers.


%%% Initialization and input handling.
c = .2997924580;  % um/fs
cinv = 1/c;
if (nargin < 7), pol = 'TE'; end
if (nargin < 8)
  dpQ = true;
elseif strcmp('nodp',approx)
  dpQ = false;
end
if (nargin < 6), theta = 0; end
M = size(ks);
if M(1) == 1
  ks = ks.';
  flipout = true;
else
  flipout = false;
end
nk = length(ks);  % number of wavelengths
M = size(ds);
if M(1) ~= 1  % d is not row vector
  ds = ds.';
  ns = ns.';
  dns = dns.';
end
ds = [ds, 0];  % create dummy layer for substrate
n = length(ds);  % number of layers
if mod(n,2) == 0  % odd number?
  oddQ = true;
else
  oddQ = false;
end

%%% Allocate memory.
pTEs = zeros(nk,4);
dpTEs = zeros(nk,4);
p0s = zeros(nk,4);
neff = zeros(n-1,1);
kneff = zeros(n-1,1);
kdneff = zeros(n-1,1);
p = zeros(n-1,1);
dp = zeros(n-1,1);
Tfor = zeros(n,2);  % T from 0 to L
Trev = zeros(n,2);  % T from N to L
dTfor = zeros(n,2);  % T from 0 to L
dTrev = zeros(n,2);  % T from N to L
Tlay = zeros(n,2);  % T at L
dTlay = zeros(n,2);  % dT at L
Tgrad = zeros(n-1,2);
dTgrad = zeros(n-1,2);
r2s = zeros(nk,1);
gds = zeros(nk,1);
r2grads = zeros(nk,n-1);
gdgrads = zeros(nk,n-1);

%%% Precalculate material variables.
% index variables: [n1, n2, nsub]
% p variables: [p01, p12, p21, p2/1sub]
% Effective indices.
neffs = sqrt(ns.^2 - (n0*sin(theta))^2);
neff0 = n0*cos(theta);
dneffs = ns.*dns./neffs;  % FIXME: not really needed for nsub
if oddQ
  lastx = 1;  % index index of last layer
else
  lastx = 2;
end
% Calculate pTE (needed regardless of polarization).
pTEs(:,1) = neff0./neffs(:,1); % n1 from air
pTEs(:,2) = neffs(:,1)./neffs(:,2); % n2 from n1
pTEs(:,3) = 1./pTEs(:,2);  % n1 from n2
pTEs(:,4) = neffs(:,lastx)./neffs(:,3);  % nsub from n2/n1
dpTEs(:,1) = -pTEs(:,1).*dneffs(:,1)./neffs(:,1);
dpTEs(:,2) = (dneffs(:,1) - pTEs(:,2).*dneffs(:,2))./neffs(:,2);
dpTEs(:,3) = (dneffs(:,2) - pTEs(:,3).*dneffs(:,1))./neffs(:,1);
dpTEs(:,4) = (dneffs(:,lastx) - pTEs(:,4).*dneffs(:,3))./neffs(:,3);
% Calculate p = pTM or p = pTE as needed.
if strcmp(pol, 'TM')
  % Calculate p = pTM from p0 and pTE.
  % p0 is just p (as above) for normal incidence.
  % p0 = ns(:,L - 1)./nL;
  % p = pTE./p0.^2;
  p0s(:,1) = n0./ns(:,1);
  p0s(:,2) = ns(:,1)./ns(:,2);
  p0s(:,3) = 1./p0s(:,2);
  p0s(:,4) = ns(:,lastx)./ns(:,3);
  dp0s(:,1) = -p0s(:,1).*dns(:,1)./ns(:,1);
  dp0s(:,2) = (dns(:,1) - p0s(:,2).*dns(:,2))./ns(:,2);
  dp0s(:,3) = (dns(:,2) - p0s(:,3).*dns(:,1))./ns(:,1);
  dp0s(:,4) = (dns(:,lastx) - p0s(:,4).*dns(:,3))./ns(:,3);
  ps = pTEs./p0s.^2;
  dps = (p0s.*dpTEs - 2*pTEs.*dp0s)./p0s.^3;
else  % pol == 'TE'
  ps = pTEs;
  dps = dpTEs;
end
if ~dpQ  % handle approximation toggle
  dps = 0*dps;
end
pps = 1 + ps;
pms = 1 - ps;

%%% Precalculations.
for kx = 1:nk,
  k = ks(kx);

  %%% Precalculate all layer matrices.
  for L = 1:(n-1),
    d = ds(L);

    % Select appropriate material parameter indices.
    if (L == n)  % last layer
      nx = 3;
      px = 4;
    elseif (L == 1)  % first layer
      nx = 1;
      px = 1;
    elseif (mod(L,2))  % odd layer, n2->n1
      nx = 1;
      px = 3;
    else  % even layer, n1->n2
      nx = 2;
      px = 2;
    end

    neff(L) = neffs(kx,nx);
    kneff(L) = k*neff(L);
    kdneff(L) = k*dneffs(kx,nx);
    p(L) = ps(kx,px);
    dp(L) = dps(kx,px);
    ephi = exp(-d*1i*kneff(L))/2;  % phasor/2 of layer
    Tlay(L,1) = ephi*pps(kx,px);  % T matrix of layer
    Tlay(L,2) = ephi*pms(kx,px);
    Dlay = -1i*d*(neffs(kx,nx) + k*dneffs(kx,nx));  % diff operator
    dTlay(L,1) = Dlay*Tlay(L,1) + ephi*dps(kx,px);
    dTlay(L,2) = Dlay*Tlay(L,2) - ephi*dps(kx,px);
  end  % precalc for over wavelengths
  % Last layer
  Tlay(n,1) = pps(kx,4)/2 + eps*1i;
  Tlay(n,2) = pms(kx,4)/2 + eps*1i;
  dTlay(n,1) = dps(kx,4)/2;
  dTlay(n,2) = -dps(kx,4)/2;

  %%% Step forward through structure, calculating forward matrices.
  Tfor(1,:) = Tlay(1,:);  % initial matrix
  dTfor(1,:) = dTlay(1,:);  % initial 1st derivative matrix
  for L = 2:n,  % changed from n-1
    Lp = L - 1;
    Tfor(L,:) = tmatrixprod(Tlay(L,:), Tfor(Lp,:));
    dTfor(L,:) = tmatrixprod(dTlay(L,:), Tfor(Lp,:)) + ...
      tmatrixprod(Tlay(L,:), dTfor(Lp,:));
  end

  %%% Step backward through structure, calculating reverse T matrices.
  Trev(n,:) = Tlay(n,:);
  dTrev(n,:) = dTlay(n,:);
  for L = (n-1):-1:2,
    Lp = L + 1;
    Trev(L,:) = tmatrixprod(Trev(Lp,:), Tlay(L,:));
    dTrev(L,:) = tmatrixprod(dTrev(Lp,:), Tlay(L,:)) + ...
      tmatrixprod(Trev(Lp,:), dTlay(L,:));
  end

  %%% Step through structure, calculating matrix gradients for T and T'.
  Tgrad(1,:) = tmatrixprod(Trev(2,:), -1i*kneff(1)*Tfor(1,:));
  
  dTdkgradmat = -(ds(1)*kneff(1) + 1i)*(neff(1) + kdneff(1)) - ...
    1i*kneff(1)*dp(1)*[1/(p(1) + 1), 1/(p(1) - 1)];
  
  dTgrad(1,:) = tmatrixprod(dTrev(2,:), -1i*kneff(1)*Tfor(1,:)) + ...
    tmatrixprod(Trev(2,:), dTdkgradmat.*Tlay(1,:));
    
  for L = 2:n-1,
    Tgrad(L,:) = tmatrixprod(Trev(L+1,:), -1i*kneff(L)*Tfor(L,:));
    
    dTdkgradmat = -(ds(L)*kneff(L) + 1i)*(neff(L) + kdneff(L)) - ...
      1i*kneff(L)*dp(L)*[1/(p(L) + 1), 1/(p(L) - 1)];
    
    dTgrad(L,:) = tmatrixprod(dTrev(L+1,:), -1i*kneff(L)*Tfor(L,:)) + ...
      tmatrixprod(Trev(L+1,:), tmatrixprod(dTdkgradmat.*Tlay(L,:), Tfor(L-1,:))) + ...
      tmatrixprod(Trev(L+1,:), tmatrixprod(-1i*kneff(L)*Tlay(L,:), dTfor(L-1,:)));
  end

  % Calculate complex reflectance and group delay scalars.
  T = Tfor(n,:);  % total transfer matrix
  T1 = T(1);
  T2 = T(2);
  T12 = T1^2;
  dT = dTfor(n,:);  % total transfer matrix derivative (wrt k)
  R = -T2/T1;  % complex reflection coefficient
  ReR = real(R);
  ImR = imag(R);
  r = abs(R);
  r2 = r^2;
  dR = (T2*dT(1) - dT(2)*T1)/T12;
  dr = (real(dR)*real(R) + imag(dR)*imag(R))/r;  % magnitude diff.
  dphi = (imag(dR)*real(R) - real(dR)*imag(R))/r2;  % phase diff.
  gd = cinv*dphi;
  
  % Calculate gradients.
  for L = 1:n-1,
    Rgrad = -(R*Tgrad(L,1) + Tgrad(L,2))/T1;
    ReRgrad = real(Rgrad);
    ImRgrad = imag(Rgrad);
    rgrad = (ReRgrad*ReR + ImRgrad*ImR)/r;
    r2grads(kx,L) = 2*r*rgrad;

    dRgrad = (Tgrad(L,2)*dT(1) + Tgrad(L,1)*dT(2) + ...
      R*(2*Tgrad(L,1)*dT(1) - T1*dTgrad(L,1)) - T1*dTgrad(L,2))/T12;
    phigrad = (ImRgrad*ReR - ReRgrad*ImR)/r2;
    dphigrad = (imag(dRgrad)*ReR - real(dRgrad)*ImR)/r2 - ...
      (phigrad*dr + rgrad*dphi)/r;
    gdgrads(kx,L) = 3.33564095198152*dphigrad;
  end
  
  % Outputs.
  r2s(kx) = r2;
  gds(kx) = gd;  
end  % foreach k

if flipout
  r2s = r2s.';
  gds = gds.';
  r2grads = r2grads.';
  gdgrads = gdgrads.';
end


function c = tmatrixprod(b, a)
% Compute the top row of the products of two T matrices given their top
% rows, as in c = b*a. Parallelized along columns.
c = [b(:,1).*a(:,1) + b(:,2).*conj(a(:,2)), ...
      b(:,1).*a(:,2) + b(:,2).*conj(a(:,1))];
