function [r2s, gds, r2grads, gdgrads] = stackgdgradfast(ks, ds, n0, ns, dns, theta, pol)
% STACKGDGRADFAST  Calculate complex reflection, group delay and GDD of
% arbitrary dielectric structure using approximate analytic methods,
% Assumes two indices.
%
% [r, gd, gdd] = stackgdgradfast(k, ds, n0, ns, dns, ddns, theta, pol)
%
% Input:
% ds: thickness vector (um)
% ns(nk,3): n1, n2, nsubstrate
% dns(nk,3): first order index
% n0: initial index, assumed to be dispersionless
% ks: wavenumbers (1/um)
% theta: incidence angle (radians)
% pol: {'TE','TM'}
%
% Output:
% r(nk): complex reflection coefficient
% gd(nk): group delay (fs)
% rgrad(nk,n)
% gdgrad(nk,n)

% Implementation:
% T matrix rows are the top row of transfer matrices with
% layer changing along columns.
% The individual layers are composed of a reflection (R matrix) followed by
% the propagation (P diagonal matrix). Since we're only dealing with the
% top row of the T matrix, the P matrix is actually just a scalar.

% TODO: Clean up initialization and add comments.


%%% Initialization and input handling.
if (nargin < 7), pol = 'TE'; end
if (nargin < 6), theta = 0; end
nk = length(ks);  % number of wavelengths
ds = [ds(:); 0];  % create dummy layer for substrate
n = length(ds);  % number of layers (including substrate)

%%% Allocate memory.
pTEs = zeros(nk,4);
p0s = zeros(nk,4);
Tfor = zeros(n,2);  % T from 0 to L
Trev = zeros(n,2);  % T from N to L
dTfor = zeros(n,2);  % T from 0 to L
dTrev = zeros(n,2);  % T from N to L
Tlay = zeros(n,2);  % T at L
Dlay = zeros(n,1);  % single layer d/dk operator
neff = zeros(n,1);
kneff = zeros(n,1);
p = zeros(n,1);
dTgradfor = zeros(n-1,2);  % d^2/(dk dd) of Tfor
Tgrad = zeros(n-1,2);
dTgrad = zeros(n-1,2);
r2s = zeros(nk,1);
gds = zeros(nk,1);
r2grads = zeros(nk,n-1);
gdgrads = zeros(nk,n-1);

%%% Precalculate material variables.
% index variables: [n1, n2, nsub]
% p variables: [p01, p12, p21, p2sub]
% Effective indices.
neffs = sqrt(ns.^2 - (n0*sin(theta))^2);
neff0 = n0*cos(theta);
dneffs = ns.*dns./neffs;  % not really needed for nsub?
% Calculate pTE (needed regardless of polarization).
pTEs(:,1) = neff0./neffs(:,1); % n1 from air
pTEs(:,2) = neffs(:,1)./neffs(:,2); % n2 from n1
pTEs(:,3) = 1./pTEs(:,2);  % n1 from n2
pTEs(:,4) = neffs(:,2)./neffs(:,3);  % nsub from n2
% Calculate p = pTM or p = pTE as needed.
if strcmp(pol, 'TM')
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
pps = 1 + ps;
pms = 1 - ps;


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
    p(L) = ps(kx,px);
    ephi = exp(-d*i*kneff(L));  % phasor/2 of layer
    Rlay(L,:) = [pps(kx,px)/2, pms(kx,px)/2];  % R matrix between L and L+1
    Play(L) = ephi;  % P matrix element for layer L
    Dlay(L) = -i*d*(neffs(kx,nx) + k*dneffs(kx,nx));
  end
  % Last layer
  Tlay(n,1) = pps(kx,4)/2 + eps*i;
  Tlay(n,2) = pms(kx,4)/2 + eps*i;
  Rlay(n,:) = [pps(kx,4)/2, pms(kx,4)/2];

  %%% Step forward through structure, calculating forward matrices.
  % Tfor(n,:) is the nth layer, up to, but not including, the interface
  % with the next layer.
  Tfor(1,:) = Play(1)*Rlay(1,:);  % initial matrix
  dTfor(1,:) = Dlay(1)*Tfor(1,:);  % initial 1st derivative matrix
  dTgradfor(1,:) = (1/ds(1) + i*kneff(1))*dTfor(1,:);  % gradient term
  for L = 2:n-1,
    Lp = L - 1;

    % Inductive T-matrix computations.
    Tlay = Play(L)*Rlay(L,:);
    Tfor(L,:) = tmatrixprod(Tlay, Tfor(Lp,:));
    TlaydTfor = tmatrixprod(Tlay, dTfor(Lp,:));
    dTlayTfor = Dlay(L)*Tfor(L,:);
    dTfor(L,:) = dTlayTfor + TlaydTfor;

    % Stuff for gradient, to be used later. (See paper.)
    dTgradfor(L,:) = (1/ds(L) - i*kneff(L))*dTlayTfor - i*kneff(L)*TlaydTfor;
  end
  % Last layer (substrate).
  Lp = n - 1;
  Tfor(n,:) = tmatrixprod(Rlay(n,:), Tfor(Lp,:));
  dTfor(n,:) = tmatrixprod(Rlay(n,:), dTfor(Lp,:));

  %%% Step backward through structure, calculating reverse T matrices.
  % Defined exactly as for the forward case, except everything
  % is reciprocal. We keep track of intermediate results without the
  % propagations (ending at interfaces) so that we can combine with the
  % forward matrices. The intermediate results are in Trev and dTrev.
  Trev(n,:) = revtmatrix(Rlay(n,:));  % (dTrev already initialized to zero)
  Trevcum = Play(n-1)*Trev(n,:);  % cumulative matrix
  dTrevcum = Dlay(n-1)*Trevcum;  % cumulative derivative matrix
  for L = (n-1):-1:2,
    Lp = L + 1;
    
    % Pick off just the part up to the interface.
    % TODO: Use variable for revtmatrix(Rlay).
    Trev(L,:) = trevprod(Rlay(L,:), Trevcum);
    dTrev(L,:) = trevprod(Rlay(L,:), dTrevcum);
    
    % Finish inductive computation with full layers.
    Trevcum = Play(L-1)*Trev(L,:);
    TlayrevdTrevcum = Play(L-1)*dTrev(L,:);
    dTlayrevTrevcum = Dlay(L-1)*Trevcum;
    dTrevcum = TlayrevdTrevcum + dTlayrevTrevcum;
  end

  %%% Step through structure, calculating matrix gradients for T and T'.
  for L = 1:n-1,
    Tgrad(L,:) = trevprod(Trev(L+1,:), -i*kneff(L)*Tfor(L,:));
    dTgrad(L,:) = trevprod(dTrev(L+1,:), -i*kneff(L)*Tfor(L,:)) + ...
      trevprod(Trev(L+1,:), dTgradfor(L,:));
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
  gd = 3.33564095198152*dphi;

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


%%%%% Subfunctions %%%%%

function c = tmatrixprod(b, a)
% Compute the top row of the products of two T matrices given their top
% rows, as in c = b*a. Parallelized along columns.
c = [b(:,1).*a(:,1) + b(:,2).*conj(a(:,2)), ...
  b(:,1).*a(:,2) + b(:,2).*conj(a(:,1))];

function c = trevprod(b, a)
% Compute matrix product with b reversed.
c = [b(:,1).*a(:,1) + conj(b(:,2).*a(:,2)), ...
  b(:,1).*a(:,2) + conj(b(:,2).*a(:,1))];
    
function c = revtmatrix(a)
% Compute reverse T matrix, using b for the determinant, if supplied.
c = [a(1), conj(a(2))];
