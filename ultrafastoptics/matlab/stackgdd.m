function [r, gd, gdd] = stackgdd(ks, ds, n0, ns, dns, ddns, theta, pol, approx)
% STACKGDD  Calculate complex reflection, group delay and GDD of
% arbitrary dielectric structure using exact analytic methods,
% including all relevent material derivatives.
%
% [r, gd, gdd] = stackgdd(ks, ds, n0, ns, dns, ddns, theta, pol)
%
% Input:
% ks: wavenumber column vector (1/um)
% ds: thickness row vector (um)
% ns(length(ks),length(ds)): index matrix, from first layer to substrate
% dns(length(ks),length(ds)): first order index derivatives at ks
% ddns(length(ks),length(ds)): second order index derivatives at ks
% n0: initial index, assumed to be dispersionless
% theta: incidence angle (radians)
% pol: {'TE','TM'}
% approx: {'none', 'nodp'}  (for testing purposes)
%
% Output:
% r(ks): complex reflection coefficient
% gd(ks): group delay (fs)
% gdd(ks): first order dispersion (fs^2)

% Implemenation Notes:
%
% This code is optimized more for simplicity and readability than speed.
%
% Each T matrix goes from the end of the last layer through the Lth layer.
% Only the first two elements of the transfer matrix are kept
% at each point due to the symmetry of T matrices.
%
% All matrices are arranged with wavenumbers in columns.


% Initialization and input handling.
if (nargin < 8), pol = 'TE'; end
if (nargin < 7), theta = 0; end
if (nargin < 9), approx = 'none'; end
k = ks(:);  % ensure k is column vector of wavenumbers
ds = ds(:).';
nk = length(ks);  % number of wavelengths
ds = [NaN, ds, 0];  % create dummy layers for outside and substrate
n = length(ds);  % number of layers
ns = [n0*ones(nk,1), ns];
if (dns == 0)
  dns = zeros(nk, n);
else
  dns = [zeros(nk,1), dns];
end
if (ddns == 0)
  ddns = zeros(nk, n);
else
  ddns = [zeros(nk,1), ddns];
end
c = .2997924580;  % um/fs

% Initialize loop variables.
n0sinth2 = (n0*sin(theta))^2;
T = [ones(nk,1), zeros(nk,1)];  % initial matrix
dT = zeros(nk,2);  % initial 1st derivative matrix
ddT = zeros(nk,2);
neffL = sqrt(ns(:,1).^2 - n0sinth2);  % neff of ambient medium
dneffL = zeros(nk,1);
ddneffL = zeros(nk,1);
%%% Step through structure, calculating all wavelengths in parallel.
for L = 2:n,
  % Initialize local variables.
  d = ds(L);

  % Calculate effective index quantities.
  neffLL = neffL; dneffLL = dneffL; ddneffLL = ddneffL;
  nL = ns(:,L); dnL = dns(:,L); ddnL = ddns(:,L);
  neffL = sqrt(nL.^2 - n0sinth2);
  dneffL = nL.*dnL./neffL;
  ddneffL = (dnL.^2 - dneffL.^2 + nL.*ddnL)./neffL.^2;

  % Calculate pTE (needed regardless of polarization).
  pTE = neffLL./neffL;
  dpTE = (dneffLL - pTE.*dneffL)./neffL;
  ddpTE = (ddneffLL - 2*dneffL.*dpTE - pTE.*ddneffL)./neffL;

  % Calculate p = pTM or p = pTE as needed.
  if (strcmp(pol,'TM'))
    % Calculate p0, which is just p (as above) for normal incidence.
    LL = L - 1;
    nLL = ns(:,LL); dnLL = dns(:,LL); ddnLL = ddns(:,LL);
    p0 = nLL./nL;
    dp0 = (dnLL - p0.*dnL)./nL;
    ddp0 = (ddnLL - 2*dnL.*dp0 - p0.*ddnL)./nL;

    % Calculate p = pTM from p0 and pTE.
    p = pTE./p0.^2;
    dp = (p0.*dpTE - 2*pTE.*dp0)./p0.^3;
    ddp = (6*pTE.*dp0.^2 - 4*p0.*dp0.*dpTE - 2*p0.*pTE.*ddp0 + ...
      p0.^2.*ddpTE)./p0.^4;
  else  % pol == 'TE'
    p = pTE;
    dp = dpTE;
    ddp = ddpTE;
  end
  if (strcmp(approx,'nodp'))
    dp = 0;
    ddp = 0;
  end
  pp = 1 + p;
  pm = 1 - p;

  % Find transfer matrices for current layer.
  ephi = exp(-1i*d*k.*neffL);  % phasor of layer
  TL = [ephi.*pp, ephi.*pm]/2;  % current layer matrix elements
  dterm = -1i*d*(neffL + k.*dneffL);
  dTL = [ephi.*(pp.*dterm + dp), ephi.*(pm.*dterm - dp)]/2;
  ddTL(:,1) = ephi.*(d*(-pp.*(2i.*dneffL + d*(neffL + k.*dneffL).^2) - ...
    2i*(neffL + k.*dneffL).*dp - 1i*k.*pp.*ddneffL) + ddp)/2;
  ddTL(:,2) = ephi.*(d*(-d*neffL.^2.*pm + ...
    2*neffL.*(-d*k.*pm.*dneffL + 1i*dp) + ...
    dneffL.*(-pm.*(2i + d*k.^2.*dneffL) + 2i*k.*dp) - ...
    1i*k.*pm.*ddneffL) - ddp)/2;

  % Calculate new matrices for full structure up to L.
  Tnew = tmatrixprod(TL, T);
  dTnew = tmatrixprod(TL, dT) + tmatrixprod(dTL, T);
  ddTnew = tmatrixprod(ddTL, T) + 2*tmatrixprod(dTL, dT) + ...
    tmatrixprod(TL, ddT);
  T = Tnew; dT = dTnew; ddT = ddTnew;
end

% Calculate complex reflectance and group delay from total T matrices.
r = -T(:,2)./T(:,1);  % complex reflection coefficient
dr = (T(:,2).*dT(:,1) - dT(:,2).*T(:,1))./T(:,1).^2;
ddr = (T(:,2).*(-2*dT(:,1).^2 + T(:,1).*ddT(:,1)) + ...
  T(:,1).*(2*dT(:,1).*dT(:,2) - T(:,1).*ddT(:,2)))./T(:,1).^3;
dphi = (imag(dr).*real(r) - real(dr).*imag(r))./abs(r).^2;
dmag = (real(dr).*real(r) + imag(dr).*imag(r))./abs(r);
ddphi = (imag(ddr).*real(r) - real(ddr).*imag(r))./abs(r).^2 - ...
  2*dmag.*dphi./abs(r);
gd = dphi/c;
gdd = ddphi/c^2;


function c = tmatrixprod(b, a)
% Compute the top row of the products of two T matrices given their top
% rows, as in c = b*a. Parallelized along columns.
c = [a(:,1).*b(:,1) + conj(a(:,2)).*b(:,2), ...
  a(:,2).*b(:,1) + conj(a(:,1)).*b(:,2)];
