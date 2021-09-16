function [r, gd] = stackgd(ks, ds, n0, ns, dns, theta, pol, approx)
% STACKGD  Calculate complex reflection and group delay of
% arbitrary dielectric structure using exact analytic methods,
% including all relevent material derivatives.
%
% [r, gd] = stackgd(ks, ds, n0, ns, dns, theta, pol)
%
% Input:
% ds: thickness vector (um)
% ns(length(ks),length(ds)): zeroth index matrix, from first layer to substrate
% dns(length(ks),length(ds)): first order index
% n0: initial index, assumed to be dispersionless
% ks: wavenumber column vector (1/um)
% theta: incidence angle (radians)
% pol: {'TE','TM'}
% approx: {'none', 'nodp'}
%
% Output:
% r(ks): complex reflection coefficient
% gd(ks): group delay (fs)

% Implemenation Notes:
%
% This code is optimized more for simplicity and readability than speed.
%
% Each T matrix goes from the end of the last layer through the Lth layer.
% Only the first two elements of the transfer matrix are kept
% at each point due to the symmetry of T matrices.
%
% All matrices are arranged with wavenumbers in a column, as MATLAB can
% access matrices quicker along columns than rows.


% Initialization and input handling.
if (nargin < 7), pol = 'TE'; end
if (nargin < 6), theta = 0; end
if (nargin < 8), approx = 'none'; end
k = ks(:);  % ensure k is column vector of wavenumbers
nk = length(k);  % number of wavelengths
M = size(ds);
if M(1) ~= 1  % d is not row vector
  ds = ds.';
  ns = ns.';
  dns = dns.';
end
ds = [NaN, ds, 0];  % create dummy layers for outside and substrate
n = length(ds);  % number of layers
ns = [n0*ones(nk,1), ns];
if (dns == 0)
    dns = zeros(nk, n);
else
    dns = [zeros(nk,1), dns];
end
c = .2997924580;  % um/fs
n0sinth2 = n0*sin(theta)^2;

% Initialize loop variables.
T = [ones(nk,1), zeros(nk,1)];  % initial matrix
dT = zeros(nk,2);  % initial 1st derivative matrix
neffL = sqrt(ns(:,1).^2 - n0sinth2);  % neff of ambient medium (dispersionless)
dneffL = zeros(nk,1);
% Step through structure, calculating all wavelengths in parallel.
for L = 2:n,
    % Initialize local variables.
    d = ds(L);
    id = 1i*d;  % layer thickness times i
    
    % Calculate effective index quantities.
    neffLL = neffL; dneffLL = dneffL;
    nL = ns(:,L); dnL = dns(:,L);
    neffL = sqrt(nL.^2 - n0sinth2);
    dneffL = nL.*dnL./neffL;
        
    % Calculate pTE (needed regardless of polarization).
    pTE = neffLL./neffL;
    dpTE = (dneffLL - pTE.*dneffL)./neffL;
    
    % Calculate p = pTM or p = pTE as needed.
    % [TM part not yet debugged.]
    if (strcmp(pol,'TM'))
        % Calculate p0, which is just p (as above) for normal incidence.
        LL = L - 1;
        nLL = ns(:,LL); dnLL = dns(:,LL);
        p0 = nLL./nL;
        dp0 = (dnLL - p0.*dnL)./nL;
        
        % Calculate p = pTM from p0 and pTE.
        p = pTE./p0.^2;
        dp = (p0.*dpTE - 2*pTE.*dp0)./p0.^3;
    else  % pol == 'TE'
        p = pTE;
        dp = dpTE;
    end
    if (strcmp(approx,'nodp'))
        dp = 0;
    end
    pp = 1 + p;
    pm = 1 - p;
    
    % Find transfer matrices for current layer.
    ephi = exp(-id*k.*neffL);  % phasor of layer
    TL = [ephi.*pp, ephi.*pm]/2;
    dterm = -id*(neffL + k.*dneffL);
    dTL = [ephi.*(pp.*dterm + dp), ephi.*(pm.*dterm - dp)]/2;
    
    % Calculate new matrices for full structure up to L.
    Tnew = tmatrixprod(TL, T);
    dTnew = tmatrixprod(TL, dT) + tmatrixprod(dTL, T);

    T = Tnew; dT = dTnew;
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
