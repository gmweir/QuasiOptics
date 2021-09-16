function [Tout, Rout] = stacktr(ks, ds, n0, ns, theta, pol)
% STACKR  Calculate complex transmission/reflection of
% arbitrary dielectric structure.
%
% [t, r] = stackt(ks, ds, n0, ns, theta, pol)
%
% Input:
% ds: thickness vector (um)
% ns(ks,ds): zeroth index matrix, from first layer to substrate
% n0: initial index, assumed to be dispersionless
% ks: wavenumber column vector (1/um)
% theta: incidence angle (radians)
% pol: {'TE','TM'}
%
% Output:
% t(ks): complex reflection coefficient

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
if (nargin < 6)
  pol = 'TE'; end
if (nargin < 5)
  theta = 0; end
k = ks(:);  % ensure k is column vector of wavenumbers
ds = ds(:).';  % ensure ds is row vector
nk = length(k);  % number of wavelengths
ds = [NaN, ds, 0];  % create dummy layers for outside and substrate
n = length(ds);  % number of layers
ns = [n0*ones(nk,1), ns];
%c = .2997924580;  % um/fs
n0sinth2 = (n0*sin(theta))^2;

% Initialize loop variables.
T = [ones(nk,1), zeros(nk,1)];  % initial matrix
neffL = sqrt(ns(:,1).^2 - n0sinth2);  % neff of ambient medium (dispersionless)
% Step through structure, calculating all wavelengths in parallel.
for L = 2:n,
    % Initialize local variables.
    d = ds(L);
    id = 1i*d;  % layer thickness times i
    
    % Calculate effective index quantities.
    neffLL = neffL;
    nL = ns(:,L);
    neffL = sqrt(nL.^2 - n0sinth2);
        
    % Calculate pTE (needed regardless of polarization).
    pTE = neffLL./neffL;
    
    % Calculate p = pTM or p = pTE as needed.
    % [TM part not yet debugged.]
    if (pol == 'TM')
        % Calculate p0, which is just p (as above) for normal incidence.
        LL = L - 1;
        nLL = ns(:,LL);
        p0 = nLL./nL;
        
        % Calculate p = pTM from p0 and pTE.
        p = pTE./p0.^2;
    else  % pol == 'TE'
        p = pTE;
    end
    pp = 1 + p;
    pm = 1 - p;
    
    % Find transfer matrices for current layer.
    ephi = exp(-id*k.*neffL);  % phasor of layer
    TL = [ephi.*pp, ephi.*pm]/2;
    
    % Calculate new matrices for full structure up to L.
    %Tnew = tmatrixprod(TL, T);
    Tnew = [T(:,1).*TL(:,1) + conj(T(:,2)).*TL(:,2), ...
      T(:,2).*TL(:,1) + conj(T(:,1)).*TL(:,2)];
    
    T = Tnew;
end

% Calculate complex reflectance and group delay from total T matrices.
Rout = -T(:,2)./T(:,1);  % complex reflection coefficient
Tout = T(:,1) - T(:,2).*conj(T(:,2))./conj(T(:,1));  % complex transmission

% function c = tmatrixprod(b, a)
% % Compute the top row of the products of two T matrices given their top
% % rows, as in c = b*a. Parallelized along columns.
% c = [a(:,1).*b(:,1) + conj(a(:,2)).*b(:,2), ...
%   a(:,2).*b(:,1) + conj(a(:,1)).*b(:,2)];
