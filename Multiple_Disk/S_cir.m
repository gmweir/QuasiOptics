function [S_m] = S_cir(er_pr, tnd, dz, rd)
% Scattering matrix of a disk inside a waveguide with radius rd and
% thickness dz.
% Complex permitivitty er = er_pr*(1-j*tnd) with tnd = loss tangent,
% Rho and Scattering matrix: Nickel and Thum
% Implemented in the function: Rho for TE11 only
% jwo, 22/02/2014
% jwo, 24.07.2017 in Git

global x_p f % mode and frequency

e0 = 8.849E-12; u0 = 4*pi*1E-7; c=1/sqrt(e0*u0);

% Complex permitivitty
er = er_pr*(1-j*tnd);

% Propagation constant gamma
k_0 = 2*pi.*f/c; % phase constant in vacuum
k_cmn = x_p/rd; % modification phase constant resulting from waveguide walls
gm_mn = sqrt((k_cmn^2-er.*k_0.^2)); % let MatLab determine real and img. part


% Reflection coefficient rho for TE
rho_mn = (sqrt(1-(k_cmn./k_0).^2) - sqrt(er_pr-(k_cmn./k_0).^2))...
       ./ (sqrt(1-(k_cmn./k_0).^2) + sqrt(er_pr-(k_cmn./k_0).^2));
   
% Scattering matrix S_m
ab0 = 1-rho_mn.^2.*exp(-2*gm_mn*dz);
a1 = rho_mn.*(1-exp(-2.*gm_mn*dz));
b1 = (1-rho_mn.^2).*exp(-1*gm_mn*dz);
b2 = (1-rho_mn.^2).*exp(-1*gm_mn*dz);
a2 = rho_mn.*(1-exp(-2*gm_mn*dz));
S_m = zeros(2);
S_m(1,1) = a1./ab0;
S_m(1,2) = b2./ab0;
S_m(2,1) = b1./ab0;
S_m(2,2) = a2./ab0;
end