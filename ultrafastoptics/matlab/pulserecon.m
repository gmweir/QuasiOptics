function [t, p, freg, spreg, phi] = pulserecon(y, sp, gd, npad)
%PULSERECON  Reconstruct pulse in time domain
%  The y vector is in nm, and need not be uniform.
%  The group delay is in femtoseconds.
%  The spectrum is proportional to power.

%  TODO: Switch to um units.


% Initialization.
nf = length(y);
if nargin > 3
  n = npad;
else
  n = nf;
end

% Convert to frequency.
f = 299./y;  % PHz
dydf = 299./f.^2;

% Fix spectrum and transform to field spectrum.
sp(sp < 0) = 0;
spf = sqrt(sp.*dydf);

% Put on uniform frequency grid.
freg = linspace(f(end), f(1), nf);
spreg = interp1(f, spf, freg, 'cubic');
gdreg = interp1(f, gd, freg, 'cubic');

% Compute reconstruction.
[maxsp, maxk] = max(spreg);
spreg = spreg/maxsp;
%gdreg = gdreg - gdreg(maxk);
phi = 2*pi*cumtrapz(freg, gdreg);
pbar = ifftshift(ifft(spreg.*exp(i*phi), n));
p = pbar/max(abs(pbar));

df = (freg(2) - freg(1));  % PHz (1/fs)
T = 1/df;
dt = T/n;
t = [-fliplr(linfill(dt, dt, floor(n/2))), linfill(0, dt, ceil(n/2))];
