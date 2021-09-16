function w = pulsefwhm(tin, uin)

n = length(tin);

% find first guesses
t = linspace(tin(1), tin(end), 10*n);
pcentered = abs(centerpulse(uin)).^2;
u = interp1(tin, pcentered, t, 'cubic');
umax = max(u);
ks = find(u/umax > 0.5);
tmin0 = t(ks(1));
tmax0 = t(ks(end));
pnorm = pcentered/umax;

% refine
x0 = [tmin0, tmax0];
optopt = optimset('LargeScale', 'off', 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'final');
x = fminunc(@fwhmerr, x0, optopt);
w = abs(x(2) - x(1));

  function z = fwhmerr(x)
    us = interp1(tin, pnorm, x, 'cubic');
    z = sum(abs(0.5 - us).^2);
  end

end