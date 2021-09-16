%% Load data

clear
load 'teststack'
s = readmatstruct;

%ds1 = ds1(1:5);

n0 = 1;  % dispersionless
nk = 256;
kx = 1;  % look at longest wavelength, which likely penetrates the most
y1 = .6; y2 = 1.2;
ks = linspace(2*pi/y2, 2*pi/y1, nk).';
dk = ks(2)-ks(1);
ys = 2*pi./ks;
nk = length(ks);
y0 = 0.65;
k0 = 2*pi/y0;


%% Data setup

nl = length(ds1);
np = nl/2;

s1 = s.SiO2_NL;
s2 = s.TiO2_NL;
ssub = s.FS;

[n1, n1d] = dndksell(ks, s1);
[n2, n2d] = dndksell(ks, s2);
[nsub, nsubd] = dndksell(ks, ssub);

ns2 = [n1 n2 nsub];
dns2 = [n1d n2d nsubd];

ksmex = ks';
ds1mex = ds1';
ns2mex = ns2.';
dns2mex = dns2.';

%% MEX-file tests

reps = 180;
tic
  for k = 1:reps
    [rmex2, gdmex2, rgradmex2, gdgradmex2] = ...
        stackgdgradfastmex(ksmex, ds1mex, n0, ns2mex, dns2mex, theta, 'TM');
  end
toc
tic
  for k = 1:reps
    [rmex3, gdmex3, rgradmex3, gdgradmex3] = ...
        stackgdgradmex(ksmex, ds1mex, n0, ns2mex, dns2mex, theta, 'TM');
  end
toc


%% M-file tests

% The gold standard...
[r2exp, gdexp, r2gradexp, gdgradexp] = ...
    stackgdgrad(ks, ds1, n0, ns2, dns2, theta, 'TM');
[r2ex, gdex, r2gradex, gdgradex] = ...
    stackgdgrad(ks, ds1, n0, ns2, dns2, theta, 'TM', 'nodp');


%% Plots

figure(4)
xlabel('layer')
plot(1:nl, [rgradmex2(:,kx).' - r2gradex(kx,:); ...
            rgradmex3(:,kx).' - r2gradex(kx,:)])
legend('stackgdgradfastmex', 'stackgdgradmex')
title('R Gradient Error')

figure(5)
xlabel('layer')
plot(1:nl, gdgradmex2(:,kx).' - gdgradex(kx,:))
title('Fast GD Gradient Error')

figure(6)
xlabel('layer')
plot(1:nl, gdgradmex3(:,kx).' - gdgradexp(kx,:))
title('gdgradmex GD Gradient Error')

figure(7)
plot(ks, gdmex3.' - gdexp)
title('gdgradmex GD error')

% figure(8)
% plot(ks, [r2exp rmex3.' rmex2.'])
