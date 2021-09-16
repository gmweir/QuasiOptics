% setup parameters and variables
n0 = 1.0;
np = 40;
ny = 1024;
y1 = 0.6; y2 = 1.2;
y0 = y1;
%ys = linspace(y1, y2, ny).'; ks = 2*pi./ys;
ks = linspace(2*pi./y2, 2*pi./y1, ny).'; ys = 2*pi./ks;

% load material indices
load demodata
[nlo, nlod, nlodd] = dndksell(ks, s.sio2);  % first layer
[nhi, nhid, nhidd] = dndksell(ks, s.tio2);  % second layer
[nsub, nsubd, nsubdd] = dndksell(ks, s.fs);  % substrate
%nsub = 1.5*ones(ny,1); nsubd = 0*nsub; nsubdd = 0*nsub;

% compute repeating matrix of indices for use in full GDD function
if mod(np,2)  % odd
  ns = [repmat([nlo, nhi], 1, np), nlo, nsub];
  dns = [repmat([nlod, nhid], 1, np), nlod, nsubd];
  ddns = [repmat([nlodd, nhidd], 1, np), nlodd, nsubdd];
else
  ns = [repmat([nlo, nhi], 1, np), nsub];
  dns = [repmat([nlod, nhid], 1, np), nsubd];
  ddns = [repmat([nlodd, nhidd], 1, np), nsubdd];
end

% compute simple two material matrix for use in fast GDD function
ns2 = [ns(:,1), ns(:,2), ns(:,end)];
dns2 = [dns(:,1), dns(:,2), dns(:,end)];

% generate contrived stack example
dhigh = y0/nhi(end)/4;
dlow = y0/nlo(end)/4;
ds0 = repmat([dlow, dhigh], 1, np);
nl = length(ds0);
ds = ds0.*linspace(1,2,nl);
theta = 9/180*pi;
pol = 'TM';

% run computations
[r, gd, gdd] = stackgdd(ks,ds,n0,ns,dns,ddns,theta,pol);
[rfast, gdfast] = stackgddfast(ks, ds, n0, ns2, dns2, theta, pol);
[rfastmex, gdfastmex, temp] = stackgddfastmexwrap(ks, ds, n0, ns2, dns2, theta, pol);

% fd
c = 0.299792458;
gdfd = [diff(unwrap(angle(r)))/(ks(2) - ks(1))/c; gd(end)];

% plot results
figure(1)
clf
plot(ys, [gdfd, gd, gdfast, gdfastmex])
legend('fd', 'stackgdd', 'new', 'fast', 'fastmex','Location','NorthWest')
title('GD Test')
axis tight

% figure(2)
% plot(ys, [gdd, gddfast, gddfastmex.'])
% legend('stackgdd', 'stackgddfast', 'gddfastmex')
% title('GDD Test')
% axis tight
