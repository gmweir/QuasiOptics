function [errsum, errgrad] = mirrorgderrmex(...
  ds, ks, n0, ns, dns, theta, pol, rgoal, gdgoal0, wr, wgd)
% Computes higher order error for reflection to approximate an inequality
% constraint. GD error is computed for an optimal offset and factor.

nord = 4;  % must be even
n = length(ds);
nk = length(ks);
wgdsum = sum(wgd);
gdgoal = gdgoal0 - sum(gdgoal0.*wgd)/wgdsum;

if (nargout == 1)
  %fprintf('mirrorgderrmex: error eval.\n');

  [r, gd] = stackgdfastmex(ks, ds.^2, n0, ns, dns, theta, pol);
  %[r2, gd2, rgrad2, gdgrad2] = stackgdgrad(ks.', (ds.^2).', n0, ns.', dns.', theta, pol);
%   [r2, gd2, rgrad2, gdgrad2] = stackgdgradfix(ks.', (ds.^2).', n0, ns.', dns.', theta, pol);
%   r = r2.';
%   gd = gd2.';
%   rgrad = rgrad2.';
%   gdgrad = gdgrad2.';

  gdmean = sum(gd.*wgd)/wgdsum;
  gdtot = gd - gdmean;
else
  %fprintf('mirrorgderrmex: gradient eval.\n');
  [r, gd, rgrad, gdgrad] = stackgdgradfastmex(ks, ds.^2, n0, ns, dns, theta, pol);
  %[r, gd, rgrad, gdgrad] = stackgradold(ks, ds.^2, n0, ns, dns, theta, pol);
  %[r2, gd2, rgrad2, gdgrad2] = stackgdgrad(ks.', (ds.^2).', n0, ns.', dns.', theta, pol);
  %[r2, gd2, rgrad2, gdgrad2] = stackgdgradfix(ks.', (ds.^2).', n0, ns.', dns.', theta, pol);
%   r = r2.';
%   gd = gd2.';
%   rgrad = rgrad2.';
%   gdgrad = gdgrad2.';
  
  gdmean = sum(gd.*wgd)/wgdsum;
  gdtot = gd - gdmean;
  
  gdmeangrad = sum(repmat(wgd,n,1).*gdgrad, 2)/wgdsum;
  gdtotgrad = gdgrad - repmat(gdmeangrad,1,nk);
  
  % compute gradient of merit function
  rerrgrad = sum(nord*repmat(wr.*(r-rgoal).^(nord-1), n, 1).*rgrad, 2);
  gderrgrad = sum(2*repmat(wgd.*(gdtot - gdgoal), n, 1).*gdtotgrad, 2);
  dserrgrad = 2.*ds.*(rerrgrad + gderrgrad);
    
  % combine into full gradient
  errgrad = dserrgrad;
end

rerr = sum(wr.*(r - rgoal).^nord);
gderr = sum(wgd.*(gdtot - gdgoal).^2);

errsum = rerr + gderr;
