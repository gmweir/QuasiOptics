function [errsum, errgrad] = mirrorfloatmex2(X, ks, n0, ns, dns, theta, pol, ...
  rgoal, gdgoalstatic, wr, wgd, gdgoalvar)
% Computes higher order error for reflection to approximate an inequality
% constraint. Zero GD point is chosen to minimize goal function. Considers
% GD error from both pair as well as individual layers. GD is set to be the
% GD of a single mirror in the pair (i.e. half the total).

nord = 4;  % must be even, and should also be input parameter.
nk = length(ks);
n = length(X)-1;
T = X(1);
ds = X(2:end);
wgdsum = sum(wgd);
gdgoal0 = T*gdgoalvar + gdgoalstatic;
gdgoal = gdgoal0 - sum(gdgoal0.*wgd)/wgdsum;

if (nargout < 2)
  [r, gd] = stackgdfastmex(ks, ds.^2, n0, ns, dns, theta, pol);
  
  gdmean = sum(gd.*wgd)/wgdsum;
  gdtot = gd - gdmean;
else
  [r, gd, r2grad, gdgrad] = stackgdgradfastmex(ks, ds.^2, n0, ns, dns, theta, pol);
  
  gdmean = sum(gd.*wgd)/wgdsum;
  gdtot = gd - gdmean;
  gdmeangrad = sum(repmat(wgd,n,1).*gdgrad, 2)/wgdsum;
  gdtotgrad = gdgrad - repmat(gdmeangrad,1,nk);
  
  wtgderr = wgd.*(gdtot - gdgoal);
  rerrgrad = sum(nord*repmat(wr.*(r-rgoal).^(nord-1), n, 1).*r2grad, 2);
  gderrgrad = sum(2*repmat(wtgderr, n, 1).*gdtotgrad, 2);
  dserrgrad = (2*ds.*(rerrgrad + gderrgrad)).';
  
  terrgrad = -sum(2*wtgderr.*(gdgoalvar - sum(wgd.*gdgoalvar)));
  
  errgrad = [terrgrad, dserrgrad];
end

rerr = sum(wr.*(r - rgoal).^nord);
gderr = sum(wgd.*(gdtot - gdgoal).^2);

errsum = rerr + gderr;
