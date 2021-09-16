function [errsum, errgrad] = mirrorpairfloatmex2(X, ks, n0, ns, dns, theta, pol, ...
  r1goal, r2goal, gdgoalstatic, wr1, wr2, wgd, gdgoalvar)
% Computes higher order error for reflection to approximate an inequality
% constraint. Zero GD point is chosen to minimize goal function. Considers
% GD error from both pair as well as individual layers. GD is set to be the
% GD of a single mirror in the pair (i.e. half the total).

indweight = 0.5e-1;  % should be input parameter!
nord = 4;  % must be even, and should also be input parameter.

nk = length(ks);
n1 = (length(X)-1)/2;
n2 = length(X) - n1 - 1;
T = X(1);
ds1 = X(2:n1+1);
ds2 = X(n1+2:end);
wgdsum = sum(wgd);
gdgoal0 = T*gdgoalvar + gdgoalstatic;
gdgoal = gdgoal0 - sum(gdgoal0.*wgd)/wgdsum;

%if (min(ds) < 0)
%  warning('negative layer.')
%end

if (nargout < 2)
  [r1, gd1] = stackgdfastmex(ks, ds1.^2, n0, ns, dns, theta, pol);
  [r2, gd2] = stackgdfastmex(ks, ds2.^2, n0, ns, dns, theta, pol);

  gd1mean = sum(gd1.*wgd)/wgdsum;
  gd2mean = sum(gd2.*wgd)/wgdsum;
  gdtot = (gd1 + gd2 - gd1mean - gd2mean)/2;
else
  [r1, gd1, r2grad1, gdgrad1] = stackgdgradfastmex(ks, ds1.^2, n0, ns, dns, theta, pol);
  [r2, gd2, r2grad2, gdgrad2] = stackgdgradfastmex(ks, ds2.^2, n0, ns, dns, theta, pol);

  gd1mean = sum(gd1.*wgd)/wgdsum;
  gd2mean = sum(gd2.*wgd)/wgdsum;
  gdtot = (gd1 + gd2 - gd1mean - gd2mean)/2;
  
  gd1meangrad = sum(repmat(wgd,n1,1).*gdgrad1, 2)/wgdsum;
  gd2meangrad = sum(repmat(wgd,n2,1).*gdgrad2, 2)/wgdsum;
  gdtotgrad = [gdgrad1 - repmat(gd1meangrad,1,nk); ...
    gdgrad2 - repmat(gd2meangrad,1,nk)]/2;

  rerrgrad = [sum(nord*repmat(wr1.*(r1-r1goal).^(nord-1), n1, 1).*r2grad1, 2); ...
    sum(nord*repmat(wr2.*(r2 - r2goal).^(nord-1), n2, 1).*r2grad2, 2)];
  wtgderr = wgd.*(gdtot - gdgoal);
  wtgd1err = wgd.*(gd1 - gd1mean - gdgoal);
  wtgd2err = wgd.*(gd2 - gd2mean - gdgoal);
  gderrgrad = sum(2*repmat(wtgderr, n1+n2, 1).*gdtotgrad, 2) + ...
    2*indweight*sum(2*[repmat(wtgd1err, n1,1); repmat(wtgd2err,n2,1)].*gdtotgrad, 2);
  dserrgrad = (2*[ds1; ds2].*(rerrgrad + gderrgrad)).';

  gdgoalvarmean = sum(wgd.*gdgoalvar);
  terrgradtot = -sum(2*wtgderr.*(gdgoalvar - gdgoalvarmean));
  terrgrad1 = -indweight*sum(2*wtgd1err.*(gdgoalvar - gdgoalvarmean));
  terrgrad2 = -indweight*sum(2*wtgd2err.*(gdgoalvar - gdgoalvarmean));

  errgrad = [terrgradtot + terrgrad1 + terrgrad2, dserrgrad];
end

rerr = sum(wr1.*(r1 - r1goal).^nord + wr2.*(r2 - r2goal).^nord);
gderr = sum(wgd.*(gdtot - gdgoal).^2);
gd1err = sum(wgd.*(gd1 - gd1mean - gdgoal).^2);
gd2err = sum(wgd.*(gd2 - gd2mean - gdgoal).^2);

errsum = rerr + gderr + indweight*(gd1err + gd2err);
