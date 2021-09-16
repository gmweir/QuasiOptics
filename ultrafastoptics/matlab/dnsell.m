function [n, dn, ddn] = dnsell(y, S)
%DNSELL  Calculate index and derivatives from Sellmeir coefficients
% [n, dn, ddn] = dnsell(y, S) is d^(k) n(y) for k = 0...2 at wavelength y, given
% the Sellmier coefficients given in matrix S, where the first column is
% the A coefficients and the second are the B coefficients.
% FIXME: Apparently the second derivatives aren't right?

% Initialization.
[ns, n] = size(S);
y2 = y.^2;

% Calculate terms over all wavelengths in parallel.
terms = 0; dterms = 0; ddterms = 0;
for k = 1:ns,
  A = S(k,1);
  if (A ~= 0)
    B = S(k,2);
    terms = terms + A*y2./(y2 - B);
    dterms = dterms - 2*A*B*y./(y2 - B).^2;
    ddterms = ddterms + 2*A*B*(3*y2 + B)./(y2 - B).^3;
  end
end

% Final expressions.
n = sqrt(1 + terms);
dn = dterms./(2*n);
ddn = (2*n.^2.*ddterms - dterms.^2)./(4*n.^3);
%ddn = (ddterms - 2*dn.^2)./(2*n);