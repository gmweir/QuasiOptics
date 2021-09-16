function n = nsell(y, S)
%NSELL  Calculate index from Sellmier coefficients
% n = nsell(y, S) is the index at wavelength y, given
% the Sellmier coefficients given in S, where the first column is
% the A coefficients (unitless) and the second are the B coefficients,
% in units of wavelength squared. S can be either a matrix representing a
% single material, or a cell array of matrices for multiple materials to
% evaluate in parallel. The wavelength vector ys MUST be a column vector.

% Initialization.
if iscell(S)
  Scell = S;
else
  Scell = {S};
end
nmat = length(Scell);
y2 = y.^2;

% Calculate index curve for each material, operating over all
% wavelengths in parallel.
n = zeros(length(y), nmat);
for sk = 1:nmat
  sellmat = Scell{sk};
  ns = size(sellmat, 1);
  terms = 0;
  for k = 1:ns,
    terms = terms + sellmat(k,1)*y2./(y2 - sellmat(k,2));
  end
  n(:,sk) = sqrt(1 + terms);
end
