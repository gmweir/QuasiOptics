function [tiac, iac] = iacsim(tin, p)
% IACSIM
% t in fs, and p is field strength.

Nt = length(tin);
dt = tin(2) - tin(1);
t2 = 0:dt:(Nt/2-1)*dt;
t = [-fliplr(t2(2:end)) t2];

iac2 = zeros(1, Nt/2);
p2 = p.^2.^2;
p2sumfor = [0 cumsum(p2)];
p2sumback = [0 cumsum(fliplr(p2))];  % this could be derived from the above
for k=1:Nt/2,
  iac2(k) = p2sumfor(k) + sum((p(1:Nt-k+1)+p(k:Nt)).^2.^2) + p2sumback(k);
end
iac2 = iac2/(2*p2sumfor(end));  % cheat with normalization
iac = [fliplr(iac2(2:Nt/2)) iac2(1:Nt/2)];
tiac = t;
