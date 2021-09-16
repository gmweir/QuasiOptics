function xout = centerpulse(xin)

n = length(xin);
x = abs(xin);
[xmax, kmax] = max(x);
dk = kmax - ceil(n/2);
xout = circshift(xin, dk);