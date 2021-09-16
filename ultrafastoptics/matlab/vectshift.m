function y = vectshift(x, dk, n)

k1 = mod(dk, n) + 1;
y = [x(k1:n), x(1:k1-1)];
