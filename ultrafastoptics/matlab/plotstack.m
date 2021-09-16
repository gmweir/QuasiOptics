function plotstack(ds, ns)

n = length(ds);

locs = [0 cumsum(ds)];

x = zeros(1, 2*n);
y = zeros(1, 2*n);
for k = 1:n;
  x(2*k-1) = locs(k);
  x(2*k) = locs(k+1);
  y(2*k-1) = ns(1,k);
  y(2*k) = ns(1,k);
end

plot(x, y)