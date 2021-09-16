function imagestack(ds)

ns = [repmat([1.6 1.9], 1, length(ds)) 1.4];
n0 = 1;

ni = 1280;
dbuf = 1;
n = length(ds);
dlocs = cumsum(ds);

x = [-dbuf, 0, dlocs, dlocs(end)+2*dbuf];
y = [n0, ns, ns(end)];

xi = linspace(min(x), max(x), ni);
yi = zeros(1,ni);

for k = 1:ni,
  kx = max(find(xi(k) >= x));
  yi(k) = y(kx);
end

imagesc(yi)
