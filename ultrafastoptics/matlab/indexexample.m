S = readmatstruct;
ys = linspace(.4, 1.2, 256);
gdd = gddsell(ys, 1e3, S.Nb2O5_AT2);
n = nsell(ys, S.Nb2O5_AT2);
plot(ys, gdd)