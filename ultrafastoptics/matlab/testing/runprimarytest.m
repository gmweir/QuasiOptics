% this runs the test case derived analytically in the Mathematica file that
% should be included in this directory, and establishes the validity of
% stackgdd.m and stackgdgrad.m

s = readmatstruct;
ds = [0.2 0.3];
ys = linspace(0.4,1.2,256).';
ys0 = 0.7;

fprintf('\nrunprimarytest.m: two layers\n')

% [Rs, gds, Rgrads, gdgrads] = stackgdgradez(ys, ds, s.Nb2O5_MLD, s.SiO2_MLD, ...
%   s.FS, 30/180*pi, 'TM');
% figure(1)
% plot(ys, Rs)
% axis([ys(1) ys(end) 0 1])
% figure(2)
% plot(ys, gds)

[R, gd, gdd] = stackgddez(ys0, ds, s.Nb2O5_MLD, s.SiO2_MLD, ...
  s.FS, 30/180*pi, 'TM');
fprintf('stackgdd: R^2 = %g, gd = %g\n', abs(R)^2, gd)

[r, gd, rgrad, gdgrad] = stackgdgradez(ys0, ds, s.Nb2O5_MLD, s.SiO2_MLD, ...
  s.FS, 30/180*pi, 'TM');
fprintf('stackgdgrad: R^2 = %g, gd = %g\n', abs(R)^2, gd)
fprintf('stackgdgrad: Rgrad = %g\n', rgrad)
fprintf('stackgdgrad: gdgrad = %g\n', gdgrad)

[r2, gd2, rgrad2, gdgrad2] = stackgdgradez(ys0, ds + [zeros(1,length(ds)-1) 1e-9], s.Nb2O5_MLD, s.SiO2_MLD, ...
  s.FS, 30/180*pi, 'TM');
drapp = (r2 - r)/1e-9
dgdapp = (gd2 - gd)/1e-9
