% Calculate multi disk response, book Hartfuss & Geist 7.3.3
% Nickel and Thum (Rho and S-matrix)
% Hartfuss & Geist (T-matrices method for cascading)
% Data on TE11 waveguide: Moreno p. 117
% REQUIRES function: S_cir -> computes progation constant and S-parameters.
% REQUIRES Excel input file: sample_data.xlsx
% jwo, 20/12/2015 read waveguide radius from Excel sheet
% jwo, 24.07.2017, put into Git (W7-X)


clear all
close all

%% Import data from spreadsheet (MatLab auto generation)
[~, ~, raw] = xlsread('sample_data.xlsx','Neha');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,5);
raw = raw(:,[1,2,3,4,6,7,8,9]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));


% Allocate imported array to column variable names
element = data(:,1);
d = data(:,2);
er_p = data(:,3);
td = data(:,4);
comment = cellVectors(:,1);
f_start = data(1,5);
f_step = data(1,6);
f_stop = data(1,7);
r_mm = data(1,8);

% Clear temporary variables
clearvars data raw cellVectors R;

%% Start code jwo

global x_p f

e0 = 8.854E-12; u0 = 4*pi*1E-7;
c = 1/sqrt(u0*e0);

% Circular waveguide
r = r_mm/1000;
x_p = 1.8412; % u prime for TEmn with m = 1, n = 1 (TE11).


%% Body

% Convert to m, Hz
d = d*1E-3;
f_start = f_start(1)*1E9;
f_step  = f_step(1)*1E6;
f_stop  = f_stop(1)*1E9;

%Pre allocate
N_points = int32((f_stop - f_start)/f_step);
S_11(1:N_points) = 0;
S_21(1:N_points) = 0;
f_v(1:N_points) = 0;


%% Compute S11 and S22 for each frequency
m = size(element); % an element is a disk or a spacing between disks
Sm = zeros(2,2,m(1)); % 2 x 2 matrices of each element
Tm = zeros(2,2,m(1));
f = f_start;
count =1;
while f < f_stop % compute cascaded S-matrix for each frequency
  for n = 1 : m(1)
    Sm(:,:,n) =  S_cir(er_p(n), td(n), d(n), r); % Scattering matrix each element n
    Tm(:,:,n) = s2t(Sm(:,:,n)); % Correponding T-matrices
  end
  
  Tm_cas = Tm(:,:,m(1)); % Compute cascaded T-matrix = T1*T2*T3 etc. (note to compute right to left)
  for n = 1 : m(1)-1
      Tm_cas = Tm(:,:,m(1)-n)*Tm_cas;
  end
  
  Sm_cas = t2s(Tm_cas); % Cascaded S-matrix for each specific frequency
  
  S_21(count) = Sm_cas(2,1);
  S_11(count) = Sm_cas(1,1);

  f_v(count) = f;
  f = f + f_step;
  count = count+1;
end


%% Process results
% Reflection, Transmission, Absorption   
T = abs(S_21).^2; % _0: parallel plates
R = abs(S_11).^2;
A = 1-T-R;

%% Approximations
A_SPE = 2*(1-exp(-2*pi.*f_v*sqrt(er_p(1))*td(1)*d(1)/c)); % careful: 2 disks
A_SP = 2*2*pi.*f_v*sqrt(er_p(1))*td(1)*d(1)/c; % careful: 2 disks
A_res = 2*pi.*f_v*(er_p(1)+1)*td(1)*d(1)/c; % careful: 2 disks (single disk is pi*.*f etc.)



%% plot
f_v=f_v./1E9;


figure(1)
plot(f_v, 10*log10(T), 'LineWidth', 1.2);
hold on
plot(f_v, 10*log10(R), 'LineWidth', 1.2);
plot(f_v, 10*log10(A), 'LineWidth', 1.2);
hold off
grid on
xlabel('Frequency [GHz]','fontsize', 12)
ylabel('Fraction of power [dB]','fontsize', 12)
legend('Transmission', 'Reflection','Absorption','Location', 'SouthEast')
[d_str, errmsg] = sprintf('%4.2f',1000*d(1));
if m(1)>1
[spac_str, errmsg] = sprintf('%3.1f',1000*d(2));
else
[spac_str, errmsg] = sprintf('%3.1f',0);
end
[er_p_str, errmsg] = sprintf('%4.2f',er_p(1));
[td_str, errmsg] = sprintf('%2.1e',td(1));
[r_str, errmsg] = sprintf('%3.1f',1000*r);
message_str = ['d = ', d_str, ' mm, spacing ',spac_str,' mm, ','er = ',er_p_str, ', tan delta = ',td_str, ', r = ',r_str, ' mm'];
title(message_str,'fontsize', 14)
set(gca,'FontSize',14,'LineWidth', 1.5)
% print(figure(1), 'DD_CTS','-depsc')
