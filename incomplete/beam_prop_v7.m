% jwo, 07/09/2014: function propagate with input q_in and output q_out
% jwo, 12/12/2015: make function propagate more transparent by accounting
% first for the new element and subsequently for the beam evolution.
function [] = beam_prop_v7()
% Brief description:
% *** Input: Excel sheet with data on system:
% - Input radius of curvature (R) and input beam radius (w)
% - Frequency f
% - Elements (mirrors) in the line with focal distance(s) and distance to next element
% *** Process: Gaussian beam parameter q_in is input to the element and gives q_pe: the
% Gaussian beam parameter immediately after the new element. Subsequently q_out is calculated
% as a function of distance over the trajectory by taking q_pe as input and
% m*stepsize as increment in distance.
% *** Output: vector S with propagation vector and horizontal and vertical w and R

close all
clear all

%% Import the data
[~, ~, raw] = xlsread('Mirror_data','notch_rig');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
raw = raw(:,[2,3,4,5]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
element = cellVectors(:,1);
fd_h = data(:,1);
fd_v = data(:,2);
L = data(:,3);

%% Extract beam parameters
f = fd_h(5)*1E9; % Frequency in Hz
win_h= fd_h(3); % Horizontal beam radius at input in m
Rin_h = fd_h(4);
win_v= fd_v(3); % Vertical beam radius at input in m
Rin_v = fd_v(4);
d0 = L(10); % distance 1-st element measured from z = 0
offset = fd_h(6);


%% Clear temporary variables
clearvars data raw cellVectors R;

%% Start of propagation code
% Constants
e0 = 8.854E-12; % permittivity of free space
u0 = 4*pi*1E-7; % permeability of free space
c = 1/sqrt(u0*e0);

% variabels
global ld S S_step n

% Derivatives
ld = c/f;

% Propagation vector
S_min = 0;
S_step = 1E-4;
% Determine S_max
fp = 11; % first pair 'element - distance' in prop. table

u = fp;
while not(strcmp(element(u), 'end')) % check for end 
    u = u+1;
end
le = u; % last element
S_max = L(u); % max of prop vector is total distance
S(:,1) = S_min:S_step:S_max; % propagation vector [m]
S(:,2) = zeros; % horizontal beam radius [m]
S(:,3) = zeros; % horizontal radius of curvature [m]
S(:,4) = zeros; % vertical beam radius [m]
S(:,5) = zeros; % vertical radius of curvature [m]


%% Compute horizontal trajectory
% calculate q_in and initialise 
q_inv = 1/Rin_h - j*ld/(pi*win_h^2); % Goldsmith e.q. 2.18
q_in = 1/q_inv;
n = 1;
S(n,2) = win_h;
S(n,3) = Rin_h;    
n = 2;

propagate(q_in, inf, d0, 2, 3) % trajectory from launch to first element
% f_c is flat in this case as only z is changed in the 1-st trajectory

u = fp;
  for i = 1:le-fp;
      propagate(q_in, fd_h(u), L(u), 2, 3) % propagate Element_n -> distance_n
      u = u+1;
  end
  
  %% Compute vertical trajectory
% calculate q_in and initialise 
q_inv = 1/Rin_v - j*ld/(pi*win_v^2); % Goldsmith e.q. 2.18
q_in = 1/q_inv;
n = 1;
S(n,4) = win_v;
S(n,5) = Rin_v;    
n = 2;

propagate(q_in, inf, d0, 4, 5) % trajectory from launch to first element
% f_c is flat in this case as only z is changed in the 1-st trajectory

u = fp;
  for i = 1:le-fp;
      propagate(q_in, fd_v(u), L(u), 4, 5) % propagate Element_n -> distance_n
      u = u+1;
  end

  
%% Function propagate
function [] = propagate(q_i, f_c, L_c, wn, Rn)
M_c = [   1     0;... % Ray transfer matrix for thin lens. If else then place outside loop
        -1/f_c  1];
q_pe = (M_c(1,1)*q_i+M_c(1,2))/(M_c(2,1)*q_i+M_c(2,2)); % Goldsmith e.q. 3.4
N = round(L_c/S_step); % number of points in this sub trajectory

for m = 1:N % Calculate beam evolution as function of distance with q_pe as input
    M_dS = [1  m*S_step;... 
            0  1];
    q_out = (M_dS(1,1)*q_pe+M_dS(1,2))/(M_dS(2,1)*q_pe+M_dS(2,2));
    S(n,wn)=(ld/(pi*imag(-1*(1/q_out))))^0.5; % Extract beam radius
    S(n,Rn) = 1/(real(1/q_out)); % Extract radius of curvature
    n = n + 1;
end

q_in = q_out; % the final q_out is the q_in of the next sub trajectory
end % function propagate

%% Plot

figure(1)
plot(S(:,1) + offset, S(:,2), S(:,1) + offset, S(:,4), 'LineWidth', 1.5)
grid on
set(gca,'FontSize',12,'LineWidth', 1.5) % paper 14 and 2 
xlabel('Propagation vector z [m]', 'fontsize', 12)
ylabel('Beam radius w [m]', 'fontsize', 12)
legend('w_{hor}','w_{ver}', 'location', 'NorthWest', 'Location', 'best')
title(['CTS gyrotron beam radii'], 'fontsize', 14)
set(gca,'FontSize',14,'LineWidth', 2)


y_max = 0.05;
ylim([0 y_max])
% xlim([0 3])

% Put in lines and text automatically
lo = 0.8*y_max; % line hight
to = 0; % horizontal text offset with respect to line
th = 0.85*y_max; % hight of text label


line([(0 + offset) (0 + offset)], [0 lo])
text((offset - to), th ,element(fp-1),'fontsize', 10, 'BackgroundColor',[1 1 1])
for i=1:u-fp
    line([(sum(L(fp-1:fp-2+i)) + offset) (sum(L(fp-1:fp-2+i)) + offset)],[0 lo]);
    text((sum(L(fp-1:fp-2+i))+ offset - to), th ,element(fp-1+i),'fontsize', 10, 'BackgroundColor',[1 1 1])
end


%{

figure(2)
plot(S(:,1) + offset, S(:,3), S(:,1) + offset, S(:,5), 'LineWidth', 1.5)
grid on
set(gca,'FontSize',12,'LineWidth', 1.5) % paper 14 and 2 
xlabel('Propagation vector z [m]', 'fontsize', 12)
ylabel('Radius of curvature [m]', 'fontsize', 12)
legend('R_{hor}','R_{ver}', 'location', 'NorthWest')
title(['Radii of curvature'], 'fontsize', 14)
set(gca,'FontSize',14,'LineWidth', 2)

%}
end