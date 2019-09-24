%%%%%%%%%%%%%
% solow_continuous.m
% Solves plain vanilla Solow model in continuous time (see example 2 in discussion 2)
% Stefano Pica
% Fall 2019
%%%%%%%%%%%%%

clear; close all; clc

% Model Parameters - annual calibration
s = 0.2;        %average share of real investment in real GDP (around 20%)
delta = 0.05;   %average ratio of depreciation to GDP (around 5%)
n = 0.02;       %population growth (around 2%)
g = 0.02;       %technology growth
alpha = 1/3;    %current level of capital share in the economy (around 33%)

% Plotting parameters
fsizenum = 14;
lwidnum = 2;

%SS quantities: analytical solution
k_ss = (s/(g+n+delta))^(1/(1-alpha)); %SS capital per efficiency units

%set up and solve differential equation
k0 = k_ss - 1; %initial condition - start not too far from SS
tspan = [0 100]; %integrates differential equations from 0 to 100
k_dot = @(t, k) s * k ^ alpha - (g + n + delta) * k; %differential equation
[time, k] = ode45(k_dot, tspan, k0); %solve it with ode45

figure % plot series of interest
plot(time, k, 'b', time, k_ss*ones(length(time)),'r', 'LineWidth', lwidnum), title('Evolution of Capital starting at k0');
xlabel('Time'), set(gca,'FontSize',fsizenum);
legend('Capital','Steady State','Location','Southeast'), legend boxoff;

%print('solow_continuous','-djpeg')
