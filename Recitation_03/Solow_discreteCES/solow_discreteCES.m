%%%%%%%%%%%%%
% solow_discreteCES.m
% Solves plain vanilla Solow model in discrete time with CES production
% function
% Stefano Pica
% Fall 2019
%%%%%%%%%%%%%

clear; close all; clc

% Model Parameters - annual calibration
s = 0.2;        %average share of real investment in real GDP (around 20%)
delta = 0.05;   %average ratio of depreciation to GDP (around 5%)
n = 0.02;       %population growth (around 2%)
gY = 0.02;      %output growth (around 2%)
sigma = 1.5;
alpha_K = 1/3;            %current level of capital share in the economy (around 33%)
alpha_L = 1-alpha_K;    %current level of labor share in the economy (around 66%)
A = (gY+delta+n) / (s * alpha_K ^ (1 / (sigma-1))); %match the 2% growth rate of output
L = 1;          %normalize labor to get capital per worker

% Plotting parameters
fsizenum = 14;
lwidnum = 2;

%time objects
T = 100;        %end of time
time = 1:T;     %grid of time periods

%preallocate vectors
k = NaN(1,T);       %capital per efficiency units
y = NaN(1,T);       %output per efficiency units
mpk = NaN(1,T);     %marginal product of capital
sharek = NaN(1,T);  %capital share

%first point in time
k(1) = 1; %initial condition on state variable
y(1) = prodfCES(k(1), L, A, sigma, alpha_K, alpha_L); %output
mpk(1) = mpkCES(k(1), L, A, sigma, alpha_K, alpha_L); %marginal product of capital
share(1) = mpk(1) * k(1) / y(1); %capital share

for i=2:T  %solve recursively law of motion
       k(i) = (1 / (1 + n)) * ((1 - delta) * k(i-1) + s * y(i-1)); 
       y(i) = prodfCES(k(i), L, A, sigma, alpha_K, alpha_L);
       mpk(i) = mpkCES(k(i), L, A, sigma, alpha_K, alpha_L); 
       sharek(i) = mpk(i) * k(i) / y(i);
end

figure % plot series of interest
subplot(2,2,1);
plot(time, sharek, 'b', 'LineWidth', lwidnum), title('Capital Share'), xlabel('Time');
set(gca,'FontSize',fsizenum);

subplot(2,2,2);
plot(time, k, 'b', 'LineWidth', lwidnum), title('Capital per worker'), xlabel('Time');
set(gca,'FontSize',fsizenum);

subplot(2,2,3);
plot(time, y, 'b', 'LineWidth', lwidnum), title('Output per worker'), xlabel('Time');
set(gca,'FontSize',fsizenum);

subplot(2,2,4);
plot(time, mpk,'b', 'LineWidth', lwidnum), title('Marginal product of capital'), xlabel('Time');
set(gca,'FontSize',fsizenum);
%print('solow_discrete','-djpeg')
