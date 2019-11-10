%% initial Parameter values (if you use shock.m, then the end values need to coincide)

mu=0; 

alpha=0.3;  % share of capital in production
            
theta =1;   % intertemporal elasticity of substitution euqal to 1. This is required for the simple code. Otherwise one needs to adjust things.

frish=0.5;

psi=1;

%% The following are YEARLY parameters. Here we interpret a period of duration 1 as corresponding to a year.
delta=0.05;         % annual depreciation rate.
nPop=0.01;          % annual rate of population growth.
gx=0.02;            % annual rate of GDP growth in advanced economies.
rho=0.04-gx/theta;  % discount rate required to match a 4% annual interest rate in steady state. 
                    % the interest rate is rho+g/theta
                    
%% parameters that only affect the level of variables:

A_neutral=2; % level of factor-neutral productivity.

ell=1;      % labor endowment per household. We can fix it at 1 if we want to think in terms of wages per person.
            % we could fix it at 1,200 if these are hours.

A_labor=1;  % initial level of A_L. Useful only if we are interested in the level of the variables per capita.
            % if I were to calibrate it, i would choose it to ensure a
            % level of wages close to the yearly wage per person. 