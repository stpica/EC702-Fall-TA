% Basic RBC Model.
% Discussion session #10 EC702B (Professor: Stephen Terry)
% Stefano Pica, TA. stpica@bu.edu
% Fall 2019

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k inv n a; % I reduce the system to 6 endogenous variables, the other 3 variables are implicit in the system
varexo e; %exogenous variables that will be shocked. Here, only productivity shock

parameters beta psi delta alpha rho;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

alpha   = 0.33;
beta    = 0.99;
delta   = 0.015;
psi     = 1.67;
rho     = 0.98;  
sigma   = 0.0072;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

%Specific functional forms assumed for the RBC framework:
%Utility function: u(c,n) = log(C) + psi*log(1-n)


%Note:In Dynare, the timing of each variable reflects when that variable
%is decided. Hence k_t becomes k(-1)
%Dynare will do linear approximation of the levels of the variables, so we
%want to consider exp(x) instead of x.

model; %Here we insert all the equations of the model
(1/exp(c)) = beta*(1/exp(c(+1)))*(1-delta+alpha*exp(a(+1))*(exp(k)^(alpha-1))*(exp(n(+1)))^(1-alpha)); %Euler equation
psi*exp(c)/(1-exp(n)) = (1-alpha)*exp(a)*(exp(k(-1))^alpha)*(exp(n)^(-alpha)); %Intratemporal optimality condition. LHS=wage
exp(c)+exp(inv) = exp(y); %resource constraint
exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(n))^(1-alpha); %production function
exp(inv) = exp(k)-(1-delta)*exp(k(-1)); %investment definition
a = rho*a(-1)+e; %productivity evolution
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

% Now, Dynare needs to know the SS of our model. We can either enter exact 
%SS values or just approximations and let Dynare find the exact SS (which 
%it will do using numerical methods based on these approximations).
%The following are approximations

initval;
  k = log(9);
  c = log(0.76);
  n = log(0.3);
  a = log(1); 
  e = 0;
end;

% Shocks are only allowed to be temporary. A permanent shock cannot be
%accommodated due to the need to stationarize the model around a SS.
%Furthermore, shocks can only hit the system today, as the expectation of
%future shocks must be zero.

shocks;
var e = sigma^2;
end;

% Adding "steady" just after your initval block will instruct Dynare to 
%consider your initial values as mere approximations and start simulations
%or impulse response functions from the exact steady state. On the
%contrary, if you don't add the command "steady", your simulations or
%impulse response functions will start from your initial values, even if
%Dynare will have calculated your model's exact steady state for the
%purpose of linearization.
% We add "steady".

steady;


% The "stoch_simul" command instructs Dynare to compute a Taylor
%approximation of the decision and transition functions for the model
%(the equations listing current values of the endogenous variables of the
%model as a function of the previous state of the model and current
%shocks), impulse response functions and various descriptive statistics
%(moments, variance decomposition, correlation and autocorrelation coefficients).

% Note that in the case of a second order approximation, Dynare will return
%the actual sample moments from the simulations. For first order
%linearizations, Dynare will instead report theoretical moments, unless
%we impose to run simulations (as we do as follows). In both cases, the
%return to steady state is asymptotic.

% We want to have a first order approximation, we want to look at
%300 periods in the IRF responses, and we want to simulate our series with T=10000.

stoch_simul(order = 1, irf=300, periods=10000);
