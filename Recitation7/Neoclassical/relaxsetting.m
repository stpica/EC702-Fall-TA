%relaxsetting.m

% number of initial boundary conditions ( = number of state variables)
n1=1;

%number of differential equations (Note: number of final boundary
%conditions equals n-n1)
n=2;

% number of static equations to be solved simultanously
n3=0;

%number of mesh points
M=100;

%number of variables altogether
N=n+n3;         %(Do not change)

%Normalization of variables
normal=[];  %no normalisation

% Steady state values (mu different from zero)
kss=ell*((delta+rho+gx/theta)/((1-mu)*alpha*A_neutral))^(1/(alpha-1));
css=A_neutral*kss^alpha*ell^(1-alpha)-(nPop+delta+gx)*kss;

% Steady state values
kss_old=ell*((delta+rho+gx/theta)/(alpha*A_neutral))^(1/(alpha-1));
css_old=A_neutral*kss_old^alpha*ell^(1-alpha)-(nPop+delta+gx)*kss_old;

%guess of final steady state values
y=ones(N,1);
y(1)=kss;    %k
y(2)=css;    %c


%In case the shock consists of a reduction of state variables, enter this
%here

%Option 1: set to zero if you are going to use initbound to feed inital
%values for states:
statev=0;

%Option 2: the initial value of capital is at 50% of its steady state
%value:
%statev=ones(N,1); statev(1)=0.5;

%Option 3: if you are interested in simulating a MIT shock using shock.m:
%statev=ones(N,1);

%--------------------------------------------------------------------------
%Specification, which differential equations are used for constructing the
%final boundary conditions. 
Endcond=0;

tol=10^-9;      %tolerance for the Newton procedure
maxit=50;       %Maximum number of iterations
nu=0.05;        %Parameter for time transformation
damp=1;         %Dampening factor of the Newton procedure. The dampening factor will be 
dampfac=2;      %multiplied by the factor dampfac in every iteration until it equals 1

