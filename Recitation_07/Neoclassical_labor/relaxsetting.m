%relaxsetting.m

% number of initial boundary conditions ( = number of state variables)
n1=1;

%number of differential equations (Note: number of final boundary
%conditions equals n-n1)
n=2;

% number of static equations to be solved simultanously
n3=1;

%number of mesh points
M=100;

%number of variables altogether
N=n+n3;         %(Do not change)

%Normalization of variables
normal=[];  %no normalisation

% Steady state values
phi=((delta+rho+gx/theta)/((1-mu)*alpha*A_neutral))^(1/(alpha-1));
lss=((1/psi)*( (1-mu)*(1-alpha)*A_neutral*phi^alpha/(A_neutral*phi^alpha-(delta+nPop+gx)*phi)))^(1/(1+frish));
kss=lss*phi;
css=A_neutral*kss^alpha*lss^(1-alpha)-(nPop+delta+gx)*kss;

%guess of final steady state values
y=ones(N,1);
y(1)=kss;    %k
y(2)=css;    %c
y(3)=lss;

%In case the shock consists of a reduction of state variables, enter this
%here

statev=0;                          %Option 1: set to zero if you are going to use initbound to feed inital
                                    %value for states
                                    
%statev=ones(N,1); statev(1)=1;      %Option 2: the initial value of capital is at 50% of its steady state value

%statev=ones(N,1);                  %Option 3: if you are interested in simulating a
                                    %MIT shock using shock.m.

%--------------------------------------------------------------------------
%Specification, which differential equations are used for constructing the
%final boundary conditions. 
Endcond=0;

tol=10^-9;      %tolerance for the Newton procedure
maxit=50;       %Maximum number of iterations
nu=0.05;        %Parameter for time transformation
damp=1;         %Dampening factor of the Newton procedure. The dampening factor will be 
dampfac=2;      %multiplied by the factor dampfac in every iteration until it equals 1

