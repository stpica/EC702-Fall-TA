%%%%%%%%%%%%%
% cake_problem.m
% Solves a cake eating problem using VFI.
% Stefano Pica, TA for EC 702
% Fall 2019
%%%%%%%%%%%%%
close all; clear; clc;

%The cake eating problem is a special case of the Ramsey model, where
%depreciation is one and production function is linear: f(k)=k.
%In this code we make use of the  analytical solution of the cake eating
%problem to construct a more clever grid.

% model to be solved
% V(k) = max_{kp<=k} [ u(c) + beta V(kprime) ]
% BC: c + kprime = k
% utility function: u(c) = c^(1-gamma) / (1-gamma)

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the cake eating problem - awesome!')
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vfmaxit = 1000; %max number of VF iterations
vftol = 1e-7; %tolerance on Bellman equation, in percentages
T=200; %final period for transition
transitspan=0:T-1; %time span for transition
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%model parameters
gamma=1.5;
beta=0.9;

%Careful choice of the state space
knum=100;
num=knum:-1:1;
K0=beta.^(num/gamma); %make use of policy: k' = beta^(1/gamma) k
%K0=linspace(0.00001,1,knum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%compute return function (utility function). We use meshgrid
[kp,k]=meshgrid(K0);
kp=0.99*kp; %this is to avoid undefined utility for the first state
Glow=zeros(knum, knum);
Gmax=k;
kp(kp>=Gmax | kp<Glow)=NaN; %restrict control variable to feasible set
utility=(k-kp).^(1-gamma)/(1-gamma);
%notice property of utility matrix. Rows represent the state (capital
%today). columns are the control (capital tomorrow).



%initialize the policy and VF arrays
Vold = zeros(knum,1); %this will store real value of old VF
V = zeros(knum,1);  %this will store real value of new VF
kprimeind=zeros(knum,1);

disp('%%%%')
disp('Starting VFI now - it is actually happening!')
disp(' ')

for vfit = 1:vfmaxit
    %form Bellman equation
    RHSMAT = utility + beta*repmat(Vold',knum,1);
   	[V,kprimeind] = max(RHSMAT,[],2);
    absdiff=abs(V-Vold);
    vferr = max(absdiff); %maximum absolute error
    if (mod(vfit,50)==1)
        disp(['VF error = ' num2str(vferr) ' on VF iteration ' num2str(vfit) '.'])
    end
       %exit if converged
    if (vferr<vftol)
        break; 
    end
    %if not converged, then Vold <- V, and redo
    Vold = V;
end

disp(' ')
disp('VF iteration complete.')
toc;
disp('%%%%')
disp(' ')


Kpol = K0(kprimeind); %real values of capital policy function
Cpol = K0-Kpol; %consumption policy function
Kanalytical = beta^(1/gamma)*K0; %exact policy function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM AN ARBITRARY INITIAL CONDITION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if you want to start in an arbitrary point that is not in the grid, can
%use the following code:
% kstart=0.5; %initial condition
% [~,kinit]=min(abs(K0-kstart)); %grid point closer to initial condition

kinit=K0(knum); %initial condition. last point in the grid
ktransit=zeros(1,T); %preallocate the capital transition vector
ctransit=zeros(1,T); %preallocate the consumption transition vector
ktransit(1)=kinit; %first entry of capital transition is initial condition
for it=2:T
    ktransit(it)=Kpol(K0==kinit); %get from the grid the initial value and apply the PF
    ctransit(it)=kinit-ktransit(it);  %consumption transition
    kinit=ktransit(it); %update initial capital value
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
subplot(2,2,1)
plot(K0,V,'b','LineWidth',lwidnum)
xlabel('Cake Today')
ylabel('Value');
title('Value Function')
set(gca,'FontSize',fsizenum)

subplot(2,2,2)
plot(K0,Kpol,'r--',K0,Kanalytical,'b-.',K0,K0,'k--','LineWidth',lwidnum)
xlabel('Cake Today')
ylabel('Cake Tomorrow')
title('Accumulation Policy')
legend('PF computational','PF analytical','45 degree line','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)
%notice that computational and analytical policy functions overlap!

subplot(2,2,3)
plot(transitspan,ktransit,'b--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Cake');
title('Transition for Cake')
set(gca,'FontSize',fsizenum)

subplot(2,2,4)
plot(transitspan(2:end),ctransit(2:end),'b--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Cake Consumption');
title('Transition for Consumption')
set(gca,'FontSize',fsizenum)
