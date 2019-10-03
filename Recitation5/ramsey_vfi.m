%%%%%%%%%%%%%
% ramsey_vfi.m
% Solves the Ramsey model using VFI.
% Stefano Pica, TA for EC 702
% Fall 2019
%%%%%%%%%%%%%
close all; clear; clc;

%reference: Stokey_Lucas_1989 5.1
% model to be solved
% V(k) = max_{kp<=f(k)+(1-delta)k} [ u(c) + beta V(kprime) ]
% BC: c + kprime = f(k) + (1-delta)*k
% prod function: f(k)=k^alpha
% utility function: u(c) = c^(1-gamma) / (1-gamma)


%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the Ramsey problem - awesome!')
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vfmaxit = 1000; %max number of VF iterations
vftol = 1e-7; %tolerance on Bellman equation, in percentages
T=20; %final period for transition
transitspan=0:T-1; %time span for transition
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%model parameters
gamma=1.5;
delta=1;
alpha=0.5;
beta=0.9;

%exact SS quantities
kss=(alpha/(1/beta-1+delta))^(1/(1-alpha));
Rss=alpha*kss^(alpha-1)+(1-delta);

%define grid
kmin=kss*0.1;
kmax=kss*2;
knum=1000;
K0=linspace(kmin,kmax,knum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%compute return function (utility function). We use meshgrid
[kp,k]=meshgrid(K0);
Glow=zeros(knum, knum);
Gmax=k.^alpha+(1-delta)*k;
kp(kp>=Gmax | kp<Glow)=NaN; %restrict control variable to feasible set
utility=(k.^alpha+(1-delta)*k-kp).^(1-gamma)/(1-gamma);
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
Cpol = K0.^alpha + (1-delta)*K0 - Kpol;
R=alpha*K0.^(alpha-1)+(1-delta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM AN ARBITRARY INITIAL CONDITION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if you want to start in an arbitrary point that is not in the grid, can
%use the following code:
% kstart=0.5; %initial condition
% [~,kinit]=min(abs(K0-kstart)); %grid point closer to initial condition

kinit=K0(1); %initial condition. first point in the grid
ktransit=zeros(T,1); %preallocate the capital transition vector
ktransit(1)=kinit; %first entry of capital transition is initial condition
for it=2:T
    ktransit(it)=Kpol(K0==kinit); %get from the grid the initial value and apply the PF
    kinit=ktransit(it); %update initial capital value
end
Rtransit=alpha*ktransit.^(alpha-1)+(1-delta); %transition for interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure 
subplot(3,2,1)
plot(K0,V,'b','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Value');
title('Value Function')
set(gca,'FontSize',fsizenum)

subplot(3,2,2)
plot(K0,Cpol,'k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Consumption');
title('Consumption Policy')
set(gca,'FontSize',fsizenum)

subplot(3,2,3)
plot(K0,Kpol,'r',K0,K0,'k--', 'LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Capital Tomorrow')
title(' Capital Policy')
legend('Capital Policy','45 degree line','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(3,2,4)
plot(K0,R,'k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Interest Rate');
title('Interest Rate')
set(gca,'FontSize',fsizenum)

subplot(3,2,5)
plot(transitspan,ktransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[kss,kss],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Capital');
title('Transition for Capital')
legend('Capital','KSS','Location','SouthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(3,2,6)
plot(transitspan,Rtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[Rss,Rss],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Interest Rate');
title('Transition for Interest Rate')
legend('Interest Rate','RSS','Location','NorthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)
