%%%%%%%%%%%%%
% ramseylabor_vfi.m
% Solves the Ramsey model with endogenous labor using VFI.
% Stefano Pica, TA for EC 702
% Fall 2019
%%%%%%%%%%%%%
close all; clear; clc;

% model to be solved
% V(k) = max_{c,l,kp} [ u(c,l) + beta V(kprime) ]
% BC: c + kp = f(k,l) + (1-delta)*k
% prod function: f(k,l)=A*k^alpha*l^(1-alpha)
% utility function: u(c) = log(c)-psi*l^(1+eps)/(1+eps)

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the Ramsey problem with endogenous labor - awesome!')
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vfmaxit = 1000; %max number of VF iterations
vftol = 1e-7; %tolerance on Bellman equation, in percentages
T=100; %final period for transition
transitspan=0:T-1; %time span for transition
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%model parameters
alpha=0.6;
beta=0.95;
delta=0.2;
psi=1;
epsilon=0.3;
A=1;

%exact SS quantities
L_ss=((1/psi)*(1-alpha)/(1-alpha*delta/(1/beta-1+delta)))^(1/(epsilon+1)); %labor
K_ss=(A*alpha/(1/beta-1+delta))^(1/(1-alpha))*L_ss; %capital
R_ss=A*alpha.*K_ss.^(alpha-1).*L_ss.^(1-alpha)+(1-delta); %interest rate
W_ss=A*(1-alpha).*K_ss.^(alpha).*L_ss.^(-alpha); %wage
C_ss=A*K_ss^alpha*L_ss^(1-alpha)+(1-delta)*K_ss-K_ss; %consumption

%capital grid
kmin=0.8*K_ss;
kmax=1.2*K_ss;
knum=300;
K0=linspace(kmin, kmax, knum);

%labor grid
lmin=0.8*L_ss;
lmax=1.5*L_ss;
lnum=250;
L0=linspace(lmin,lmax,lnum);

[kp, k, l]=meshgrid(K0, K0, L0); %current k is row, tomorrow k column, employment is l var (third dimension)
% cutkp=kp(:,:,2);
% cutl=squeeze(l(2,:,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%feasible set for kprime:
Glow=zeros(knum, knum, lnum);
Gmax=A*k.^alpha.*l.^(1-alpha)+(1-delta)*k;
%kp(kp>Gmax | kp<Glow)=NaN; %restrict kp to feasible set
%cutkp=kp(:,:,1);

%Objective function:
matrix_uxy=log(A*k.^alpha.*l.^(1-alpha)+(1-delta)*k-kp)-psi*l.^(1+epsilon)/(1+epsilon);
matrix_uxy(kp>Gmax | kp<Glow)=NaN; %restrict to feasible set
cututils=matrix_uxy(:,:,1); %cut the utility to look at it

%initialize the policy and VF arrays
Vold = zeros(knum,1); %this will store real value of old VF
V = zeros(knum,1);  %this will store real value of new VF
kprimeind=zeros(knum,1);

disp('%%%%')
disp('Starting VFI now - it is actually happening!')
disp(' ')

%Useful to illustrate role of repmat:
% Vtest=[1:1:knum]'; %continuation values
% test=repmat(Vtest',[knum,1,lnum]);
% test(:,:,5) %along any z, we get copies of Vtest with same continuation values

%Now, there are two ways we can solve Ramsey with endogenous labor. This
%depends on the fact that the labor choice is a static one, which means we
%can maximize labor out before even starting the VFI. A second approach
%would be to maximize labor at each step of the iteration.


% METHOD 1: solve out optimal labor (max over z---dim 3) so to get the
% optimal utility matrix (matrix_uxy_Lsolved). Then use this matrix to
% solve as usual the VFI.
[matrix_uxy_Lsolved,lprimeind] = max(matrix_uxy,[],3);

for vfit = 1:vfmaxit
    
    %METHOD 2: uncomment the following lines.
    %form Bellman equation (this solves for L at each step. Might be useful in some applications, but not here)   
    %RHSMAT = matrix_uxy + beta*repmat(Vold',[knum,1,lnum]); % the 1 is the dimension corresponding to the kp variable
   	%[RHSMAT_Lsolved,lprimeind] = max(RHSMAT,[],3); %solve out optimal labor (max over z---dim 3)
    %[V,kprimeind] = max(RHSMAT_Lsolved,[],2); % solve out capital next period (max over column---dim 2)
    
    RHSMAT = matrix_uxy_Lsolved + beta*repmat(Vold',knum,1);
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

ValueFunction = V; %value function
Kpol = K0(kprimeind); %capital policy function
Lpol = NaN(1,knum); %preallocate labor policy function
for i=1:knum %labor policy function
Lpol(i) = L0(lprimeind(i,kprimeind(i)));
end
Cpol = A*K0.^alpha.*Lpol.^(1-alpha)+(1-delta)*K0-Kpol;
R = A*alpha.*K0.^(alpha-1).*Lpol.^(1-alpha)+(1-delta); %interest rate
W = A*(1-alpha).*K0.^(alpha).*Lpol.^(-alpha); %wage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM AN ARBITRARY INITIAL CONDITION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if you want to start in an arbitrary point that is not in the grid, can
%use the following code:
% kstart=0.5; %initial condition
% [~,kinit]=min(abs(K0-kstart)); %grid point closer to initial condition

kinit=K0(1); %initial condition. first point in the grid
ktransit=zeros(1,T); %preallocate the capital transition vector
ltransit=zeros(1,T); %preallocate the labor transition vector
ctransit=zeros(1,T); %preallocate the consumption transition vector
ktransit(1)=kinit; %first entry of capital transition is initial condition
ltransit(1)=Lpol(K0==kinit); %first entry of capital transition is initial condition
ctransit(1)=Cpol(K0==kinit);
for it=2:T
    ktransit(it)=Kpol(K0==kinit); %get from the grid the initial value and apply the PF
    ltransit(it)=Lpol(K0==ktransit(it));
    ctransit(it)=Cpol(K0==ktransit(it));
    kinit=ktransit(it); %update initial capital value
end
Rtransit = A*alpha.*ktransit.^(alpha-1).*ltransit.^(1-alpha)+(1-delta); %transition for interest rate
Wtransit = A*(1-alpha).*ktransit.^(alpha).*ltransit.^(-alpha); %transition for wage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure % value and policy functions
subplot(2,3,1)
plot(K0,ValueFunction,'b','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Value')
title('Value Function')
set(gca,'FontSize',fsizenum)

subplot(2,3,2)
plot(K0,Cpol,'k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Consumption');
title('Consumption Policy')
set(gca,'FontSize',fsizenum)

subplot(2,3,3)
plot(K0,Kpol,'r',K0,K0,'k--', 'LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Capital Tomorrow')
title(' Capital Policy')
legend('Capital Policy','45 degree line','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,4)
plot(K0,Lpol,'r',K0,ones(knum,1)*L_ss,'b--','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Labor')
title('Labor Policy')
legend('Labor Policy','SS Labor','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,5)
plot(K0,R,'k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Interest Rate')
title('Interest Rate')
set(gca,'FontSize',fsizenum)

subplot(2,3,6)
plot(K0,W,'k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Wage')
title('Wage')
set(gca,'FontSize',fsizenum)

figure %transition from arbitrary initial point
subplot(2,3,1)
plot(transitspan,ktransit,'b--', transitspan, ones(1,T)*K_ss, 'r--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Capital')
title('Transition for Capital')
legend('Capital','SS Capital','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,2)
plot(transitspan,ltransit,'b--',transitspan,ones(1,T)*L_ss,'r--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Labor')
title('Transition for Labor')
legend('Labor','SS Labor','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,3)
plot(transitspan,ctransit,'b--',transitspan,ones(1,T)*C_ss,'r--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Consumption')
title('Transition for Consumption')
legend('Consumption','SS Consumption','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,4)
plot(transitspan,Wtransit,'b--',transitspan,ones(1,T)*W_ss,'r--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Labor')
title('Transition for Wage')
legend('Wage','SS Wage','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,5)
plot(transitspan,Rtransit,'b--',transitspan,ones(1,T)*R_ss,'r--','LineWidth',lwidnum)
xlabel('Time')
ylabel('Interest Rate')
title('Transition for Interest Rate')
legend('Interest Rate','RSS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)
