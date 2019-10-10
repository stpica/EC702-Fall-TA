%%%%%%%%%%%%%
% ramseylabor_vfi.m
% Solves the Ramsey model using VFI.
% Stefano Pica, TA for EC 702
% 10/11/2018
%%%%%%%%%%%%%
close all; clear; clc;

% model to be solved
% V(k) = max_{c,l,kp} [ u(c,l) + beta V(kprime) ]
% BC: c + kp = f(k,l) + (1-delta)*k
% prod function: f(k,l)=A]k^alpha*l^(1-alpha)
% utility function: u(c) = log(c)-psi*l^(1+eps)/(1+eps)

%here we solve twice the Ramsey model with endogenous labor of choice, each
%time with a different alpha (the share of capital). Then we look at the 
%transition from a steady state to another.

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the Ramsey problem with endogenous labor - awesome!')
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vfmaxit = 2000; %max number of VF iterations
vftol = 1e-7; %tolerance on Bellman equation, in percentages
T=200; %final period for transition
transitspan=0:T-1; %time span for transition
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%model parameters
alpha=0.6;
alphap=0.7; %go from 0.6 to 0.7
beta=0.95;
delta=0.2;
psi=1;
epsilon=0.3;
A=1;

%exact SS quantities - lower alpha
L_ss1=((1/psi)*(1-alpha)/(1-alpha*delta/(1/beta-1+delta)))^(1/(epsilon+1)); %labor
K_ss1=(A*alpha/(1/beta-1+delta))^(1/(1-alpha))*L_ss1; %capital
R_ss1=A*alpha.*K_ss1.^(alpha-1).*L_ss1.^(1-alpha)+(1-delta); %interest rate
W_ss1=A*(1-alpha).*K_ss1.^(alpha).*L_ss1.^(-alpha); %wage
Y_ss1=A*K_ss1^alpha*L_ss1^(1-alpha); %output
C_ss1=Y_ss1+(1-delta)*K_ss1-K_ss1; %consumption

%exact SS quantities - bigger alpha
L_ss2=((1/psi)*(1-alphap)/(1-alphap*delta/(1/beta-1+delta)))^(1/(epsilon+1)); %labor
K_ss2=(A*alphap/(1/beta-1+delta))^(1/(1-alphap))*L_ss2; %capital
R_ss2=A*alphap.*K_ss2.^(alphap-1).*L_ss2.^(1-alphap)+(1-delta); %interest rate
W_ss2=A*(1-alphap).*K_ss2.^(alphap).*L_ss2.^(-alphap); %wage
Y_ss2=A*K_ss2^alphap*L_ss2^(1-alphap); %output
C_ss2=Y_ss2+(1-delta)*K_ss2-K_ss2; %consumption

%capital grid
kmin=0.8*K_ss1;
kmax=1.2*K_ss2;
knum=300;
K0=linspace(kmin, kmax, knum);

%labor grid
lmin=0.6*L_ss1;
lmax=1.4*L_ss2;
lnum=300;
L0=linspace(lmin,lmax,lnum)';
[kp, k, l]=meshgrid(K0, K0, L0); %current k is row, tomorrow k column, employment is l var (third dimension)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%feasible set for kprime:
Glow=zeros(knum, knum, lnum);
Gmax1=A*k.^alpha.*l.^(1-alpha)+(1-delta)*k;
Gmax2=A*k.^alphap.*l.^(1-alphap)+(1-delta)*k;
%kp(kp>Gmax | kp<Glow)=NaN; %restrict kp to feasible set
%cutkp=kp(:,:,1);

%Objective function:
matrix_uxy1=log(A*k.^alpha.*l.^(1-alpha)+(1-delta)*k-kp)-psi*l.^(1+epsilon)/(1+epsilon);
matrix_uxy2=log(A*k.^alphap.*l.^(1-alphap)+(1-delta)*k-kp)-psi*l.^(1+epsilon)/(1+epsilon);
matrix_uxy1(kp>Gmax1 | kp<Glow)=NaN; %restrict to feasible set
matrix_uxy2(kp>Gmax2 | kp<Glow)=NaN; %restrict to feasible set


%cututils=matrix_uxy(:,:,.5*lnum); %cut the utility to look at it

%initialize the policy and VF arrays
Vold1 = zeros(knum,1); %this will store real value of old VF
V1 = zeros(knum,1);  %this will store real value of new VF
kprimeind1=zeros(knum,1);
Vold2 = zeros(knum,1); %this will store real value of old VF
V2 = zeros(knum,1);  %this will store real value of new VF
kprimeind2=zeros(knum,1);

disp('%%%%')
disp('Starting VFI now - it is actually happening!')
disp(' ')

%Now, there are two ways we can solve Ramsey with endogenous labor. This
%depends on the fact that the labor choice is a static one, which means we
%can maximize labor out before even starting the VFI. A second approach
%would be to maximize labor at each step of the iteration.


% METHOD 1: solve out optimal labor (max over z---dim 3) so to get the
% optimal utility matrix (matrix_uxy_Lsolved). Then use this matrix to
% solve as usual the VFI.
[matrix_uxy_Lsolved1,lprimeind1] = max(matrix_uxy1,[],3);
[matrix_uxy_Lsolved2,lprimeind2] = max(matrix_uxy2,[],3);


for vfit = 1:vfmaxit
    
    RHSMAT1 = matrix_uxy_Lsolved1 + beta*repmat(Vold1',knum,1);
   	[V1,kprimeind1] = max(RHSMAT1,[],2);
    
    absdiff=abs(V1-Vold1);
    vferr = max(absdiff); %maximum absolute error
    
    if (mod(vfit,50)==1)
        disp(['VF error = ' num2str(vferr) ' on VF iteration ' num2str(vfit) '.'])
    end
       %exit if converged
    if (vferr<vftol)
        break; 
    end
    %if not converged, then Vold <- V, and redo
    Vold1 = V1;
end

for vfit = 1:vfmaxit
    
    RHSMAT2 = matrix_uxy_Lsolved2 + beta*repmat(Vold2',knum,1);
   	[V2,kprimeind2] = max(RHSMAT2,[],2);
    
    absdiff=abs(V2-Vold2);
    vferr = max(absdiff); %maximum absolute error
    
    if (mod(vfit,50)==1)
        disp(['VF error = ' num2str(vferr) ' on VF iteration ' num2str(vfit) '.'])
    end
       %exit if converged
    if (vferr<vftol)
        break; 
    end
    %if not converged, then Vold <- V, and redo
    Vold2 = V2;
end

disp(' ')
disp('VF iterations complete.')
toc;
disp('%%%%')
disp(' ')

ValueFunction1 = V1; %value function
Kpol1 = K0(kprimeind1); %capital policy function
Lpol1 = NaN(1,knum);
for i=1:knum %labor policy function
Lpol1(i) = L0(lprimeind1(i,kprimeind1(i)));
end
Cpol1 = A*K0.^alpha.*Lpol1.^(1-alpha)+(1-delta)*K0-Kpol1;
R1 = A*alpha.*K0.^(alpha-1).*Lpol1.^(1-alpha)+(1-delta); %interest rate
W1 = A*(1-alpha).*K0.^(alpha).*Lpol1.^(-alpha); %wage
Y1 = A*K0.^alpha.*Lpol1.^(1-alpha); %output

ValueFunction2 = V2; %value function
Kpol2 = K0(kprimeind2); %capital policy function
Lpol2 = NaN(1,knum);
for i=1:knum %labor policy function
Lpol2(i) = L0(lprimeind2(i,kprimeind2(i)));
end
Cpol2 = A*K0.^alphap.*Lpol2.^(1-alphap)+(1-delta)*K0-Kpol2;
R2 = A*alphap.*K0.^(alphap-1).*Lpol2.^(1-alphap)+(1-delta); %interest rate
W2 = A*(1-alphap).*K0.^(alphap).*Lpol2.^(-alphap); %wage
Y2 = A*K0.^alphap.*Lpol2.^(1-alphap); %output


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
ytransit=zeros(1,T); %preallocate the consumption transition vector
Rtransit=zeros(1,T); %preallocate the consumption transition vector
Wtransit=zeros(1,T); %preallocate the consumption transition vector
ktransit(1)=kinit; %first entry of capital transition is initial condition
ltransit(1)=Lpol1(K0==kinit); %first entry of capital transition is initial condition
ctransit(1)=Cpol1(K0==kinit);
for it=2:T/2 %from initial condition to high alpha SS
    ktransit(it)=Kpol1(K0==kinit); %get from the grid the initial value and apply the PF
    ltransit(it)=Lpol1(K0==ktransit(it));
    ctransit(it)=Cpol1(K0==ktransit(it));
    ytransit(it)=Y1(K0==ktransit(it));    
    Rtransit(it) = R1(K0==ktransit(it));
    Wtransit(it) = W1(K0==ktransit(it));
    kinit=ktransit(it); %update initial capital value
end
for it=T/2+1:T %from high alpha SS to low alpha SS
    ktransit(it)=Kpol2(K0==kinit); %get from the grid the initial value and apply the PF
    ltransit(it)=Lpol2(K0==ktransit(it));
    ctransit(it)=Cpol2(K0==ktransit(it));
    ytransit(it)=Y2(K0==ktransit(it));    
    Rtransit(it) = R2(K0==ktransit(it));
    Wtransit(it) = W2(K0==ktransit(it));
    kinit=ktransit(it); %update initial capital value
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure % value and policy functions
subplot(2,3,1)
plot(K0,ValueFunction1,'b',K0,ValueFunction2,'r','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Value')
title('Value Functions')
legend('Low alpha','High alpha','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,2)
plot(K0,Cpol1,'b',K0,Cpol2,'r','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Consumption');
title('Consumption Policy')
legend('Low alpha','High alpha','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,3)
plot(K0,Kpol1,'b',K0,Kpol2,'r',K0,K0,'k--', 'LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Capital Tomorrow')
title(' Capital Policy')
legend('Low alpha','High alpha','45 degree line','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,4)
plot(K0,Lpol1,'b',K0,Lpol2,'r','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Labor')
title('Labor Policy')
legend('Low alpha','High alpha','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,5)
plot(K0,R1,'b',K0,R2,'r','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Interest Rate')
title('Interest Rate')
legend('Low alpha','High alpha','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,6)
plot(K0,W1,'b',K0,W2,'r','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Wage')
title('Wage')
legend('Low alpha','High alpha','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

figure %transition from arbitrary initial point, then change in SS
subplot(2,3,1)
plot(transitspan,ktransit,'k--',transitspan,ones(1,T)*K_ss1,'b',transitspan,ones(1,T)*K_ss2,'r','LineWidth',lwidnum)
xlabel('Time')
ylabel('Capital')
title('Transition for Capital')
legend('Capital','Low alpha SS','High alpha SS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,2)
plot(transitspan,ltransit,'k--',transitspan,ones(1,T)*L_ss1,'b',transitspan,ones(1,T)*L_ss2,'r','LineWidth',lwidnum)
xlabel('Time')
ylabel('Labor')
title('Transition for Labor')
legend('Labor','Low alpha SS','High alpha SS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,3)
plot(transitspan,ctransit,'k--',transitspan,ones(1,T)*C_ss1,'b',transitspan,ones(1,T)*C_ss2,'r','LineWidth',lwidnum)
xlabel('Time')
ylabel('Consumption')
title('Transition for Consumption')
legend('Consumption','Low alpha SS','High alpha SS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,4)
plot(transitspan,Wtransit,'k--',transitspan,ones(1,T)*W_ss1,'b',transitspan,ones(1,T)*W_ss2,'r','LineWidth',lwidnum)
xlabel('Time')
ylabel('Wage')
title('Transition for Wage')
legend('Wage','Low alpha SS','High alpha SS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,5)
plot(transitspan,Rtransit,'k--',transitspan,ones(1,T)*R_ss1,'b',transitspan,ones(1,T)*R_ss2,'r','LineWidth',lwidnum)
xlabel('Time')
ylabel('Interest Rate')
title('Transition for Interest Rate')
legend('Interest Rate','Low alpha SS','High alpha SS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,3,6)
plot(transitspan,ytransit,'k--',transitspan,ones(1,T)*Y_ss1,'b',transitspan,ones(1,T)*Y_ss2,'r','LineWidth',lwidnum)
xlabel('Time')
ylabel('Output')
title('Transition for Output')
legend('Output','Low alpha SS','High alpha SS','Location','best')
legend boxoff
set(gca,'FontSize',fsizenum)
