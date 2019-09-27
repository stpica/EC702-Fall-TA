%%%%%%%%%%%%%
% reaching_the_center.m
% Solves the "Reaching the Center" problem using VFI.
% Stefano Pica, TA for EC 702
% Fall 2019
%%%%%%%%%%%%%

close all; clear; clc;

% model to be solved by value function iteration
% V(x) = min_{y in R} [ (x-y)^2 + y^2 + beta V(y) ]


%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the "Reaching the Center" problem - awesome!')
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vfmaxit = 1000; %max number of VF iterations
pfmaxit = 20; %max number of PF iterations
vftol = 1e-7; %tolerance on Bellman equation, in percentages
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%utility parameters
beta = 0.95; %discount factor

%grid
dmin = -1;
dmax = 1;
dnum = 1000;
d0 = linspace(dmin,dmax,dnum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%compute return matrix. Can do it outside loop of value function iteration
%as it is indipendent of it.

Rmat = zeros(dnum,dnum);

for dct = 1:dnum
    for dprimect = 1:dnum
        dval = d0(dct);
        dprimeval = d0(dprimect);
       	Rmat(dct,dprimect) = (dval - dprimeval) ^ 2 + dprimeval ^ 2;
    end
end

%this matrix tells me: for each state (position I am in the line - fix the
%row to see this), what is my utility if I jump to reach distance y (look
%at columns here)

disp('Finished setting up the return matrix.')
toc;
disp(' ')


%%%%%%%%%%%%%%%%%%% OPTIONAL-begin %%%%%%%%%%%%%%
%alternative method to computing the return matrix via loop

[xvar,yvar] = meshgrid(d0, d0);
%you can restrict the set of "feasible jumps" (u can jump within 0 and your
%current x)
%Glow=zeros(knum, knum);
%Gmax=xvar;
%yvar(yvar>=Gmax | yvar<Glow)=NaN;
Rmat2 = (xvar - yvar) .^ 2 + yvar .^ 2;
Rmat2 = Rmat2';

%%%%%%%%%%%%%%%%%%% OPTIONAL-end %%%%%%%%%%%%%%


%initialize the policy and VF arrays
Vold = zeros(dnum,1); %this will store real value of old VF
V = zeros(dnum,1);  %this will store real value of new VF
dprimeind = V;

disp('%%%%')
disp('Starting VFI now - it is actually happening!')
disp(' ')


for vfit = 1:vfmaxit
    %form Bellman equation
    RHSMAT = Rmat + beta * repmat(Vold',dnum,1);
   	[V, dprimeind] = min(RHSMAT,[],2);
    absdiff = abs(V-Vold);
    vferr = max(absdiff); %maximum absolute error
    if (mod(vfit,25)==1) %show print every 25 steps of the loop
        disp(['VF error = ' num2str(vferr) ' on VF iteration ' num2str(vfit) '.'])
    end
       
    if (vferr < vftol) %exit if converged
        break; 
    end
    
    Vold = V; %if not converged, then Vold <- V, and redo
end

disp(' ')
disp('VF iteration complete.')
toc;
disp('%%%%')
disp(' ')

dprime = d0(dprimeind); %real values of policy function. Each point tells me where I jump


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM AN ARBITRARY INITIAL CONDITION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x0 = 1; %initial condition. let's start jumping from x0
transit = zeros(pfmaxit,1); %preallocate the transition vector
transit(1) = x0;
for pfit = 2:pfmaxit
    transit(pfit) = dprime(d0==x0); %get from the grid the initial value and apply the PF
    x0 = transit(pfit); %the new initial value is where we just jumped
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
subplot(2,2,1)
plot(d0,dprime,'r',d0,d0,'k','LineWidth',lwidnum)
xlabel('Line')
ylabel('Jump')
title('Policy Function')
legend('How much do I jump?','45 Line','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,2)
plot(d0,V,'b','LineWidth',lwidnum)
xlabel('Line')
ylabel('Value'); title('Value Function')
set(gca,'FontSize',fsizenum)

subplot(2,2,3)
plot(0:pfmaxit-1,transit,'b','LineWidth',lwidnum)
xlabel('Time')
ylabel('Jumps'); title('Transition from x0 at t=0')
set(gca,'FontSize',fsizenum)
