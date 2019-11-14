%%%%%%%%%%%%%
% rbc_vfi.m
% Solves a simple RBC model using VFI.
% Stephen Terry, EC 702
% 11/13/17
%%%%%%%%%%%%%
close all; clear; clc;

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving an RBC model globally - woohoo, EC702 is awesome!')
disp(' ')

%%%
% Boom-bust RBC model without labor
% V(A,K) = max_{K'} { C^(1-gamma)/(1-gamma) + beta E V(A',K')}
% C = A K^alpha - K' + (1-delta)K
% A in {1-sigma,1+sigma}, Pi_A = [1-p p; p 1-p];
%%%%

%%%%%%%%%
%block of code where you can change some stuff to control the model
%%%%%%%%%

%set some useful parameters relating to the model
beta = 1/1.04; %4% real interest rate, annual
delta = 0.1; %10% depreciation rate, annual
alpha = 0.5; %capital share
gamma = 2; %CRRA parameter (if this is 1, the code switches to log C utility)
sigma = 0.05; %TFP volatility parameter
p = 0.075; %TFP persistence parameter (inversely indexing persistence)

%set up some solution parameters and display parameters
Knum = 750; %size of grid for endogenous capital (exogenous TFP dimensions pre-set, see below)
vfmaxit = 500; %max number of VF iterations
vftol = 1e-5; %tolerance on Bellman equation, in percentages
Tsim = 2500; %how many periods to simulate
Terg = 500; %how many periods to discard to remove influence of initialization
Ttot = Tsim + Terg; %implied total number of periods
plotrange = (500:750); %periods of the simulation to plot
Ainit = 1; %starting point for productivity simulation
Kinit = floor(Knum/2); %starting point for capital simulation
randseed = 2501; %random number generator seed, for reproducibility
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%report time
disp('Parameter input complete.')
toc;
disp(' ')

%%%%%%%%%
%this block actually starts the numerical work in earnest, by setting up
%some grids and pre-computing some useful arrays
%%%%%%%%%

%set up the exogenous productivity process
A0 = [1-sigma; 1+sigma]; %TFP grid
Anum = length(A0); 
PrmatA = [1-p , p; p, 1-p];

%set up the capital grid
Kmin = ((alpha*(1-sigma))/((1/beta)-1+delta))^(1/(1-alpha)); %guess lower bound equal to SS with A = A_{bust} forever
Kmax = ((alpha*(1+sigma))/((1/beta)-1+delta))^(1/(1-alpha)); %guess upper bound equal to SS with A = A_{bust} forever
Kmin = Kmin*0.75; %adjust b/c these are guesses, not analytical solutions
Kmax = Kmax*1.25; %adjust b/c these are guesses, not analytical solutions
K0 = linspace(Kmin,Kmax,Knum)'; %use linear grid for simplicity

%we'll want to treat the state (A,K) as a single vectorized state, so we set up a composite grid
Statenum = Anum*Knum;

Endog0 = zeros(Statenum,2); %store values of states on composite grid
Endog0(:,1) = kron(A0,ones(Knum,1));
Endog0(:,2) = kron(ones(Anum,1),K0);

EndogInd = zeros(Statenum,2); %store indexes of states on composite grid
EndogInd(:,1) = kron((1:Anum)',ones(Knum,1));
EndogInd(:,2) = kron(ones(Anum,1),(1:Knum)');


%report time
disp('Grid setup complete.')
toc;
disp(' ')

%this block computes the current period return matrix outside of the VFI
%loop
Rmat = zeros(Statenum,Knum); %(i,j) = the current period return to choosing K'(j) given (A,K)(i)

Statect = 0;
for Act=1:Anum %loop over prod values
    for Kct=1:Knum % loop over today's capital values
        Statect = Statect + 1; % iterate composite state index
        for Kprimect=1:Knum % iterate over tomorrow's capital policy
            
            %extract states and policies
            Aval = A0(Act);
            Kval = K0(Kct);
            Kprimeval = K0(Kprimect);
            
            %compute implied output, investment, and consumption values
            Yval = Aval * (Kval^alpha);
            Ival = Kprimeval - (1-delta)*Kval;
            Cval = Yval - Ival;
            
            %the period return is CRRA, so choose a large negative number
            %if Cval<0
            
            %do you use the CRRA formulas?
            if (abs(gamma-1)>1e-3)
                if (Cval>0)
                    Rmat(Statect,Kprimect) = (1/(1-gamma))*Cval^(1-gamma);
                else 
                    Rmat(Statect,Kprimect) = - 10^8;
                end
            %or do you use log(C)?
            else
               if (Cval>0)
                    Rmat(Statect,Kprimect) = log(Cval);
                else 
                    Rmat(Statect,Kprimect) = - 10^8;
                end 
            end
        end
    end
end

%report time
disp('Return matrix calculation complete.')
toc;
disp(' ')


%%%%%%%%%
%%%%this block actually does the VFI
%%%%%%%%%

%initialize the policy and VF arrays
Vold = log(Endog0(:,2)); %this will store real value of old VF
V = 0*Vold; %this will store real value of new VF

disp('%%%%')
disp('Starting VFI now - its actually happening!')
disp(' ')

%do the VF iterations
for vfct=1:vfmaxit
    
    
    %given tomorrow's value, or the last VF iteration's function Vold,
    %construct an array giving continuation values E(V(A',K'))(i,j) as a
    %function of (A,K)(i), K'(j)
    Voldreshape = reshape(Vold,[Knum Anum]);
    EVmat = zeros(Statenum,Knum);
    EVmat(1:Knum,:) = repmat((PrmatA(1,1)*Voldreshape(:,1)+PrmatA(1,2)*Voldreshape(:,2) )',Knum,1);
    EVmat((Knum+1):Statenum,:) = repmat((PrmatA(2,1)*Voldreshape(:,1)+PrmatA(2,2)*Voldreshape(:,2) )',Knum,1);
    
    %now, construct Bellman equation
    RHSV = Rmat + beta*EVmat;
    [V,Kprimeind] = max(RHSV,[],2);
    
    %compute VF error & display
    vferr = max(100*abs(log(V./Vold)));
    if (mod(vfct,25)==1)
        disp(['VF error = ' num2str(vferr) ' on VF iteration ' num2str(vfct) '.'])
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


%%%%%%%%%
%%%%this block does some checks, plots policies, etc...
%%%%%%%%%

%now that you've converged, check to make sure the grid isn't too small
minKprime = K0(min(Kprimeind));
maxKprime = K0(max(Kprimeind));
disp('Grid Diagnostics')
disp('K0(1), min(Kprime), max(Kprime), K0(Knum)')
disp(num2str([K0(1) minKprime maxKprime K0(Knum)]))
if ((max(Kprimeind)==Knum)||(min(Kprimeind)==1))
    disp('Uh-oh, hit grid endpoints! Expand and retry.')
else
    disp('Phew, the grid is ok.')
end
disp(' ')

%plot policies and values
Kprime = K0(Kprimeind); %real values of policy
Kprime = reshape(Kprime,[Knum Anum]);
V = reshape(V,[Knum Anum]);

figure; 
subplot(1,2,1)
plot(K0,Kprime(:,1),'r',K0,Kprime(:,2),'b',K0,K0,'k','LineWidth',lwidnum)
xlabel('Capital Today'); ylabel('Capital Tomorrow'); title('Policy')
axis([K0(1) K0(Knum) min(min(Kprime))*0.9 max(max(Kprime))*1.1])
legend('Low TFP','High TFP','45-Degree Line','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(1,2,2)
plot(K0,V(:,1),'r',K0,V(:,2),'b','LineWidth',lwidnum)
axis([K0(1) K0(Knum) min(min(V))-1 max(max(V))+1])
xlabel('Capital Today K'); ylabel('Value'); title('Value Function')
legend('Low TFP','High TFP','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)

disp('Finished some basic diagnostics.')
toc;
disp(' ')


%%%%%%%%%
%%%%this block computes the model's ergodic distribution
%%%%%%%%%

%Note that the solution to the model induces a Markov chain mapping states
%(A,K) --> (A',K') with some probability, so compute implied transition 
%matrix from state to state in equilibrium. The ergodic or stationary
%distribution of this induced Markov chain is the ergodic distribution of
%the model, i.e., the distribution over states that you would get after
%many, many simulations of the model.

Prmat_big = zeros(Statenum,Statenum);

for Statect=1:Statenum

        %extract state
        Act = EndogInd(Statect,1);
        Kprimect = Kprimeind(Statect);
        Prmat_big(Statect,Kprimect) = PrmatA(Act,1);
        Prmat_big(Statect,Kprimect+Knum) = PrmatA(Act,2);

end

%extract ergodic distribution from the Markov chain on the composite grid
Prlim = (Prmat_big')^1000;
Dist = Prlim(:,1);

%plot ergodic distribution
Dist = reshape(Dist,[Knum Anum]);

figure; 
plot(K0,Dist(:,1)/sum(Dist(:,1)),'r',K0,Dist(:,2)/sum(Dist(:,2)),'b','LineWidth',lwidnum)
xlabel('Capital'); ylabel('Density'); title('Ergodic Distribution')
axis([K0(1) K0(Knum) 0 max(max(Dist))*2.5])
legend('Low TFP','High TFP','Location','North');
legend boxoff;
set(gca,'FontSize',fsizenum)

disp('Finished computing ergodic distribution.')
toc;
disp(' ')

%%%%%%%%%
%%%%this block simulates the model
%%%%%%%%%

Asimind = zeros(Ttot,1);
Ksimind = zeros(Ttot,1);
Ysim = zeros(Ttot,1);
Csim = zeros(Ttot,1);
Isim = zeros(Ttot,1);

%first, need to draw some shocks
s = RandStream('mt19937ar','Seed',randseed);
Ashocks = rand(s,Ttot,1);

%then, need to simulate the exogenous process, which you can do by
%comparing uniform shocks to intervals on the transition matrix
PrsumA = [PrmatA(:,1) ones(2,1)];
Asimind(1) = Ainit;

for t=2:Ttot
    
    %extract yesterday's exog state
    Act = Asimind(t-1);
    
    %compare today's uniform shock to transition matrix intervals
    if (Ashocks(t)<PrsumA(Act,1))
        Aprimect = 1;
    else
        Aprimect = 2;
    end
    
    %store simulated value
    Asimind(t) = Aprimect;
    
end
Asim = A0(Asimind); %store real values, not indexes

%then, need to simulation the endogenous process
Ksimind(1) = Kinit;

for t=2:Ttot
    
    %extract states
    Act = Asimind(t-1);
    Kct = Ksimind(t-1);
    Statect = Kct + Knum*(Act-1);
    
    %extract policy
    Kprimect = Kprimeind(Statect);
    
    %store simulated policy
    Ksimind(t) = Kprimect;
end
Ksim = K0(Ksimind); %store real values, not indexes

%now, simulate other endogenous variables
for t=1:(Ttot-1)
    
    Aval = Asim(t);
    Kval = Ksim(t);
    Kprimeval = Ksim(t+1);
    Ival = Kprimeval - (1-delta)*Kval;
    Yval = Aval * (Kval^alpha);
    Cval = Yval - Ival;
    
    Ysim(t) = Yval;
    Isim(t) = Ival;
    Csim(t) = Cval;
    
    
end

%discard initial periods
Asim = Asim((Terg:(Ttot-1)));
Ksim = Ksim((Terg:(Ttot-1)));
Ysim = Ysim((Terg:(Ttot-1)));
Isim = Isim((Terg:(Ttot-1)));
Csim = Csim((Terg:(Ttot-1)));


%plot a snippet of the simulation of the model

%first, form a matrix with the simulation data, in 100*logs
Xsim = 100*log([Ysim Isim Csim Asim]);
Xsim = Xsim - repmat(mean(Xsim,1),Tsim,1);

%then, actually plot the data
figure; 
subplot(2,2,1)
plot(plotrange,Xsim(plotrange,1),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% from Mean'); title('Output')

subplot(2,2,2)
plot(plotrange,Xsim(plotrange,2),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% from Mean'); title('Investment')

subplot(2,2,3)
plot(plotrange,Xsim(plotrange,3),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% from Mean'); title('Consumption')

subplot(2,2,4)
plot(plotrange,Xsim(plotrange,4),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% from Mean'); title('TFP')

%check that we're not crazy in our calculations of the ergodic
%distribution, by plotting simulating capital histograms over computed
%ergodic distributions, an overlay which should become more accurate with
%more simulated periods
figure;
subplot(1,2,1)
histogram(Ksim(Asim<1),Knum,'Normalization','probability'); 
xlabel('Capital'); 
title('Low TFP Periods'); 
ylabel('Ergodic Prob. over Sim. Prob.')
hold on;
plot(K0,Dist(:,1)/sum(Dist(:,1)),'r'); 
set(gca,'FontSize',fsizenum)
axis([K0(1) K0(Knum) 0 max(max(Dist))*3])

subplot(1,2,2)
histogram(Ksim(Asim>1),Knum,'Normalization','probability'); 
xlabel('Capital'); 
title('High TFP Periods'); 
ylabel('Ergodic Prob. over Sim. Prob.')
hold on;
plot(K0,Dist(:,2)/sum(Dist(:,2)),'b'); 
set(gca,'FontSize',fsizenum)
axis([K0(1) K0(Knum) 0 max(max(Dist))*3])


disp('Finished simulating the model.')
toc;
disp(' ')

%%%%%%%%%
%%%%this block computes moments from the model simulation
%%%%%%%%%

%compute covariances
CovX = cov(Xsim);

%extract standard deviations, relative standard deviations, correlations
StDevX = diag(CovX).^0.5;
RelStDevX = StDevX/StDevX(1);
CorrX = corr(Xsim);
CorrXY = CorrX(:,1);

%form BC moment table
BC_MOMS = [StDevX RelStDevX CorrXY];
disp('Simulated Business Cycle Moments')
disp('      StDev(X)    StDev(X)/StDev(Y)      Corr(X,Y)')
disp(['Y     ' num2str(BC_MOMS(1,:))])
disp(['I     ' num2str(BC_MOMS(2,:))])
disp(['C     ' num2str(BC_MOMS(3,:))])
disp(['TFP   ' num2str(BC_MOMS(4,:))])
disp(' ')

toc;
disp('Done with solving, simulating, and analyzing the model - nice!')
disp('%%%%%%%%%%%%%%%%%%%')
