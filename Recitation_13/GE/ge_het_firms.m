%%%%%%%%%%%%%
% ge_het_firms.m
% Solves a GE heterogeneous firms investment model
% Stephen Terry, EC 702
% 12/13/17
%%%%%%%%%%%%%
close all; clear; clc;

%%%%%%%%model to be solved
%V(z,k) = max_{k',n} [ d + (1/R) E (V(z',k')|z) ]
%d = y - i - Wn - AC(k,k')
%i = k' - (1-delta)k
%AC(k,k) = gamma/2 * i^2
%z 2pt Markov chain
%clearing conditions R = 1/beta, W = phi C
% C = Y - I
%%%%%%%%

%%%%initial program setup

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the heterogeneous firms model in GE - woohoo, EC702 is awesome!')
disp(' ')

%set parameters
beta = 1/1.04; %household subjective discount factor
R = 1/beta; %implies the real interest rate in GE
phi = 2; %household disutility of labor
delta = 0.1; %depreciation rate
alpha = 0.25; %capital share
nu = 0.5; %labor share
zl = 0.85; %low income realization
zh = 1.15; %high income realization
plh = 0.25; %P(z' = zh | z = zl)
phh = 0.75; %P(z' = zh | z = zh)
gamma = 2; %adjustment cost parameter

%set grid dimensions for both states (z,k)
knum = 300; %capital grid
znum = 2; %productivity grid

%set capital grid boundaries
kmin = 0.6;
kmax = 0.9;



%set up solution parameters & some plotting parameters
maxWit = 50; %max number of iterations on the wage
maxWerr = 1e-4; %max error on the wage consistency condition
maxvfit = 1000; %max number of iterations on the value function
maxvferr = 1e-7; %max abs error or change in the value functions
maxdistit = 1000; %max number of distributional iterations
maxdisterr = 1e-7; %max abs error of change in stationary dist
lwidnum = 2; %line width for graphs
fsizenum = 12; %font size for graphs
Tsim = 5000; %how many periods to simulate
Terg = 500; %how many periods to discard to remove influence of initialization
Ttot = Tsim + Terg; %implied total number of periods
plotrange = (700:750); %periods of the simulation to plot
kinit = floor(knum/2); %starting point for capital simulation
zinit = 2; %starting point for productivity simulation
randseed = 2501; %random number generator seed, for reproducibility

%set up capital grid
k0 = linspace(kmin,kmax,knum)';

%set up productivity grid
z0 = [zl; zh];
Piz = [ 1-plh plh; 1-phh phh];


disp('Finished setting up the model.')
toc;
disp(' ')

%%%%%%%%%%%%
%%%%This block begins iteration on the wage clearing loop
%%%%%%%%%%%%

%clear markets using bisection
a = 0.01; %lower bound for bisection on wage
c = 5; %upper bound for bisection on wage



for Wit = 1:maxWit

%INSERT GUESS FOR THE WAGE
b = (a+c)/2;
W = b;


%%%%%%%%%%%%
%%%%This block sets up the return function, i.e., it creates a matrix
%%%%R(z,k,k') such that V(z,k) = max_{k'} R(z,k,k') + (1/R) E (V(z',k')|z).
%%%%Since the return matrix doesn't depend on V, pre-computation saves
%%%%time during the later VFI process
%%%%%%%%%%%%

%initialize return matrix
Rmat = zeros(znum,knum,knum);

%loop over potential states and policies
for zct=1:znum
for kct=1:knum
for kprimect=1:knum
    
    %extract values
    zval = z0(zct);
    kval = k0(kct);
    kprimeval = k0(kprimect);
    
    %extract value of the dividend flow d = y - Wn - i - AC
    dval = dividend(zval,kval,kprimeval,alpha,nu,W,delta,gamma);
    
    %flow return in the Bellman equation is just the dividend
     Rmat(zct,kct,kprimect) = dval;

    
end
end
end


disp('Finished setting up the return matrix.')
toc;
disp(' ')


%%%%%%%%%%%%
%%%%This block performs the VFI process to solve for optimal policies
%%%%%%%%%%%%

%start with a guess for the value function, and track policies
Vold = zeros(znum,knum); %will track guessed value
kprimeindold = Vold; %will track indexes of last guess of capital policy
kprimeold = Vold; %will track real values of last guess of capital policy

V = Vold; %will track firm value
kprimeind = V; %will track indexes of capital policy
kprime = V; %will track real values of capital policy

%now, actually do the value function iteration
for vfit = 1:maxvfit
    
    %now, loop over productivity values today
    for zct=1:znum
        
        %NOTE: to save on time-consuming loops, this block of code is 
        %optimizing the Bellman equation in matrix form for
        %all capital values simultaneously given that z=z(zct)
        
        %form return matrix 
        RHSMAT = squeeze(Rmat(zct,:,:));
        
        %form expected value matrix
        EVMAT = zeros(knum,knum);
        for zprimect=1:znum
            EVMAT = EVMAT +  Piz(zct,zprimect)*repmat(V(zprimect,:),knum,1);
        end
        
        %form Bellman equation
        RHSMAT = RHSMAT + (1/R)*EVMAT;
        
        %optimize Bellman equation
        [Vdum,kprimeinddum] = max(RHSMAT,[],2);
        V(zct,:) = Vdum;
        kprimeind(zct,:) = kprimeinddum;
        kprime(zct,:) = k0(kprimeind(zct,:));
        
    end
    
    %now, check error in VF and policies
    vferr = max(abs(V(:)-Vold(:)));
    polerr = max(abs(kprime(:)-kprimeold(:)));
    
    %display diagnostics periodically
    if (mod(vfit,50)==1)
        disp(['On VF iter ' num2str(vfit) ' |V-Vold| = ' num2str(vferr) ' |k - kold| = ' num2str(polerr)])   
    end
    
    %exit VF loop if converged
    if (vferr<maxvferr) 
        break
    end
    
    %if haven't converged, update and move on
    kprimeold = kprime;
    kprimeindold = kprimeind;
    Vold = V;
    
end

%now, store the labor policy as well
npol = zeros(znum,knum);
for zct=1:znum
    for kct=1:knum
        npol(zct,kct) =  labor(z0(zct),k0(kct),alpha,nu,W);
    end
end


disp(['Finished solving the Bellman equation for policies and firm value. Max err = ' num2str(vferr) '.'])
toc;
disp(' ')


%%%%%%%%%%%%
%%%%this block computes the ergodic distribution of the stationary Markov
%%%%chain for (z',k')|(z,k) induced by the optimal policies
%%%%%%%%%%%%

dist = zeros(znum,knum);
distold = dist;
distold(:,:) = 1.0;
distold = distold/sum(distold(:));


for distit = 1:maxdistit
   
    %loop over states today
    for zct=1:znum
    for kct=1:knum
       
        %what is the optimal capital policy, in index form?
        kprimect = kprimeind(zct,kct);
        
        %push forward weight today
        for zprimect=1:znum
            dist(zprimect,kprimect) = dist(zprimect,kprimect) + ...
                Piz(zct,zprimect)*distold(zct,kct);
        end
        
    end
    end
    
    %round to ensure it's a distribution
    %dist = dist/sum(dist(:));
    
    %check distributional stationary error
    disterr = max(abs(dist(:)-distold(:)));
    
    %display diagnostics periodically
    if (mod(distit,5)==1)
        disp(['On dist iter ' num2str(distit) ' |dist - distold| = ' num2str(disterr)])   
    end
    
    %exit distributional loop if converged
    if (disterr<maxdisterr) 
        break
    end
    
    %if haven't converged, update and move on
    distold = dist;
    dist(:) = 0.0;
    
    
end

%if you've got any weight on the edges, output a warning
endwgt = dist(1,1)+dist(1,knum)+dist(2,1) +dist(2,knum);
if (endwgt>0)
    disp('%%%%%')
    disp('GRID HITS ENDPTS W/POS WEIGHT!!! EXPAND AND RETRY')
    disp('%%%%%')
end

disp(['Finished solving for the stationary distribution. Max err = ' num2str(disterr) '.'])
toc;
disp(' ')


%%%%%%%%%%%%
%%%%this block computes the aggregates needed to determine whether market
%%%%clearing has occurred or not
%%%%%%%%%%%%

Yagg = 0;
Cagg = 0;
Nagg = 0;
Kagg = 0;
Iagg = 0;
Dagg = 0;
ACagg = 0;

for zct=1:znum
   for kct=1:knum
      
       %extract states & policies
       zval = z0(zct);
       kval = k0(kct);
       nval = npol(zct,kct);
       kprimeval = kprime(zct,kct);
       ival = kprimeval - (1-delta)*kval;
       yval = output(zval,kval,nval,alpha,nu);
       dval = dividend(zval,kval,kprimeval,alpha,nu,W,delta,gamma);
       ACval = AC(kval,kprimeval,delta,gamma);
       
       %increment aggregates
       Yagg = Yagg + yval*dist(zct,kct);
       Kagg = Kagg + kval*dist(zct,kct);
       Iagg = Iagg + ival*dist(zct,kct);
       Nagg = Nagg + nval*dist(zct,kct);
       Dagg = Dagg + dval*dist(zct,kct);
       ACagg = ACagg + ACval*dist(zct,kct);
       
   end
end

%some of the aggregates are implied
Cagg = Yagg - Iagg - ACagg;

%now, check for market clearing
Werr = W - phi * Cagg;

disp('%%%%')
disp(['On W loop ' num2str(Wit) ' W - phi * C = ' num2str(Werr)])  
disp(['W guess = ' num2str(W)])
disp(['C = ' num2str(Cagg) ', so phi*C = ' num2str(phi*Cagg) ])
disp(['N = ' num2str(Nagg)])
disp(['Y = ' num2str(Yagg)])
disp(['I = ' num2str(Iagg)])
disp(['K = ' num2str(Kagg)])
disp(['AC = ' num2str(ACagg)])
disp('%%%%')

%if converged, then done
if (abs(Werr)<maxWerr)||((c-a)<maxWerr)
    break
end

%otherwise, update bisection bounds and return
fb = Werr;
if (fb>0) 
    c = b;
elseif (fb<0)
    a = b;
end 
disp(' ') 

end
%%%THIS ENDS THE ITERATION ON THE WAGE

%at this point, the model is solved, and simulation, plots, etc... can
%continue as before for the PE version of the model
disp(' ')
disp('%%%%')
disp('Model is solved')
disp('Final aggregate values')
disp(['W = ' num2str(W)])
disp(['C = ' num2str(Cagg)])
disp(['N = ' num2str(Nagg)])
disp(['Y = ' num2str(Yagg)])
disp(['I = ' num2str(Iagg)])
disp(['K = ' num2str(Kagg)])
disp(['AC = ' num2str(ACagg)])
disp(['U(C,1-N) = ' num2str(log(Cagg)-phi*Nagg)])
TFPagg = Yagg/((Kagg^alpha)*(Nagg^nu));
disp(['TFP = Y/K^alpha * N^nu = ' num2str(TFPagg)])
disp('%%%%')

%%%%%%%%%%%%
%%%%this block simulates the model, allowing for continuous choice of 
%%%%assets by the household
%%%%%%%%%%%%
zsimind = zeros(Ttot,1);
ksimind = zeros(Ttot,1);
isim = zeros(Ttot,1);
dsim = zeros(Ttot,1);
nsim = zeros(Ttot,1);
ysim = zeros(Ttot,1);

%first, need to draw some shocks
s = RandStream('mt19937ar','Seed',randseed);
zshocks = rand(s,Ttot,1);

%then, need to simulate the exogenous process, which you can do by
%comparing uniform shocks to intervals on the transition matrix
Pisumz = [Piz(:,1) ones(2,1)];
zsimind(1) = zinit;

for t=2:Ttot
    
    %extract yesterday's exog state
    zct = zsimind(t-1);
    
    %compare today's uniform shock to transition matrix intervals
    if (zshocks(t)<Pisumz(zct,1))
        zprimect = 1;
    else
        zprimect = 2;
    end
    
    %store simulated value
    zsimind(t) = zprimect;
    
end
zsim = z0(zsimind); %store real values, not indexes

%now, need to simulate the endogenous process
ksimind(1) = kinit;

for t=2:Ttot
    
    %extract last period's states
    kct = ksimind(t-1);
    zct = zsimind(t-1);
    
    %what is the optimal value of capital chosen for today?
    kprimect = kprimeind(zct,kct);
    ksimind(t) = kprimect;
    
end
ksim=k0(ksimind);


%now, need to obtain the implied investment, labor, dividend values
for t=2:Ttot-1
   
    %extract policies and states
    zval = zsim(t);
    kval = ksim(t);
    kprimeval = ksim(t+1);
    ival = kprimeval - (1-delta)*kval;
    nval = labor(zval,kval,alpha,nu,W);
    dval = dividend(zval,kval,kprimeval,alpha,nu,W,delta,gamma);
    yval = output(zval,kval,nval,alpha,nu);
    
    %insert results
    isim(t) = ival;
    nsim(t) = nval;
    dsim(t) = dval;
    ysim(t) = yval;
    
end


%remove the influence of initial conditions
zsim = zsim((Terg:(Ttot-1)));
zsimind = zsimind((Terg:(Ttot-1)));
ksim = ksim((Terg:(Ttot-1)));
ksimind = ksimind((Terg:(Ttot-1)));
isim = isim((Terg:(Ttot-1)));
ysim = ysim((Terg:(Ttot-1)));
dsim = dsim((Terg:(Ttot-1)));
nsim = nsim((Terg:(Ttot-1)));
iratesim = isim./ksim;
dratesim = dsim./ksim;

%plot a snippet of the simulation of the model

%first, form a matrix with the simulation data, in 100*logs
Xsim = 100*[log(zsim)  log(ysim) iratesim log(nsim) dratesim log(ksim)];

%then, actually plot the data
figure; 
subplot(2,3,1)
plot(plotrange,Xsim(plotrange,1),'b','LineWidth',lwidnum);
xlabel('Period'); ylabel('100 x Log'); title('Productivity')


subplot(2,3,2)
plot(plotrange,Xsim(plotrange,2),'b','LineWidth',lwidnum);
xlabel('Period'); ylabel('100 x Log'); title('Output')


subplot(2,3,3)
plot(plotrange,Xsim(plotrange,3),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% of Capital'); title('Investment Rate')


subplot(2,3,4)
plot(plotrange,Xsim(plotrange,4),'b','LineWidth',lwidnum);
xlabel('Period'); ylabel('100 x Log'); title('Labor')

subplot(2,3,5)
plot(plotrange,Xsim(plotrange,5),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% of Capital'); title('Dividends')

subplot(2,3,6)
plot(plotrange,Xsim(plotrange,6),'b','LineWidth',lwidnum);
xlabel('Period'); ylabel('100 x Log'); title('Capital')


%saveas(gcf,'UNCOND_SIM.pdf')

disp('Finished simulating the model.')
toc;
disp(' ')

%%%%%%%%%
%%%%this block computes moments from the model simulation
%%%%%%%%%

CovX = cov(Xsim);
StDevX = diag(CovX).^0.5; StDevX=StDevX';

disp('Simulated Volatilies')
disp(' Productivity ,  Output   ,  Investment Rate , Labor , Dividends, Capital')
disp(['StDev ' num2str(StDevX)])
disp(' ')

CorrX = corr(Xsim);
disp('Simulated Correlation Matrix')
disp(' Productivity ,  Output   ,  Investment Rate , Labor , Dividends, Capital')
disp(num2str(CorrX))
disp(' ')

%%%%%%%%%
%%%%this block plots policies, distributions, etc...
%%%%%%%%%

figure;
plot(k0,dist(1,:)/sum(dist(1,:)),'r',k0,dist(2,:)/sum(dist(2,:)),'b','LineWidth',lwidnum)
axis([kmin kmax 0 max(dist(:))*4])
title('Stationary Distributions Conditional on Productivity')
xlabel('Capital Today')
legend('Low Prod','High Prod','Location','NorthWest')
legend boxoff;
set(gca,'FontSize',fsizenum)
%saveas(gcf,'STATIONARY_DIST.pdf')

figure;
subplot(2,2,1)
plot(k0,k0,'k',k0,kprime(1,:),'r',k0,kprime(2,:),'b','LineWidth',lwidnum)
axis([kmin kmax kmin kmax*1.2])
title('Capital Choice')
xlabel('Capital Today')
ylabel('Capital Tomorrow')
legend('45-Degree Line','Low Prod','High Prod','Location','NorthWest')
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,2)
plot(k0,npol(1,:),'r',k0,npol(2,:),'b','LineWidth',lwidnum)
axis([kmin kmax 0 max(npol(:))*1.1])
title('Labor Choice')
xlabel('Capital Today')
ylabel('Labor Today')
set(gca,'FontSize',fsizenum)

iratemat = kprime;
iratemat(1,:) = iratemat(1,:)./(k0') - (1-delta);
iratemat(2,:) = iratemat(2,:)./(k0') - (1-delta);
iratemat = 100*iratemat;

subplot(2,2,3)
plot(k0,100*delta*ones(1,knum),'k',k0,iratemat(1,:),'r',k0,iratemat(2,:),'b','LineWidth',lwidnum)
axis([kmin kmax 0 max(iratemat(:))*1.1])
title('Investment Rate')
xlabel('Capital Today')
ylabel('i/k, Percent')
legend('Dep. Rate','Low Prod','High Prod','Location','NorthEast')
legend boxoff;
set(gca,'FontSize',fsizenum)


subplot(2,2,4)
plot(k0,V(1,:),'r',k0,V(2,:),'b','LineWidth',lwidnum)
axis([kmin kmax min(V(:))*0.9 max(V(:))*1.1])
title('Value')
xlabel('Capital Today')
ylabel('V')
set(gca,'FontSize',fsizenum)


%saveas(gcf,'POLICIES.pdf')



toc;
disp('Done with solving, simulating, and analyzing the model - nice!')
disp('%%%%%%%%%%%%%%%%%%%')

