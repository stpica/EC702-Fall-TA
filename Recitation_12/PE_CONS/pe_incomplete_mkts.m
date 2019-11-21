%%%%%%%%%%%%%
% pe_incomplete_mkts.m
% Solves an incomplete mkts model using the EGM in partial equilibrium.
% Stephen Terry, EC 702
% 11/25/17
%%%%%%%%%%%%%
close all; clear; clc;

%%%%initial program setup

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the incomplete markets model in PE - woohoo, EC702 is awesome!')
disp(' ')

%set parameters
R = 1.02; %gross real interest rate
beta = 1/1.04; %HH subjective discount rate
gamma = 2; %CRRA parameter
yl = 0.25; %low income realization
yh = 1.0; %high income realization
plh = 0.7; %P(y' = yh | y = yl)
phh = 0.95; %P(y' = yh | y = yh)
abar = 0; %borrowing constradense0nt a' >= abar 

%set grid dimensions and bounds
anum = 300; %assets grid size over whiche to solve the Euler equations
adensenum = 750; %assets grid size over which to plot policies, compute ergodic distribution
amin = abar;
amax = 1.75;
ynum = 2;

%set up solution parameters & some plotting parameters
maxpolit = 500; %max number of time iterations on the value function
maxpolerr = 1e-7; %max abs error or change in the asset mapping
maxdistit = 1000; %max number of distributional iterations
maxdisterr = 1e-7; %max abs error of change in stationary dist
lwidnum = 2; %line width for graphs
fsizenum = 12; %font size for graphs
Tsim = 5000; %how many periods to simulate
Terg = 500; %how many periods to discard to remove influence of initialization
Ttot = Tsim + Terg; %implied total number of periods
plotrange = (700:800); %periods of the simulation to plot
ainit = (amin+amax)/2; %starting point for asset simulation
yinit = 2; %starting point for income simulation
randseed = 2501; %random number generator seed, for reproducibility

%set up assets grid
aprime0 = linspace(amin,amax,anum)';
adense0 = linspace(amin,amax,adensenum)';

%set up income grid & transition matrix
y0 = [yl; yh];
Piy = [ 1-plh plh; 1-phh phh];


disp('Finished setting up the model.')
toc;
disp(' ')

%%%%%%%%%%%%
%%%%this block solves for optimal policies by employing time iteration on
%%%%the Euler equation with an implementation of the EGM
%%%%%%%%%%%%

%this analysis is based on the endogenous grid points method of Carroll
%(2006)

%start with a guess for (aold -> aprime0) mapping, i.e.
%a matrix aold(aprimenum,ynum) s.t. a'(aold(i),y0(j)) = aprime0(i);
aold =[aprime0 aprime0]; %initial guess
anew = 0*aold;

%now, actually do the time iteration
for polit = 1:maxpolit
    
    %loop over income values today
    for yct=1:ynum
    
        %extract today's state
        yval = y0(yct);
        
        %loop over policies tomorrow
        for aprimect=1:anum

            %extract tomorrow's assets
            aprimeval = aprime0(aprimect);
            
            %compute the expected marginal utility tomorrow, given the
            %policy rule induced by aold --> aprime0
            expecmu = 0;
            %loop over possible income levels
            for yprimect=1:ynum

                %extract tomorrow's income
                yprimeval = y0(yprimect);

                %interpolate asset choice tomorrow
                aprimeprimeval =  aprimeinterp(aold,aprime0,aprimeval,yprimect);
                
                %tomorrow's consumption
                cprimeval = aprimeval + yprimeval - (1/R)*aprimeprimeval;

                %create expected value tomorrow of marginal utility
                expecmu = expecmu + Piy(yct,yprimect)*muofc(cprimeval,gamma);

            end

            %construct RHS of time iteration formula
            RHS = invmuofc(beta*R*expecmu,gamma);
            RHS = (1/R)*aprimeval - yval + RHS;
            
            %store new implied asset value, which must have been held today
            anew(aprimect,yct) = RHS;

        end
    end

    %now, check for error in the interpolation mapping
    polerr = max(abs(anew(:)-aold(:)));
    
    %output diagnostics, but not too frequently
    if (mod(polit,25)==1)
    disp(['On iteration ' num2str(polit) ' error is ' num2str(polerr)])
    end
    
    %exit time iteration loop if error is small enough
    if (polerr<maxpolerr)
        break;
    end
    
    %if error isn't small enough, update and move on
    aold = anew;
    
end

disp(['Finished solving for policies. Max err = ' num2str(polerr) '.'])
toc;
disp(' ')

%%%%%%%%%%%%
%%%%this block creates a dense interpolation of the policy, then plots
%%%%assets/consumption on a fine grid
%%%%%%%%%%%%

%store converged policy matrix
apol = anew;

%now, plot the savings policies, by interpolating them on a fine grid
aprimedense = zeros(adensenum,ynum);
for yct=1:ynum
    aprimedense(:,yct) = interp1q(apol(:,yct),aprime0,adense0); 
    aprimedense(adense0<apol(1,yct),yct)=amin;
    aprimedense(adense0>apol(anum,yct),yct)=amax;
end

figure;
plot(adense0,adense0,'k',adense0,aprimedense(:,1),'r',adense0,aprimedense(:,2),'b','LineWidth',lwidnum)
axis([amin amax amin max(aprimedense(:))*1.1])
title('Savings Functions')
xlabel('Assets Today')
ylabel('Assets Tomorrow')
legend('45-Degree Line','Low Income','High Income','Location','NorthWest')
legend boxoff;
set(gca,'FontSize',fsizenum)
% saveas(gcf,'SAVE_POL.pdf')


%now, plot the consumption policies, by interpolating them on a fine grid
cdense = 0*aprimedense;
for yct=1:ynum
    cdense(:,yct) = adense0 + y0(yct) - (1/R)*aprimedense(:,yct);
end

figure;
plot(adense0,cdense(:,1),'r',adense0,cdense(:,2),'b',adense0,(adense0+y0(1)),'k--','LineWidth',lwidnum)
axis([amin amax 0 max(cdense(:))*1.1])
title('Consumption Functions')
xlabel('Assets Today')
ylabel('Consumption Today')
legend('Low Income','High Income','Cash on Hand = a + y','Location','NorthWest')
legend boxoff;
set(gca,'FontSize',fsizenum)
% saveas(gcf,'CONS_POL.pdf')


disp(['Finished plotting policies. Max err = ' num2str(polerr) '.'])
toc;
disp(' ')

%%%%%%%%%%%%
%%%%this block computes and plot the stationary distribution of the model
%%%%over the finely discretized grid for assets
%%%%%%%%%%%%

%this analysis is based on the non-stochastic simulation method of Young
%(2010)

%do some initializations
adistold = ones(adensenum,2);
adistold = adistold/sum(adistold(:));
adistnew = 0*adistold;

%actually loop over distributional pushforwards
for distct=1:maxdistit
   
    %loop over states today
    for act=1:adensenum
       for yct=1:ynum
          
           wgtnow = adistold(act,yct);
           aprimeval = aprimedense(act,yct);
           
           aprimect = sum(aprimeval>adense0);
           if (aprimect>0)&&(aprimect<adensenum)
                aprimewgt = (adense0(aprimect+1)-aprimeval)/...
                    (adense0(aprimect+1)-adense0(aprimect));
           elseif (aprimect==0)
               aprimect=1;
               aprimewgt = 1.0; %should be 1.0
           elseif (aprimect==adensenum)
               aprimect=adensenum-1;
               aprimewgt = 0.0;  %should be 0.0             
           end
           
           for yprimect=1:ynum
                adistnew(aprimect,yprimect) = adistnew(aprimect,yprimect) + ...
                    Piy(yct,yprimect)*wgtnow*aprimewgt;
                adistnew(aprimect+1,yprimect) = adistnew(aprimect+1,yprimect) + ...
                    Piy(yct,yprimect)*wgtnow*(1.0-aprimewgt);
           end
           
       end
    end
    
    
    %now, check for error in stationarity
    disterr = max(abs(adistold(:)-adistnew(:)));
    
    if (mod(distct,25)==1)
    disp(['On dist iteration ' num2str(distct) ' error is ' num2str(disterr)])
    end
    
%    exit if error is small enough
    if (disterr<maxdisterr)
        break;
    end

    %if error isn't small enough, update and move on
    adistold = adistnew;
    adistnew = 0*adistold;
    
end

disp(['Finished computing stationary distribution. Max err = ' num2str(disterr) '.'])
toc;
disp(' ')

%check grid upper bound, and output warning if you hit the upper bound
if (max(adistnew(end,:))>0.0)
    disp('%%WARNING!')
    disp('Positive weight on upper bound of grid, expand and retry.')
    disp('%%WARNING!')
    disp(' ')
end


%%%now, plot the conditional distributions for each state
adist = adistnew;
adistcond = adist;
adistcond(:,1) = adistcond(:,1)/sum(adistcond(:,1));
adistcond(:,2) = adistcond(:,2)/sum(adistcond(:,2));

figure;
plot(adense0,adistcond(:,1),'r',adense0,adistcond(:,2),'b','LineWidth',lwidnum)
axis([amin amax 0 max(adistcond(:))*1.1])
title('Stationary Distributions Conditional on Income')
xlabel('Assets Today')
legend('Low Income','High Income','Location','NorthWest')
legend boxoff;
set(gca,'FontSize',fsizenum)
% saveas(gcf,'STATIONARY_DIST.pdf')


%%%%%%%%%%%%
%%%%this block simulates the model, allowing for continuous choice of 
%%%%assets by the household
%%%%%%%%%%%%
ysimind = zeros(Ttot,1);
asim = zeros(Ttot,1);
csim = zeros(Ttot,1);

%first, need to draw some shocks
s = RandStream('mt19937ar','Seed',randseed);
yshocks = rand(s,Ttot,1);

%then, need to simulate the exogenous process, which you can do by
%comparing uniform shocks to intervals on the transition matrix
Pisumy = [Piy(:,1) ones(2,1)];
ysimind(1) = yinit;

for t=2:Ttot
    
    %extract yesterday's exog state
    yct = ysimind(t-1);
    
    %compare today's uniform shock to transition matrix intervals
    if (yshocks(t)<Pisumy(yct,1))
        yprimect = 1;
    else
        yprimect = 2;
    end
    
    %store simulated value
    ysimind(t) = yprimect;
    
end
ysim = y0(ysimind); %store real values, not indexes

%now, need to simulate the endogenous process
asim(1) = ainit;

for t=2:Ttot
    
    %extract last period's asset holdings & income level
    aval = asim(t-1);
    yct = ysimind(t-1);
    
    %what is the interpolated level of assets today?
    aprimeval = interp1q(apol(:,yct),aprime0,aval); 
    aprimeval = min([max([amin aprimeval]) amax]);
    
    asim(t) = aprimeval;

end

%now, need to obtain the implied consumption value
csim(:) = 0.0;

for t=1:(Ttot-1)
    
   %extract today's assets, income, and savings
   aval = asim(t);
   yval = ysim(t);
   aprimeval = asim(t+1);
   cval = aval + yval - (1/R) * aprimeval;
   
   csim(t) = cval;
    
    
end

%remove the influence of initial conditions
asim = asim((Terg:(Ttot-1)));
ysim = ysim((Terg:(Ttot-1)));
csim = csim((Terg:(Ttot-1)));
ysimind = ysimind((Terg:(Ttot-1)));
sratesim = (ysim - csim)./asim;
sratesim(asim==0) = 0;


%plot a snippet of the simulation of the model

%first, form a matrix with the simulation data, in 100*logs
Xsim = [ysim  csim asim 100*sratesim ];

%then, actually plot the data
figure; 
subplot(2,2,1)
plot(plotrange,Xsim(plotrange,1),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('y'); title('Income')
axis([plotrange(1) plotrange(end) 0 1.1])

subplot(2,2,2)
plot(plotrange,Xsim(plotrange,2),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('c'); title('Consumption')
axis([plotrange(1) plotrange(end) 0 1.1])

subplot(2,2,3)
plot(plotrange,Xsim(plotrange,3),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('a'); title('Assets')
axis([plotrange(1) plotrange(end) 0 amax*1.1])

subplot(2,2,4)
plot(plotrange,Xsim(plotrange,4),'b',plotrange,0*plotrange,'k','LineWidth',lwidnum);
xlabel('Period'); ylabel('% of Assets'); title('Savings Rate')
axis([plotrange(1) plotrange(end) -125 250])
% saveas(gcf,'UNCOND_SIM.pdf')

disp('Finished simulating the model.')
toc;
disp(' ')

%%%%%%%%%
%%%%this block computes moments from the model simulation
%%%%%%%%%

CorrX = corr(Xsim);

disp('Simulated Correlation Matrix')
disp(' Income ,    Cons   ,  Savings Rate  , Assets')
disp(num2str(CorrX))
disp(' ')

toc;
disp('Done with solving, simulating, and analyzing the model - nice!')
disp('%%%%%%%%%%%%%%%%%%%')

