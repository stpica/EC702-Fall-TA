%%%%%%%%%%%%%
% ge_incomplete_mkts.m
% Solves the incomplete mkts model using the EGM in general equilibrium.
% Stephen Terry, EC 702
% 11/26/17
%%%%%%%%%%%%%
close all; clear; clc;

%%%%initial program setup

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the incomplete markets model in GE - woohoo, EC702 is awesome!')
disp(' ')

%set parameters
beta = 0.96; %HH subjective discount factor
alpha = 0.33; %capital share
delta = 0.1; %depreciation rate
gamma = 2; %CRRA risk aversion
eh = 1.0891; %high labor efficiency level 
el = 0.1980; %low labor efficiency level 
phh = 0.95; %high efficiency persistence
plh = 0.45; %low efficiency (inverse) persistence
%NOTE: code adjusts (el,eh) to ensure the stationary lab supply is 1
abar = 0; %lower bound on assets

%set grid dimensions and bounds
anum = 500; %assets grid size over which to solve the Euler equations
adensenum = 1000; %assets grid size over which to plot policies, compute ergodic distribution
amin = abar;
amax = 15;
enum = 2;

%set up solution parameters & some plotting parameters
maxpolit = 500; %max number of time iterations on the value function
maxpolerr = 1e-5; %max abs error or change in the asset mapping
maxdistit = 1500; %max number of distributional iterations
maxdisterr = 1e-5; %max abs error of change in stationary dist
maxRit = 50; %max number of iterations to find equilibrium real interest rate
maxKerr = 1e-5; %max abs error of asset market clearing
maxRerr = 1e-5; %max abs error of R bounds
lwidnum = 2; %line width for graphs
fsizenum = 12; %font size for graphs
Tsim = 5000; %how many periods to simulate
Terg = 500; %how many periods to discard to remove influence of initialization
Ttot = Tsim + Terg; %implied total number of periods
plotrange = (700:800); %periods of the simulation to plot
ainit = (amin+amax)/2; %starting point for asset simulation
einit = 2; %starting point for labor shock simulation
randseed = 2501; %random number generator seed, for reproducibility

%set up assets grid
aprime0 = linspace(amin,amax,anum)';
adense0 = linspace(amin,amax,adensenum)';

%set up labor shock grid & transition matrix
e0 = [el; eh];
Pie = [ 1-plh plh; 1-phh phh];
edist = (Pie')^1000;
edist = edist(:,1);
e0 = e0/sum(e0.*edist);

disp('Finished setting up the model.')
toc;
disp(' ')


%%%%%%%%%%%%
%%%%this block loops over different real interest guesses
%%%%%%%%%%%%


%set up bisection algorithm in the real interest rate
a = 0.9;
c = (1/beta)-(1e-4);
b = (a+c)/2; %initial guess for the real interest rate

for Rit = 1:maxRit
    
    %choose R from the bisection algorithm
    R = b;
    
    disp('%%%%%%%%%%')
    disp(['%%% Entering iteration ' num2str(Rit) ' for guesses of R.'])
    disp(['%%% Rguess = ' num2str(R) '.'])
    disp(' ')

    %determine the aggregate capital level and wage implied by the guess for R
    K = capdemand(R,alpha,delta);
    W = wage(K,alpha);

    %%%%%%%%%%%%
    %%%%given the real interest rate, this block solves for optimal policies by
    %%%%employing time iteration on the Euler equation with an implementation
    %%%%of the EGM
    %%%%%%%%%%%%

    %this analysis is based on the endogenous grid points method of Carroll
    %(2006)
    
    disp('Solving for optimal savings policies of the HHs.')

    
    %start with a guess for (aold -> aprime0) mapping, i.e.
    %a matrix aold(aprimenum,enum) s.t. a'(aold(i),e0(j)) = aprime0(i);
    if (Rit==1)
        aold =[aprime0 aprime0]; %initial guess
        anew = 0*aold;
    end
    
    %now, actually do the time iteration
    for polit = 1:maxpolit

        %loop over labor shock values today
        for ect=1:enum

            %extract today's state
            eval = e0(ect);

            %loop over policies tomorrow
            for aprimect=1:anum

                %extract tomorrow's assets
                aprimeval = aprime0(aprimect);

                %compute the expected marginal utility tomorrow, given the
                %policy rule induced by aold --> aprime0
                expecmu = 0;
                %loop over possible labor levels
                for eprimect=1:enum

                    %extract tomorrow's labor shock
                    eprimeval = e0(eprimect);

                    %interpolate asset choice tomorrow
                    aprimeprimeval =  aprimeinterp(aold,aprime0,aprimeval,eprimect);

                    %tomorrow's consumption
                    cprimeval = R*aprimeval + W*eprimeval - aprimeprimeval;

                    %create expected value tomorrow of marginal utility
                    expecmu = expecmu + Pie(ect,eprimect)*muofc(cprimeval,gamma);

                end

                %construct RHS of time iteration formula
                RHS = invmuofc(beta*R*expecmu,gamma);
                RHS = (1/R)*(aprimeval - W*eval + RHS);

                %store new implied asset value, which must have been held today
                anew(aprimect,ect) = RHS;

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

    %now, interpolate the savings policies on a fine grid
    aprimedense = zeros(adensenum,enum);
    for ect=1:enum
        aprimedense(:,ect) = interp1q(apol(:,ect),aprime0,adense0); 
        aprimedense(adense0<apol(1,ect),ect)=amin;
        aprimedense(adense0>apol(anum,ect),ect)=amax;        
    end

    %%%%%%%%%%%%
    %%%%this block computes and plot the stationary distribution of the model
    %%%%over the finely discretized grid for assets
    %%%%%%%%%%%%

    %this analysis is based on the non-stochastic simulation method of Young
    %(2010)

    disp('Solving for stationary distribution from HH policies.')

    %do some initializations
    if (Rit==1)
        adistold = zeros(adensenum,2);
        adistold(1,:) = 0.5;
        adistnew = 0*adistold;
    end
    
    %actually loop over distributional pushforwards
    for distct=1:maxdistit

        %loop over states today
        for act=1:adensenum
           for ect=1:enum

               wgtnow = adistold(act,ect);
               aprimeval = aprimedense(act,ect);

               aprimect = sum(aprimeval>adense0);
               if (aprimect>0)&&(aprimect<adensenum)
                    aprimewgt = (adense0(aprimect+1)-aprimeval)/...
                        (adense0(aprimect+1)-adense0(aprimect));
               elseif (aprimect==0)
                   aprimect=1;
                   aprimewgt = 0.0;
               elseif (aprimect==adensenum)
                   aprimect=adensenum-1;
                   aprimewgt = 1.0;               
               end

               for eprimect=1:enum
                    adistnew(aprimect,eprimect) = adistnew(aprimect,eprimect) + ...
                        Pie(ect,eprimect)*wgtnow*(1.0-aprimewgt);
                    adistnew(aprimect+1,eprimect) = adistnew(aprimect+1,eprimect) + ...
                        Pie(ect,eprimect)*wgtnow*aprimewgt;
               end

           end
        end

        adistnew = adistnew/sum(adistnew(:));
        
        %now, check for error in stationarity
        disterr = max(abs(adistold(:)-adistnew(:)));

        if (mod(distct,50)==1)
        disp(['On dist iteration ' num2str(distct) ' error is ' num2str(disterr)])
        end

        %exit if error is small enough
        if (disterr<maxdisterr)
            break;
        end

        %if error isn't small enough, update and move on
        adistold = adistnew;
        adistnew = 0*adistold;

    end
    adist= adistold;
    disp(['Finished computing stationary distribution. Max err = ' num2str(disterr) '.'])
    toc;
    disp(' ')

    %check grid upper bound, and output warning if you hit the upper bound
    if (max(adist(end,:))>1e-7)
        disp('%%WARNING!')
        disp('Positive weight on upper bound of grid, expand and retry.')
        disp('%%WARNING!')
        disp(' ')
    end

    %%%%%%%%%%%%
    %%%%this block computes the implied macro asset supply implied by the R
    %%%%guess
    %%%%%%%%%%%%
    
    disp('Computing implied household capital supply.')

    %simply add up capital holdings over the stationary distribution
    Ksupply = sum(adense0.*adist(:,1)+adense0.*adist(:,2));
    
    %what is the interest rate implied by the HH savings behavior?
    Rimplied = intrate(Ksupply,alpha,delta);
    
    %compute errors
    Kerr = Ksupply-K;
    
    disp(['On iteration ' num2str(Rit) ' we have: '])
    disp(['K implied by firm demand = ' num2str(K)])
    disp(['K supplied by HHs = ' num2str(Ksupply)])
    disp(['K error = K supply - K demand = ' num2str(Kerr)])
    
   %exit if GE error is small enough or if interval has converged
   if (abs(Kerr)<maxKerr)||(abs(c-a)<maxRerr)
        break;
   end
   
   %if not converged, update R guess and move on
   fb = Kerr;
   if (fb>0)
      c = b;
      b = (a+c)/2;
   elseif (fb<0)
      a = b;
      b = (a+c)/2;
   end
   
    
end
%this ends the loop over real interest rate guesses


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
for ect=1:enum
    cdense(:,ect) = R*adense0 + e0(ect) - aprimedense(:,ect);
end

figure;
plot(adense0,cdense(:,1),'r',adense0,cdense(:,2),'b',adense0,(R*adense0+e0(1)),'k--','LineWidth',lwidnum)
axis([amin amax 0 max(cdense(:))*1.1])
title('Consumption Functions')
xlabel('Assets Today')
ylabel('Consumption Today')
legend('Low Income','High Income','Cash on Hand = Ra + y','Location','NorthWest')
legend boxoff;
set(gca,'FontSize',fsizenum)
% saveas(gcf,'CONS_POL.pdf')

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

disp(' ')
toc;
disp('Solved for GE, so we are done with solving and analyzing the model - nice!')
disp('%%%%%%%%%%%%%%%%%%%')
