%%%%%%%%%%%%%
% SolveRBC.m
% Solves an RBC model using Dynare. We then do plots and compute moments
% Stefano Pica, TA for EC702. stpica@bu.edu
% Fall 2019
%%%%%%%%%%%%%
close all; clear; clc;

addpath '/Applications/MathWorks/dynare/4.5.7/matlab' %add path for system files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK RUNS DYNARE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dynare rbc %run rbc.mod. The model follow the slides "Lecture_RBC"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK COMPUTES MOMENTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%there are two ways of computing first and second year moments. First,
%Dyanre computes for us theoretical moments if we ask for a first order
%linearization (as we have done in the mod file). Second, we can compute
%moments from the simulation we ask Dynare to run for us once the model is
%solved. These two should obviously coincide.

% 1. Get business cycle second moments using Dynare moments
varX=diag(oo_.var); %Variables in order: X=y,c,k,i,n,a
varX(3)=[]; %drop capital
stdX=(varX).^(0.5);
relstdX=stdX./stdX(1); %first vector of interest
covX=oo_.var(:,1);
covX(3)=[]; %drop capital
corrXY=(covX./stdX)./stdX(1); %second vector of interest
autocorrXX_1=diag(oo_.autocorr{1,1}); %third vector of interest
autocorrXX_1(3)=[]; %drop capital
tabXDyn=[relstdX corrXY autocorrXX_1]; %matrix of moments, using Dynare moments

% 2. Get business cycle second moments using the simulated series.
autocorrX = zeros(1,5); %preallocate autocorrelation vector
X = [y c inv n a];
stdevX = var(X,1).^0.5;
stdevX = stdevX/stdevX(1);
corrX = corr(X,X(:,1));
for ct=1:5
dum = corrcoef(X(2:end,ct),X(1:(end-1),ct));
autocorrX(ct)=dum(1,2);
end

%autocorrX = diag(corr(X(2:end,:),X(1:(end-1),:)));
tabX = [stdevX' corrX autocorrX']; %matrix of moments, using simulated series

% Export results table in Latex - this is optional
% rowLabels = {'Y','C','I','N','A'};
% columnLabels = {'$\sigma(X)/\sigma(Y)$', '$corr(X,Y)$','$corr(X,X(-1))$'};
% matrix2latex(tabX, 'disc3_output/momentsX.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'large');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare objects for plots
nsim=10000; % #periods in simulation
timesim=linspace(1,nsim,nsim)';
nirf=300; % #periods for IRF
timeirf=linspace(1,nirf,nirf)';
zerline=zeros(length(y_e),1); %zero line
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

%unpack steady state values
yss=oo_.steady_state(1);
css=oo_.steady_state(2);
kss=oo_.steady_state(3);
iss=oo_.steady_state(4);
nss=oo_.steady_state(5);
ass=oo_.steady_state(6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PLOTS IRFs %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure %plot IRFs
subplot(2,2,1)
plot(timeirf,100*y_e,'--b',timeirf,100*a_e,'r',timeirf,zerline,'k','LineWidth',lwidnum)
ylabel('Output'), xlabel('time')
axis tight
title('Labor')
legend('Output','Productivity','Zero Line')
legend boxoff;
set(gca,'FontSize',fsizenum)
title('Output')

subplot(2,2,2)
plot(timeirf,100*inv_e,'--b',timeirf,100*a_e,'r',timeirf,zerline,'k','LineWidth',lwidnum)
ylabel('Investment'), xlabel('time')
axis tight
title('Labor')
legend('Investment','Productivity','Zero Line')
legend boxoff;
set(gca,'FontSize',fsizenum)
title('Investment')

subplot(2,2,3)
plot(timeirf,100*c_e,'--b',timeirf,100*a_e,'r',timeirf,zerline,'k','LineWidth',lwidnum)
ylabel('Consumption'), xlabel('time')
axis tight
title('Labor')
legend('Consumption','Productivity','Zero Line')
legend boxoff;
set(gca,'FontSize',fsizenum)
title('Consumption')

subplot(2,2,4)
plot(timeirf,100*n_e,'--b',timeirf,100*a_e,'r',timeirf,zerline,'k','LineWidth',lwidnum)
ylabel('Labor'), xlabel('time')
axis tight
title('Labor')
legend('Labor','Productivity','Zero Line')
legend boxoff;
set(gca,'FontSize',fsizenum)
% print(gcf,'output_relevant/IRF.png','-dpng','-r400')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PLOTS SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure %plot unconditional simulations
subplot(2,2,1)
plot(timesim(nsim-99:end),100*(y(nsim-99:end)-yss),'b',timesim(nsim-99:end),100*(a(nsim-99:end)-ass),'r','LineWidth',lwidnum)
ylabel('Output'), xlabel('time')
axis tight
title('Output')
legend('Output','Productivity')
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,2)
plot(timesim(nsim-99:end),100*(inv(nsim-99:end)-iss),'b',timesim(nsim-99:end),100*(a(nsim-99:end)-ass),'r','LineWidth',lwidnum)
ylabel('Investment'), xlabel('time')
axis tight
title('Investment')
legend('Investment','Productivity')
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,3)
plot(timesim(nsim-99:end),100*(c(nsim-99:end)-css),'b',timesim(nsim-99:end),100*(a(nsim-99:end)-ass),'r','LineWidth',lwidnum)
ylabel('Consumption'), xlabel('time')
axis tight
title('Consumption')
legend('Consumption','Productivity')
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,4)
plot(timesim(nsim-99:end),100*(n(nsim-99:end)-nss),'b',timesim(nsim-99:end),100*(a(nsim-99:end)-ass),'r','LineWidth',lwidnum)
ylabel('Labor'), xlabel('time')
axis tight
title('Labor')
legend('Labor','Productivity')
legend boxoff;
set(gca,'FontSize',fsizenum)
% print(gcf,'output_relevant/sims.png','-dpng','-r400')






