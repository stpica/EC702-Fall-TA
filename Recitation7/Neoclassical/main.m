% Version 3.1
% RELAXATION algorithm to solve infinite-horizon continuous time models.
%
% Description of procedure:
% Trimborn, Koch, Steger (2007) Multidimensional Transitional Dynamics: A
% Simple Numerical Procedure, forthcoming in Macroeconomic Dynamics
% 
% Copyright by Trimborn, Koch, Steger, 2008
%
% For further information contact Timo Trimborn, University of Hannover
% or visit http://www.relaxation.uni-siegen.de


% This main file solves the Neoclassical growth model in continuous time
% that has been shown in class.


clear; close all; clc
disp('Initialize Relaxation algorithm');

oldpath = path; %define path for system files
path('/Applications/MathWorks/relaxsys',path)


tic

% Initializes the global parameter:
globalpar

% Loads the Parameter Values:
parini

% Loads the settings for the Relaxation algorithm, i.e. dimensions,
% boundary conditions etc.:
relaxsetting

% Converts settings to a form suitable for relax.m:
[guess, start, errorcode] = initrelax(@funcODE, @funcSTAT, n, n1, n3, nu, y, M, statev);       

% Executes the relaxation algorithm if no error occured during
% initilization:
if errorcode==0
    
    [t, x]=relax(@funcODE, @funcSTAT, @funcINI, @funcfinal, n, n1, n3, nu, guess, M, start, Endcond, maxit, tol, damp, dampfac);    
    
    %Normalization of specified variables
    for i=1:M
        x(normal,i)=x(normal,i)./x(normal,end);
    end
    
    % Extracts the variables and stores them in the memory:
    varex
 
%    If you want to calculate and display the eigenvalues at the steady
%    state, remove the comment of the two subsequent lines.
%    [EVa EVe Jac]=eigDAS(@funcODE, @funcSTAT, x(:,end));
%    disp(['Eigenvalues: ',num2str(EVa')]);disp([' ']);
    
end


%Calculation time
time=toc;
timesec=mod(time,60);
timemin=floor(time/60);
disp(['Calculation time: ',num2str(time),' seconds (',num2str(timemin),' min ',num2str(timesec),' sec)']);

%To get a first impression of the results, remove the comment of the
%subsequent line
%plotrelax(t, x, n1, 100)


lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs


%Start plots
if statev==ones(N,1) %transition after MIT shock

%This is to plot the new SS and the transition from the old SS and the new
%SS (after the MIT shock). If there is no MIT shock, then this plot
%delivers the same result if the initial condition is exactly the SS value.
figure % Plot functions of interest
subplot(2,1,1)
plot(t(1:30),c(1:30), 'b', t(1:30),css*ones(1,30),'-r', t(1:30),css_old*ones(1,30),'--r','LineWidth',lwidnum)
ylabel('Consumption'), xlabel('Time')
legend('Transition', 'New SS', 'Old SS','Location','best')
%axis([0 30 -inf inf]);
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,1,2)
plot(t(1:30),k(1:30), 'b', t(1:30),kss*ones(1,30),'-r', t(1:30),kss_old*ones(1,30),'--r','LineWidth',lwidnum)
ylabel('Capital'), xlabel('Time')
set(gca,'FontSize',fsizenum)

%print('transitioncont','-djpeg')

elseif statev==0 %transition from initial condition

figure % Plot functions of interest
subplot(2,1,1)
plot(t(1:30),c(1:30), 'b', t(1:30),css*ones(1,30),'-r','LineWidth',lwidnum)
ylabel('Consumption'), xlabel('Time')
legend('Transition', 'New SS','Location','best')
%axis([0 30 -inf inf]);
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(2,1,2)
plot(t(1:30),k(1:30), 'b', t(1:30),kss*ones(1,30),'-r','LineWidth',lwidnum)
ylabel('Capital'), xlabel('Time')
set(gca,'FontSize',fsizenum)

end
