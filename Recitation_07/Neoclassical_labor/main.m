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
path('/Applications/MathWorks/relaxsys',oldpath)

tic

globalpar                           % Initializes the global parameter

parini                              % Loads the Parameter Values

relaxsetting                        % Loads the settings for the Relaxation algorithm, i.e. dimensions, boundary conditions etc. 

% Converts settings to a form suitable for relax.m
[guess, start, errorcode]=initrelax(@funcODE, @funcSTAT, n, n1, n3, nu, y, M, statev);       


if errorcode==0            % Executes the relaxation algorithm if no error occured during initilization
    
    [t, x]=relax(@funcODE, @funcSTAT, @funcINI, @funcfinal, n, n1, n3, nu, guess, M, start, Endcond, maxit, tol, damp, dampfac);    
    
    %Normalization of specified variables
    for i=1:M
        x(normal,i)=x(normal,i)./x(normal,end);
    end
    
    varex                               % Extracts the variables and stores them in the memory
 
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


%for plots, I only prepare them for the case of solution with initial
%conditions


lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs

figure % Plot functions of interest
subplot(3,1,1);  
plot(t(1:50),c(1:50), 'b', t(1:50),css*ones(1,50),'-r','LineWidth',lwidnum)
ylabel('Consumption'), xlabel('Time')
legend('Transition', 'SS','Location','best')
%axis([0 30 -inf inf]);
legend boxoff
set(gca,'FontSize',fsizenum)

subplot(3,1,2)
plot(t(1:50),k(1:50), 'b', t(1:50),kss*ones(1,50),'-r','LineWidth',lwidnum)
ylabel('Capital'), xlabel('Time')
set(gca,'FontSize',fsizenum)

subplot(3,1,3)
plot(t(1:50),l(1:50), 'b', t(1:50),lss*ones(1,50),'-r','LineWidth',lwidnum)
ylabel('Labor'), xlabel('Time')
set(gca,'FontSize',fsizenum)
 