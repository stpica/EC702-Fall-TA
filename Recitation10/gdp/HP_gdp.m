%%%%%%%%%%%%%
% HP_gdp.m
% HP filter of GDP, plot trend and cyclical component
% Stefano Pica, TA for EC702. stpica@bu.edu
% Fall 2019
%%%%%%%%%%%%%
close all; clear; clc;

%read data from xls
data = xlsread('data.xls','gdp'); % Quarterly real GDP, Billions of Chained 2009 Dollars, Seasonally Adjusted Annual Rate. From January 1947 to July 2017
gdp=data(:,2); 
lgdp=log(gdp); % log gdp
[trend,cyclical] = hpfilter(lgdp,1600); %use HP filter with lambda=1600 (appropriate for quarterly data). By construction, gdp=trend+cyclical


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PREPARES PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

date = datetime(1947,01,01):calmonths(3):datetime(2017,07,01); %create corresponding dates (quarterly). From January 1947 to July 2017
zerline=zeros(length(cyclical),1); %zero line
lwidnum = 2; %line width on graphs
fsizenum = 14; %font size on graphs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK PLOTS GDP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure %plot raw GDP data
plot(date,gdp,'b','LineWidth',lwidnum)
recessionplot
ylabel('gdp'), xlabel('time')
axis tight
title('Raw GDP Data')
set(gca,'FontSize',fsizenum)
% print(gcf,'output_relevant/rawgdp.png','-dpng','-r400')


figure %plot gdp, trend, and cyclical component
subplot(1,2,1)
plot(date,lgdp,'b',date,trend,'r','LineWidth',lwidnum)
recessionplot
ylabel('gdp'), xlabel('time')
axis tight
legend('Log GDP','Trend')
legend boxoff
set(gca,'FontSize',fsizenum)
title('Log GDP and Trend')

subplot(1,2,2)
plot(date,100*cyclical, 'b',date,zerline,'k','LineWidth',lwidnum)
recessionplot
ylabel('gdp'), xlabel('time')
axis tight
title('Cyclical Component of GDP')
set(gca,'FontSize',fsizenum)
% print(gcf,'output_relevant/hp_gdp.png','-dpng','-r400')

