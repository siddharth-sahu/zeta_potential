% Code for Zeta Potential calculation from streaming current, version 17 - by Siddharth S. Sahu

%% Data Import

clc;
clear;
close all;
tic; % to measure the total processing time (ended with toc)
warning('off','all');


% input file name along with extension
input_file = 'wet-current-24hrwater-1-500mbar.nwI';
delimiter = '\t';

% store the x and y data into the vectors time and current respectively
startRow = 7;       % as the first row is not data

% specifiy that the number of rows in the data file
% formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]'; % use this for the capillary setup
formatSpec = '%f%f%[^\n\r]';  % use this for the nanowire setup

filename=[input_file];      % specify the input file name

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);  
time_in = dataArray{:, 1};     % stores the time data
current_in = dataArray{:, 2};   %stores the current data
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% Plotting the raw data
figure(1);
plot(time_in/60,current_in/(10^-9),'o-');
title("Raw Data") ;
xlabel("Time(min)");
ylabel("Current (nA)");
set(gca,'fontsize',12);
grid on;
file1=[input_path,'1 raw_data'];

%% Data Selection

% Select chunks of data by specifying the starting and ending points
timeSelect = ([time_in(1876:end)]);
currentSelect = ([current_in(1876:end)]);


% Eliminating outliers in each square pulse
intraPulseNoise= 2*(0.5592-0.558)/0.558;
interPulseNoise= 3*(0.5592-0.558)/0.558;
j=1;

for i=1+1:length(currentSelect)-1
    relativeForwardChange = (currentSelect(i)-currentSelect(i-1))/currentSelect(i);
    relativeBackwardChange = (currentSelect(i)-currentSelect(i+1))/currentSelect(i);
    if ((relativeForwardChange*relativeBackwardChange < intraPulseNoise^2) && ...
            (relativeForwardChange*relativeBackwardChange > -1*interPulseNoise^2))
        % the data point is tested for noisiness by comparing it to the two adjacent data points
        timeRefine(j)=timeSelect(i);
        currentRefine(j)=currentSelect(i);
        j=j+1;
    end
end

dataReduction = 100*(1-length(timeRefine)/length(timeSelect));

% Plotting the selected and refined data
figure(2);
plot(timeRefine/60,currentRefine/(10^-9),'o','color','[0.2 0.5 0.1]');
title("Refined Data");
xlabel("Time(min)");
ylabel("Current (nA)");
set(gca,'fontsize',12);
h = line(nan, nan, 'Color', 'none');
legend(h, ['Data reduction = ' num2str(dataReduction) '%'], 'Location', 'best')
grid on;
file2=[input_path,'2 selected_data'];

%% Data Segregation 

%Specify if the data begins from the lower (0) or higher (1) flow
switch_s = 0;

% indices and temporary vaiables
startCheck=switch_s;

% Specify how many initial points of the pulse to skip, for both the lower
% and higher flow rate portions of the data
lowerInitialSkip=0;
upperInitialSkip=1;

% Select a threshold of separation and a 2nd purifier for the data
thresh=(0.1147-0.1112)*1e-9*0.4;

timeLow = {};
currentLow = {};
timeHigh = {};
currentHigh = {};

% separating the data 
for i=1+lowerInitialSkip:length(timeRefine)-upperInitialSkip-1
    if(switch_s == 0)
        i = i + lowerInitialSkip;
        timeLow = [timeLow, timeRefine(i)];
        currentLow = [currentLow, currentRefine(i)];
    elseif (switch_s == 1)
        i = i + upperInitialSkip;
        timeHigh = [timeHigh, timeRefine(i)];
        currentHigh = [currentHigh, currentRefine(i)];
    end
    if( (currentRefine(i+1)-currentRefine(i)) > thresh)  % switching to the stepped up data
        switch_s=1;
    elseif ( (currentRefine(i+1)-currentRefine(i)) < (-1*thresh) ) % switching to the stepped down data
        switch_s=0;
    end
end

timeLow = cell2mat(timeLow);
currentLow = cell2mat(currentLow);
timeHigh = cell2mat(timeHigh);
currentHigh = cell2mat(currentHigh);

% Plotting the segregated data
figure(3);
scatter(timeLow/60, currentLow/(10^-9));
hold on;
scatter(timeHigh/60, currentHigh/(10^-9));
legend('low', 'high', 'Location', 'best');
title("Segregated Data");
xlabel("Time(min)");
ylabel("Current (nA)");
set(gca,'fontsize',12);
grid on;
file3=[input_path,'3 segregated data'];
hold off;

%% Averaging
% allocating initial values to the summing and loop variables
time_sum=0;
current_sum=0;
j=0;

timeAvLow={};
currentAvLow={};
timeAvHigh={};
currentAvHigh={};

% averaging over each set of stepped down data
for i=1:length(timeLow)-1
    if (timeLow(i+1)-timeLow(i) < 10)
        time_sum = time_sum + timeLow(i);
        current_sum = current_sum + currentLow(i);
        j=j+1;
    elseif (j==0)
        timeAvLow = [timeAvLow, timeLow(i)] ;
        currentAvLow = [currentAvLow, currentLow(i)];
    else 
        timeAvLow = [timeAvLow, time_sum/j] ;
        currentAvLow = [currentAvLow, current_sum/j];
        j=0;
        time_sum=0;
        current_sum=0;
    end
end 


time_sum=0;
current_sum=0;
i=1;
j=0;

% averaging over each set of stepped up data
for i=1:length(timeHigh)-1
    if (timeHigh(i+1)-timeHigh(i) < 10)
        time_sum = time_sum + timeHigh(i);
        current_sum = current_sum + currentHigh(i);
        j=j+1;
    elseif (j==0)
        timeAvHigh = [timeAvHigh, timeHigh(i)];
        currentAvHigh = [currentAvHigh, currentHigh(i)];
    else
        timeAvHigh = [timeAvHigh, time_sum/j] ;
        currentAvHigh = [currentAvHigh, current_sum/j];
        j=0;
        time_sum=0;
        current_sum=0;
    end
end 

timeAvLow = cell2mat(timeAvLow);
currentAvLow = cell2mat(currentAvLow);
timeAvHigh = cell2mat(timeAvHigh);
currentAvHigh = cell2mat(currentAvHigh);

% Plotting the average data
figure(4);
plot(timeAvLow/60, currentAvLow/(10^-9),'s','color','[0.1 0.5 0.7]','linewidth',3.5);
hold on;
plot(timeAvHigh/60, currentAvHigh/(10^-9),'s','color','[0.7 0.5 0.1]','linewidth',3.5);
legend('average low', 'average high','Location', 'best');
title("Average over each cycle");
xlabel("Time(min)");
ylabel("Current (nA)");
set(gca,'fontsize',12);
grid on;
file4=[input_path,'4 average'];
hold off;

%% Piecewise non-linear fit of the averaged data

fitOrder=3;
fitWindow=4;
% fitWindow must be greater than fitOrder!

timeFit={};
currentFitLow={};
currentFitHigh={};

if (length(timeAvLow) < length(timeAvHigh))
    tmp1=length(timeAvLow);
else
    tmp1=length(timeAvHigh);
end


if(startCheck==0)
    for k=2:tmp1-fitWindow-1
        timeFit = [timeFit, timeAvLow(k+2)];
        currentFitLow = [currentFitLow, currentAvLow(k+2)];
        
        for m=1:fitWindow
            x_range(m)=timeAvHigh(m+k-1);
            y_range(m)=currentAvHigh(m+k-1);
        end
        
        P = polyfit(x_range,y_range,fitOrder);
        currentFitHigh = [currentFitHigh, polyval(P, timeAvLow(k+2))];
        currentFitHigh = [currentFitHigh, currentAvHigh(k+2)];
        timeFit = [timeFit, timeAvHigh(k+2)];
    
        for m=1:fitWindow
            x_range(m)=timeAvLow(m+k-1);
            y_range(m)=currentAvLow(m+k-1);
        end
    
        P=polyfit(x_range, y_range, fitOrder);
        currentFitLow = [currentFitLow, polyval(P, timeAvHigh(k+2))];
    end
else    
    for k=2:tmp1-fitWindow-1
        timeFit = [timeFit, timeAvHigh(k+2)];
        currentFitHigh = [currentFitHigh, currentAvHigh(k+2)];
        
        for m=1:fitWindow
            x_range(m)=timeAvLow(m+k-1);
            y_range(m)=currentAvLow(m+k-1);
        end
        
        P = polyfit(x_range,y_range,fitOrder);
        currentFitLow = [currentFitLow, polyval(P, timeAvHigh(k+2))];
        currentFitLow = [currentFitLow, currentAvLow(k+2)];
        timeFit =[timeFit, timeAvLow(k+2)];
    
        for m=1:fitWindow
            x_range(m)=timeAvHigh(m+k-1);
            y_range(m)=currentAvHigh(m+k-1);
        end
    
        P = polyfit(x_range,y_range,fitOrder);
        currentFitHigh = [currentFitHigh, polyval(P, timeAvLow(k+2))];
    
    end    
end

timeFit = cell2mat(timeFit);
currentFitLow = cell2mat(currentFitLow);
currentFitHigh = cell2mat(currentFitHigh);

clearvars fitOrder fitWindow;

% Plotting the fits. Each intervals is marked by *
figure(5);
plot (timeFit/60, currentFitLow/(10^-9), '*-','color','m');
hold on;
plot (timeFit/60, currentFitHigh/(10^-9), '*-','color','b');
title("Fitted Data");
xlabel("Time(min)");
ylabel("Current (nA)");
set(gca,'fontsize',12);
legend('low-fit','high-fit','Location', 'best');
grid on;
file5=[input_path,'5 fitted_data'];
hold off;

%% Delta I

for i=3:length(currentFitHigh)-4
    delta_I_x(i-2)=timeFit(i);
    delta_I_y(i-2) = currentFitHigh(i) - currentFitLow(i);
end


figure(6);
plot (delta_I_x/60, delta_I_y/(10^-12),'color','[0.2 0.65 0.9]', 'linewidth', 1.2);
% Uncomment this to add the noise level to the plot
% AddText2Plot(['Noise: ',num2str(noise),' pA'],'Location','best');
title("\Delta I");
xlabel("Time(min)");
ylabel("Current (pA)");
set(gca,'fontsize',12);
grid on;
file6=[input_path,'6 delta_I'];
hold off;

%% Zeta

visc = 0.8 * 10^-3;
L = 4 * 10^-2 ;
A = 490 * 10^-12;
perm = 681.85 * 10^-12;
delta_P = (3.0-1.5) * 10^5;
zeta =1000*delta_I_y*L*visc/(A*perm*delta_P);

% %Use this to calculate the noise (std. dev.) after specifying the range of
% %points
% i = 40; % starting index of noise window
% j=58;   % ending index of noise window
% noise=nanstd(zeta(i:j));
% 
% % Drifted adjusted noise (std. dev.)
% i = 40; % starting index of noise window
% j=58;   % ending index of noise window
% noise_fit = polyfit(delta_I_x(i:j), zeta(i:j), 1);
% for k=i:j
%     zeta_adjusted(k-i+1) = zeta(k) - noise_fit(1)*delta_I_x(k);
% end
% noise_adjusted = nanstd(zeta_adjusted);


figure(7);
plot(delta_I_x/60, zeta,'color','k', 'Linewidth',1.15);
% % Uncomment this to show the fitted portion used for estimating drift
% % adjusted noise
% hold on;
% plot (delta_I_x(i:j)/60, polyval(noise_fit, delta_I_x(i:j)), 'color','r');
% % Uncomment this to add the noise level to the plot
% AddText2Plot(['Drift adjusted noise: ',num2str(noise_adjusted),' mV'],'Location','best');
title("Zeta Potential");
xlabel("Time(min)");
ylabel("Zeta Potential (mV)");
set(gca,'fontsize',12);
grid on;
file7=[input_path,'7 Zeta'];
hold off;

figure(8);
% Plotting Zeta with with Savitzky-Golay filter
plot(delta_I_x/60, sgolayfilt(zeta,2,21),'color','k', 'Linewidth',1.15);
title("Filtered Zeta Potential");
xlabel("Time(min)");
ylabel("Zeta Potential (mV)");
set(gca,'fontsize',12);
grid on;
file8=[input_path,'8 Filtered Zeta'];
hold off;

toc;

%% Export graphs

% Run this section after activating the lines corresponding to the figures
% you want to save in high resolution. Turned off by default for faster code execution.

% print(figure(1), file1,'-dpng','-r1000');  % "raw data" in PNG format 
% print(figure(2),file2,'-dpng','-r1000');   % "selected data" in PNG format
% print(figure(3),file3,'-dpng','-r1000');   % "segregated data in PNG format
% print(figure(4),file4,'-dpng','-r1000');   % "averaged data" in PNG format
% print(figure(5),file5,'-dpng','-r1000');   % "fitted data" in PNG format
% print(figure(6),file6,'-dpng','-r1000');   % "delta I" in PNG format
% saveas(figure(6),file6);                   % "delta I" FIG format
% print(figure(7),file7,'-dpng','-r1000');     % "original zeta" in PNG format
% print(figure(8),file8,'-dpng','-r1000');     % "filtered zeta" in PNG format
% saveas(figure(8),file8);                   % "filtered zeta" in FIG format
