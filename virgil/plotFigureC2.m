%% THIS SCRIPT PRODUCES CHART C2 OF UNCERTAIN POLICY PROMISES PAPER
% This script draws plot of pre and post announcement yield curves for
% FOMC threshold-based guidance announcement.

%% HOUSEKEEPING
close all;                                                                  
%restoredefaultpath;
clear variables;
delete('*.asv');
clc; 

%% ADD DIRECTORIES
fullPathNameForThisFile = mfilename('fullpath');
fullPathFileNameSplitByFinalBackSlash = regexp(...
    fullPathNameForThisFile,'[^\\]+$','split');
thisDir = fullPathFileNameSplitByFinalBackSlash{1};
addpath(genpath(thisDir));

%% SPECIFY FIGURE SAVING & NAMING OPTION
saveFigure = true;
figName = 'TBFGannouncementCurve';

%% DEFINE TEX FIGURES DIRECTORY
TeXfolderName = [thisDir,'Figures\'];

%% SETTINGS
xAxis = 0:0.25:4.5;
H = length(xAxis);

%% DATA
TBFGdata = [  0.1425    0.1392
            0.1265    0.1235
            0.1198    0.1232
            0.1307    0.1332
            0.1500    0.1491
            0.1604    0.1604
            0.1709    0.1764
            0.1957    0.2093
            0.2445    0.2665
            0.3181    0.3476
            0.4103    0.4461
            0.5142    0.5549
            0.6230    0.6671
            0.7331    0.7791
            0.8439    0.8906
            0.9553    1.0015
            1.0671    1.1117
            1.1789    1.2212
            1.2906    1.3299
            1.4019    1.4379
            1.5127    1.5450];

%% PLOT
figHandle = figure;
hold on;
plot(xAxis,TBFGdata(1:H,1)','b','linewidth',2);
plot(xAxis,TBFGdata(1:H,2)','co','linewidth',2);
% Legends
legend('Pre-announcement','Post-announcement','location','NorthWest');
legend boxoff;
xlabel('Years');
title('December 2012 threshold-based guidance');

%% SAVE
if saveFigure
    figDimensions = [0 0 10 5];
    set(figHandle,'PaperUnits','inches', 'PaperPosition',figDimensions);
    figSaveName = [TeXfolderName,figName,'.eps'];
    print(figHandle,'-depsc',figSaveName);
end