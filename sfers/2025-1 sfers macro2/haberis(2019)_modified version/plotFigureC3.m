%% THIS SCRIPT PRODUCES CHART C2 OF UNCERTAIN POLICY PROMISES PAPER
% This script examines data from NY Fed primary dealer surveys.

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

%% SAVE FIGURE OPTION
saveFigure = true;

%% DEFINE TEX FIGURES DIRECTORY
TeXfolderName = [thisDir,'Figures\'];

%% DATA AND DATES: HARD-WIRED BECAUSE THEY ARE LOW-DIMENSIONAL
DecRawSurvey = [1; 2; 5; 10; 22; 28; 18; 8; 7];
JanRawSurvey = [0; 1; 5; 10; 21; 29; 18; 8; 4 ; 4];
datesLabels = {'13H1';'13H2';'14H1';'14H2';'15H1';'15H2';'16H1';'16H2';...
    '17H1';'17H2'};

%% RE-NORMALISE RAW DATA TO ENSURE IT SUMS TO 100
DecSurvey = 100*DecRawSurvey/sum(DecRawSurvey);
JanSurvey = 100*JanRawSurvey/sum(JanRawSurvey);

%% COMPUTE BASIC INFORMATION ABOUT SERIES
nDates = size(datesLabels,1);
nDecObs = size(DecSurvey,1);
nJanObs = size(JanSurvey,1);

%% A BASIC PLOT OF THE SERIES
DecSurveyPadded = [DecSurvey; nan(nDates-nDecObs,1)];
JanSurveyPadded = [JanSurvey; nan(nDates-nJanObs,1)];
dataToPlot = [DecSurveyPadded JanSurveyPadded];
rawDataHandle = figure;
barHandle = bar(1:nDates,dataToPlot);
set(barHandle(1),'EdgeColor','black','FaceColor',[0.95 0.95 0.95]);
set(barHandle(2),'EdgeColor','black','FaceColor',[0.15 0.15 0.15]);
set(barHandle,'BarWidth',1)
set(gca,'XTickLabel',datesLabels,'XLim',[0.5 nDates+0.5]);
legend('December 2012','January 2013','Location','NorthWest');
legend boxoff;
x1 = nDates-1.1;
y1 = JanSurvey(end)-5;
x2 = nDates-3.35;
y2 = DecSurvey(end);

%% SAVE THE FIGURES IF REQUIRED
if saveFigure
    figSaveName = [TeXfolderName,'rawPDSprobsTBFG.eps'];
    set(rawDataHandle,'PaperUnits','inches','PaperPosition',[0 0 10 4]);    
    print(rawDataHandle,'-depsc',figSaveName);
end