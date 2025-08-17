%% THIS SCRIPT PRODUCES CHART B2 OF UNCERTAIN POLICY PROMISES PAPER
% The chart compares 3 alternative calibrations of the function mapping net
% benefit of renege to the probability of doing so.

%% HOUSEKEEPING
close all;
restoredefaultpath;
clear variables;
delete('*.asv');
clc;

%% DEFINE ALTERNATIVE PARAMETERS TO CHART & ASSOCIATED LINE COLOURS ETC
alfa1s = [0.075;0.075;0.075];
alfa2s = [3; 2; 1];
lineStyles = {'-rd';'-k';'-bs';':c';'-m*'};

%% DEFINE GRID OF NET BENEFITS OVER WHICH TO COMPUTE PROBABILITIES
netBenefits = (0:0.005:0.2);

%% FIGURE SAVE OPTIONS
saveTeXfigure = true;
saveFiguresInPaperFormat = true;

%% EXTRACT DIRECTORY FROM THIS FILE NAME
fullPathNameForThisFile = mfilename('fullpath');
fullPathFileNameSplitByFinalBackSlash = regexp(...
    fullPathNameForThisFile,'[^\\]+$','split');
thisDir = fullPathFileNameSplitByFinalBackSlash{1};

%% DEFINE FIGURES DIRECTORY
TeXfiguresDir = [thisDir,'Figures'];

%% COMPUTE PROBABILITIES
nParams = size(alfa1s,1);
nNetBenefits = size(netBenefits,2);
renegeProbs = NaN*ones(nParams,nNetBenefits);
for iParam = 1:nParams
    iAlpha1 = alfa1s(iParam);
    iAlpha2 = alfa2s(iParam);
    renegeProbs(iParam,:) = ...  
        1- exp(-(netBenefits/iAlpha1).^iAlpha2);
end

%% FIGURE FORMATTING OPTIONS
if saveFiguresInPaperFormat
    figPrintVec = [0 0 10 6];
    figFontSize = 12;
else
    figPrintVec = [0 0 10 6];                                               %#ok<UNRCH>
    figFontSize = 14;
end

%% CREATE LEGEND ENTRY
legStrs = cell(nParams,1);
for iParam = 1:nParams
    legStrs{iParam} = ['\alpha_{1}=',num2str(alfa1s(iParam)),...
        '; \alpha_{2}=',num2str(alfa2s(iParam))];
end

%% CREATE FIGURE
figHandle = figure;
for iParam = 1:nParams
    hold on;
    plot(netBenefits,renegeProbs(iParam,:),...
        lineStyles{iParam},'linewidth',2);
    hold off;
end
xlabel('Net benefit of reneging','FontSize',figFontSize);
ylabel('Probability','FontSize',figFontSize)
yLimVec = ylim;
xLimVec = xlim;
AxisHandle = get(figHandle,'CurrentAxes');
set(AxisHandle,'FontSize',figFontSize);
set(AxisHandle,'xLim',xLimVec);
set(AxisHandle,'yLim',yLimVec);
legend(legStrs,'location',[0.74 0.33 0.1 0.1]);
legend boxoff;

%% PRINT FIGURE FOR TEX
if saveTeXfigure
    figFullPathName = [TeXfiguresDir,'\reversionProbFunc.eps'];
    set(figHandle,'PaperUnits','inches','PaperPosition',figPrintVec);
    print(figHandle,'-depsc',figFullPathName);
end