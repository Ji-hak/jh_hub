%% THIS SCRIPT PRODUCES CHART C1 OF UNCERTAIN POLICY PROMISES PAPER
% The chart shows a monetary policy shock impulse response using the 
% FRBUS model.

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
experimentName = 'OneOffMPshock';

%% SPECIFY WHAT TO PLOT
simHorizon = 35;
LFLsimPeriods = 1;
rFixValue = 1;
policyShockMnemonic = 'etar';
policyRateMnemonic = 'rffe';
varsToPlot = {'rffe';'xgap2';'lur';'inflation'};
varNames = {'Federal funds rate';'Output gap';...
    'Unemployment rate';'Annual PCE inflation'};

%% DEFINE TEX FIGURES DIRECTORY
TeXfolderName = [thisDir,'Figures\'];

%% CREATE THE MODEL
Model = create_model('FRBUS.maps');

%% CREATE SIMULATION JOURNEY BASE
simBase = create_simulation_journey_base(Model,simHorizon);

%% RUN THROUGH THE SIMULATIONS
polRatePath = rFixValue;
instrumentUsage = [1 zeros(1,simHorizon-1)];
Judgements.AnticipatedFixes.modelVariables = ...
    {policyRateMnemonic, polRatePath};
Judgements.AnticipatedFixes.shockUsages = ...
    {policyShockMnemonic, instrumentUsage};

irfSim = ...
    impose_judgement_using_LSS_model(Model,simBase,Judgements);

%% PLOT THE RESULTS
varInfo = lookup_LSS_mnemonic_type_and_position(Model,varsToPlot);
simFigHandle = figure;
nVars = length(varsToPlot);
nRows = ceil(nVars/2);
for iVar = 1:nVars
    subplot(nRows,2,iVar);
    hold on;
    iLine = irfSim.Forecast.(varInfo{iVar,2})(varInfo{iVar,3},:);
    plot(iLine,'r','linewidth',2);
    title(varNames{iVar});
end

%% FORMAT AND SAVE FIGURES
if saveFigure
    figName = [TeXfolderName 'FRBUSexperiment' experimentName '.eps'];
    set(simFigHandle,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
    print(simFigHandle,'-depsc',figName);
end