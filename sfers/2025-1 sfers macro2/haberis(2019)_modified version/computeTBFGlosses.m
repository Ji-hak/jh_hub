%% THIS SCRIPT PRODUCES LOSS CALCULATIONS FOR TBFG POLICIES IN SECTION 4
% It uses a Yellen (2012) loss function to compute losses in the baseline
% forecast and in the two alternative TBFG policies.

%% HOUSEKEEPING
close all;
restoredefaultpath;
clear variables;
delete('*.asv');
clc; 

%% ADD DIRECTORIES
fullPathNameForThisFile = mfilename('fullpath');
fullPathFileNameSplitByFinalBackSlash = regexp(...
    fullPathNameForThisFile,'[^\\]+$','split');
thisDir = fullPathFileNameSplitByFinalBackSlash{1};
addpath(genpath(thisDir));

%% SPECIFY RESULTS FILE NAMES
firstPolicyResultsFileName = 'TBFGusingFRBUSuThreshold5point75trigger';
firstPolicyName = '''trigger'' 5.75%% unemployment threshold';
secondPolicyResultsFileName = 'TBFGusingFRBUSuThreshold6openEnded';
secondPolicyName = '''open-ended'' 6%% unemployment threshold';

%% SPECIFY LOSS FUNCTION INFO
objVarMnems = {'PCEinflation';'lur';'drffe'};
objVarWeights = [1;1;1];
discountFactor = 0.99;

%% LOAD BASELINE FORECAST & EXPECTED PATHS FOR BOTH POLICIES
resultsDir = [thisDir,'Results\'];
firstPolicyResultsFullPathFileName = [...
    resultsDir,firstPolicyResultsFileName,'.mat'];
secondPolicyResultsFullPathFileName = [...
    resultsDir,secondPolicyResultsFileName,'.mat'];
[BaselineForecast,Model,FirstPolicyPaths,firstPolicyPathProbabilities] = ...
    load_variables_from_mat_file(firstPolicyResultsFullPathFileName,...
    'baselineForecast','Model','TBFGpaths','TBFGpathProbs');
[SecondPolicyPaths,secondPolicyPathProbabilities] = ...
    load_variables_from_mat_file(secondPolicyResultsFullPathFileName,...
    'TBFGpaths','TBFGpathProbs');

%% COMPUTE OUTPUT WEIGHT VIA OKUN'S LAW
% Okun's law from Yellen's speech, footnote 17: Ygap{t} = 2.3*(U* - U{t})
% Implies: U{t}-U* = -(1/2.3)*Ygap{t}
outputWeight =  objVarWeights(2,1)*(-1/2.3)^2;

%% BASELINE LOSS
[nx,H,nForecasts] = size(FirstPolicyPaths.modelVariables);
xBase = nan(nx,H,1);
xBase(:,:,1) = BaselineForecast.Forecast.modelVariables;
baseLosses = compute_yellen_losses(Model,xBase,objVarMnems,...
    objVarWeights,discountFactor);
% CONVERT TO CONSUMPTION EQUIVALENCE
baseLossesConsumptionEquiv = ...
    ((1-discountFactor)/outputWeight)^0.5*baseLosses^0.5;

%% CREATE LOSS FUNCITON HANDLE
lossFuncHandle = @(Model,xPaths) ...
    compute_yellen_losses(Model,xPaths,objVarMnems,...
    objVarWeights,discountFactor);

%% COMPUTE LOSSES ALONG PATHS FOR TWO POLICIES
firstPolicyLossesAlongPaths = ...
    compute_losses_in_sim_with_exogenous_imperfect_credibility(...
    Model,FirstPolicyPaths.modelVariables,lossFuncHandle);
secondPolicyLossesAlongPaths = ...
    compute_losses_in_sim_with_exogenous_imperfect_credibility(...
    Model,SecondPolicyPaths.modelVariables,lossFuncHandle);

%% COMPUTE EX ANTE LOSS
[~,firstPolicyLossConsumptionEquiv] = compute_ex_ante_TBFG_loss(...
    firstPolicyLossesAlongPaths,firstPolicyPathProbabilities,...
    discountFactor,outputWeight);
[~,secondPolicyLossConsumptionEquiv] = compute_ex_ante_TBFG_loss(...
    secondPolicyLossesAlongPaths,secondPolicyPathProbabilities,...
    discountFactor,outputWeight);

%% LOSS DIFFERENCES
firstPolicyLossDiff = ...
    baseLossesConsumptionEquiv-firstPolicyLossConsumptionEquiv;
secondPolicyLossDiff = ...
    baseLossesConsumptionEquiv-secondPolicyLossConsumptionEquiv;

%% PRINT RESULTS
fprintf(['Loss difference for ',firstPolicyName,' policy is: %3.2f\n'],...
    firstPolicyLossDiff);
fprintf(['Loss difference for ',secondPolicyName,' policy is: %3.2f\n'],...
    secondPolicyLossDiff);