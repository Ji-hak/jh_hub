%% THIS SCRIPT PRODUCES FIGURES FOR SECTION 3 OF UNCERTAIN POLICY PROMISES
% This script loads results for the simple model application section 3
% of the paper, plots charts and saves figures for TeX. The script also
% extracts and prints some of the numerical results referred to in the 
% text.

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
addpath(genpath('Models'));
addpath(genpath('Functions'));
addpath(genpath('Figures'));
addpath(genpath('Data'));
addpath(genpath('MAPSlite'));
addpath(genpath('Results'));

%% SPECIFY RESULTS FILE NAME
resultsFileName = 'E:\github_jh\jh_hub\sfers\2025-1 sfers macro2\haberis(2019)_modified_real\Results\bggresult2';

%% SPECIFY SAVING & FORMAT OPTIONS
saveTeXfigures = true;
useSlideFormat = false;

%% SPECIFY PLOT HORIZON
plotHorizon = 14;

%% DEFINE TEX FIGURES DIRECTORY
%% DEFINE TEX FIGURES DIRECTORY
thisDir       = fileparts(mfilename('fullpath'));   % 현재 파일 위치
TeXfiguresDir = fullfile(thisDir, 'Figures');       % 경로 결합

% 없으면 생성
if ~exist(TeXfiguresDir, 'dir')
    mkdir(TeXfiguresDir);
end


%% LOAD RESULTS
load('E:\github_jh\jh_hub\sfers\2025-1 sfers macro2\haberis(2019)_modified_real\Results\bggresult2');

%% PLOT THE RESULTS AND SAVE IF REQUIRED
if saveTeXfigures
    plot_small_model_app_results(...
        results,plotHorizon,TeXfiguresDir,useSlideFormat);
else
    plot_small_model_app_results(results,plotHorizon);                      %#ok<UNRCH>
end

%% PRINT NUMERICAL RESULTS USED IN MAIN TEXT
stars = '******************************************************';
% Difference in expected yields for baseline function setup.
rAnnIndex = ...
    unpack_model_metadata_and_lookup_index_numbers(results.Model,...
    'Ymnems','rAnn');
baseYC = results.Baseline.Paths.modelObservables(rAnnIndex,:,end);
diffYCk = ...
    results.LFLk.ImpCred.ExpectedPaths.modelObservables(rAnnIndex,:,1) - ...
    baseYC;
diffYCkPlus1 = ...
    results.LFLkPlus1.ImpCred.ExpectedPaths.modelObservables(rAnnIndex,...
    :,1) - baseYC;
displayHorizon = 12;
resultsToShow = cell(displayHorizon+1,3);
resultsToShow(1,:) = {'Date','LFL K','LFL K+1'}; 
resultsToShow(2:end,:) = ...
    num2cell([(1:displayHorizon)', diffYCk(1:displayHorizon)', ...
    diffYCkPlus1(1:displayHorizon)']);
disp(resultsToShow);
% Chance of implementing 10-quarter promise in high cred variant
if isfield(results.LFLk,'HighCred')
    disp(['Chance of implementing 10-quarter promise in high cred case = ',...
        num2str(results.LFLk.HighCred.ExAnteExitProbs(end))]);
end
% Full credibility losses
disp(stars);
disp('Losses under full credibility');
disp(['Baseline = ',num2str(results.Baseline.Losses(1,1,end))]);
disp(['K-period = ',num2str(results.LFLk.FullCred.Losses(1,1,end))]);
disp(['K+1-period = ',num2str(results.LFLkPlus1.FullCred.Losses(1,1,end))]);
disp(stars);