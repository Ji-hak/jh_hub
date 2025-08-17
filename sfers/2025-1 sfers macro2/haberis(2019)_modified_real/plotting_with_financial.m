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
resultsFileName = 'E:\github_jh\jh_hub\sfers\2025-1 sfers macro2\haberis(2019)_modified_real\Results\withfinancial';

%% SPECIFY SAVING & FORMAT OPTIONS
saveTeXfigures = true;
useSlideFormat = false; %#ok<NASGU> % (미사용)

%% SPECIFY PLOT HORIZON
plotHorizon = 14;

%% DEFINE TEX FIGURES DIRECTORY
thisDir       = fileparts(mfilename('fullpath'));   % 현재 파일 위치
TeXfiguresDir = fullfile(thisDir, 'Figures');       % 경로 결합
if ~exist(TeXfiguresDir, 'dir')
    mkdir(TeXfiguresDir);
end

%% LOAD RESULTS
load(resultsFileName,'results');

%% ========================= PLOTTING (MODIFIED) ==========================
% [변경 사항]
% - 기존에는 plot_small_model_app_results()를 호출했으나,
%   "빨강=Step 2 (Full cred, liftoff in period 11)" &
%   "파랑=Step 3 (Endog cred, liftoff in period 12)"가 되도록
%   본 스크립트에서 직접 3패널 그림을 그리도록 수정함.
%   (회색은 Step 1 Baseline)

% ---- 인덱스 가져오기 ----
rAnnIndex = unpack_model_metadata_and_lookup_index_numbers(results.Model,'Ymnems','rAnn');
piAnnIndex = unpack_model_metadata_and_lookup_index_numbers(results.Model,'Ymnems','piAnn');
yIndex    = unpack_model_metadata_and_lookup_index_numbers(results.Model,'xMnems','y');

H = plotHorizon;

% ---- Step 1: Baseline ----
baseObs = results.Baseline.Paths.modelObservables;        % [ny x T x S]
baseX   = results.Baseline.Paths.modelVariables;          % [nx x T x S]

% ---- Step 2: Full credibility, liftoff in period 11 (k)  --> 빨강 ----
fullCredK_Obs = results.LFLkPlus1.FullCred.ExpectedPaths.modelObservables;   % [ny x T x 1]
fullCredK_X   = results.LFLkPlus1.FullCred.ExpectedPaths.modelVariables;     % [nx x T x 1]

% ---- Step 3: Imperfect credibility, liftoff in period 12 (k+1)  --> 파랑 ----
impCredK1_Obs = results.LFLkPlus1.ImpCred.ExpectedPaths.modelObservables; % [ny x T x 1]
impCredK1_X   = results.LFLkPlus1.ImpCred.ExpectedPaths.modelVariables;   % [nx x T x 1]

% ---- 시계열 추출 ----
t = 1:H;

% 연율 인플레이션(%), 연율 정책금리(%), 산출갭(%)
pi_base = squeeze(baseObs(piAnnIndex,t,end));
pi_full = squeeze(fullCredK_Obs(piAnnIndex,t,1));
pi_imp1 = squeeze(impCredK1_Obs(piAnnIndex,t,1));

y_base  = squeeze(baseX(yIndex,t,end));
y_full  = squeeze(fullCredK_X(yIndex,t,1));
y_imp1  = squeeze(impCredK1_X(yIndex,t,1));

r_base = squeeze(baseObs(rAnnIndex,t,end));
r_full = squeeze(fullCredK_Obs(rAnnIndex,t,1));
r_imp1 = squeeze(impCredK1_Obs(rAnnIndex,t,1));

% ---- Plot ----
fig = figure('Color','w','Position',[100 100 980 360]);

coGrey = [0.35 0.35 0.35];
coBlue = [0.00 0.45 0.74];   % 파랑 (Step 3)
coRed  = [0.85 0.33 0.10];   % 빨강 (Step 2)

lwBase = 2.5; lwAlt = 1.8; ms = 5; % 선/마커 설정

% (1) Inflation
subplot(1,3,1); hold on; box on;
plot(t, pi_base, 'Color', coGrey, 'LineWidth', lwBase);
plot(t, pi_imp1, 's-', 'Color', coBlue, 'MarkerSize', ms, 'LineWidth', lwAlt);
plot(t, pi_full, 'o-', 'Color', coRed,  'MarkerSize', ms, 'LineWidth', lwAlt);
title('Quarterly inflation (%)'); xlabel('');
yline(0,'k-','LineWidth',0.5);

% (2) Output gap
subplot(1,3,2); hold on; box on;
plot(t, y_base, 'Color', coGrey, 'LineWidth', lwBase);
plot(t, y_imp1, 's-', 'Color', coBlue, 'MarkerSize', ms, 'LineWidth', lwAlt);
plot(t, y_full, 'o-', 'Color', coRed,  'MarkerSize', ms, 'LineWidth', lwAlt);
title('Output gap (%)'); yline(0,'k-','LineWidth',0.5);

% (3) Policy rate (annualised)
subplot(1,3,3); hold on; box on;
plot(t, r_base, 'Color', coGrey, 'LineWidth', lwBase);
plot(t, r_imp1, 's-', 'Color', coBlue, 'MarkerSize', ms, 'LineWidth', lwAlt);
plot(t, r_full, 'o-', 'Color', coRed,  'MarkerSize', ms, 'LineWidth', lwAlt);
title('Policy rate (annualised, %)');
legend({'Baseline', 'Endog cred: liftoff in period 11', 'Full cred: liftoff in period 11'}, ...
       'Location','southeast');

% 저장
if saveTeXfigures
    outName = fullfile(TeXfiguresDir, 'section3_baseline_fullcred_impcred.pdf');
    exportgraphics(fig, outName, 'ContentType','vector');
end

%% =================== PRINT NUMERICAL RESULTS (MODIFIED) ==================
stars = '******************************************************';
% [변경 사항]
% - 아래 표는 Baseline 대비 (i) Step 2: Full cred, k  (ii) Step 3: Imp cred, k+1
%   의 연율 정책금리(rAnn) 차이를 보여줌.

baseYC = results.Baseline.Paths.modelObservables(rAnnIndex,:,end);
diffYC_fullK = results.LFLk.FullCred.ExpectedPaths.modelObservables(rAnnIndex,:,1) - baseYC;   % Step 2

diffYC_impK1 = results.LFLkPlus1.ImpCred.ExpectedPaths.modelObservables(rAnnIndex,:,1) - baseYC; % Step 3

displayHorizon = 12;
resultsToShow = cell(displayHorizon+1,3);
resultsToShow(1,:) = {'Date','Full cred: K','Endog cred: K+1'}; 
resultsToShow(2:end,:) = num2cell([(1:displayHorizon)', diffYC_fullK(1:displayHorizon)', diffYC_impK1(1:displayHorizon)']);
disp(resultsToShow);

% High-cred 변형이 있으면 (원래 코드 유지)
if isfield(results.LFLk,'HighCred')
    disp(['Chance of implementing 10-quarter promise in high cred case = ',...
        num2str(results.LFLk.HighCred.ExAnteExitProbs(end))]);
end

% Full credibility losses (원문 출력 유지)
disp(stars);
disp('Losses under full credibility');
disp(['Baseline = ',num2str(results.Baseline.Losses(1,1,end))]);
disp(['K-period = ',num2str(results.LFLk.FullCred.Losses(1,1,end))]);
disp(['K+1-period = ',num2str(results.LFLkPlus1.FullCred.Losses(1,1,end))]);
disp(stars);


% 1) K의 내생신뢰 조기 종료 확률(사전확률)
results.LFLk.ImpCred.ExAnteExitProbs    % 길이 K+1 벡터
results.LFLkPlus1.ImpCred.ExAnteExitProbs

% 2) 각 시점 유인(continuation - renege)
T_K   = squeeze(results.LFLk.ImpCred.ContinuationLosses ...
               - results.LFLk.ImpCred.RenegeLosses);
T_K1  = squeeze(results.LFLkPlus1.ImpCred.ContinuationLosses ...
               - results.LFLkPlus1.ImpCred.RenegeLosses);

