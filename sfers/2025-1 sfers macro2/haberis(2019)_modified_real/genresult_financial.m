%% THIS SCRIPT PRODUCES RESULTS FOR SECTION 3 OF UNCERTAIN POLICY PROMISES
% The following results are computed:
% -- Fully credible lower for longer policy (as in Figure 2)
% -- Endogenously imperfectly credible versions of the same "lower for
% longer" policies (as in Figure 3).
% -- Sensitivity analysis for the function mapping renege incentives to
% renege probabilities (as in the two equivalents to Figure 3 in the 
% appendix).

%% HOUSEKEEPING
close all;
%%restoredefaultpath;
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

%% HIGH LEVEL OPTIONS
saveResults = true;
doSensitivityAnalysis = false;

%% DEFINE MODEL
modelFileName = 'E:\github_jh\jh_hub\sfers\2025-1 sfers macro2\haberis(2019)_modified_real\Models\bgg_ver2.maps';

%% INFO ABOUT LOSS FUNCTION
objVarMnems = {'pie';'y';'lev'};    
objVarWeights = [1; 0.025/10; 0.025/10;];
beta = 0.99;
% Loss function handle
lossFuncHandle = @(Model,xPaths) ...
    compute_scaled_square_root_losses(...
    Model,xPaths,objVarMnems{1},objVarMnems{2},objVarWeights(2),beta);

%% SET INFORMATION ABOUT POLICY VARIABLE, EQUATIONS, LOWER BOUND & SHOCKS
policyEqName = 'Monetary policy rule';
policyVarMnem = 'rn';
policyShockMnem = 'etamp';
policyLBval = 0.25*100*(1-(1/beta)^4); % NB: quarterly units
% policyLBval = 0;

%% SET OPTIONS FOR POLICY EXPERIMENTS
H = 15; % Simulation horizon
shockMnem = 'etarn';
% shockSize = -12/4;
shockSize = -16/4;
nLFLperiods = 10;                       % Liftoff in period 11 case

%% SET CALIBRATION OF PROBABILITY MAPPING FUNCTION FOR MAIN & SENSITIVITY
alfa1 = 0.075;
alfa2 = 2
if doSensitivityAnalysis
    alfa1HighCred = alfa1;
    alfa2HighCred = 3;
    alfa1LowCred = alfa1;
    alfa2LowCred = 1;
end

%% SET RESULTS FILE NAME
saveFileName = [thisDir,'Results\withfinancial.mat'];

%% STEP 0: SET-UP
% Create model and unpack
Model = create_model(modelFileName);
[xMnems,zMnems,B,PHI,F] = unpack_model(Model,...
    {'xMnems';'zMnems';'B';'PHI';'F'});
shockIndex = lookup_model_index_numbers(zMnems,shockMnem);
[rVarType,rVarIndex] = ...
    lookup_LSS_model_variable_type_and_index_number(Model,policyVarMnem);
% Create steady state base
ssBase = create_simulation_journey_base(Model,H);
% %% STEP 1: BASELINE SIMULATION (1기 ELB Fallback 포함)
% simBase = ssBase;
% simBase.Forecast.Shocks.anticipated(shockIndex,1) = shockSize;
% 
% pSimBase    = 0;
% rPegSimBase = policyLBval;
% 
% try
%     basePaths = execute_policy_sim_with_imperfect_credibility(...
%         Model, simBase, pSimBase, rPegSimBase, policyLBval, ...
%         policyEqName, policyVarMnem, policyShockMnem);
% catch ME
%     if contains(ME.message,'Quadratic programming problem')
%         % ─── QP 실패 시 수동 바인딩 ───
%         % 1) Past 모델 변수를 [nx×nsim] → [nx×1×nsim] 으로 reshape
%         tmpX = ssBase.Past.modelVariables;        % size = [nx × nsim]
%         [nx, nsim] = size(tmpX);
%         initialX = reshape(tmpX, [nx, 1, nsim]);  % size = [nx × 1 × nsim]
% 
%         % 2) 빈 경로 할당
%         basePaths = struct();
%         basePaths.modelVariables = nan(nx, H, nsim);
% 
%         % 3) t = 1: steady state + ZLB binding
%         basePaths.modelVariables(:,1,:) = initialX;
%         basePaths.modelVariables(rVarIndex,1,:) = policyLBval;
% 
%         % 4) t = 2…H: 상태전이 + Taylor rule
%         for t = 2:H
%             Xprev = squeeze(basePaths.modelVariables(:,t-1,:));   % [nx×nsim]
%             z_t   = simBase.Forecast.Shocks.anticipated(:,t);     % [nz×1]
%             basePaths.modelVariables(:,t,:) = B*Xprev + PHI*z_t;
%             basePaths.modelVariables(rVarIndex,t,:) = ...
%                 F(rVarIndex,:)*squeeze(basePaths.modelVariables(:,t,:));
%         end
%     else
%         rethrow(ME);
%     end
% end
% 
% % 5) 손실 함수 계산
% baseLosses = nan(1,nLFLperiods+1);
% for t = 1:nLFLperiods+1
%     baseLosses(t) = lossFuncHandle(...
%         Model, basePaths.modelVariables(:,t:H,end) );
% end



%% STEP 1: CREATE THE BASELINE SIMULATION
% In this simulation, monetary policy is set according to the conventional
% reaction function throughout, though subject to the ZLB.
% First add shock to steady state base
simBase = ssBase;
simBase.Forecast.Shocks.anticipated(shockIndex,1) = shockSize;
% Simulation base solved using fully credible promise to hold rates at ZLB
% for one period (code will compute duration for which ZLB binds (which is 
% longer than one period in this example!).
pSimBase = 0;
rPegSimBase = policyLBval;
basePaths = execute_policy_sim_with_imperfect_credibility(...
    Model,simBase,pSimBase,rPegSimBase,policyLBval,policyEqName,...
    policyVarMnem,policyShockMnem);
% Compute losses along base path
baseLosses = nan(1,nLFLperiods+1);
for t = 1:nLFLperiods+1
    baseLosses(1,t) = lossFuncHandle(...
        Model,basePaths.modelVariables(:,t:H,end));
end



%% STEP 2: LOWER FOR LONGER WITH FULL CREDIBILITY
% These simulations are conducted under full credibility so that we can
% build intuition about the incentive to renege. This is done for two
% experiments one in which the policy rate is held at the ZLB for
% nLFLperiods and the other when it is held there for an additional period.
% Simulate LFL announcement under full credibility
rLFL = policyLBval*ones(1,nLFLperiods);
pLFL = zeros(1,nLFLperiods);   % zero probability of reneging
[fullCredLFLpaths,fullCredLFLpathProbs,fullCredLFLexpectedPaths,...
    fullCredLFLexitProbs] = ...
    execute_policy_sim_with_imperfect_credibility(...
    Model,simBase,pLFL,rLFL,policyLBval,policyEqName,...
    policyVarMnem,policyShockMnem);
lossesAlongPathsFullCred = ...
        compute_losses_in_sim_with_exogenous_imperfect_credibility(...
        Model,fullCredLFLpaths.modelVariables,lossFuncHandle);   
% Extend LFL horizon by one quarter
rLFLkPlus1 = policyLBval*ones(1,nLFLperiods+1);
pLFLkPlus1 = zeros(1,nLFLperiods+1);
[fullCredLFLpathsKplus1,fullCredLFLpathProbsKplus1,...
    fullCredLFLexpectedPathsKplus1,fullCredLFLexitProbsKplus1] = ...
    execute_policy_sim_with_imperfect_credibility(...
    Model,simBase,pLFLkPlus1,rLFLkPlus1,policyLBval,policyEqName,...
    policyVarMnem,policyShockMnem);
lossesAlongPathsFullCredKplus1 = ...
        compute_losses_in_sim_with_exogenous_imperfect_credibility(...
        Model,fullCredLFLpathsKplus1.modelVariables,lossFuncHandle);
    
%% SETP 3: IMPERFECT CREDIBILITY MAIN RESULTS
% Set probability mapping function handle
probFuncHandle = @(temptationToRenege) ...
    compute_reversion_probabilities_using_Weibull_distribution(...
    temptationToRenege,alfa1,alfa2);
% Define function handle to update probabilities
probUpdateFuncHandle = @(simResults) ...
    update_renege_probabilities_using_loss_comparison(simResults,...
    lossFuncHandle,probFuncHandle);
% (i) Case with liftoff in period 11
rLFL = policyLBval*ones(1,nLFLperiods);
% Make a guess for the equilibrium probabilities -- this is a hardcoded 
% informed guess ... an arbitrary guess would work in most cases
pLFL = [0,0,0,0.0214,0.1246,0.0880,0.0334,0.0086,0.0009,0];  
[LFLpaths,LFLpathProbs,LFLexpectedPaths,LFLexitProbs] = ...
    execute_policy_sim_with_imperfect_credibility(...
    Model,simBase,pLFL,rLFL,policyLBval,policyEqName,...
    policyVarMnem,policyShockMnem,probUpdateFuncHandle);
% Extract ex ante exit probabilities
exanteLFLexitProbs = ...
    reshape(LFLpathProbs(1,1,:),1,nLFLperiods+1);
% Compute liftoff date probabilities
sumLiftOffDateProbs = ...
    compute_liftoff_probabilities(LFLpaths,rVarType,rVarIndex,...
    policyLBval,exanteLFLexitProbs);
% Compute renege and continuation losses
lossesAlongPaths = ...
    compute_losses_in_sim_with_exogenous_imperfect_credibility(...
    Model,LFLpaths.modelVariables,lossFuncHandle);
[renegeLosses,continuationLosses,netLosses] = ...
    compute_renege_and_continuation_losses(lossesAlongPaths,...
    LFLpathProbs);
% (ii) Case iwith liftoff in period 12 ... same steps as (i)
pLFLkPlus1 = [0,0,0,0.0531,0.2878,0.3449,0.1686,0.0596,0.0142,0.0014,0];
[LFLpathsKplus1,LFLpathProbsKplus1,LFLexpectedPathsKplus1,...
    LFLexitProbsKplus1] = ...
    execute_policy_sim_with_imperfect_credibility(...
    Model,simBase,pLFLkPlus1,rLFLkPlus1,policyLBval,policyEqName,...
    policyVarMnem,policyShockMnem,probUpdateFuncHandle);
exanteLFLexitProbsKplus1 = ...
    reshape(LFLpathProbsKplus1(1,1,:),1,(nLFLperiods+1)+1);
sumLiftOffDateProbsKplus1 = ...
    compute_liftoff_probabilities(LFLpathsKplus1,rVarType,rVarIndex,...
    policyLBval,exanteLFLexitProbsKplus1);
lossesAlongPathsKplus1 = ...
    compute_losses_in_sim_with_exogenous_imperfect_credibility(...
    Model,LFLpathsKplus1.modelVariables,lossFuncHandle);
[renegeLossesKplus1,continuationLossesKplus1,netLossesKplus1] = ...
    compute_renege_and_continuation_losses(lossesAlongPathsKplus1,...
    LFLpathProbsKplus1);

%% SETP 3a: IMPERFECT CREDIBILITY - HIGH CREDIBILITY VARIANT
if doSensitivityAnalysis
    % Follows the same recipe as per main results from step 3
    probFuncHandleHighCred = @(temptationToRenege) ...
        compute_reversion_probabilities_using_Weibull_distribution(...
        temptationToRenege,alfa1HighCred,alfa2HighCred);
    probUpdateFuncHandleHighCred = @(simResults) ...
        update_renege_probabilities_using_loss_comparison(simResults,...
        lossFuncHandle,probFuncHandleHighCred);    
    % (i) Case with liftoff in period 11
    pLFLhighCred = zeros(1,nLFLperiods);
    [LFLpathsHighCred,LFLpathProbsHighCred,LFLexpectedPathsHighCred,...
        LFLexitProbsHighCred] = ...
        execute_policy_sim_with_imperfect_credibility(...
        Model,simBase,pLFLhighCred,rLFL,policyLBval,policyEqName,...
        policyVarMnem,policyShockMnem,probUpdateFuncHandleHighCred);
    exanteLFLexitProbsHighCred = ...
        reshape(LFLpathProbsHighCred(1,1,:),1,nLFLperiods+1);
    sumLiftOffDateProbsHighCred = ...
        compute_liftoff_probabilities(LFLpathsHighCred,rVarType,...
        rVarIndex,policyLBval,exanteLFLexitProbsHighCred);
    lossesAlongPathsHighCred = ...
        compute_losses_in_sim_with_exogenous_imperfect_credibility(...
        Model,LFLpathsHighCred.modelVariables,lossFuncHandle);
    [renegeLossesHighCred,continuationLossesHighCred,...
        netLossesHighCred] = compute_renege_and_continuation_losses(...
        lossesAlongPathsHighCred,LFLpathProbsHighCred);
    % (ii) Case with liftoff in period 12
    pLFLkPlus1HighCred = zeros(1,nLFLperiods+1);
    [LFLpathsKplus1HighCred,LFLpathProbsKplus1HighCred,...
        LFLexpectedPathsKplus1HighCred,LFLexitProbsKplus1HighCred] = ...
        execute_policy_sim_with_imperfect_credibility(...
        Model,simBase,pLFLkPlus1HighCred,rLFLkPlus1,policyLBval,...
        policyEqName,policyVarMnem,policyShockMnem,...
        probUpdateFuncHandleHighCred);
    exanteLFLexitProbsKplus1HighCred = ...
        reshape(LFLpathProbsKplus1HighCred(1,1,:),1,(nLFLperiods+1)+1);
    sumLiftOffDateProbsKplus1HighCred = ...
        compute_liftoff_probabilities(LFLpathsKplus1HighCred,rVarType,...
        rVarIndex,policyLBval,exanteLFLexitProbsKplus1HighCred);
    lossesAlongPathsKplus1HighCred = ...
        compute_losses_in_sim_with_exogenous_imperfect_credibility(...
        Model,LFLpathsKplus1HighCred.modelVariables,lossFuncHandle);
    [renegeLossesKplus1HighCred,continuationLossesKplus1HighCred,...
        netLossesKplus1HighCred] = ...
        compute_renege_and_continuation_losses(...
        lossesAlongPathsKplus1HighCred,LFLpathProbsKplus1HighCred);
end

%% SETP 3b: IMPERFECT CREDIBILITY - LOW CREDIBILITY VARIANT
if doSensitivityAnalysis
    % Follows the same recipe as per main results from step 3
    probFuncHandleLowCred = @(temptationToRenege) ...
        compute_reversion_probabilities_using_Weibull_distribution(...
        temptationToRenege,alfa1LowCred,alfa2LowCred);
    probUpdateFuncHandleLowCred = @(simResults) ...
        update_renege_probabilities_using_loss_comparison(simResults,...
        lossFuncHandle,probFuncHandleLowCred);    
    % (i) Case with liftoff in period 11
    pLFLlowCred = zeros(1,nLFLperiods);
    [LFLpathsLowCred,LFLpathProbsLowCred,LFLexpectedPathsLowCred,...
        LFLexitProbsLowCred] = ...
        execute_policy_sim_with_imperfect_credibility(...
        Model,simBase,pLFLlowCred,rLFL,policyLBval,policyEqName,...
        policyVarMnem,policyShockMnem,probUpdateFuncHandleLowCred);
    exanteLFLexitProbsLowCred = ...
        reshape(LFLpathProbsLowCred(1,1,:),1,nLFLperiods+1);
    sumLiftOffDateProbsLowCred = ...
        compute_liftoff_probabilities(LFLpathsLowCred,rVarType,...
        rVarIndex,policyLBval,exanteLFLexitProbsLowCred);
    lossesAlongPathsLowCred = ...
        compute_losses_in_sim_with_exogenous_imperfect_credibility(...
        Model,LFLpathsLowCred.modelVariables,lossFuncHandle);
    [renegeLossesLowCred,continuationLossesLowCred,netLossesLowCred] = ...
        compute_renege_and_continuation_losses(lossesAlongPathsLowCred,...
        LFLpathProbsLowCred);
    % (ii) Case with liftoff in period 12
    pLFLkPlus1LowCred = zeros(1,nLFLperiods+1);
    [LFLpathsKplus1LowCred,LFLpathProbsKplus1LowCred,...
        LFLexpectedPathsKplus1LowCred,LFLexitProbsKplus1LowCred] = ...
        execute_policy_sim_with_imperfect_credibility(...
        Model,simBase,pLFLkPlus1LowCred,rLFLkPlus1,policyLBval,...
        policyEqName,policyVarMnem,policyShockMnem,...
        probUpdateFuncHandleLowCred);
    exanteLFLexitProbsKplus1LowCred = ...
        reshape(LFLpathProbsKplus1LowCred(1,1,:),1,(nLFLperiods+1)+1);
    sumLiftOffDateProbsKplus1LowCred = ...
        compute_liftoff_probabilities(LFLpathsKplus1LowCred,rVarType,...
        rVarIndex,policyLBval,exanteLFLexitProbsKplus1LowCred);
    lossesAlongPathsKplus1LowCred = ...
        compute_losses_in_sim_with_exogenous_imperfect_credibility(...
        Model,LFLpathsKplus1LowCred.modelVariables,lossFuncHandle);
    [renegeLossesKplus1LowCred,continuationLossesKplus1LowCred,...
        netLossesKplus1LowCred] = ...
        compute_renege_and_continuation_losses(...
        lossesAlongPathsKplus1LowCred,LFLpathProbsKplus1LowCred);
end

%% PACK RESULTS & SAVE
results.Model = Model;
results.nLFLperiods = nLFLperiods;
% Key settings
results.settings.shockSize = shockSize;
results.settings.objVarMnems = objVarMnems;
results.settings.objVarWeights = objVarWeights;
% Save baseline results
results.Baseline.Paths = basePaths;
results.Baseline.Losses = baseLosses;
% Save the results under perfect credibility
results.LFLk.FullCred.Paths = fullCredLFLpaths;
results.LFLk.FullCred.PathProbs = fullCredLFLpathProbs;
results.LFLk.FullCred.ExpectedPaths = fullCredLFLexpectedPaths;
results.LFLk.FullCred.ExitProbs = fullCredLFLexitProbs;
results.LFLk.FullCred.Losses = lossesAlongPathsFullCred;
results.LFLkPlus1.FullCred.Paths = fullCredLFLpathsKplus1;
results.LFLkPlus1.FullCred.PathProbs = fullCredLFLpathProbsKplus1;
results.LFLkPlus1.FullCred.ExpectedPaths = fullCredLFLexpectedPathsKplus1;
results.LFLkPlus1.FullCred.ExitProbs = fullCredLFLexitProbsKplus1;
results.LFLkPlus1.FullCred.Losses = lossesAlongPathsFullCredKplus1;
% Save the results under imperfect credibility
results.LFLk.ImpCred.Paths = LFLpaths;
results.LFLk.ImpCred.PathProbs = LFLpathProbs;
results.LFLk.ImpCred.ExpectedPaths = LFLexpectedPaths;
results.LFLk.ImpCred.ExitProbs = LFLexitProbs;
results.LFLk.ImpCred.Losses = lossesAlongPaths;
results.LFLk.ImpCred.RenegeLosses = renegeLosses;
results.LFLk.ImpCred.ContinuationLosses = continuationLosses;
results.LFLk.ImpCred.NetLosses = netLosses;
results.LFLk.ImpCred.ExAnteExitProbs = exanteLFLexitProbs;
results.LFLk.ImpCred.SumLiftOffDateProbs = sumLiftOffDateProbs;
results.LFLkPlus1.ImpCred.Paths = LFLpathsKplus1;
results.LFLkPlus1.ImpCred.PathProbs = LFLpathProbsKplus1;
results.LFLkPlus1.ImpCred.ExpectedPaths = LFLexpectedPathsKplus1;
results.LFLkPlus1.ImpCred.ExitProbs = LFLexitProbsKplus1;
results.LFLkPlus1.ImpCred.Losses = lossesAlongPathsKplus1;
results.LFLkPlus1.ImpCred.RenegeLosses = renegeLossesKplus1;
results.LFLkPlus1.ImpCred.ContinuationLosses = continuationLossesKplus1;
results.LFLkPlus1.ImpCred.NetLosses = netLossesKplus1;
results.LFLkPlus1.ImpCred.ExAnteExitProbs = exanteLFLexitProbsKplus1;
results.LFLkPlus1.ImpCred.SumLiftOffDateProbs = sumLiftOffDateProbsKplus1;
if doSensitivityAnalysis
    results.LFLk.HighCred.Paths = LFLpathsHighCred;
    results.LFLk.HighCred.PathProbs = LFLpathProbsHighCred;
    results.LFLk.HighCred.ExpectedPaths = LFLexpectedPathsHighCred;
    results.LFLk.HighCred.ExitProbs = LFLexitProbsHighCred;
    results.LFLk.HighCred.Losses = lossesAlongPathsHighCred;
    results.LFLk.HighCred.RenegeLosses = renegeLossesHighCred;
    results.LFLk.HighCred.ContinuationLosses = continuationLossesHighCred;
    results.LFLk.HighCred.NetLosses = netLossesHighCred;
    results.LFLk.HighCred.ExAnteExitProbs = exanteLFLexitProbsHighCred;
    results.LFLk.HighCred.SumLiftOffDateProbs = ...
        sumLiftOffDateProbsHighCred;
    results.LFLkPlus1.HighCred.Paths = LFLpathsKplus1HighCred;
    results.LFLkPlus1.HighCred.PathProbs = LFLpathProbsKplus1HighCred;
    results.LFLkPlus1.HighCred.ExpectedPaths = ...
        LFLexpectedPathsKplus1HighCred;
    results.LFLkPlus1.HighCred.ExitProbs = LFLexitProbsKplus1HighCred;
    results.LFLkPlus1.HighCred.Losses = lossesAlongPathsKplus1HighCred;
    results.LFLkPlus1.HighCred.RenegeLosses = ...
        renegeLossesKplus1HighCred;
    results.LFLkPlus1.HighCred.ContinuationLosses = ...
        continuationLossesKplus1HighCred;
    results.LFLkPlus1.HighCred.NetLosses = netLossesKplus1HighCred;
    results.LFLkPlus1.HighCred.ExAnteExitProbs = ...
        exanteLFLexitProbsKplus1HighCred;
    results.LFLkPlus1.HighCred.SumLiftOffDateProbs = ...
        sumLiftOffDateProbsKplus1HighCred;    
    results.LFLk.LowCred.Paths = LFLpathsLowCred;
    results.LFLk.LowCred.PathProbs = LFLpathProbsLowCred;
    results.LFLk.LowCred.ExpectedPaths = LFLexpectedPathsLowCred;
    results.LFLk.LowCred.ExitProbs = LFLexitProbsLowCred;
    results.LFLk.LowCred.Losses = lossesAlongPathsLowCred;
    results.LFLk.LowCred.RenegeLosses = renegeLossesLowCred;
    results.LFLk.LowCred.ContinuationLosses = ...
        continuationLossesLowCred;
    results.LFLk.LowCred.NetLosses = netLossesLowCred;
    results.LFLk.LowCred.ExAnteExitProbs = exanteLFLexitProbsLowCred;
    results.LFLk.LowCred.SumLiftOffDateProbs = ...
        sumLiftOffDateProbsLowCred;    
    results.LFLkPlus1.LowCred.Paths = LFLpathsKplus1LowCred;
    results.LFLkPlus1.LowCred.PathProbs = LFLpathProbsKplus1LowCred;
    results.LFLkPlus1.LowCred.ExpectedPaths = ...
        LFLexpectedPathsKplus1LowCred;
    results.LFLkPlus1.LowCred.ExitProbs = LFLexitProbsKplus1LowCred;
    results.LFLkPlus1.LowCred.Losses = lossesAlongPathsKplus1LowCred;
    results.LFLkPlus1.LowCred.RenegeLosses = ...
        renegeLossesKplus1LowCred;
    results.LFLkPlus1.LowCred.ContinuationLosses = ...
        continuationLossesKplus1LowCred;
    results.LFLkPlus1.LowCred.NetLosses = netLossesKplus1LowCred;
    results.LFLkPlus1.LowCred.ExAnteExitProbs = ...
        exanteLFLexitProbsKplus1LowCred;
    results.LFLkPlus1.LowCred.SumLiftOffDateProbs = ...
        sumLiftOffDateProbsKplus1LowCred;
end
if saveResults
    save(saveFileName,'results');
end