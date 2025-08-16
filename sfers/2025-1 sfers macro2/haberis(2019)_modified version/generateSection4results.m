%% THIS SCRIPT PRODUCES RESULTS FOR SECTION 4 OF UNCERTAIN POLICY PROMISES
% The following threshold-based forward guidance experiments using the FRB-
% US model are computed:
% -- 5.75% unemployment threshold with a trigger-like liftoff function
% -- 5.75% unemployment threshold with a more-open-ended liftoff function
% -- 6% unemployment threshold with a trigger-like liftoff function
% -- 6% unemployment threshold with a more-open-ended liftoff function
% NB. FRBS is a large model and so this script will take a while to run --
% there are some default initialisations for the reversion probabilities
% taken from the results in the paper which ensures that those variants run
% quickly.

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

%% CHOOSE WHETHER TO SAVE RESULTS OR NOT
saveResults = true;

%% DEFINE ALTERNATIVE TBFG EXPERIMENTS
% The combinations will be computed below -- i.e. TBFG experiments with the
% two alternative unemployment thresholds both computed using the
% alternative probability mapping function calibrations.
uThresholdPct1 = 6.0;
uThresholdPct2 = 5.75;
alfa1empirical = -0.25/log(0.50); 
alfa1trigger = -0.01/log(0.001); 
alfa2 = 1; 
% Naming convention for the saving of results below
probFuncMapEmpiricalName = 'openEnded';
probFuncMapTriggerName = 'trigger';
uThreshold1Name = strrep(num2str(uThresholdPct1),'.','point');
uThreshold2Name = strrep(num2str(uThresholdPct2),'.','point');

%% CHOOSE GENERAL EXPERIMENT OPTIONS
H = 32;                         % Forecast horizon
elb = 0.125;                    % ELB for FFR
maxFGhorizon = 24;
maxFGhorizon = min(maxFGhorizon,H);
% Options that depend on guidance announcement date
LFLhorizon = 8;
baseDataBook = 'Data\FRBUSdataTBFG.xlsx';
baseDataSheet = 'Projection';
FGsheet = 'FGprobabilities';
% Options for building the base
maxTriesBase = 2;
baseDeviationTolerance = 1e-04;
dampingFactors = [0;0;0;0;0];

%% MODEL & DATA OPTIONS
modelName = 'FRBUS.maps';
% Info about ploicy rule, instrument and shock
policyEqName = 'Equation 2';
policyVarMnem = 'rffe';
policyShockMnem = 'etar';
rssMnem = 'FFRbar'; 
% Data source
dataBookName = 'Data\FRBUSdataTBFG.xlsx';
dataSheetName = 'DataForExport';
% Variables to match in building baseline forecast
baseVarMnems = {'annualPCEinf';'core';'annGDPgrowth';'fedFunds';'urate'};

%% BUILD MODEL 
Model = create_model(modelName);
% Specify a (near) full invert to deliver the base
residualsToExclude = {policyShockMnem;'poilr_'};
zMnems = unpack_model(Model,'zMnems');
baseShockMnems = zMnems;
resExcludeInds = lookup_model_index_numbers(zMnems,residualsToExclude);
baseShockMnems(resExcludeInds) = [];

%% LOAD DATA AND PRODUCE PLAIN VANILLA PROJECTION
Options.isStraightEdge = false;
[YtildeAll,runDataDates] = load_model_time_series_data_from_excel(...
    Model,{'YtildeMnems'},dataBookName,dataSheetName,...
    Options);
% COMPUTE DIMENSIONS OF DATA
isNaNYtilde = all(isnan(YtildeAll));
T = sum(~isNaNYtilde);
Ytilde = YtildeAll(:,1:T);
% CREATE FORECAST RUN DATASET CONTAINING PLAIN VANILLA PROJECTION OF DATA
SimBase = create_simulation_journey_base(Model,H,T-1);
DataToFilter.rawObservables = Ytilde;
plainVanillaForecast = filter_data_and_project_using_LSS_model(...
    Model,SimBase,DataToFilter);

%% BUILD BASELINE FORECAST BY IMPOSING JUDGEMENT ON KEY VARIABLES
% Load the forecasts to be imposed in the baseline
Options.isStraightEdge = false;
Options.isComplete = false;
[baseData,~,baseForecastDates] = ...
    load_specified_time_series_data_from_excel(...
    baseVarMnems,baseDataBook,baseDataSheet,Options);
% Build the judgements on the endogenous variables
baseVarInfo = lookup_LSS_mnemonic_type_and_position(Model,baseVarMnems);
nBaseVars = size(baseVarMnems,1);
baseJudgements = struct;
ROcount = 1;
MOcount = 1;
MVcount = 1;
baseTargets = cell(nBaseVars,2);
for iBaseVar = 1:nBaseVars
    iBaseVarMnem = baseVarMnems{iBaseVar};
    iBaseVarType = baseVarInfo{iBaseVar,2};
    iBaseData = baseData(iBaseVar,:);
    iBaseDataLength = size(iBaseData,2);
    % Build judgement structure
    switch iBaseVarType
        case 'rawObservables'
            baseJudgements.AnticipatedFixes.rawObservables(ROcount,:) = ...
                {iBaseVarMnem, iBaseData};
            ROcount = ROcount + 1;
        case 'modelObservables'
            baseJudgements.AnticipatedFixes.modelObservables(MOcount,:) = ...
                {iBaseVarMnem, iBaseData};
            MOcount = MOcount+1;
        case 'modelVariables'
            baseJudgements.AnticipatedFixes.modelVariables(MVcount,:) = ...
                {iBaseVarMnem, iBaseData};
            MVcount = MVcount+1;
    end
    baseTargets(iBaseVar,:) = {iBaseVarMnem,iBaseData};
end
% Incorporate the information about which shocks to use for inversion
nBaseShocks = size(baseShockMnems,1);
for iBaseShock = 1:nBaseShocks
    baseJudgements.AnticipatedFixes.shockUsages(iBaseShock,:) = ...
        {baseShockMnems{iBaseShock}, ones(1,H)};
end
% BUILD THE BASE
baselineForecast = impose_judgement_using_LSS_model(Model,...
    plainVanillaForecast,baseJudgements);

%% SET UP INPUTS FOR EXPERIMENTS
% Compute interest rate measurement equation intercept
[thetaMnems,thetaValues] = unpack_model(Model,{'thetaMnems';'theta'});
intRateAdjustment = ...
    thetaValues(lookup_model_index_numbers(thetaMnems,rssMnem));
elbModelUnits = elb - intRateAdjustment;
FGelbPath = elbModelUnits*ones(1,maxFGhorizon);

%% ADJUST THE UNEMPLOYMENT THRESHOLDS TO BE IN MODEL UNITS
theta = unpack_model(Model,'theta');
ustarIndex = ...
    unpack_model_metadata_and_lookup_index_numbers(Model,'thetaMnems',...
    'ustar');
uThesholdAdj = theta(ustarIndex);
uThreshold1 = uThresholdPct1 - uThesholdAdj;
uThreshold2 = uThresholdPct2 - uThesholdAdj;
% Compute index number of unemployment rate in model
uVarIndex = unpack_model_metadata_and_lookup_index_numbers(...
    Model,'xMnems','lur');

%% DEFINE ALL EXPERIMENTS TO BE COMPUTED
uThresholdSet = [uThreshold1;   uThreshold1; uThreshold2;   uThreshold2];
alfa1Set =      [alfa1empirical;alfa1trigger;alfa1empirical;alfa1trigger];
experimentNames = {
    [uThreshold1Name,probFuncMapEmpiricalName];
    [uThreshold1Name,probFuncMapTriggerName];
    [uThreshold2Name,probFuncMapEmpiricalName];
    [uThreshold2Name,probFuncMapTriggerName]
    };
% Unemployment thresholds in pct space
uThresholdPctSet = [...
    uThresholdPct1;uThresholdPct1;uThresholdPct2;uThresholdPct2];

%% SPECIFY DEFAULT REVERSION PROBABILITY INITIALISATIONS FOR THE 4 CASES
% These are hardwired in from earlier computation 
FGprobsInitSet = {
    [0	0	0	0	0	0	0	0	0	0	0.0992546434233766	...
    0.284013049393979	0.422092976971615	0.476983721979927	...
    0.538540900337388	0.606543686168726	0.675865328960447	...
    0.740772855798428	0.797560693762419	0.844726258198955	...
    0.874517227353670	0.899754039676317	0.920479241934032	...
    0.937141263737413];
    [0	0	0	0	0	0	0	0	0	0	0	0.0001	...
    1	1	1   1	1	1	1	1	1	1	1	1];
    [0  0   0   0   0   0   0   0   0   0   0.0459	...
    0.2261  0.3555  0.3938  0.4425  0.5050  0.5765  ...
    0.6493  0.7174 0.7770 0.8150 0.8486 0.8773 0.9011];
    [0	0	0	0	0	0	0	0	0	0	0	0	...
    0.585562596357987	0.999904384448783	0.999991457936047   ...
    0.999999629042455	0.999999993997651	0.999999999866402	...
    0.999999999992854	0.999999999999145	0.999999999999806	...
    0.999999999999963	0.999999999999996 	0.999999999999996];
    };
    
%% COMPUTE EXPERIMENTS & SAVE RESULT
nEXperiments = size(uThresholdSet,1);
for iExperiment = 1:nEXperiments
    % Set function handles
    thresholdFuncHandle = ...
        @(simResults) compute_unemp_threshold_distance(...
        simResults,uVarIndex,uThresholdSet(iExperiment));
    probFuncHandle = @(thresholdGaps) ...
        compute_reversion_probabilities_using_Weibull_distribution(...
        thresholdGaps,alfa1Set(iExperiment),alfa2);
    probUpdateFuncHandle = @(xPaths) ...
        update_reversion_probabilities_using_threshold_gap(xPaths,...
        thresholdFuncHandle,probFuncHandle);
    % Initialise reversion probabilities
    FGprobsInit = FGprobsInitSet{iExperiment};
    % Run the simulation
    [TBFGpaths,TBFGpathProbs,TBFGexpectedPaths,TBFGrenegeProbs] = ...
        execute_policy_sim_with_imperfect_credibility(...
        Model,baselineForecast,FGprobsInit,FGelbPath,elbModelUnits,...
        policyEqName,policyVarMnem,policyShockMnem,probUpdateFuncHandle);
    % Save the results
    if saveResults
        Output.Model = Model;
        Output.baselineForecast = baselineForecast;
        Output.uThresholdPct = uThresholdPctSet(iExperiment);
        Output.alfa1 = alfa1Set(iExperiment);
        Output.alfa2 = alfa2;
        Output.TBFGpaths = TBFGpaths;
        Output.TBFGpathProbs = TBFGpathProbs;
        Output.TBFGexpectedPaths = TBFGexpectedPaths;
        Output.TBFGrenegeProbs = TBFGrenegeProbs;
        resultsFileName = [thisDir,'Results\TBFGusingFRBUSuThreshold',...
            experimentNames{iExperiment},'.mat'];
        save_content_of_structure_to_mat_file(resultsFileName,Output);
    end
end