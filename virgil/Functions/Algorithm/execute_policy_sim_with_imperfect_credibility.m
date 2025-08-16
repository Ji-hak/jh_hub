function [ForecastPaths,xPathProbs,ExpectedForecastPaths,renegeProbs] = ...
    execute_policy_sim_with_imperfect_credibility(...
    Model,BaseForecastRun,renegeProbs,pegVals,policyLBval,policyEqName,...
    policyVarMnem,policyShockMnem,probUpdateFuncHandle)
% Computes a set of peg sims with imperfect credibility from base forecast.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> BaseForecastRun: Base forecast run data structure
%   -> renegeProbs: 1*K row vector of renege probabilities (where K is
%      maximum number of periods over which the peg policy might be active)
%   -> pegVals: 1*K row vector of peg values (where K is max number of
%      periods over which the peg might be active)
%   -> policyLBval: scalar value for the policy ELB (in model space)
%   -> policyEqName: name of the policy equation in the model
%   -> policyVarMnem: mnemonic of the policy variable in the model
%   -> policyShockMnem: mnemonic of the policy variable in the model
%   -> probUpdateFuncHandle (mode dependent): function handle that maps
%   simulation outcomes to probabilities
%
% OUTPUTS:  
%   -> ForecastPaths: structure with alternative forecast paths for model
%      variables, observables and raw observables
%   -> xPathProbs: 1*(K+1)*(K+1) matrix of probabilities to attach to each
%      path conditional on no renge up to that point (1,1,:) is the vector
%      of ex-ante probabilities
%   -> ExpectedForecastPaths: structure with expected forecast paths for 
%      model variables, observables and raw observables
%   -> renegeProbs: renege probabilities
%
% DETAILS:  
%   -> See the main function called below for details of the computation
%      and the format of the outputs.
%   -> This function checks and prepares the inputs and then transforms the
%      endogenous model variable outputs into model observables and raw
%      observables.
%
% NOTES:
%   -> The switching between exogenous and endogenous probabilities is a
%      bit inelegant (esp as it doesn't allow the separate options to be
%      passed in around the policy shock and probability iterations).
%
% This version: 05/11/2015
% Author(s): Matt & Rich

%% CHECK OF INPUTS
% This could be upgraded at some point to include more checks if we start
% to have trouble with particular error cases.
if nargin<8 
   error([mfilename,' requires either 8 inputs (for exogenous ',...
       'imperfect credibility sims) or 9 inputs for endogenous ',...
       'imperfect credibility sims']) 
end

%% UNPACK BASE FORECAST RUN & COMPUTE FORECAST HORIZON
% Could do more checks here - validity of forecast run with model input;
% model must be foreward looking etc.
z = BaseForecastRun.Forecast.Shocks.anticipated;
x0 = BaseForecastRun.Constraint.modelVariables;

%% DETERMINE WHETHER THE SIM IS FOR EXOGENOUS OR ENDOGENOUS PROBS
if nargin > 8
    runEndogenousImperfectCredSims = true;
else
    runEndogenousImperfectCredSims = false;
end

%% CALL LOWER LEVEL SIMULATION FUNCTION
if runEndogenousImperfectCredSims
    [xPaths,xPathProbs,ExPaths,renegeProbs] = ...
        compute_policy_sim_with_endogenous_imperfect_credibility(...
        Model,x0,z,renegeProbs,pegVals,policyLBval,policyEqName,...
        policyVarMnem,policyShockMnem,probUpdateFuncHandle);
else
    [xPaths,xPathProbs,ExPaths] = ...
        compute_policy_sim_with_exogenous_imperfect_credibility(...
        Model,x0,z,renegeProbs,pegVals,policyLBval,policyEqName,...
        policyVarMnem,policyShockMnem);
end

%% COMPUTE ALTERNATIVE OBSERVABLE PATHS (IF APPLICABLE)
ForecastPaths.modelVariables = xPaths;
K = size(renegeProbs,2);
H = size(z,2);
[modelHasMeasurementEqs,modelHasDataTransformationEqs] = unpack_model(...
    Model,{'modelHasMeasurementEqs','modelHasDataTransformationEqs'});
if modelHasMeasurementEqs
    [D,G] = unpack_model(Model,{'D','G'});
    nY = size(D,1);
    Ypaths = NaN*ones(nY,H,K+1);
    for i = 1:K+1
        Ypaths(:,:,i) = compute_observables_from_model_variables(...
            xPaths(:,:,i),D,G);
    end
    ForecastPaths.modelObservables = Ypaths;
end
if modelHasDataTransformationEqs
    [RTfunHandle,modelHasTimeVaryingTrends] = unpack_model(...
        Model,{'RTfunHandle','modelHasTimeVaryingTrends'});
    Ytilde0 = BaseForecastRun.Constraint.rawObservables;
    if modelHasTimeVaryingTrends
        etatf = BaseForecastRun.Forecast.timeVaryingTrends;
        etat0 = BaseForecastRun.Constraint.timeVaryingTrends;
    end
    YtildePaths = NaN*ones(nY,H,K+1);
    for i = 1:K+1
        if modelHasTimeVaryingTrends
            YtildePaths(:,:,i) = ...
                transform_observables_from_model_to_raw_space(...
                RTfunHandle,Ypaths(:,:,i),Ytilde0,etatf,etat0);
        else
            YtildePaths(:,:,i) = ...
                transform_observables_from_model_to_raw_space(...
                RTfunHandle,Ypaths(:,:,i),Ytilde0);
        end
    end
    ForecastPaths.rawObservables = YtildePaths;    
end

%% COMPUTE ALTERNATIVE EXPECTED OBSERVABLE PATHS (IF APPLICABLE)
ExpectedForecastPaths.modelVariables = ExPaths;
if modelHasMeasurementEqs
    EYpaths = NaN*ones(nY,H,K+1);
    for i = 1:K
        EYpaths(:,i:H,i) = compute_observables_from_model_variables(...
            ExPaths(:,i:H,i),D,G);
    end
    ExpectedForecastPaths.modelObservables = EYpaths;
end
if modelHasDataTransformationEqs
    EYtildePaths = NaN*ones(nY,H,K+1);
    for i = 1:K
        if i == 1
            Ytilde0 = BaseForecastRun.Constraint.rawObservables;
        else
            Ytilde0 = YtildePaths(:,i-1,i);
        end
        if modelHasTimeVaryingTrends
            if i == 1
                etat0 = BaseForecastRun.Constraint.timeVaryingTrends;
            else
                etat0 = etatf(:,i-1);
            end
            EYtildePaths(:,i:H,i) = ...
                transform_observables_from_model_to_raw_space(...
                RTfunHandle,EYpaths(:,i:H,i),Ytilde0,etatf(:,i:H),etat0);
        else
            EYtildePaths(:,i:H,i) = ...
                transform_observables_from_model_to_raw_space(...
                RTfunHandle,EYpaths(:,i:H,i),Ytilde0);
        end
    end
    ExpectedForecastPaths.rawObservables = EYtildePaths;    
end

end