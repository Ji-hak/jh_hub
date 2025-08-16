function updatedProbabilities = ...
    update_renege_probabilities_using_loss_comparison(simResults,...
    lossFuncHandle,probFuncHandle)

% This function updated renege probabilities using a function mapping the
% net benefits of reneging to the probabilities of reneging.

% INPUTS:
% -> simResults: structure containing all of the results from an imperfect
% credibility policy simulation
% -> lossFuncHandle: handle of function used to compute losses incurred by
% alternative policy choices (renege dates)
% -> probFuncHandle: handle of function that maps net benefits of reneging
% to the probability of reneging

% OUTPUTS:
% -> updatedProbabilities: probabilities implied by the net benefits of
% reneging.

% Author: Rich
% This version: 05/11/2015

%% COMPUTE LOSSES ALONG PATHS CORRESPONDING TO RENEGE AT DATES 1,...,K+1
lossesAlongPaths = ...
    compute_losses_in_sim_with_exogenous_imperfect_credibility(...
    simResults.Model,simResults.xPaths,lossFuncHandle);

%% COMPUTE RENEGE LOSSES & CONTINUATION LOSSES
pathProbs = simResults.xPathProbs;
[renegeLosses,continuationLosses] = ...
    compute_renege_and_continuation_losses(lossesAlongPaths,pathProbs);

%% MAP NET BENEFITS OF RENEGE INTO PROBABILITIES OF RENEGE
temptationToRenege = continuationLosses - renegeLosses;
updatedProbabilities = probFuncHandle(temptationToRenege);