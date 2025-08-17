function updatedProbabilities = ...
    update_reversion_probabilities_using_threshold_gap(xPaths,...
    thresholdFuncHandle,probFuncHandle)
%   This function updates probabilities using a function mapping the
%   results of a simulation with uncertainty over policy behaviour into the
%   probabilities that the policymaker reverts to usual behaviour.

%% COMPUTE THRESHOLD GAPS
thresholdGaps = thresholdFuncHandle(xPaths);
thresholdGaps = thresholdGaps(1,1:end-1);

%% UPDATE PROBABILITIES
updatedProbabilities = probFuncHandle(thresholdGaps);

end