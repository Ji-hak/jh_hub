function [renegeLosses,continuationLosses,netLosses] = ...
    compute_renege_and_continuation_losses(lossesAlongPaths,pathProbs)
%compute_renege_and_continuation_losses 
%   This function computes the losses experienced when reneging from the
%   announced policy plan or continuing with it, by weighting the losses
%   attached to each path by the probabilities of those paths.


K = size(lossesAlongPaths,2)-1;
%% COMPUTE RENEGE LOSSES & CONTINUATION LOSSES
renegeLosses = zeros(1,K);
continuationLosses = zeros(1,K);
for i = 1:K
    renegeLosses(i) = lossesAlongPaths(1,i,i);
    for s = i+1:K+1
        continuationLosses(i) = continuationLosses(i)+...
            (pathProbs(1,i,s)/(1-pathProbs(1,i,i)))*...
            lossesAlongPaths(1,i,s);
    end
end

%% COMPUTE NET LOSSES
netLosses  = continuationLosses - renegeLosses;
netLosses = max(zeros(1,K),netLosses);


end

