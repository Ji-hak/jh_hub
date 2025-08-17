function uThresholdDistance = ...
    compute_unemp_threshold_distance(simResults,varIndex,uBar)

%% EXTRACT THE FORECASTS FOR THE VARIABLES OF RELEVANCE
xPaths = simResults.xPaths;
uPaths = xPaths(varIndex,:,:);

%% COMPUTE FORECAST HORIZON & NUMBER OF FORECAST PATHS
nPaths = size(xPaths,3);

%% POPULATE VECTOR OF THRESHOLD DISTANCES
% Note that the convention in the code is that threshold gaps are positive,
% so we reverse the sign of the unemployment gap accordingly.
uThresholdDistance = nan(1,nPaths);
for iPath = 1:nPaths
    iThresholdDistance = uPaths(1,iPath,iPath) - uBar;
    iThresholdDistance = min(0,iThresholdDistance);
    uThresholdDistance(iPath) = -1*iThresholdDistance;
end

end