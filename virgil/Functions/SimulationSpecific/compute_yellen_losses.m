function losses = compute_yellen_losses(Model,xPaths,objVarMnems,...
    objVarWeights,discountFactor)
% This function computes losses for Yellen's loss function.

%% SET UP WEIGHTS MATRIX ON OBJECTIVE VARIABLES IN LOSS FUNCTION
W = diag(objVarWeights);

%% UNPACK XMNEMS
xMnems = unpack_model(Model,'xMnems');

%% EXTRACT THE FORECASTS FOR THE VARIABLES OF RELEVANCE
objVarInd = lookup_model_index_numbers(xMnems,objVarMnems);
objVarForecasts = xPaths(objVarInd,:,:);

%% COMPUTE FORECAST HORIZON & NUMBER OF FORECAST PATHS
[~,H,nForecasts] = size(xPaths);

%% INITIALISE OUTPUT
losses = NaN*ones(1,1,nForecasts);

%% COMPUTE LOSSES
for iForecast = 1:nForecasts
     L = 0;
    for t = 1:H
        x_t = objVarForecasts(:,t,iForecast);
        L_t = x_t'*W*x_t;
        L = L + discountFactor^(t-1)*L_t;
    end
    losses(1,1,iForecast) = L;
end

end