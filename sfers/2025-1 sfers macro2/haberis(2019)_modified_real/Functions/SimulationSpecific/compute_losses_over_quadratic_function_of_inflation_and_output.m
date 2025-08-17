function losses = ...
    compute_losses_over_quadratic_function_of_inflation_and_output(...
    Model,xPaths,piMnem,yMnem,outputWeight,discountFactor,lossScalar)
% This function computes quadratic losses over inflation and output gap.
% More precisely, it computes the discounted sum of squared deviations for
% annual inflation from target and the output gap from zero (with a choice
% of weight to apply to output) across a set of forecasts with a choice of
% horizon up to which to compute the losses.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> xPaths: nx*H*nForecasts matrix of forecasts for the model variables
%   -> piMnem: mnemonic for the inflation variable to use 
%   -> yMnem: mnemonic for the output variable to use
%   -> outputWeight (optional): relative weight to apply to output
%   -> discountFactor (optional): discount factor to apply
%   -> lossScalar (otpional): scaling factor to apply
%
% OUTPUTS:  
%   -> losses: 1*1*nForecasts vector of losses for each forecast in the set
%
% DETAILS:  
%   -> This function computes a quadratic loss over deviations of inflation
%      and output from steady state (eg inflation from target & output gap
%      from 0).
%   -> It does so over a set of forecast paths, which could represent the
%      outcomes based on alternative policies or based on a stochastic
%      simulation conditional on a particular policy.
%   -> The loss will depend on the loss-specific parameters passed into
%      this function, as well as the forecast paths themselves.
%
% NOTES:
%   -> The error handling in this function is designed to protect users
%      from accidentally passing the parameters in the wrong order or 
%      providing something that doesn't make any sense.
%
% This version: 09/10/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 4
    error([mfilename,'cannot proceed because too few inputs passed in'])
end

%% CHECK ON OPTIONAL INPUTS
if nargin < 5
    outputWeight = 1;
end
if nargin < 6
    discountFactor = 1;
end
if nargin < 7
    lossScalar = 1;
end

%% EXTRACT THE FORECASTS FOR THE VARIABLES OF RELEVANCE
xMnems = unpack_model(Model,'xMnems');
yInd = lookup_model_index_numbers(xMnems,yMnem);
piInd = lookup_model_index_numbers(xMnems,piMnem);
yForecasts = xPaths(yInd,:,:);
piForecasts = xPaths(piInd,:,:);

%% COMPUTE FORECAST HORIZON & NUMBER OF FORECAST PATHS
[~,H,nForecasts] = size(xPaths);

%% INITIALISE OUTPUT
losses = NaN*ones(1,1,nForecasts);

%% COMPUTE LOSSES
% This just computes the discounted sum of squared deviations of annual
% inflation from target (which is 0 in model variable space) and the output
% gap from zero with some weight, all divided by an abitrary scalar.
discountVector = discountFactor.^(0:H-1);
for iForecast = 1:nForecasts
    losses(1,1,iForecast) = discountVector*...
        (piForecasts(1,1:H,iForecast).^2'+...
        outputWeight*yForecasts(1,1:H,iForecast).^2')/lossScalar;
end

end