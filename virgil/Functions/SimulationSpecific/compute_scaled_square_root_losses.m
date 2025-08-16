function squareRootLosses = compute_scaled_square_root_losses(...
    Model,xPaths,piMnem,yMnem,outputWeight,discountFactor)
% This function computes scaled quadratic losses over inflation and output gap.
% More precisely, it computes the discounted sum of squared deviations for
% annual inflation from target and the output gap from zero (with a choice
% of weight to apply to output) across a set of forecasts. Scaling factors
% and a square root operator are applied to map the losses into units that
% are equivalent to a permanent change in the inflation target.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> xPaths: nx*H*nForecasts matrix of forecasts for the model variables
%   -> piMnem: mnemonic for the inflation variable to use 
%   -> yMnem: mnemonic for the output variable to use
%   -> outputWeight: relative weight to apply to output
%   -> discountFactor: discount factor to apply

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
% This version: 05/10/2015
% Author(s): Richard Harrison

%% CHECK INPUTS
if nargin < 6
    error([mfilename,'cannot proceed because too few inputs passed in'])
end


%% COMPUTE LOSSES WITH RELEVANT SCALING FACTOR
% lossScalar = 16*(1-discountFactor); % For (annualised) infln equivalent 
lossScalar = (1-discountFactor)/outputWeight; % For output equivalent

losses = ...
    compute_losses_over_quadratic_function_of_inflation_and_output(...
    Model,xPaths,piMnem,yMnem,outputWeight,discountFactor,lossScalar);

%% COMPUTE SQUARE ROOT OF LOSSES
squareRootLosses = losses.^0.5;


end