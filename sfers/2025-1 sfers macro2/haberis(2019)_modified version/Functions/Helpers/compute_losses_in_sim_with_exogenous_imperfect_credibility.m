function lossesAlongPaths = ...
    compute_losses_in_sim_with_exogenous_imperfect_credibility(...
    Model,xPaths,lossFuncHandle)
% Computes recursive losses in a deterministic LFL sim with imperfect cred.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> xPaths: nx*H*(K+1) matrix of simulation paths
%   -> lossFuncHandle: function handle for the loss function to evaluate
%
% OUTPUTS:  
%   -> lossesAlongPaths: 1*(K+1)*(K+1) set of loss vector.
%
% DETAILS:  
%   -> In the output to this function, the 2nd dimension references time
%      (i.e. a standard forecast horizon) and the 3rd dimension references
%      states of the world. For example, the 2nd slice of this 3rd
%      dimension, lossesAlongPaths(:,2,2:K+1), would contain losses the 
%      policy-maker could achieve in period 2 conditional on all available
%      choices (i.e. conditional on reneging in any of the K available
%      periods).
%
% NOTES:
%   -> Error handling is minimal. Please take care.
%
% DEV NOTES:
%   -> This could be generalised to take in a structure that includes
%      model/raw observables rather than just the model variables. The loss
%      function would then have to deal with that (using the lookup
%      variable type and index number helper).
%
% This version: 09/10/2013
% Author(s): Matt

%% CHECK INPUTS
% This could be upgraded at some point to include more checks if we start
% to have trouble with particular error cases.
if nargin < 3
   error([mfilename,' cannot proceed because too few inputs passed in']);
end

%% COMPUTE FORECAST HORIZON & NUMBER OF STATES FOR POLICY
[~,H,nPaths] = size(xPaths);
K = nPaths-1;

%% INITIALISE LOSSES OUTPUT
lossesAlongPaths = NaN*ones(1,K+1,K+1);

%% COMPUTE LOSSES RECURSIVELY
for t = 1:K+1
    lossesAlongPaths(1,t,t:K+1) = lossFuncHandle(...
        Model,xPaths(:,t:H,t:K+1));
end

end