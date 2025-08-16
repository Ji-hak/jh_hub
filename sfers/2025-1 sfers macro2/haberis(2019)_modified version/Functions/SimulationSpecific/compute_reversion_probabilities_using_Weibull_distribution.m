function p = compute_reversion_probabilities_using_Weibull_distribution(...
    reversionGap,alfa1,alfa2)
% This function computes reversion probabilities as a function of the 
% reversion gap. It does so using the Weibull distribution with parameters 
% alpha1 and alpha2.
%
% INPUTS:   
%   -> reversionGap: if the reversion gap is positive, there is a positive
%      probability that the policymaker reverts to the rule. See Details
%      below for examples.
%   -> alfa1, alfa1: cdf parameters
%
% OUTPUTS:  
%   -> p: 1*K vector of probabilities
%
% DETAILS:  
%   -> This function computes the probabilities of early exit of an
%      announced LFL policy path as an exponential function of the
%      reversion gap. 
%      Example 1, imperfect credibility: the reversion gap is the
%      temptation to renege. There is a temptation to renege if the
%      continuation loss is greater than the renege loss.
%      Example 2, threshold-based forward guidance: the reversion gap is
%      the gap between the threshold variable and the threshold level. To
%      the extent that the TBFG is aiming to introduce stimulus, care needs
%      to be taken in how this gap is specified depending on the threshold
%      variable.
%         > For unemployment, this is ubar - u{t}. This means if u{t} is 
%           below ubar (so that ubar-u{t}>0), there is a positive 
%           probability of reverting to the rule.
%         > For the output, it would be y{t}-ybar. That is, if output is
%           above the threshold level, there is a positive probability of
%           reverting to the rule.
%
% NOTES:
%   -> The error handling in this function is minimal.
%
% This version: 20/01/2017
% Author(s): Richard Harrison and Alex Haberis

%% CHECK INPUTS
if nargin < 3
    error([mfilename,' cannot proceed because too few inputs passed in']);
end

%% INITIALISE OUTPUT
K = size(reversionGap,2);
p = NaN*ones(1,K);

%% COMPUTE UPDATED PROBABILITIES
for i = 1:K
    if reversionGap(i) > 0
        p(i) = evaluate_weibull_cdf(reversionGap(i),alfa1,alfa2);
    else
        p(i) = 0;
    end
end

end