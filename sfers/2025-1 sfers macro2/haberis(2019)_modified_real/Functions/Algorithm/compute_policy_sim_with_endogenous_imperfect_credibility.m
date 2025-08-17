function [xPaths,xPathProbs,eStarPaths,p] = ...
    compute_policy_sim_with_endogenous_imperfect_credibility(...
    Model,x0,z,p,rPeg,b,policyEqName,policyVarMnem,policyShockMnem,...
    updateProbFuncHandle,ProbabilityIterOpts,PolicyShockIterOpts)
% Computes a deterministic peg sim with endogenous imperfect credibility.
% This function computes a complete set of simulation paths, probabilities
% associated with each of those paths and a set of consistent expectations
% given an interest rate peg announcement and an endogenous (and time-
% varying) probabiity that the policy maker will abandon the peg & revert 
% to the policy rule.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> x0: nx*1 vector of initial conditions for the endogneous vars
%   -> z: nz*H matrix of anticipated shocks
%   -> p: 1*K row vector initialisation of renege probabilities (where K is
%      maximum number of periods over which the LFL policy might be active)
%   -> rPeg: 1*K row vector of instrument peg values (in model space)
%   -> b: scalar value for the policy ELB (in model space)
%   -> policyEqName: name of the policy equation in the model
%   -> policyVarMnem: mnemonic of the policy variable in the model
%   -> policyShockMnem: mnemonic of the policy shock in the model
%   -> updateProbFuncHandle: function handle to update per-period renege
%   probabilities (this if the function F(\widetilde X,P,B) in the paper
%   -> ProbabilityIterOpts (optional): structure of options for iteration
%      over probabilities:
%       - pChangeTol: amount be which changes in probabilities are treated
%         as identical
%       - maxIter: maximum number of iterations allowed
%       - bifurcationFrac: numeric scalar controlling speed of update
%         between iterations (eg value of 0 is no update, value of 1 is
%         "complete" update)
%   -> PolicyShockIterOpts (optional): structure of options for iteration
%      over policy shocks necessary to impose the ELB in event that it 
%      binds on renege:
%       - bViolationTol: amount be which ELB is allowed to be violated (to
%         allow for some numerical imprecision in quadprog routine)
%       - maxIter: maximum number of iterations allowed
%       - bifurcationFrac: numeric scalar controlling speed of update
%         between iterations (eg value of 0 is no update, value of 1 is
%         "complete" update)
%       - epsrChangeTol: exit criteria tolerance for change in shock values
%         between iterations (i.e. if max abs change is smaller than this
%         value the search will be treated as successful)
%
% OUTPUTS:  
%   -> xPaths: nx*H*(K+1) matrix of simulation paths
%   -> xPathProbs: 1*(K+1) vector of probabilities to attach to each path
%   -> eStarPaths: nx*H*K matrix of expectations conditional on policy 
%      maker not having reneged
%   -> p: 1*K row vector of endogenous renege probabilities (where K is
%      maximum number of periods over which the LFL policy might be active)
%
% DETAILS:  
%   -> In the output to this function, the 2nd dimension references time
%      (i.e. a standard forecast horizon) and the 3rd dimension references
%      states of the world. For example, the 2nd slice of this 3rd
%      dimension, xPaths(:,:,2), would contain the paths for all the
%      endogenous variables over the full simulation horizon conditional on
%      the policy-maker reneging from the LFL policy in period 2. The
%      probability associated with being on that path is given by
%      xPathsProbs(1,2) and the policy shocks to impose the ELB after 
%      renege has happened are given by epsrPaths(1,:,2) (which will be 
%      zero if the ELB is non-binding after renege).
%   -> The expectations paths are computed conditional on the policy-maker
%      not already having reneged (on the grounds that the expectations 
%      paths on reneging are not interesting because they are just given by
%      the rational expectations solution. The dimensions can be 
%      interpreted in a similar way to the simulation paths in that the 2nd
%      dimension refers to time (in the standard MAPS way) and the third 
%      dimension refers to the period in which the expectation is taken 
%      conditional on the LFL policy still being in play. Using the 
%      notation of the paper, exPaths(:,:,2) would be agents' expectations
%      of the endogenous variables across the forecast horizon given that 
%      the policy maker did not renge in period 1, 
%      \mathbb{E}^{\ast}_{2}x_{t} for t = 2:H. Note that exPaths(:,1,2)
%      would contain NaNs as it is simply given by the non-renege outturn
%      in period 1 and so is uninteresting/undefined.
%   -> The function does one thing: iterates over the renege probabilities
%      using compute_LFL_sim_with_exogenous_imperfect_credibility and 
%      compute_losses_in_LFL_sim_with_exogenous_imperfect_credibility.
%      Given a guess for the vector of probabilities, it computes the
%      simulation paths and the losses associated with those paths. It then
%      updates the probabilities and repeats.
%   -> The bulk of the work is done in the function which computes the
%      simulation given a set of exogenous probabilities.
%
% NOTES:
%   -> Latest version incorporates very general function handle control for
%   updating the probabilities.
%
% This version: 05/11/2015
% Author(s): Matt & Rich

%% CHECK INPUTS
if nargin < 10
    error([mfilename,' cannot proceed because too few inputs passed in']);
end

%% HANDLE OPTIONAL POLICY ITERATION OPTIONS INPUT
if nargin < 12
    PolicyShockIterOpts = struct;
end

%% HANDLE OPTIONAL PROBABILITY ITERATION OPTIONS INPUT
DefaultProbabilityIterOpts = get_default_probability_iteration_options();
%if nargin < 11
if nargin < 12
    ProbabilityIterOpts = DefaultProbabilityIterOpts;
else
    ProbabilityIterOpts = overlay_default_structure(...
        ProbabilityIterOpts,DefaultProbabilityIterOpts);
end
maxIter = ProbabilityIterOpts.maxIter;
bifurcationFrac = ProbabilityIterOpts.bifurcationFrac;
pChangeTol = ProbabilityIterOpts.pChangeTol;

%% ITERATE OVER PROBABILITIES
iter = 0;
pMaxChange = inf;

% Populate iteration results structure with elements that do not change
lastIterResults = struct;
lastIterResults.Model = Model;
lastIterResults.rPeg = rPeg;
lastIterResults.b = b;

while (pMaxChange>pChangeTol) && iter<maxIter
    pOld = p;
    iter = iter+1;
    fprintf(['PROBABILITIES ITERATION ',num2str(iter),'\n']);

    [xPaths,xPathProbs,eStarPaths] = ...
        compute_policy_sim_with_exogenous_imperfect_credibility(...
        Model,x0,z,p,rPeg,b,policyEqName,policyVarMnem,policyShockMnem,...
        PolicyShockIterOpts);

    % Updating of probabilities is taken care of by a function that is
    % passed all of the arguments from the latest iteration, packed into a
    % structure.
    lastIterResults.xPaths = xPaths;
    lastIterResults.xPathProbs = xPathProbs;
    lastIterResults.eStarPaths = eStarPaths;
    lastIterResults.p = p;
    pFullUpdate = updateProbFuncHandle(lastIterResults);
        
    p = bifurcationFrac*pFullUpdate+(1-bifurcationFrac)*pOld;
    pMaxChange = max(abs(pFullUpdate-pOld));
    fprintf(['MAX CHANGE IN PROBS: ',num2str(pMaxChange),'\n']);
end

end

%% FUNCTION TO DEFINE DEFAULT OPTIONS USED IN ITERATION OVER PROBABILITIES
function DefaultProbabilityIterOpts = ...
    get_default_probability_iteration_options()
% These defaults seem OK and should not need to be changed.
% They can be overridden by the user on input to this function.
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> DefaultPolicyShockIterOpts: structure of default options

%% DEFINE DEFAULT OPTIONS
DefaultProbabilityIterOpts.pChangeTol = 0.0001;
DefaultProbabilityIterOpts.maxIter = 1e7;
DefaultProbabilityIterOpts.bifurcationFrac = 1/6;

end