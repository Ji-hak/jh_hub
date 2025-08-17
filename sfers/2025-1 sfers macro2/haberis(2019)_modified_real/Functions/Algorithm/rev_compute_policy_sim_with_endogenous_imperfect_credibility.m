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
%% ITERATE OVER PROBABILITIES
%% ITERATE OVER PROBABILITIES
iter = 0;
pMaxChange = inf;

% Populate iteration results structure with elements that do not change
lastIterResults = struct;
lastIterResults.Model = Model;
lastIterResults.rPeg  = rPeg;
lastIterResults.b     = b;

% --- 안정화 파라미터 (동적 스텝/감쇠) ---------------------------
% 기본값: 옵션 구조체에 없으면 디폴트 사용
if ~isfield(ProbabilityIterOpts,'bifurcationFrac') || ...
        isempty(ProbabilityIterOpts.bifurcationFrac) || ...
        ProbabilityIterOpts.bifurcationFrac <= 0 || ...
        ProbabilityIterOpts.bifurcationFrac > 1
    probRelax = 0.3;           % 0.2~0.5 권장
else
    probRelax = ProbabilityIterOpts.bifurcationFrac;
end
if ~isfield(ProbabilityIterOpts,'stepCap') || isempty(ProbabilityIterOpts.stepCap)
    stepCap = 0.15;            % 각 원소별 1회 최대 이동폭
else
    stepCap = ProbabilityIterOpts.stepCap;
end
if ~isfield(ProbabilityIterOpts,'stepCapMax') || isempty(ProbabilityIterOpts.stepCapMax)
    stepCapMax = 0.6;          % 스텝캡 상한
else
    stepCapMax = ProbabilityIterOpts.stepCapMax;
end
if ~isfield(ProbabilityIterOpts,'relaxMin') || isempty(ProbabilityIterOpts.relaxMin)
    relaxMin = 0.05;
else
    relaxMin = ProbabilityIterOpts.relaxMin;
end
if ~isfield(ProbabilityIterOpts,'relaxMax') || isempty(ProbabilityIterOpts.relaxMax)
    relaxMax = 0.6;
else
    relaxMax = ProbabilityIterOpts.relaxMax;
end
if ~isfield(ProbabilityIterOpts,'plateauWindow') || isempty(ProbabilityIterOpts.plateauWindow)
    plateauWindow = 5;         % plateau 감지 횟수
else
    plateauWindow = ProbabilityIterOpts.plateauWindow;
end
if ~isfield(ProbabilityIterOpts,'plateauTol') || isempty(ProbabilityIterOpts.plateauTol)
    plateauTol = 0.02;         % 포화스텝 대비 2% 이내면 plateau로 간주
else
    plateauTol = ProbabilityIterOpts.plateauTol;
end

pMin = 1e-6; pMax = 1-1e-6;

prevMaxChange = inf;
plateauCount  = 0;

while (pMaxChange > pChangeTol) && iter < maxIter
    pOld = p;
    iter = iter + 1;
    fprintf('PROBABILITIES ITERATION %d\n', iter);

    [xPaths,xPathProbs,eStarPaths] = ...
        compute_policy_sim_with_exogenous_imperfect_credibility( ...
            Model,x0,z,p,rPeg,b,policyEqName,policyVarMnem,policyShockMnem, ...
            PolicyShockIterOpts);

    % 업데이트용 입력
    lastIterResults.xPaths     = xPaths;
    lastIterResults.xPathProbs = xPathProbs;
    lastIterResults.eStarPaths = eStarPaths;
    lastIterResults.p          = p;

    % 1) 제안 확률
    pProp = updateProbFuncHandle(lastIterResults);

    % 2) (0,1) 범위로 클리핑
    pProp = max(pMin, min(pMax, pProp));

    % 3) 스텝 캡(각 원소별 변화 폭 제한)
    dp = pProp - pOld;
    dp = max(-stepCap, min(stepCap, dp));

    % 4) 감쇠 적용
    pNew = pOld + probRelax * dp;

    % 5) 최종 클리핑
    pNew = max(pMin, min(pMax, pNew));

    % 변화량 계산
    pMaxChange = max(abs(pNew - pOld));
    fprintf('MAX CHANGE IN PROBS: %g\n', pMaxChange);

    % ----- 동적 적응 로직 --------------------------------------
    % (a) plateau 감지: 변화량이 probRelax*stepCap 근처에서 고정되면 (=포화)
    satStep = probRelax * stepCap;
    if abs(pMaxChange - satStep) <= plateauTol * max(1e-12, satStep)
        plateauCount = plateauCount + 1;
    else
        plateauCount = 0;
    end
    % plateau가 연속 발생하면 과감히 전진하기 위해 스텝/감쇠 키움
    if plateauCount >= plateauWindow
        stepCap  = min(stepCap * 1.5, stepCapMax);
        probRelax = min(relaxMax, probRelax * 1.25);
        plateauCount = 0;
    end

    % (b) 악화/진동 감지: 변화량이 커지면 보수적으로 축소
    if pMaxChange > prevMaxChange * 1.25
        probRelax = max(relaxMin, probRelax * 0.5);
        stepCap   = max(0.05, stepCap * 0.7);
        % 축소 적용 후 재계산
        dp   = pProp - pOld;
        dp   = max(-stepCap, min(stepCap, dp));
        pNew = pOld + probRelax * dp;
        pNew = max(pMin, min(pMax, pNew));
        pMaxChange = max(abs(pNew - pOld));
    end
    prevMaxChange = pMaxChange;
    % ----------------------------------------------------------

    % 6) 업데이트
    p = pNew;
end

if iter >= maxIter && pMaxChange > pChangeTol
    warning('Probability iteration hit maxIter (%d) with maxChange=%g', maxIter, pMaxChange);
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
DefaultProbabilityIterOpts.pChangeTol = 1e-3;
DefaultProbabilityIterOpts.maxIter = 1000;
DefaultProbabilityIterOpts.bifurcationFrac = 0.3;

end