function [xPaths,xPathProbs,eStarPaths,epsrPaths] = ...
    compute_policy_sim_with_exogenous_imperfect_credibility(...
    Model,x0,z,p,rbar,b,policyEqName,policyVarMnem,policyShockMnem,...
    PolicyShockIterOpts)
% Computes a deterministic peg sim with exogenous imperfect credibility.
% This function computes a complete set of simulation paths, probabilities
% associated with each of those paths and a set of consistent expectations
% given a peg policy (announced path for interest rates) and an exogenous 
% (but time-varying) probabiity that the policy maker will abandon the peg 
% & revert to the policy rule..
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> x0: nx*1 vector of initial conditions for the endogneous vars
%   -> z: nz*H matrix of anticipated shocks
%   -> p: 1*K row vector of renege probabilities (where K is
%      maximum number of periods over which the LFL policy might be active)
%   -> rbar: 1*K row vector of instrument peg values (in model space)
%   -> b: Scalar value for the ELB/ZLB
%   -> policyEqName: name of the policy equation in the model
%   -> policyVarMnem: mnemonic of the policy variable in the model
%   -> policyShockMnem: mnemonic of the policy shock in the model
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
%   -> xPathProbs: 1*(K+1)*(K+1) matrix of probabilities to attach to each
%      path conditional on no renge up to that point (1,1,:) is the vector
%      of ex-ante probabilities
%   -> eStarPaths: nx*H*K matrix of expectations conditional on policy 
%      maker not having reneged
%   -> epsrPaths: 1*H*(K+1) matrix of values for the policy shock necessary
%      to impose an occasionally binding ELB on renege
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
%      the rational expectations solution). The dimensions can be 
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
%   -> The implementation of this function tries to follow the notation in
%      the paper as closely as possible.
%   -> The function does three things: a) computes the simulation paths 
%      taking into account that reversion to the policy rule may violate 
%      the ELB either on exit or afterwards in which case anticipated 
%      policy shocks are used to impose it; b) computes the probabilities 
%      of those paths; c) computes the expected paths.
%   -> The bulk of the work is done in computing the simulation paths,
%      which is done in a sub-function. Within that there is a sub-function
%      which solves the stacked-time inversion and a function that computes
%      the scalar value for the interest rate on reneging. Within the
%      former, there is a function that constructs the matrices in the
%      stack-time, which summarises most of the algebra in the paper. The
%      implementation of the occasionally binding constraint on renege is
%      done iteratively using the Holden and Paetz approach.
%
% NOTES:
%   -> The algorithm in this function (and accompanying paper) works under
%      perfect foresight only and so this function does not permit
%      unanticipated shocks.
%
% DEV NOTES:
%   -> The stack-time inversion used in this function could be upgraded to
%      use sparse in the event that multi-period simulations crash with
%      out-of-memory errors - tests so far have revealed that this is
%      unlikely to be necessary.
%   -> We have also identified that the implementation of the occasionally
%      binding ELB after renege could be upgraded to solve for the paths of
%      the anticipated policy shocks necessary to impose the ELB in all
%      states of the world simultaneously (rather than working
%      state-by-state).
%
% This version: 27/10/2015
% Author(s): Matt

%% CHECK INPUTS
% This could be upgraded at some point to include more checks if we start
% to have trouble with particular error cases.
if nargin < 9
   error([mfilename,' cannot proceed because too few inputs passed in']);
end

%% HANDLE OPTIONAL OPTIONS INPUT
DefaultPolicyShockIterOpts = get_default_policy_shock_iteration_options();
if nargin < 10
    PolicyShockIterOpts = DefaultPolicyShockIterOpts;
else
    PolicyShockIterOpts = overlay_default_structure(...
        PolicyShockIterOpts,DefaultPolicyShockIterOpts);
end

%% COMPUTE MAX NUMBER OF PERIODS OF LFL & FORECAST HORIZON
% In principle, the code could be generalised to allow for H<K+1 (on the
% assumption that all shocks after period H are 0). The trickiness in
% implementing that would be around the anticipated policy shocks
% necessary to impose any occasionally binding ELB on renege.
K = size(p,2);
H = size(z,2);
if H < K+1
    error([mfilename,' cannot proceed because the logic of the ',...
        'implementation requires that the forecast horizon (as ',...
        'measured by the dimension of the shocks passed in) be larger ',...
        'than the maximum possible number of periods in the LFL ',...
        'policy (as measured by the dimesnion of the renege ',...
        'probabilities passed in)']); 
elseif size(rbar,2) ~= K
    error([mfilename,' cannot proceed because the logic of the ',...
        'implementation requires that the vector of renege ',...
        'probabilities be of the same length as the vector of peg ',...
        'values']);
end

%% UNPACK MODEL CONTENT REQUIRED FOR SIMULATION
[HB,HC,HF,PSI,B,F,PHI,xEqNames,xMnems,zMnems] = unpack_model(Model,...
    {'HB','HC','HF','PSI','B','F','PHI','xEqNames','xMnems','zMnems'});

%% COMPUTE ALTERNATIVE PATHS
[xPaths,epsrPaths] = ...
    compute_simulation_paths_with_occasionally_binding_ELB(...
    HB,HC,HF,PSI,B,F,PHI,xEqNames,xMnems,zMnems,x0,z,p,rbar,b,...
    policyEqName,policyVarMnem,policyShockMnem,K,H,PolicyShockIterOpts);

%% COMPUTE PROBABILITIES OF PATHS
% Implemented with "if-else" statements to make it clear of the assumptions
% at either end of the sim: in period 1, the vector of probabilities
% represents the ex-ante outcome of each state. In all other periods, the
% probabilities are conditional on not already having reneged.
xPathProbs = NaN*ones(1,K+1,K+1);
for t = 1:K+1
    for i = t:K+1
        if t == 1
            if i == 1
                xPathProbs(1,1,1) = p(1);
            elseif i < K+1
                xPathProbs(1,1,i) = prod(1-p(1:i-1))*p(i);
            else
                xPathProbs(1,1,K+1) = prod(1-p);
            end
        else
            % Adjusted by RH to deal with issue that 
            % (xPathProbs(1,1,1:t-1)) may equal 1 to machine precision, for
            % finite t
            xPathProbs(1,t,i) = min(1,xPathProbs(1,1,i)/...
                (1-sum(xPathProbs(1,1,1:t-1))));
%             xPathProbs(1,t,i) = xPathProbs(1,1,i)/...
%                 (1-sum(xPathProbs(1,1,1:t-1)));
        end
    end
end

%% COMPUTE EXPECTATIONS
% These are expectations conditional on the policy-maker not already having
% reneged up to that point in time. For instance, estarPaths(:,5,1),
% represents agents' expectations for the endogenous variables in period 5, 
% conditional on information available at the beginning of period 1. And 
% estarPaths(:,5,4) represents agents' expectations for the endogenous 
% variables in period 5, conditional on information available at the 
% beginning of period 4 in the state of the world where the policy maker 
% has not already reneged.
nx = size(B,1);
eStarPaths = NaN*ones(nx,H,K);
for s = 1:K
    for t = s:H
        eStarPaths(:,t,s) = zeros(nx,1);
        for i = s:K+1
            eStarPaths(:,t,s) = eStarPaths(:,t,s)+...
                xPathProbs(1,s,i)*xPaths(:,t,i);
        end
    end
end

end

%% FUNCTION TO COMPUTE COMPLETE SET OF SIMULATION PATHS
function [xPaths,epsrPaths] = ...
    compute_simulation_paths_with_occasionally_binding_ELB(...
    HB,HC,HF,PSI,B,F,PHI,xEqNames,xMnems,zMnems,x0,z,p,rbar,b,...
    policyEqName,policyVarMnem,policyShockMnem,K,H,PolicyShockIterOpts)
% This function contains the logic for computation of all simulation paths.
% It makes allowance for the fact that the policy rule may not respect the
% ELB on renege.
%
% INPUTS:   
%   -> HB: loadings on lags of endogenous vars from structural equations 
%   -> HC: loadings on time t endogenous vars from structural equations 
%   -> HF: loadings on leads of endogenous vars from structural equations 
%   -> PSI: loadings on shocks from structural equations 
%   -> B: loadings on lagged endogenous variables from model solution 
%   -> F: forward loadings on shock loadings from model solution 
%   -> PHI: loadings on shocks from model solution 
%   -> xEqNames: names of the model equations
%   -> xMnems: mnemonics of the endogenous variables
%   -> x0: nx*1 vector of initial conditions for the endogneous variables
%   -> z: nz*H matrix of anticipated shocks
%   -> p: 1*K row vector of renege probabilities
%   -> rbar: 1*K row vector of instrument peg values (in model space)
%   -> b: scalar value for the ELB on interest rates
%   -> policyEqName: name of the policy rule equation
%   -> policyVarMnem: mnemonic of the policy variable
%   -> policyShockMnem: mnemonic of the policy shock in the model
%   -> K: maximum number of periods for LFL policy
%   -> H: forecast/simulation horizon 
%   -> PolicyShockIterOpts: structure of options for iteration over policy
%      shocks necessary to impose the ELB in event that it binds on renege:
%       - bViolationTol: amount be which ELB is allowed to be violated
%       - maxIter: maximum number of iterations allowed
%       - bifurcationFrac: numeric scalar controlling speed of update
%       - epsrChangeTol: exit criteria tolerance for change in shock values
%         between iterations
%       - displayProgress: true/false
%
% OUTPUTS:  
%   -> xPaths: nx*H*(K+1) matrix of simulation paths
%   -> epsrPaths: 1*H*(K+1) matrix of values for the policy shock

%% LOOKUP INDEX NUMBER OF POLICY RULE, VARIABLE & SHOCK IN MODEL
% Note that these functions will throw an error if the equation name,
% policy shock or policy variables mnemonics are not part of the model
% passed in (i.e. this also acts as a compatibility check).
policyEqInd = lookup_index_numbers_in_string_array(...
    xEqNames,policyEqName,true);
policyVarInd = lookup_index_numbers_in_string_array(...
    xMnems,policyVarMnem,true);
policyShockInd = lookup_index_numbers_in_string_array(...
    zMnems,policyShockMnem,true);

%% COMPUTE DIMENSION OF ENDOGENOUS VARIABLES & SHOCKS
[nx,nz] = size(PHI);

%% SET LOGICAL VECTORS FOR POLICY & NON-POLICY VARIABLES, SHOCKS & EQS
nonPolicyVarLogics = true(nx,1);
nonPolicyShockLogics = true(nz,1);
nonPolicyEqLogics = true(nx,1);
nonPolicyVarLogics(policyVarInd) = false;
nonPolicyShockLogics(policyShockInd) = false;
nonPolicyEqLogics(policyEqInd) = false;
policyVarLogics = ~nonPolicyVarLogics;
policyShockLogics = ~nonPolicyShockLogics;
policyEqLogics = ~nonPolicyEqLogics;

%% PARTITION MATRICES
HtildeBxtilde = HB(nonPolicyEqLogics,nonPolicyVarLogics);
HtildeBr = HB(nonPolicyEqLogics,policyVarLogics);
HtildeCxtilde = HC(nonPolicyEqLogics,nonPolicyVarLogics);
HtildeCr = HC(nonPolicyEqLogics,policyVarLogics);
HtildeFxtilde = HF(nonPolicyEqLogics,nonPolicyVarLogics);
HtildeFr = HF(nonPolicyEqLogics,policyVarLogics);
PSItildeztilde = PSI(nonPolicyEqLogics,nonPolicyShockLogics);
HhatBxtilde = HB(policyEqLogics,nonPolicyVarLogics);
HhatBr = HB(policyEqLogics,policyVarLogics);
HhatCxtilde = HC(policyEqLogics,nonPolicyVarLogics);
HhatCr = HC(policyEqLogics,policyVarLogics);
HhatFxtilde = HF(policyEqLogics,nonPolicyVarLogics);
HhatFr = HF(policyEqLogics,policyVarLogics);
PSIhatepsr = PSI(policyEqLogics,policyShockLogics);
PSIhatztilde = PSI(policyEqLogics,nonPolicyShockLogics);
PSItildeepsr = PSI(nonPolicyEqLogics,policyShockLogics);
Bxtildextilde = B(nonPolicyVarLogics,nonPolicyVarLogics);
Bxtilder = B(nonPolicyVarLogics,policyVarLogics);
Brxtilde = B(policyVarLogics,nonPolicyVarLogics);
Brr = B(policyVarLogics,policyVarLogics);
PHIztilde = PHI(:,nonPolicyShockLogics);
PHIepsr = PHI(:,policyShockLogics);

%% CHECK THAT MODEL IS VALID GIVEN LOGIC OF IMPLEMENTATION
% This is also a check that the mnemonic and names passed in were in fact
% correct (given the model being used).
if HhatCr == 0
    error([mfilename,' cannot proceed because the interest rate,''',...
        policyVarMnem,''', does not appear contemporaneously in the ',...
        'policy rule equation, ''',policyEqName,'''']); 
elseif PSIhatepsr == 0
    error([mfilename,' cannot proceed because the policy shock,''',...
        policyShockMnem,''', does not appear in the policy rule ',...
        'equation, ''',policyEqName,'''']); 
elseif any(PSItildeepsr)
    error([mfilename,' cannot proceed because the policy shock,''',...
        policyShockMnem,''', appears in at least one non-policy rule ',...
        'equation']);  
elseif any(PSIhatztilde)
    error([mfilename,' cannot proceed because at least one non-policy ',...
        'shock appears in the policy rule equation, ''',...
        policyEqName,'''']); 
end

%% PARTITION INITIAL CONDITIONS FOR THE ENDOGENOUS VARIABLES
xtilde0 = x0(nonPolicyVarLogics);
r0 = x0(policyVarLogics);

%% PARTITION THE SHOCKS
ztilde = z(nonPolicyShockLogics,:);

%% CREATE SELECTOR MATRICES FOR SHOCK IMPACTS
S = eye(nx,nx);
Sxtilde = S(nonPolicyVarLogics,:);
Sr = S(policyVarLogics,:); 

%% COMPUTE NON-POLICY SHOCK IMPACT MATRIX
% The matrix here is an nx*H matrix of non-policy (anticipated) shock 
% impacts, \widetilde{\mathbb{Z}} in the paper. Note that it is a not a
% function of the policy shocks, so we only need to compute it once.
ZZtilde = compute_non_policy_anticipated_shocks_impact_matrix(...
    PHIztilde,ztilde,F,nx,H);

%% COMPUTE COEFFICIENTS MATRIX IN STACK-TIME INVERSION
% Note that as in the non-policy shocks impact matrix, this is not a
% function of the policy shocks so need only be computed once.
JJ = construct_stack_time_coefficients_for_inversion(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    HhatBxtilde,HhatCxtilde,HhatCr,HhatFxtilde,HhatFr,Bxtildextilde,...
    Bxtilder,Brxtilde,Brr,xtilde0,p,K);

%% SET GUESS FOR POLICY SHOCKS REQUIRED TO MAINTAIN ELB IN CASE OF RENEGE
% The policy shock necessary to impose an oaccasionally binding ELB is not
% defined in periods before the renege has happened (i.e. for t<i), so the
% values in these cases are set to NaNs. This is also useful as a very
% rough check that the code below is working properly - if the answer comes
% back as NaN, then one reason could be that the code is incorrectly 
% picking out invalid shock values.
epsrPaths = NaN*ones(1,H,K+1);
for i = 1:K+1
    for t = i:H
        epsrPaths(1,t,i) = 0;
    end
end

%% COMPUTE SIMULATION PATHS CONDITIONAL ON GUESS AT POLICY SHOCKS ON RENEGE
xPaths = compute_simulation_paths(...
    JJ,HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,HhatBxtilde,HhatBr,HhatCxtilde,HhatCr,HhatFxtilde,...
    HhatFr,PSIhatepsr,B,Bxtildextilde,Bxtilder,Brxtilde,Brr,F,PHIepsr,...
    xtilde0,r0,ztilde,epsrPaths,ZZtilde,Sxtilde,Sr,nonPolicyVarLogics,...
    policyVarLogics,p,rbar,K,H,nx);

%% CHECK THAT THE ELB IS SATISFIED ON EXIT (WITHIN TOLERANCE) 
% If the ELB constraint is satisfied, then we are done.
bViolationTol = PolicyShockIterOpts.bViolationTol;
rPaths = xPaths(policyVarLogics,:,:);
if all(all(rPaths>b-bViolationTol))
    return
end

%% PREPARE FOR ITERATION OVER ANTICIPATED POLICY SHOCKS
% Unpack iteration options, print out message to user and get anticipated
% policy shocks multiplier matrix for implementation of the quadratic
% programming problem solved on each iteration.
maxIter = PolicyShockIterOpts.maxIter;
bifurcationFrac = PolicyShockIterOpts.bifurcationFrac;
epsrChangeTol = PolicyShockIterOpts.epsrChangeTol;
displayProgress = PolicyShockIterOpts.displayProgress;
if displayProgress
    fprintf('The interest rate violated the ELB in at least one state\n');
end
aPolicyShockMultMatForQuadProg = ...
    get_anticipated_policy_shock_multiplier_matrix(HhatCr,...
    HhatFxtilde,HhatFr,PSIhatepsr,B,Bxtilder,Brr,F,PHIepsr,Sxtilde,Sr,K,H);

%% ITERATE OVER POLICY SHOCKS
%% ITERATE OVER POLICY SHOCKS
isELBviolatedOnExit = true;
iter = 0;
epsrMaxChange = inf;

% >>> 추가: 직전 경로 보관 변수 초기화 (에러 원인 해결)
rPaths = xPaths(policyVarLogics,:,:);   % 위에서 이미 계산된 xPaths 기반
rPathsOld = rPaths;                     % 이전 금리 경로
epsrPathsOld = epsrPaths;               % 이전 정책충격 경로
rMaxChange = Inf;

while (isELBviolatedOnExit && epsrMaxChange > epsrChangeTol) && iter < maxIter
    iter = iter+1;
    if displayProgress
        fprintf('\tIteration %d\n', iter);
    end

    % 새 QP 해 구하기 (소프트 제약/폴백 버전이면 그 코드 사용)
    [epsrPathsNew,isSolutionAdmissable] = ...
        compute_policy_shocks_to_enforce_occasionally_binding_ELB( ...
            epsrPaths, aPolicyShockMultMatForQuadProg, rPaths, b, K, H);

    % 언더릴랙스(이완) 업데이트
    epsrPaths = bifurcationFrac*epsrPathsNew + (1-bifurcationFrac)*epsrPaths;

    % 새 경로 재계산
    xPaths = compute_simulation_paths( ...
        JJ,HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,HtildeFr, ...
        PSItildeztilde,HhatBxtilde,HhatBr,HhatCxtilde,HhatCr,HhatFxtilde, ...
        HhatFr,PSIhatepsr,B,Bxtildextilde,Bxtilder,Brxtilde,Brr,F,PHIepsr, ...
        xtilde0,r0,ztilde,epsrPaths,ZZtilde,Sxtilde,Sr, ...
        nonPolicyVarLogics,policyVarLogics,p,rbar,K,H,nx);

    rPaths = xPaths(policyVarLogics,:,:);

    % ELB 위반 체크
    if ~any(any(rPaths < b - bViolationTol))
        isELBviolatedOnExit = false;
        if displayProgress
            fprintf('\t\tELB was respected in all periods\n');
        end
    else
        maxViolation = max(max(b - rPaths)) - bViolationTol;
        if displayProgress
            fprintf('\t\tMax violation of ELB: %g\n', maxViolation);
        end
    end

    % >>> 변경량 계산 (이전값과 비교) — 여기서 rPathsOld가 이미 정의됨
    rMaxChange    = max(max(abs(rPaths - rPathsOld)));
    epsrMaxChange = max(max(abs(epsrPaths - epsrPathsOld)));

    if displayProgress
        fprintf('\t\tMax change in rates is: %g\n', rMaxChange);
    end

    % >>> 다음 반복을 위해 "이전값" 갱신
    rPathsOld    = rPaths;
    epsrPathsOld = epsrPaths;

    % 조기 종료(수렴 감지)
    if iter >= 3 && rMaxChange < 1e-6 && epsrMaxChange < epsrChangeTol
        break;
    end
end


    
end

%% FUNCTION TO COMPUTE SIMULATION PATHS CONDITIONAL ON GUESS AT ELB BINDING
function xPaths = compute_simulation_paths(...
    JJ,HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,HhatBxtilde,HhatBr,HhatCxtilde,HhatCr,HhatFxtilde,...
    HhatFr,PSIhatepsr,B,Bxtildextilde,Bxtilder,Brxtilde,Brr,F,PHIepsr,...
    xtilde0,r0,ztilde,epsrPaths,ZZtilde,Sxtilde,Sr,nonPolicyVarLogics,...
    policyVarLogics,p,rbar,K,H,nx)
% This function contains the logic for computation of the non-renege paths.
%
% INPUTS:   
%   -> JJ: (nxtilde*K)*(nxtilde*K) loadings on stacked non-renege paths
%   -> HtildeBxtilde: loadings on lags of non-policy vars in non-policy eqs
%   -> HtildeBr: loadings on lag of policy var in non-policy eqs 
%   -> HtildeCr: loadings on time t policy var in non-policy eqs
%   -> HtildeFxtilde:loadings on leads of non-policy vars in non-policy eqs
%   -> HtildeFr: loadings on lead of policy var in non-policy eqs
%   -> PSItildeztilde: loadings on non-policy shocks in non-policy eqs
%   -> HhatBxtilde: loadings on lags of non-policy vars in policy eq
%   -> HhatBr: loadings on lag of policy var in policy eq
%   -> HhatCxtilde: loadings on time t non-policy vars in policy eq
%   -> HhatCr: loadings on time t policy var in policy eq
%   -> HhatFxtilde: loadings on leads of non-policy vars in policy eq
%   -> HhatFr: loadings on lead of policy var in policy eq
%   -> PSIhatepsr: loading on policy shock in policy eq
%   -> B: loadings on lagged endogenous variables from model solution 
%   -> Bxtildextilde: loadings of non-policy vars on non-policy vars
%   -> Bxtilder: loading of policy var on non-policy vars
%   -> Brxtilde: loadings of non-policy vars on policy var
%   -> Brr: loading of policy var on policy var
%   -> F: forward loadings on shock loadings from model solution 
%   -> PHIepsr: loadings on policy shock from model solution
%   -> xtilde0: initial condition for non-policy vars
%   -> r0: initial condition for policy var
%   -> ztilde: (nz-1)*H matrix of non-policy shocks
%   -> epsrPaths: 1*H*(K+1) matrix of shocks to impose ZLB after renege
%   -> ZZtilde: nx*H matrix of non-policy shock impacts
%   -> Sxtilde: (nx-1)*nx selector matrix for non-policy variables
%   -> Sr: 1*nx selector matrix for policy variable
%   -> nonPolicyVarLogics: nx*1 logical vector for the non-policy variables
%   -> policyVarLogics: nx*1 logical vector for the non-policy variables
%   -> p: 1*K row vector of renege probabilities
%   -> rbar: 1*K row vector of instrument peg values (in model space)
%   -> K: maximum number of periods for LFL policy
%   -> H: forecast horizon
%   -> nx: number of endogenous variables
%
% OUTPUTS:  
%   -> xPaths: nx*H*(K+1) matrix of simulation paths

%% COMPUTE POLICY SHOCK IMPACT MATRICES
% The matrix here is an nx*H*K+1 matrix of impacts of the policy shocks 
% conditional on reversion to the rule in each of the K+1 states and 
% necessary if the ELB is violated on reversion to the rule (this is \Xi^r 
% in the paper).
XIr = compute_anticipated_policy_shock_impact_matrices(...
    PHIepsr,epsrPaths,F,nx,H,K);

%% COMPUTE NON-RENEGE PATH USING STACK TIME
% Compute the non-renege path, \widetilde{X}^{\ast}, from periods 1 to K 
% (this is a stacked (nx-1)*K column vector).
XtildeStar = compute_non_renege_path_using_stack_time(...
    JJ,HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,HhatBxtilde,HhatBr,HhatCr,HhatFxtilde,HhatFr,...
    PSIhatepsr,Bxtilder,Brr,xtilde0,r0,ztilde,epsrPaths,ZZtilde,XIr,...
    Sxtilde,Sr,p,rbar,K);

%% COMPUTE ENDOGENOUS VARIABLES CONDITIONAL ON NO RENEGE & RENEGE
% These are the vectors, x^{\ast}_{i} & x^{<i>}_{i} for i=1...K. In the
% case of no renege, the vector is the relevant segment of the stack-time
% solution augmented with the ELB for the policy rate. In the renege case,
% the vector still cntains the stack-time solution (because the private
% sector does not know about renege within period because the policy maker
% reneges at the end of the period) augmented with a scalar for the
% interest rate on renege, which is computed in a sub-function.
xstar = NaN*ones(nx,K);
xOnRenege = NaN*ones(nx,K);
for t = 1:K
    xstar(nonPolicyVarLogics,t) = XtildeStar((t-1)*(nx-1)+1:t*(nx-1));
    xstar(policyVarLogics,t) = rbar(t);
    xOnRenege(nonPolicyVarLogics,t) = xstar(nonPolicyVarLogics,t);
    if t == 1
        xtildeLag = xtilde0;
        rLag = r0;
    else
        xtildeLag = xstar(nonPolicyVarLogics,t-1);
        rLag = rbar(t-1);
    end
    xOnRenege(policyVarLogics,t) = ...
        compute_policy_rate_within_period_on_reneging(...
        HhatBxtilde,HhatBr,HhatCxtilde,HhatCr,HhatFxtilde,HhatFr,...
        PSIhatepsr,Bxtildextilde,Bxtilder,Brxtilde,Brr,xtildeLag,rLag,...
        xstar(nonPolicyVarLogics,t),epsrPaths(1,t,t),ZZtilde(:,t+1),...
        XIr(:,t+1,t),Sxtilde,Sr);
end

%% COMPUTE THE COMPLETE SET OF SIMULATION PATHS
% The logic of this is as follows: a) all paths are the same and are given
% by the non-renege path up until renege occurs; b) when renege occurs, the
% vector is given by the vector on reneging computed above (if i<=K); 
% c) after renege the path is just the RE solution; d) in K+1 reversion to
% the policy rule happens with certainty and the solution is just the RE 
% solution.
xPaths = NaN*ones(nx,H,K+1);
for i = 1:K+1
    for t = 1:H
        if t < i
            xPaths(:,t,i) = xstar(:,t);
        elseif i<K+1 && t==i
            xPaths(:,t,i) = xOnRenege(:,t);
        else
            xPaths(:,t,i) = B*xPaths(:,t-1,i)+ZZtilde(:,t)+XIr(:,t,i);
        end
    end
end

end

%% FUNCTION TO COMPUTE ANTICIPATED POLICY SHOCKS MULTIPLIERS MATRIX
function aPolicyShockMultMatForQuadProg = ...
    get_anticipated_policy_shock_multiplier_matrix(HhatCr,HhatFxtilde,...
    HhatFr,PSIhatepsr,B,Bxtilder,Brr,F,PHIepsr,Sxtilde,Sr,K,H)
% This function creates the policy shock multiplier mat for quadprog.
%
% INPUTS:
%   -> HhatCr: loadings on time t policy var in policy eq
%   -> HhatFxtilde: loadings on leads of non-policy vars in policy eq
%   -> HhatFr: loadings on lead of policy var in policy eq
%   -> PSIhatepsr: loading on policy shock in policy eq
%   -> B: loadings on lagged endogenous variables from model solution 
%   -> Bxtilder: loading of policy var on non-policy vars
%   -> Brr: loading of policy var on policy var
%   -> F: forward loadings on shock loadings from model solution 
%   -> PHIepsr: loadings on policy shock from model solution
%   -> Sxtilde: (nx-1)*nx selector matrix for non-policy variables
%   -> Sr: 1*nx selector matrix for policy variable
%   -> K: maximum number of periods for LFL policy
%   -> H: forecast horizon
%
% OUTPUTS:  
%   -> aPolicyShockMultMatForQuadProg: H*H*(K+1) impact matrix for
%      anticipated policy shocks on renege

%% BUILD HELPER MATRICES
% Here, we build 5 helper vectors (2 of which are scalars) to avoid
% repetition of the matrix calculations below.  This follows from noting
% that the impact of a unit anticipated policy shock one period after 
% renege has happened is the regardless of when renege happens.
HHhatr = HhatCr+HhatFr*Brr+HhatFxtilde*Bxtilder;
cLoadingAfterRenege = Sr*PHIepsr;
cLoadingOnRenege = PSIhatepsr/HHhatr;
bLoadingsAfterRenege = NaN*ones(H-1,1);
fLoadingsAfterRenege = NaN*ones(H-1,1);
fLoadingsOnRenege = NaN*ones(H-1,1);
bLoadingsCumulator = PHIepsr;
fLoadingsCumulator = PHIepsr;
for s = 1:H-1
    fLoadingsOnRenege(s) = ...
        -((HhatFxtilde*Sxtilde+HhatFr*Sr)/HHhatr)*fLoadingsCumulator;
    bLoadingsCumulator = B*bLoadingsCumulator;
    fLoadingsCumulator = F*fLoadingsCumulator;        
    bLoadingsAfterRenege(s) = Sr*bLoadingsCumulator;  
    fLoadingsAfterRenege(s) = Sr*fLoadingsCumulator;
end

%% INITIALISE MATRIX FOR IMPACT OF EACH POLICY SHOCK ON POLICY RATE
% Note that this matrix stores the impact of a unit policy shock on the
% policy rate in each period separately for each possible state of the 
% world. Indexing is as follows: aPolicyShockMultMatForQuadProg(t,s,i)
% stores the impact of an anticipated policy shock in period s on the
% policy rate in period t conditional on renege in period i. For examples:
% aPolicyShockMultMatForQuadProg(:,1,1) describes the impact of a policy 
% shock realised in period 1 conditional on renege in period 1 on the 
% policy rate in all H periods; aPolicyShockMultMatForQuadProg(:,2,1) 
% describes the impact of an anticipated policy shock realised in period 2
% conditional on renege in period 1 on the policy rate in all H periods; 
% aPolicyShockMultMatForQuadProg(:,1,2) is equal NaNs because conditional
% on renege in period 2, there cannot exist anticipated policy shocks in
% period 1.
aPolicyShockMultMatForQuadProg = NaN*ones(H,H,K+1);

%% COMPUTE MULTIPLIERS MATRIX FOR EACH SEPARATE POLICY SHOCK REALISATION
% The anticipated policy shocks matrix is built from 5 parts labelled
% below. This is necessary because: a) the impact of past, contemporaneous 
% and future anticipated shocks is different; b) the impact of
% contemporaneous and future shocks on the policy rate within the period of
% renege is different reflecting the timing assumption that the private 
% sector makes decisions before renege happens. However, note that the
% timing assumption does not apply in renege period K+1 because renege
% happens in that period with certainty and so the RE solution applies
% there.
for i = 1:K+1
    % i) Contemporaneous impact of shock within renege period
    if i < K+1
        aPolicyShockMultMatForQuadProg(i,i,i) = cLoadingOnRenege;
    else
        aPolicyShockMultMatForQuadProg(i,i,i) = cLoadingAfterRenege;
    end
    % ii) Impact of future shocks within renege period
    if i < K+1
        for s = i+1:H
            aPolicyShockMultMatForQuadProg(i,s,i) = fLoadingsOnRenege(s-i);
        end
    else
        for s = i+1:H
            aPolicyShockMultMatForQuadProg(i,s,i) = ...
                fLoadingsAfterRenege(s-i);
        end      
    end  
    % iii) Impact of contemporaneous shocks after renege period
    for s = i+1:H
        aPolicyShockMultMatForQuadProg(s,s,i) = cLoadingAfterRenege;
    end
    % iv) Impacts of future shocks after renege period
    for s = i+2:H
        for t = i+1:s-1
            aPolicyShockMultMatForQuadProg(t,s,i) = ...
                fLoadingsAfterRenege(s-t);
        end
    end 
    % iv) Impact of past shocks after renege period
    for s = i:H-1
        for t = s+1:H
            aPolicyShockMultMatForQuadProg(t,s,i) = ...
                bLoadingsAfterRenege(t-s);
        end
    end
end

end

%% FUNCTION TO COMPUTE POLICY SHOCKS TO ENFORCE OCCASIONALLY BINDING ELB
function [epsrPathsNew,isSolutionAdmissable] = ...
    compute_policy_shocks_to_enforce_occasionally_binding_ELB(...
    epsrPaths,aPolicyShockMultMatForQuadProg,rPaths,b,K,H)
% This function uses quadprog to compute policy shocks to impose ELB.
%
% INPUTS:
%   -> epsrPaths: 1*H*(K+1) matrix of existing guesses for shocks to impose 
%      ZLB after renege
%   -> aPolicyShockMultMatForQuadProg: H*H*(K+1) impact matrix for
%      anticipated policy shocks on renege
%   -> rPaths: interest rate paths consistent with existing guesses for
%      policy shocks
%   -> b: scalar value for the ELB on interest rates
%   -> K: maximum number of periods for LFL policy
%   -> H: forecast horizon
%
% OUTPUTS:  
%   -> epsrPathsNew: policy shocks updated with new best guesses
%   -> isSolutionAdmissable: true/false

%% RECOMPUTE PATHS FOR POLICY RATE WITH POLICY SHOCKS STRIPPED OUT
% This cell recomputes the paths for the policy rate in all renege states
% to exclude the effect of the existing guesses for anticipated shocks 
% necessary to impose the ZLB.  This wipes the slate clean for the quadprog
% routine below to avoid situations where the ELB becomes "over-binding" in
% which case the quadprog routine would return zeros in a marginal update.
rPathsNoPolicyShocks = rPaths;
for i = 1:K+1
    for t = i:H
        rPathsNoPolicyShocks(1,t,i) = rPaths(1,t,i)-...
            aPolicyShockMultMatForQuadProg(t,i:H,i)*epsrPaths(1,i:H,i)';
    end
end

%% INITIALISE UPDATED POLICY SHOCKS OUTPUT
epsrPathsNew = epsrPaths;

%% SETUP OPTIONS FOR THE QUADRATIC PROGRAMMING OPTIMISATION
% If necessary, this could be made to be a user-defined input.
options = optimset(...
    'MaxIter',10^4,...
    'Display','off',...
    'Algorithm','interior-point-convex'); 

%% COMPUTE POLICY SHOCKS TO IMPLEMENT THE ELB AFTER EXIT IN EACH STATE
% The quadratic programming problem setup below returns the vector of
% shocks that minimises 0.5*x'*H*x + f'*x s.t. A*x<=b & x>=0. In words, it
% finds the minimum size of shocks (which must be positive) that deliver
% the interest rate at or above the ELB. Notice that that is a statement
% that describes a greater than or equal inequality constraint, which is
% why the sign of the policy shock impact multipliers is flipped. Note also
% that the certain exit of the policy-maker in period K+1 creates a
% discontinuity because the RE solution applies in every period from K+1
% onwards whereas, in all periods prior to K+1, the RE solution only
% applies in each period after exit has happened (because exit is a
% surprise to agents and they make their decisions before they observe the
% outcome for policy). Note that in some circumstances quadprag may fail to
% find a sequence of shocks which satisfy the ELB in which case an error is
% thrown. Alternatively, it may find a sequence that satisfies the ELB, but
% which is not a minimum. In those circumstances, a solution admissability
% flag is set to true on output to this function (which does not
% provide any information as to the period in which the failure ocurred).
isSolutionAdmissable = true;
for i = 1:K+1
    % --- 0) ELB 위반 없으면 스킵 ---
    fquad = rPathsNoPolicyShocks(1,i:H,i)' - b;   % (H-i+1)x1
    if all(fquad >= -1e-12)
        epsrPathsNew(1,i:H,i) = epsrPaths(1,i:H,i);
        continue
    end

    % --- 1) 준비 ---
    M = aPolicyShockMultMatForQuadProg(i:H,i:H,i);  % m x n
    [m,n] = size(M);

    % NaN/Inf 방어: 있으면 폴백
    if any(~isfinite(fquad)) || any(~isfinite(M(:)))
        epsrPathsNew(1,i:H,i) = epsrPaths(1,i:H,i);
        isSolutionAdmissable = false;
        continue
    end

    % 컬럼 스케일링(조건수 개선): x = Dinv*x_sc  →  M x = Msc x_sc
    colscale = sqrt(sum(M.^2,1))'; colscale(colscale==0) = 1;
    D    = diag(colscale);
    Dinv = diag(1./colscale);
    Msc  = M * Dinv;

    % --- 2) 소프트 ELB: r_no + M x + s ≥ b (s ≥ 0) ---
    % 목적: ||x_sc||^2 + mu ||s||^2 최소화  (항상 볼록/feasible)
    mu   = 1e4;                             % slack 패널티 (필요시 1e5~1e6)
    Hqp  = blkdiag(eye(n), mu*eye(m));
    Hqp  = Hqp + 1e-9*eye(n+m);             % 릿지로 수치 안정화
    f_lin = zeros(n+m,1);
    Aineq = [-Msc, -eye(m)];
    bineq = -fquad;
    lb = [-inf(n,1); zeros(m,1)];           % x 자유, s ≥ 0
    ub =  inf(n+m,1);

    % 웜스타트: 기존 epsr 이용
    x0_sc = D * (epsrPaths(1,i:H,i)');      % n x 1
    s0    = max(0, -fquad);                 % 위반분만큼 s 초기화
    x0    = [x0_sc; s0];

    % --- 3) Solver 1: interior-point-convex ---
    opts1 = optimset('MaxIter',1e3,'Display','off','Algorithm','interior-point-convex', ...
                     'TolFun',1e-8,'TolX',1e-8,'TolCon',1e-8);
    [xs,~,exitflag] = quadprog(Hqp,f_lin,Aineq,bineq,[],[],lb,ub,x0,opts1);

    % --- 4) Solver 2: active-set 폴백 ---
    if exitflag ~= 1 || any(~isfinite(xs))
        opts2 = optimset('MaxIter',1e3,'Display','off','Algorithm','active-set', ...
                         'TolFun',1e-8,'TolX',1e-8,'TolCon',1e-8);
        [xs,~,exitflag] = quadprog(Hqp,f_lin,Aineq,bineq,[],[],lb,ub,x0,opts2);
    end

    % --- 5) 최종 폴백: slack만으로 강제 (절대 에러 던지지 않음) ---
    if exitflag ~= 1 || any(~isfinite(xs))
        x_sc = zeros(n,1);
        s    = max(0, -fquad);
        xs   = [x_sc; s];
        isSolutionAdmissable = false;    % 상위 루틴에서 경고만 출력
    end

    % --- 6) 해 복원/적용 ---
    x_sc = xs(1:n);  s = xs(n+1:end);
    x    = Dinv * x_sc;                   % 원 스케일로 환산
    epsrPathsNew(1,i:H,i) = x';

    % slack 남아있으면 '완전한' 해는 아님(경고용 플래그)
    if any(s > 1e-6)
        isSolutionAdmissable = false;
    end
end


end

%% FUNCTION TO COMPUTE THE NON-RENEGE PATH BY STACK-TIME
function XtildeStar = compute_non_renege_path_using_stack_time(...
    JJ,HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,HhatBxtilde,HhatBr,HhatCr,HhatFxtilde,HhatFr,...
    PSIhatepsr,Bxtilder,Brr,xtilde0,r0,ztilde,epsrPaths,ZZtilde,XIr,...
    Sxtilde,Sr,p,rbar,K)
% This function contains the stack-time to compute the non-renege paths.
%
% INPUTS:   
%   -> JJ: (nxtilde*K)*(nxtilde*K) loadings on stacked non-renege paths
%   -> HtildeBxtilde: loadings on lags of non-policy vars in non-policy eqs
%   -> HtildeBr: loadings on lag of policy var in non-policy eqs
%   -> HtildeCr: loadings on time t policy var in non-policy eqs
%   -> HtildeFxtilde:loadings on leads of non-policy vars in non-policy eqs
%   -> HtildeFr: loadings on lead of policy var in non-policy eqs
%   -> PSItildeztilde: loadings on non-policy shocks in non-policy eqs
%   -> HhatBxtilde: loadings on lags of non-policy vars in policy eq
%   -> HhatBr: loadings on lag of policy var in policy eq
%   -> HhatCr: loadings on time t policy var in policy eq
%   -> HhatFxtilde: loadings on leads of non-policy vars in policy eq
%   -> HhatFr: loadings on lead of policy var in policy eq
%   -> PSIhatepsr: loading on policy shock in policy eq
%   -> Bxtilder: loading of policy var on non-policy vars
%   -> Brr: loading of policy var on policy var
%   -> xtilde0: initial condition for non-policy vars
%   -> r0: initial condition for policy var
%   -> ztilde: (nz-1)*H matrix of non-policy shocks
%   -> epsrPaths: 1*H*(K+1) matrix of shocks to impose ZLB on renege
%   -> ZZtilde: nx*H matrix of non-policy shock impacts
%   -> XIr: nx*H*(K+1) matrix of state-dependent policy shock impacts
%   -> Sxtilde: (nx-1)*nx selector matrix for non-policy variables
%   -> Sr: 1*nx selector matrix for policy variable
%   -> p: 1*K row vector of renege probabilities
%   -> rbar: 1*K row vector of instrument peg values (in model space)
%   -> K: maximum number of periods for LFL policy
%
% OUTPUTS:  
%   -> XtildeStar: (nx*K)*1 stacked vector of non-renege paths

%% CONSTRUCT CONSTANTS MATRIX FOR STACK-TIME INVERSION
CC = construct_stack_time_constants_for_inversion(...
    HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,HhatBxtilde,HhatBr,HhatCr,HhatFxtilde,HhatFr,...
    PSIhatepsr,Bxtilder,Brr,xtilde0,r0,ztilde,epsrPaths,ZZtilde,XIr,...
    Sxtilde,Sr,p,rbar,K);

%% INVERT
try
    XtildeStar = JJ\CC;
catch InversionE
   BigProblemE = MException(mfilename,'Sadly, the stack-time ',...
       'inversion failed: this is not good. Here is what was the ',...
       'MATLAB error on attempted inversion. You''re going to need ',...
       'to debug and work out what''s gone wrong. Good luck.');
   BigProblemE = addCause(BigProblemE,InversionE);
   throw(BigProblemE);
end

end

%% FUNCTION TO COMPUTE POLICY RATE ON RENEGEING
function rOnRenege = compute_policy_rate_within_period_on_reneging(...
    HhatBxtilde,HhatBr,HhatCxtilde,HhatCr,HhatFxtilde,HhatFr,PSIhatepsr,...
    Bxtildextilde,Bxtilder,Brxtilde,Brr,xtildeLag,rLag,xtilde,...
    epsrOnRenege,ZZtildeLead,XIrLead,Sxtilde,Sr)
% This function computes the policy rate within period on reneging.
%
% INPUTS:   
%   -> HhatBxtilde: loadings on lags of non-policy vars in policy eq
%   -> HhatBr: loadings on lag of policy var in policy eq
%   -> HhatCxtilde: loadings on time t non-policy vars in policy eq
%   -> HhatCr: loadings on time t policy var in policy eq
%   -> HhatFxtilde: loadings on leads of non-policy vars in policy eq
%   -> HhatFr: loadings on lead of policy var in policy eq
%   -> PSIhatepsr: loading on policy shock in policy eq
%   -> Bxtildextilde: loadings of non-policy vars on non-policy vars
%   -> Bxtilder: loading of policy var on non-policy vars
%   -> Brxtilde: loadings of non-policy vars on policy var
%   -> Brr: loading of policy var on policy var
%   -> xtildeLag: vector of non-policy vars in previous period
%   -> rLag: lagged policy rate
%   -> xtilde: vector of non-policy variables in current period
%   -> epsrOnRenege: scalar existing guess at policy shock on renege
%   -> ZZtildeLead: nx*1 matrix of non-policy shock impacts in next period 
%   -> XIrLead: nx*1 matrix of policy shock impacts in next period
%      (conditional on renege this period)
%   -> Sxtilde: (nx-1)*nx selector matrix for non-policy variables
%   -> Sr: 1*nx selector matrix for policy variable
%
% OUTPUTS:  
%   -> rOnRenege: scalar value for the interest rate on reneging

%% EQUATION FOR INTEREST RATE ON RENEGING
% This uses the structural equation for the policy rate together with the
% fact that the endogenous variables follow the RE solution after renege.
rOnRenege = (...
    PSIhatepsr*epsrOnRenege-HhatBr*rLag-HhatBxtilde*xtildeLag-...
    (HhatCxtilde+HhatFr*Brxtilde+HhatFxtilde*Bxtildextilde)*xtilde-...
    (HhatFxtilde*Sxtilde+HhatFr*Sr)*(XIrLead+ZZtildeLead))/...
    (HhatCr+HhatFr*Brr+HhatFxtilde*Bxtilder);

end

%% FUNCTION TO COMPUTE IMPACT OF NON-POLICY SHOCKS
function ZZtilde = compute_non_policy_anticipated_shocks_impact_matrix(...
    PHIztilde,ztilde,F,nx,H)
% Contains the logic for impact matrices of non-policy shocks.
%
% INPUTS:   
%   -> ztilde: (nz-1)*H matrix of non-policy shocks
%   -> PHIztilde: loadings on non-policy shocks from model solution
%   -> F: forward loadings on shock loadings from model solution 
%   -> nx: number of model variables
%   -> H: forecast horizon
%
% OUTPUTS:  
%   -> ZZtilde: nx*H matrix of non-policy shock impacts

%% COMPUTE NON-POLICY SHOCK IMPACT MATRIX
% The matrix here is the nx*H matrix of non-policy (anticipated) shock 
% impacts, \widetilde{\mathbb{Z}} in the paper.
ZZtilde = NaN*ones(nx,H);
for t = H:-1:1
    if t == H
        ZZtilde(:,H) = PHIztilde*ztilde(:,H);
    else
        ZZtilde(:,t) = F*ZZtilde(:,t+1)+PHIztilde*ztilde(:,t);
    end
end

end

%% FUNCTION TO COMPUTE ANTICIPATED POLICY SHOCK IMPACT MATRICES
function XIr = compute_anticipated_policy_shock_impact_matrices(...
    PHIepsr,epsrPaths,F,nx,H,K)
% Contains the logic for impact marices of policy shocks in each state.
% Non-zero anticipated policy shocks will occur in the event that the ELB 
% is binding on exit from the LFL policy in any of the K+1 states.
%
% INPUTS:   
%   -> PHIepsr: loadings on policy shock from model solution
%   -> epsrPaths: 1*H*(K+1) matrix of shocks to impose ZLB after renege
%   -> F: forward loadings on shock loadings from model solution 
%   -> nx: number of model variables
%   -> H: forecast horizon
%   -> K: maximum number of periods for LFL policy
%
% OUTPUTS:  
%   -> XIr: nx*H*(K+1) matrix of state-dependent policy shock impacts

%% COMPUTE POLICY SHOCK IMPACT MATRICES
% The matrix here is an nx*H*K+1 matrix of impacts of the policy shocks 
% conditional on reversion to the rule in each of the K+1 states and 
% necessary if the ELB is violated on reversion to the rule (this is \Xi^r 
% in the paper).
XIr = NaN*ones(nx,H,K+1);
for t = H:-1:1
    % Changed by RH as min-max operation appears redundant:
    for i = 1:K+1 %min(max(t,K+1),K+1)
        if t == H
            XIr(:,H,i) = PHIepsr*epsrPaths(:,H,i);
        else
            XIr(:,t,i) = F*XIr(:,t+1,i)+PHIepsr*epsrPaths(:,t,i);
        end
    end
end

end

%% FUNCTION TO CONSTRUCT THE COEFFICIENTS FOR STACK-TIME INVERSION
function JJ = construct_stack_time_coefficients_for_inversion(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    HhatBxtilde,HhatCxtilde,HhatCr,HhatFxtilde,HhatFr,Bxtildextilde,...
    Bxtilder,Brxtilde,Brr,xtilde0,p,K)
% Contains the logic for construction of the stack-time coefficients.
%
% INPUTS:   
%   -> HtildeBxtilde: loadings on lags of non-policy vars in non-policy eqs
%   -> HtildeCxtilde: loadings on time t non-policy vars in non-policy eqs 
%   -> HtildeCr: loadings on time t policy var in non-policy eqs
%   -> HtildeFxtilde:loadings on leads of non-policy vars in non-policy eqs
%   -> HtildeFr: loadings on lead of policy var in non-policy eqs
%   -> HhatBxtilde: loadings on lags of non-policy vars in policy eq
%   -> HhatCxtilde: loadings on time t non-policy vars in policy eq
%   -> HhatCr: loadings on time t policy var in policy eq
%   -> HhatFxtilde: loadings on leads of non-policy vars in policy eq
%   -> HhatFr: loadings on lead of policy var in policy eq
%   -> Bxtildextilde: loadings of non-policy vars on non-policy vars
%   -> Bxtilder: loading of policy var on non-policy vars
%   -> Brxtilde: loadings of non-policy vars on policy var
%   -> Brr: loading of policy var on policy var
%   -> xtilde0: initial condition for non-policy vars
%   -> p: 1*K row vector of renege probabilities
%   -> K: maximum number of periods for LFL policy
%
% OUTPUTS:  
%   -> JJ: (nxtilde*K)*(nxtilde*K) loadings on stacked non-renege paths

%% INITIALISE
% Note that JJ is \mathbb{J} in the paper.
nxtilde = size(xtilde0,1);
JJ = zeros(nxtilde*K,nxtilde*K);

%% COMPUTE HELPER MATRICES IN THE STACK-TIME MATRIX CONSTRUCTION
% These are the row vector \widehat{\mathbb{H}}_{\tilde x}, the constant 
% \widehat{\mathbb{H}}_{r}, and the column vector 
% \widetilde{\mathbb{H}}_{r} from the paper.
HHhatxtilde = HhatCxtilde+HhatFr*Brxtilde+HhatFxtilde*Bxtildextilde;
HHhatr = HhatCr+HhatFr*Brr+HhatFxtilde*Bxtilder;
HHtilder = HtildeCr+HtildeFr*Brr+HtildeFxtilde*Bxtilder;

%% CONSTRUCT THE MATRIX JJ IN A LOOP
% This follows the logic of the paper. The stacked matrix is constructed
% on a block-row basis starting with the loadings on the vector of non-
% policy endogenous variables conditional on no renege in period 1 and 
% finishing with those on the non-renege outcome in period K of the 
% simulation. Note the following: a) there are no loadings on the lagged 
% endogenous variables in the stack-time in period 1 because the initial 
% conditions are predetermined and so appear in the vector of constants; 
% b) similarly, there are no loadings on the leads of the endogenous 
% variables in the stack-time in period K because the policy-maker reverts
% to the rule with certainty in period K+1 and so the path for the 
% endogenous variables is determined by the RE solution; c) that also means
% that the loadings on the contemporaenous endogenous variables in the 
% stack-time are different in period K as compared to any other arbitrary 
% period.
for t = 1:K
    tRowInds = ((t-1)*nxtilde+1:t*nxtilde)';
    tColInds = tRowInds;
    tLagColInds = ((t-2)*nxtilde+1:(t-1)*nxtilde)';
    tLeadColInds = (t*nxtilde+1:(t+1)*nxtilde)';
    if t > 1
       JJ(tRowInds,tLagColInds) = ...
           HtildeBxtilde-p(t)*(1/HHhatr)*HHtilder*HhatBxtilde; 
    end 
    if t < K
        JJ(tRowInds,tColInds) = HtildeCxtilde+...
            p(t)*(HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde-...
            (1/HHhatr)*HHtilder*HHhatxtilde)-...
            (1-p(t))*p(t+1)*(1/HHhatr)*HtildeFr*HhatBxtilde;
    else
        JJ(tRowInds,tColInds) = HtildeCxtilde+...
            HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde-...
            p(K)*(1/HHhatr)*HHtilder*HHhatxtilde;  
    end
    if t < K
        JJ(tRowInds,tLeadColInds) = (1-p(t))*(HtildeFxtilde-...
            p(t+1)*(1/HHhatr)*HtildeFr*HHhatxtilde);
    end
end

end

%% FUNCTION TO CONSTRUCT THE CONSTANTS FOR STACK-TIME INVERSION
function CC = construct_stack_time_constants_for_inversion(...
    HtildeBxtilde,HtildeBr,HtildeCr,HtildeFxtilde,...
    HtildeFr,PSItildeztilde,HhatBxtilde,HhatBr,HhatCr,...
    HhatFxtilde,HhatFr,PSIhatepsr,Bxtilder,Brr,...
    xtilde0,r0,ztilde,epsrPaths,ZZtilde,XIr,Sxtilde,Sr,p,rbar,K)
% Contains the logic for construction of the stack-time constants.
%
% INPUTS:   
%   -> HtildeBxtilde: loadings on lags of non-policy vars in non-policy eqs
%   -> HtildeBr: loadings on lag of policy var in non-policy eqs 
%   -> HtildeCr: loadings on time t policy var in non-policy eqs
%   -> HtildeFxtilde:loadings on leads of non-policy vars in non-policy eqs
%   -> HtildeFr: loadings on lead of policy var in non-policy eqs
%   -> PSItildeztilde: loadings on non-policy shocks in non-policy eqs
%   -> HhatBxtilde: loadings on lags of non-policy vars in policy eq
%   -> HhatBr: loadings on lag of policy var in policy eq
%   -> HhatCr: loadings on time t policy var in policy eq
%   -> HhatFxtilde: loadings on leads of non-policy vars in policy eq
%   -> HhatFr: loadings on lead of policy var in policy eq
%   -> PSIhatepsr: loading on policy shock in policy eq
%   -> Bxtilder: loading of policy var on non-policy vars
%   -> Brr: loading of policy var on policy var
%   -> xtilde0: initial condition for non-policy vars
%   -> r0: initial condition for policy var
%   -> ztilde: (nz-1)*H matrix of non-policy shocks
%   -> epsrPaths: 1*H*(K+1) matrix of shocks to impose ZLB on renege
%   -> ZZtilde: nx*H matrix of non-policy shock impacts
%   -> XIr: nx*H*(K+1) matrix of state-dependent policy shock impacts
%   -> Sxtilde: (nx-1)*nx selector matrix for non-policy variables
%   -> Sr: 1*nx selector matrix for policy variable
%   -> p: 1*K row vector of renege probabilities
%   -> rbar: 1*K row vector of instrument peg values (in model space)
%   -> K: maximum number of periods for LFL policy
%
% OUTPUTS:  
%   -> CC: (nxtilde*K)*1 vector of constants

%% INITIALISE
% Note that CC is \mathbb{C} in the paper.
nxtilde = size(xtilde0,1);
CC = zeros(nxtilde*K,1);

%% COMPUTE HELPER MATRICES IN THE STACK-TIME MATRIX CONSTRUCTION
% These are the constant \widehat{\mathbb{H}}_{r}, and the column vector 
% \widetilde{\mathbb{H}}_{r} from the paper.
HHhatr = HhatCr+HhatFr*Brr+HhatFxtilde*Bxtilder;
HHtilder = HtildeCr+HtildeFr*Brr+HtildeFxtilde*Bxtilder;

%% CONSTRUCT THE MATRIX CC IN A LOOP
% This follows the logic of the paper. The stacked vector is constructed
% on a block-row basis starting with the block row corresponding to period
% 1 and finishing with period K of the simulation. Note the following:
% a) the vector of constants includes the pre-determined initial conditions
% in period 1; b) in period, K, it must also take into account the known 
% reversion to the RE solution in period K+1; c) there is an edge case when
% t==K==1, which must account for both (a) and (b).
for t = 1:K
    tRowInds = ((t-1)*nxtilde+1:t*nxtilde)';
    if t == 1
        if t < K
            CC(tRowInds,1) = PSItildeztilde*ztilde(:,1)+...
                p(1)*(HHtilder*(1/HHhatr)*...
                (HhatFxtilde*Sxtilde+HhatFr*Sr)-...
                (HtildeFxtilde*Sxtilde+HtildeFr*Sr))*...
                (XIr(:,2,1)+ZZtilde(:,2))+...
                (1-p(1))*p(2)*HtildeFr*(1/HHhatr)*...
                (HhatFxtilde*Sxtilde+HhatFr*Sr)*...
                (XIr(:,3,2)+ZZtilde(:,3))-...
                p(1)*HHtilder*(1/HHhatr)*PSIhatepsr*epsrPaths(1,1,1)-...
                (1-p(1))*p(2)*HtildeFr*(1/HHhatr)*PSIhatepsr*...
                epsrPaths(1,2,2)+...
                (p(1)*(1/HHhatr)*HHtilder*HhatBxtilde-HtildeBxtilde)*...
                xtilde0+...
                (p(1)*(1/HHhatr)*HHtilder*HhatBr-HtildeBr)*r0+...
                (1-p(1))*(p(2)*(1/HHhatr)*HtildeFr*HhatBr-HtildeCr)*...
                rbar(1)-...
                (1-p(1))*(1-p(2))*HtildeFr*rbar(2);
        else
            CC(tRowInds,1) = PSItildeztilde*ztilde(:,1)+...
                p(1)*(HHtilder*(1/HHhatr)*...
                (HhatFxtilde*Sxtilde+HhatFr*Sr)-...
                (HtildeFxtilde*Sxtilde+HtildeFr*Sr))*...
                (XIr(:,2,1)+ZZtilde(:,2))-...
                p(1)*HHtilder*(1/HHhatr)*PSIhatepsr*epsrPaths(1,1,1)-...
                (1-p(1))*(HtildeFxtilde*Sxtilde+HtildeFr*Sr)*...
                (XIr(:,2,2)+ZZtilde(:,2))+...
                (p(1)*(1/HHhatr)*HHtilder*HhatBxtilde-HtildeBxtilde)*...
                xtilde0+...
                (p(1)*(1/HHhatr)*HHtilder*HhatBr-HtildeBr)*r0-...
                (1-p(1))*HHtilder*rbar(1);
        end
    elseif t < K
        CC(tRowInds,1) = PSItildeztilde*ztilde(:,t)+...
            p(t)*(HHtilder*(1/HHhatr)*(HhatFxtilde*Sxtilde+HhatFr*Sr)-...
            (HtildeFxtilde*Sxtilde+HtildeFr*Sr))*...
            (XIr(:,t+1,t)+ZZtilde(:,t+1))+...
            (1-p(t))*p(t+1)*HtildeFr*(1/HHhatr)*...
            (HhatFxtilde*Sxtilde+HhatFr*Sr)*...
            (XIr(:,t+2,t+1)+ZZtilde(:,t+2))-...
            p(t)*HHtilder*(1/HHhatr)*PSIhatepsr*epsrPaths(1,t,t)-...
            (1-p(t))*p(t+1)*HtildeFr*(1/HHhatr)*PSIhatepsr*...
            epsrPaths(1,t+1,t+1)+...
            (p(t)*(1/HHhatr)*HHtilder*HhatBr-HtildeBr)*rbar(t-1)+...
            (1-p(t))*(p(t+1)*(1/HHhatr)*HtildeFr*HhatBr-HtildeCr)*...
            rbar(t)-...
            (1-p(t))*(1-p(t+1))*HtildeFr*rbar(t+1);
    else
        CC(tRowInds,1) = PSItildeztilde*ztilde(:,K)+...
            p(K)*(HHtilder*(1/HHhatr)*(HhatFxtilde*Sxtilde+HhatFr*Sr)-...
            (HtildeFxtilde*Sxtilde+HtildeFr*Sr))*...
            (XIr(:,K+1,K)+ZZtilde(:,K+1))-...
            p(K)*HHtilder*(1/HHhatr)*PSIhatepsr*epsrPaths(1,K,K)-...
            (1-p(K))*(HtildeFxtilde*Sxtilde+HtildeFr*Sr)*...
            (XIr(:,K+1,K+1)+ZZtilde(:,K+1))+...
            (p(K)*(1/HHhatr)*HHtilder*HhatBr-HtildeBr)*rbar(K-1)-...
            (1-p(K))*HHtilder*rbar(K);                
    end
end

end

%% FUNCTION TO DEFINE DEFAULT OPTIONS USED IN ITERATION OVER POLICY SHOCKS
function DefaultPolicyShockIterOpts = ...
    get_default_policy_shock_iteration_options()
% These defaults seem OK and should not need to be changed.
% They can be overridden by the user on input to this function.
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> DefaultPolicyShockIterOpts: structure of default options

%% DEFINE DEFAULT OPTIONS
% Note that the tolerance of the ELB being violated is set to the quarterly 
% equivalent of one-tenth of a basis point on the annualised interest rate:
% 0.1*0.01/4.
DefaultPolicyShockIterOpts.bViolationTol = 0.00025;
DefaultPolicyShockIterOpts.maxIter = 100;
DefaultPolicyShockIterOpts.bifurcationFrac = 0.5;
DefaultPolicyShockIterOpts.epsrChangeTol = 0.01;
DefaultPolicyShockIterOpts.displayProgress = false;

end