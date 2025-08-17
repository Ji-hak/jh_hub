function [exAnteTBFGloss,exAnteTBFGloss_cEquiv] = ...
    compute_ex_ante_TBFG_loss(lossesAlongPaths,TBFGpathProbs,...
    discountFactor,outputWeight)

%% COMPUTE EX ANTE LOSS
K = size(lossesAlongPaths,2)-1;
exAnteTBFGloss = 0;
for s = 2:K+1
        exAnteTBFGloss = exAnteTBFGloss+...
            (TBFGpathProbs(1,1,s)/(1-TBFGpathProbs(1,1,1)))*...
            lossesAlongPaths(1,1,s);
end

%% CONVERT TO CONSUMPTION EQUIVALENCE
exAnteTBFGloss_cEquiv = ...
    ((1-discountFactor)/outputWeight)^0.5*exAnteTBFGloss^0.5;

end