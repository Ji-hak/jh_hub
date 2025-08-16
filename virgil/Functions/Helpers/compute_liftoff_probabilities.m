function sumLiftOffDateProbsExogLow = ...
    compute_liftoff_probabilities(LFLpaths,rVarType,rVarIndex,...
    policyLBval,exanteLFLexitProbs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Kplus1 = size(exanteLFLexitProbs,2);

liftOffDateProbsExogLow = zeros(Kplus1,Kplus1);
for t=1:Kplus1
    rPath = LFLpaths.(rVarType)(rVarIndex,:,t);
    tLiftOffDate = find(rPath>policyLBval,1);
    liftOffDateProbsExogLow(t,tLiftOffDate) = exanteLFLexitProbs(t);
end
sumLiftOffDateProbsExogLow = sum(liftOffDateProbsExogLow,1);


end
