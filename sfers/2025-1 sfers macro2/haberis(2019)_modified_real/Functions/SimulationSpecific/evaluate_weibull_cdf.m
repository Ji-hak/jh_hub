function p = evaluate_weibull_cdf(x,alfa1,alfa2)
% EVALUATES CDF FOR WEIBULL DISTRIBUTION
p = 1-exp(-(x/alfa1)^alfa2);
end

