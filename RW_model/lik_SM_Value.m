function [LH, probSpike] = lik_SM_Value(startValues, spikeCounts, Value)

% cued experiment
slope = startValues(1);
intercept = startValues(2);
rateParam = exp(slope*Value + intercept);
probSpike = poisspdf(spikeCounts, rateParam);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -sum(log(probSpike));
end