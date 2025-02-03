function [LH, probSpike] = lik_SM_mean(startValues, spikeCounts)

probSpike = poisspdf(spikeCounts, mean(spikeCounts));

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -sum(log(probSpike));
end