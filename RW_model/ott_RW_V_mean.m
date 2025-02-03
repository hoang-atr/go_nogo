function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_mean(startValues, spikeCounts, rewards, timeLocked)

probSpike = poisspdf(spikeCounts, mean(spikeCounts)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = mean(spikeCounts);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;