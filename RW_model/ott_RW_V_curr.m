function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_curr(startValues, spikeCounts, rewards, timeLocked)

% cued experiment
slope = startValues(1);
intercept = startValues(2);
alphaLearn = 0; % only current reward can be encoded
Vinit = 0.5;

trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1); % RPE for trial-by-trial learning

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    
    V(t + 1) = V(t) + alphaLearn*RPE(t);
end

rateParam = exp(slope*rewards + intercept);

probSpike = poisspdf(spikeCounts(timeLocked), rateParam(timeLocked)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(timeLocked);
V = V(1:trials);
V = V(timeLocked);
RPE = RPE(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end