function [LH, probSpike, V, mean_predictedSpikes, RPE] = ott_RW_V_base(startValues, spikeCounts, rewards, timeLocked)

% cued experiment

alphaLearn = startValues(1);
slope = startValues(2);
intercept = startValues(3);
Vinit = 0.5;


trials = length(rewards);
V = zeros(trials + 1, 1);
RPE = zeros(trials, 1);

V(1) = Vinit;
% Call learning rule
for t = 1:trials

    RPE(t) = rewards(t) - V(t);
    V(t + 1) = V(t) + alphaLearn*RPE(t);
end
rateParam = exp(slope*V(1:trials) + intercept);

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