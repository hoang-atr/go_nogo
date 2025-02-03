function [LH, probSpike, Q, mean_predictedSpikes, RPE] = ott_RW_V_RPE2(startValues, spikeCounts, s,a,rewards, timeLocked)

% cued experiment
alpha = startValues(1);
slope = startValues(2);
intercept = startValues(3);
Vinit = 0.5;

trials = length(rewards);
Q = zeros(2,2);
Q(1,1) = Vinit;
Q(2,1) = Vinit;

RPE = zeros(trials, 1);

% Call learning rule
for t = 1:trials       
    % update values
    RPE(t) = rewards(t) - Q(s(t),a(t));
    Q(s(t),a(t)) = Q(s(t),a(t)) + alpha * RPE(t);       
end
rateParam = exp(slope*RPE + intercept);
probSpike = poisspdf(spikeCounts(timeLocked), rateParam(timeLocked));

mean_predictedSpikes = rateParam(timeLocked);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end