function models = fit_RW_spike_models(models_of_interest, RPE, R, Q, spikeCount)

% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

slope_range = [0 20];
intercept_range = [-20 20];

runs = 20;

Nmodels = numel(models_of_interest);

for m = 1:Nmodels
    currMod = models_of_interest{m};
    
    if strcmp(currMod, 'RPE')
        paramNames = {'slope','intercept'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);
           -slope_range(1); -intercept_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_SM_Value, startValues(rr, :), A, b, [], [], [], [], [], options, spikeCount, RPE);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).probSpike] = ...
            lik_SM_Value(models.(currMod).bestParams, spikeCount, RPE);
     
    elseif strcmp(currMod, 'Q')
        paramNames = {'slope','intercept'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);
           -slope_range(1); -intercept_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_SM_Value, startValues(rr, :), A, b, [], [], [], [], [], options, spikeCount, Q);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).probSpike] = ...
            lik_SM_Value(models.(currMod).bestParams, spikeCount, RPE);
        
    elseif strcmp(currMod, 'outcome')
        paramNames = {'slope','intercept'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);
           -slope_range(1); -intercept_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_SM_Value, startValues(rr, :), A, b, [], [], [], [], [], options, spikeCount, R);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).probSpike] = ...
            lik_SM_Value(models.(currMod).bestParams, spikeCount, R);
    elseif strcmp(currMod, 'mean')
        paramNames = {''};
        numParam = 0;
        [NegLL, ms.(currMod).probSpike] = lik_SM_mean([], spikeCount);
        bestFit = 1;
        ms.(currMod).bestParams = [];    
    else 
        error('RW model: Model name not found')
    end
    models.(currMod).paramNames = paramNames;
    models.(currMod).LL = -1 * NegLL(bestFit, :);
    models.(currMod).BIC = log(length(spikeCount))*numParam  - 2*models.(currMod).LL;    
end

