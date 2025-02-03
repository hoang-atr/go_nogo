function models = fit_RW_models(models_of_interest, s, a, r)

% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alpha_range = [0.001 0.1];
qinit_range = [0 1];
temperature_range = [0.01 0.5];
penalty_range = [0 1];

runs = 50;

Nmodels = numel(models_of_interest);

for m = 1:Nmodels
    currMod = models_of_interest{m};
    
    if strcmp(currMod, 'base')
        paramNames = {'alpha','q-init1','q-init2','temperature','penalty'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1) ...
                       rand(runs, 1)*diff(temperature_range) + temperature_range(1)...
                       rand(runs, 1)*diff(penalty_range) + penalty_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);  qinit_range(2);  temperature_range(2);  penalty_range(2);
           -alpha_range(1); -qinit_range(1); -qinit_range(1); -temperature_range(1); -penalty_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Mbase, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Mbase(models.(currMod).bestParams, s, a, r);
        
    elseif strcmp(currMod, 'optimized1')
        paramNames = {'alpha','q-init1','q-init2'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);  qinit_range(2);
           -alpha_range(1); -qinit_range(1); -qinit_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Moptimized1, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Moptimized1(models.(currMod).bestParams, s, a, r);
   
    elseif strcmp(currMod, 'optimized2')
        paramNames = {'alpha','q-init1'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);
           -alpha_range(1); -qinit_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Moptimized2, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Moptimized2(models.(currMod).bestParams, s, a, r); 
        
    elseif strcmp(currMod, 'optimized3')
        paramNames = {'alpha','q-init1'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);
           -alpha_range(1); -qinit_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Moptimized3, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Moptimized3(models.(currMod).bestParams, s, a, r); 
        
    elseif strcmp(currMod, 'fixed0')
        paramNames = {'alpha','q-init','temperature','penalty'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1) ...
                       rand(runs, 1)*diff(temperature_range) + temperature_range(1)...
                       rand(runs, 1)*diff(penalty_range) + penalty_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);  temperature_range(2);  penalty_range(2);
           -alpha_range(1); -qinit_range(1); -temperature_range(1); -penalty_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Mfix0, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Mfix0(models.(currMod).bestParams, s, a, r);    
        
    elseif strcmp(currMod, 'fixed1')
        paramNames = {'alpha','q-init1','q-init2'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);  qinit_range(2);
           -alpha_range(1); -qinit_range(1); -qinit_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Mfix1, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Mfix1(models.(currMod).bestParams, s, a, r);
        
    elseif strcmp(currMod, 'fixed2')
        paramNames = {'alpha','q-init1'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1) ...
                       rand(runs, 1)*diff(qinit_range) + qinit_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);  qinit_range(2);
           -alpha_range(1); -qinit_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Mfix2, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Mfix2(models.(currMod).bestParams, s, a, r);
        
    elseif strcmp(currMod, 'fixed3')
        paramNames = {'alpha'};
        startValues = [rand(runs, 1)*diff(alpha_range) + alpha_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        NegLL = zeros(runs,1);
        A=[eye(numParam); -eye(numParam)];
        b=[ alpha_range(2);
           -alpha_range(1)];
        parfor rr = 1:runs
            [allParams(rr, :), NegLL(rr, :)] = ...
                fmincon(@lik_Mfix3, startValues(rr, :), A, b, [], [], [], [], [], options, s, a, r);
        end
        [~, bestFit] = min(NegLL);
        models.(currMod).bestParams = allParams(bestFit, :);
        [~, models.(currMod).choiceProb, models.(currMod).Q, models.(currMod).RPE, models.(currMod).R, models.(currMod).V] = ...
            lik_Mfix3(models.(currMod).bestParams, s, a, r);
        
    else 
        error('RW model: Model name not found')
    end
    models.(currMod).paramNames = paramNames;
    models.(currMod).numParam = numParam;
    models.(currMod).LL = -1 * NegLL(bestFit, :);
    models.(currMod).BIC = log(length(a))*numParam  - 2*models.(currMod).LL;    
end

