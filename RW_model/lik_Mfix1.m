function [NegLL, choiceProb, Q, delta, r, V]= lik_Mfix1(startValues, s, a, r)

alpha = startValues(1);
Qinit1 = startValues(2);
Qinit2 = startValues(3);
temperature = 0.1;
penalty = 0.5;

Q = zeros(2,2);
Q(1,1) = Qinit1;
Q(2,1) = Qinit2;

T = length(a);

choiceProb = zeros(T,1);
delta = zeros(T,1);

V = zeros(T+1,1);
if a(1)==1 
    if s(1)==1, V(1) = Qinit1; else, V(1) = Qinit2; end
end

ix_fa = (s == 2) & (a == 1);
r(ix_fa) = -penalty;

for t = 1:T
    p(1,:) = exp(Q(1,:)/temperature) / sum(exp(Q(1,:)/temperature));
    p(2,:) = exp(Q(2,:)/temperature) / sum(exp(Q(2,:)/temperature));
    
    % compute choice probability for actual choice
    choiceProb(t) = p(s(t),a(t));
    
    % update values
    delta(t) = r(t) - Q(s(t),a(t));
    Q(s(t),a(t)) = Q(s(t),a(t)) + alpha * delta(t);
    
    % just store the state values
    V(t+1) = Q(s(t),a(t));
end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));