function plot_metric(metric, trial)
stage = {[1:100], [400:500]};
Nstage = numel(stage);

mm = zeros(2,Nstage);
sd = zeros(2,Nstage);

p = zeros(2,1);
for i = 1:2    
    xx = cell(Nstage,1);
    for s = 1:Nstage
        ts = trial(stage{s});
        ms = metric(stage{s});
        
        mm(i,s)=mean(ms(ts==i));
        sd(i,s)=std(ms(ts==i));                
        
        xx{s} = ms(ts==i);
    end  
    [~,p(i)] = ttest2(xx{1},xx{2});
end


for i = 1:2
    fprintf('first = %.2f(%.2f), last = %.2f(%.2f), p = %.5f\n',mm(i,1),sd(i,1),mm(i,2),sd(i,2),p(i))
end




figure
barwitherr(sd,mm)
