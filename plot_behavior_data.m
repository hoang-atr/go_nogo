clear all
close all

raw_folder = 'Fulldata/';
files = dir([raw_folder 'Ald*.mat']);
Ndata = numel(files);

for a = 1:17
    Animal(a) = struct('session',[]);    
    Animal(a).session = [];
end

perf_range = [0 0.6 0.8 1.1];
Nperf = numel(perf_range)-1;

session_count = zeros(17,1);
lick_window = 0.5;

count_trial = zeros(4,1);

for ndata = 1:Ndata
            
    fname = files(ndata).name;
    sp_char = strfind(fname, '_');
    
    training_no = str2double(fname(sp_char(1)+1));
    
    mouse_name = fname(1:sp_char(1)-1);
    mouse_id = mouse_name_to_id(mouse_name);
    
    % load data
    load([raw_folder fname])
    
    Nt = numel(datst.L);
    t = (1:Nt)*(1/datst.rate);
    t_lick = t(datst.L>0);
                
    % trial
    if isfield(datst.task,'trial_def')        
        Ntrial = numel(datst.task.trial_def);               
                
        session_count(mouse_id) = session_count(mouse_id)+1;    
        
        hit_perf =  sum(datst.task.trial_def==1)  / (sum(datst.task.trial_def==1) + sum(datst.task.trial_def==4));
        fa_perf =  sum(datst.task.trial_def==2)  / (sum(datst.task.trial_def==2) + sum(datst.task.trial_def==3));
        
        perf = ( sum(datst.task.trial_def==1) + sum(datst.task.trial_def==3) ) / numel(datst.task.trial_def);
        ss = perf_to_learning_stage(perf, perf_range);
                
        Nt = numel(datst.L);
        t = (1:Nt)*(1/datst.rate);
        t_lick = t(datst.L>0);
        Nlick = numel(t_lick);  
        
        trial_onset = datst.task.trial_onset;
        
        lick_latency = zeros(Ntrial,1);
        lick_rate = zeros(Ntrial,1);
        trial_def = zeros(Ntrial,1);
        
        for trial = 1:Ntrial
            t_trial_onset = t(datst.task.trial_onset(trial));
            if trial<Ntrial
                iti = t(datst.task.trial_onset(trial+1)) - t(datst.task.trial_onset(trial));
            else
                iti = Inf;
            end
            
            td = t_lick - t_trial_onset ;
            td(td<0) = Inf;     
            
            if min(td)<iti
                lick_latency(trial) = min(td);
            else
                lick_latency(trial) = nan;
            end
            
            lick_rate(trial) = sum(t_lick>t_trial_onset & t_lick<=t_trial_onset+lick_window);
            trial_def(trial) = datst.task.trial_def(trial);
        end      
        
        [tl,yl] = spike_time_to_spike_signal({t_lick}, 0.1);
        lick_psth = zeros(4,50);
        for event_type = 1:4
            t_event = get_event_timing(datst, event_type);
            lick_psth(event_type,:)=event_trigger_response(tl,yl,t_event,5,50);
            
            count_trial(event_type) = count_trial(event_type) + numel(t_event);
        end
                
                
        Animal(mouse_id).session(session_count(mouse_id)).hit_perf = hit_perf;
        Animal(mouse_id).session(session_count(mouse_id)).fa_perf = fa_perf;
        Animal(mouse_id).session(session_count(mouse_id)).lick_latency = lick_latency;
        Animal(mouse_id).session(session_count(mouse_id)).early_lick_rate = lick_rate;
        Animal(mouse_id).session(session_count(mouse_id)).trial_def = trial_def;
        Animal(mouse_id).session(session_count(mouse_id)).lick_psth = lick_psth;
        Animal(mouse_id).session(session_count(mouse_id)).correct_fraction = perf;
        Animal(mouse_id).session(session_count(mouse_id)).learning_stage = ss;
        Animal(mouse_id).session(session_count(mouse_id)).trial_onset = trial_onset;
    end       
end


%% plot the results
colors = [[0 0 1]; [0 0.8 0]; [1 0 0]];

figure
for event_type = 1:4
    Y = cell(3,1);
    for a = 1:17
        data = Animal(a);
        Nsession = numel(data.session);
        
        for i = 1:Nsession
            s = data.session(i).learning_stage;
            x = data.session(i).lick_psth(event_type,:);
            Y{s} = [Y{s}; x];
        end
    end
    
    subplot(2,2,event_type)
    hold on
    for s = 1:3
        X = Y{s};
        X(any(isnan(X), 2), :) = [];
        plot(mean(X),'color',colors(s,:),'LineWidth',3)
    end
    ylim([0 0.8])
end

lick_latency = cell(17,4);
lick_freq = cell(17,4);

hit_perf_all = zeros(7,17);
fa_perf_all = zeros(7,17);
lick_rate_all = zeros(7,17);

figure
for a = 1:17
    data = Animal(a);
    Nsession = numel(data.session);
    
    hit_perf = zeros(Nsession,1);
    fa_perf = zeros(Nsession,1);
    
    lick_rate_mean = zeros(Nsession,1);
    lick_rate_sem = zeros(Nsession,1);
    
    for s = 1:Nsession
        hit_perf(s) = data.session(s).hit_perf;
        fa_perf(s) = data.session(s).fa_perf;
        
        trial_onset = data.session(s).trial_onset;
        lick_rate = data.session(s).early_lick_rate;        
        
        ix_hit = data.session(s).trial_def==1;
        ix_fa = data.session(s).trial_def==2;         
        ix_cr = data.session(s).trial_def==3; 
        ix_miss = data.session(s).trial_def==4;
        
        ix_go = (data.session(s).trial_def==1 | data.session(s).trial_def==4); 
        ix_nogo = (data.session(s).trial_def==2 | data.session(s).trial_def==3); 
        
        lick_rate_mean(s) = mean(lick_rate(ix_nogo));
        lick_rate_sem(s) = std(lick_rate(ix_nogo))/sqrt(sum(ix_nogo));        
        
        for i = 1:4
            ix = data.session(s).trial_def==i;            
            lick_freq{a,i} = [lick_freq{a,i}; [lick_rate(ix) zeros(sum(ix),1)+s trial_onset(ix)]];
            latency_session = data.session(s).lick_latency(ix);
            lick_latency{a,i} = [lick_latency{a,i}; [latency_session zeros(sum(ix),1)+s trial_onset(ix)] zeros(sum(ix),1)+median(latency_session,'omitnan')];
        end    
    end
    
    hit_perf_all(1:Nsession,a) = hit_perf;
    fa_perf_all(1:Nsession,a) = fa_perf;
    lick_rate_all(1:Nsession,a) = lick_rate_mean;
    
    subplot(2,2,1)
    hold on
    plot(1:Nsession, hit_perf, 'o-', 'color', [0.7 0.7 0.7], 'LineWidth',2);
    axis([0.5 7.5 0.5 1])
    
    subplot(2,2,2)
    hold on
    plot(1:Nsession, fa_perf, 'o-', 'color', [0.7 0.7 0.7], 'LineWidth',2);
    axis([0.5 7.5 0 1])
    
    subplot(2,2,4)
    hold on
    errorbar(1:Nsession, lick_rate_mean, lick_rate_sem, 'o-', 'color', [0.7 0.7 0.7], 'LineWidth',2)
    axis([0.5 7.5 0 3])
end

[ii,~,v] = find(hit_perf_all);
out = accumarray(ii,v,[],@mean);
out2 = accumarray(ii,v,[],@std);
subplot(2,2,1)
plot(1:7,out,'k','LineWidth',4)

[ii,~,v] = find(fa_perf_all);
out = accumarray(ii,v,[],@mean);
out2 = accumarray(ii,v,[],@std);
subplot(2,2,2)
plot(1:7,out,'k','LineWidth',4)

[ii,~,v] = find(lick_rate_all);
out = accumarray(ii,v,[],@mean);
out2 = accumarray(ii,v,[],@std);
subplot(2,2,4)
plot(1:7,out,'k','LineWidth',4)


% AIC = []; BIC = [];
% for event_type = 1
%     for a = 1:17   
%         y = lick_latency{a,event_type}(:,1);
%         x = (1:numel(y))';
%         tbl = table(x,y);
%         
%         m0 = fitlm(tbl, 'y~1');
%         m1 = fitlm(tbl, 'y~1+x');
%         m2 = fitlm(tbl, 'y~1+x+x^2');
%         m3 = fitlm(tbl, 'y~1+x+x^2+x^3');
%         m4 = fitlm(tbl, 'y~1+x+x^2+x^3+x^4');
%         m5 = fitlm(tbl, 'y~1+x+x^2+x^3+x^4+x^5');
%         
%         aic0 = 2*m0.NumCoefficients-2*m0.LogLikelihood;
%         aic1 = 2*m1.NumCoefficients-2*m1.LogLikelihood;
%         aic2 = 2*m2.NumCoefficients-2*m2.LogLikelihood;
%         aic3 = 2*m3.NumCoefficients-2*m3.LogLikelihood;
%         aic4 = 2*m4.NumCoefficients-2*m4.LogLikelihood;
%         aic5 = 2*m5.NumCoefficients-2*m5.LogLikelihood;
%         AIC = [AIC; [aic0 aic1 aic2 aic3 aic4 aic5]];
%         
%         bic0 = log(m0.NumObservations)*m0.NumCoefficients-2*m0.LogLikelihood;
%         bic1 = log(m1.NumObservations)*m1.NumCoefficients-2*m1.LogLikelihood;
%         bic2 = log(m2.NumObservations)*m2.NumCoefficients-2*m2.LogLikelihood;
%         bic3 = log(m3.NumObservations)*m3.NumCoefficients-2*m3.LogLikelihood;
%         bic4 = log(m4.NumObservations)*m4.NumCoefficients-2*m4.LogLikelihood;
%         bic5 = log(m5.NumObservations)*m5.NumCoefficients-2*m5.LogLikelihood;
%         BIC = [BIC; [bic0 bic1 bic2 bic3 bic4 bic5]];
%     end
% end
% mean(AIC)
% mean(BIC)

% lick residual
for event_type = 1:4
    for a = 1:17   
        y = lick_latency{a,event_type}(:,1);
        x = (1:numel(y))';
        
        p = polyfit(x,y,4);
        yfit = polyval(p,x);
        
        residual = abs(yfit-y)/mean(y,'omitnan'); 
        
%         median_y = lick_latency{a,event_type}(:,4);        
%         residual = abs(median_y-y)/mean(y,'omitnan');
                
        lick_latency{a,event_type} = [lick_latency{a,event_type} yfit residual];       
    end
end

% plot lick residual for HIT
lick_residual_all = zeros(7,17);
for a = 1:17
    lick_residual_mean = zeros(7,1);
    lick_residual_sem = zeros(7,1);
    xx = [];

    for s = 1:7
        ix = lick_latency{a,1}(:,2) == s;
        if sum(ix)==0, continue; end;
        xx = [xx; s];
        lick_residual_mean(s) = mean(lick_latency{a,1}(ix,end));
        lick_residual_sem(s) = std(lick_latency{a,1}(ix,end))/sqrt(sum(ix));
    end
    
    subplot(2,2,3)
    hold on
    errorbar(xx,lick_residual_mean(xx), lick_residual_sem(xx), 'o-', 'color', [0.7 0.7 0.7], 'LineWidth',2)
    axis([0.5 7.5 0 1])
    
    lick_residual_all(:,a) = lick_residual_mean;   
end
        
[ii,~,v] = find(lick_residual_all(:,1:16));
out = accumarray(ii,v,[],@mean);
out2 = accumarray(ii,v,[],@std);
subplot(2,2,3)
plot(1:7,out,'k','LineWidth',4)


% lick latency for HIT
figure
for a = 1:17
    y = lick_latency{a,1}(:,1);
    
    y_pred = lick_latency{a,1}(:,end-1);
    y_pred2 = lick_latency{a,1}(:,end-2);
    x_pred = (1:numel(y_pred))';
    
    subplot(3,6,a)
    plot(y,'ko','LineWidth',1)
    hold on
    plot(x_pred, y_pred,'r','LineWidth',2)
    %plot(x_pred, y_pred2,'g','LineWidth',2)
    axis([0 1100 0 1])
    box off
end

coeff_hit = zeros(17,1);
figure
for a = 1:17
    y = lick_latency{a,1}(:,end);
    x = (1:numel(y))';
    
    [m,x_pred,y_pred] = regression_analysis(x,y);    
    coeff_hit(a) = m.Coefficients{2,2};
    
    subplot(3,6,a)
    plot(y,'ko','LineWidth',1)
    hold on
    plot(x_pred, y_pred,'r','LineWidth',2)
    axis([0 1100 0 1])
    box off
end

% lick freq coeffient for No-go
coeff_fa = zeros(17,1);
figure
for a = 1:17
    lick_fa = lick_freq{a,2};
    lick_cr = lick_freq{a,3};
    
    max_training = max(max(lick_fa(:,2)),max(lick_cr(:,2)));
    y = [];
    for t = 1:max_training
        ix_fa = lick_fa(:,2)==t;
        ix_cr = lick_cr(:,2)==t;
        trial_onset = [lick_fa(ix_fa,3); lick_cr(ix_cr,3)];
        tmp = [lick_fa(ix_fa,1); lick_cr(ix_cr,1)];
        [~,sort_ix] = sort(trial_onset);
        y = [y; tmp(sort_ix)];
    end    
    
    x = (1:numel(y))';    
    
	[m,x_pred,y_pred] = regression_analysis(x,y);    
    coeff_fa(a) = m.Coefficients{2,2};
    
    subplot(3,6,a)
    plot(y,'ko','LineWidth',1)
    hold on
    plot(x_pred, y_pred,'r','LineWidth',2)
    axis([0 1200 0 4.2])
    box off
end

figure
plot(coeff_hit,coeff_fa,'ko','MarkerFaceColor','k','MarkerSize',8)
hold on
[m,x_pred,y_pred] = regression_analysis(coeff_hit,coeff_fa);    
m
%plot(x_pred, y_pred, 'r-', 'LineWidth',4)

corrcoef(coeff_hit,coeff_fa)
save('learning_index.mat','coeff_hit','coeff_fa','lick_freq','lick_latency')