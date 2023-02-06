clear all
close all

addpath('tensor_toolbox-v3.2.1/')

hatime_folder = 'ml2ha_data/';
raw_folder = 'raw_data/';
files = dir([raw_folder 'Ald*.mat']);
Ndata = numel(files);

load tensor_data_3trials.mat

%load nn_tensor_3trials_select_idx_top300.mat
load nn_tensor_3trials_select_idx_top20p.mat
Ncls = 4;

% load select_idx_zone_3trials_top100.mat
% Ncls = 2;

Nperf = 3;
colors = [[0 0 1]; [0 0.8 0]; [1 0 0]];

for c = 1:Ncls
    tt = (-9:40)*0.05;
    
    for s = 1:Nperf
%         ix = select_idx{c,s};
        
        ix_stage = find(stage==s);        
        ix = intersect(select_idx{c}, ix_stage);        
        
        nstage = numel(ix);
        
        select_psth = x_b0(ix,:,:);                
        
        for trial = 1:4
            trial_psth = squeeze(select_psth(:,trial,:));
                        
            subplot(Ncls,4,(c-1)*4+trial)
            plot(tt,mean(trial_psth),'LineWidth',2,'color',colors(s,:))                        
            hold on            
            axis([-.5 4 -.05 0.5])
        end
        
    end   
end