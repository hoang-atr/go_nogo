clear all
close all

load tensor_data_3trials.mat
load optimal_TCA.mat
Ncom = 4;
M = optimalM{Ncom};

Nstage = 3;

colors = [[1 0 0];[0 0 1];[1 0 0];[0 0 1];[1 0 0];[0 0 1];[1 0 0];[0 0 1]];

tt = (-9:40)*0.05;

figure
for c = 1:Ncom
    weight = M.u{1}(:,c);
    cue_condition = M.u{2}(:,c);
    temporal_profile = M.u{3}(:,c);
    
    MEAN = zeros(Nstage,8);
    SD = zeros(Nstage,8);
    for s = 1:Nstage
        for z = 1:8
            ix = stage == s & zone == z;
            ww = weight(ix);
            
            MEAN(s,z) = mean(ww);
            SD(s,z) = std(ww)/sqrt(sum(ix));
        end                
    end    
    
    subplot(Ncom,3,(c-1)*3+1)
    barwitherr(SD,MEAN)
    ylim([0 0.02])
%     colormap(hot)
%     imagesc(MEAN')
    
    
    subplot(Ncom,3,(c-1)*3+2)
    bar(cue_condition)
    ylim([0 1])
    
    subplot(Ncom,3,(c-1)*3+3)
    plot(tt,temporal_profile)
    axis([-0.2 1 0 0.7])
end