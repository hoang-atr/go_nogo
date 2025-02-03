clear all
close all

load tensor_data_3trials.mat
Nperf = 3;

tt = (-9:1:40)*0.05;

for s = 1:Nperf
    for t = 1:4
        psth = squeeze(zonal_x_b0(:,s,t,:));

        figure
        colormap hot
        imagesc(tt,1:8,psth)
        xlim([-.5 2])
        caxis([0 0.25])
        box off
        set(gca,'visible','off')
    end
end
