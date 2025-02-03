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

return
str_idx = 10;
end_idx = 15;

X = [];
for t = 1:2
    psth1 = squeeze(zonal_x_b0(:,1,t,str_idx:end_idx));
    psth3 = squeeze(zonal_x_b0(:,3,t,str_idx:end_idx));
    diff = (psth3-psth1);
    cum_psth = sum( diff , 2);
        
    %plot(cum_psth,'-o','LineWidth',4)    
    X = [X, cum_psth];
end

for t = 1:4
    psth = squeeze(zonal_x_b0(:,:,t,str_idx:end_idx));
    cum_psth = sum( squeeze(mean(psth,2)) , 2);
        
    %plot(cum_psth,'-o','LineWidth',4)
    X = [X, cum_psth];
end

tmp1=X(:,1);
tmp2=X(:,2);
X(:,1)=X(:,3);
X(:,2)=tmp1;
X(:,3)=X(:,4);
X(:,4)=tmp2;
X(:,4)=-X(:,4);
figure
bar(X)

% colors = [[0 0 1]; [0 0.8 0]; [1 0 0]];
% for z = [2 5]    
%     for s = 1:Nperf
%         for t = 1:2
%             psth = squeeze(zonal_x(z,s,t,:));
%             
%             subplot(2,2,(fix(z/2)-1)*2+t)
%             hold on
%             plot(tt, psth, 'LineWidth',4,'color',colors(s,:))
%             axis([-.2 2 -.05 0.35])
%         end
%     end
% end

% load select_idx_zone_3trials_top100.mat
% 
% figure
% 
% for c = 1:2
%     for s = 1:Nperf
%         ix2 = select_idx{c,s};
%         
%         for t = 1:2            
%             psth = squeeze(x(ix2,t,:));
%             
%             subplot(2,2,(c-1)*2+t)
%             hold on
%             plot(tt, mean(psth), 'LineWidth',4,'color',colors(s,:))
%             axis([-.2 1 -.05 0.5])
%         end
%     end
% end