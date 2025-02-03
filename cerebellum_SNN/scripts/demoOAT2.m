clear 
close all

% add Offline Analysis Toolbox to path
initOAT;

trial = load ('../Release/trial.txt');
go=find(trial==1);
nogo=find(trial==2);

rslt_folder = '../Release/results/';

cc = {'b','r'};

trial_duration = 500;
Nrun = 10;
Ntrial = 500;

MAX_LICK_RATE = 6;

%% IO couplings
all_coup1 = zeros(Ntrial, Nrun);
all_coup2 = zeros(Ntrial, Nrun);

for run = 1:Nrun    
    filewt = [rslt_folder sprintf('IOcoupling1_r%d.dat',run)];
    CR = ConnectionReader(filewt);
    [~, allWeights] = CR.readWeights();
    weight = mean(allWeights,2,'omitnan');
    all_coup1(:,run) = weight(2:end);
    
    filewt = [rslt_folder sprintf('IOcoupling2_r%d.dat',run)];
    CR = ConnectionReader(filewt);
    [~, allWeights] = CR.readWeights();
    weight = mean(allWeights,2,'omitnan');
    all_coup2(:,run) = weight(2:end);
end
mcoup1 = mean(all_coup1,2);
mcoup2 = mean(all_coup2,2);

figure
subplot(1,2,1)
title('IO coupling 1')
plot(go,mcoup1(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
hold on
plot(nogo,mcoup1(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})
ylim([0 10])

subplot(1,2,2)
title('IO coupling 2')
plot(go,mcoup2(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
hold on
plot(nogo,mcoup2(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})
ylim([0 10])


%% Lick
all_lick = zeros(Ntrial,Nrun);

for run = 1:Nrun
    file = [rslt_folder sprintf('spk_lick_r%d.dat',run)];
    fid=fopen(file);
    lick = fread(fid,'float');
    fclose(fid);
    
    all_lick(:,run) = lick;
end
mlick = mean(all_lick,2);

figure
plot(go,mlick(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
hold on
plot(nogo,mlick(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})

%plot_metric(mlick,trial);

for tc = 1:2
    all_rateIO = zeros(Ntrial, Nrun);
    all_rateCN = zeros(Ntrial, Nrun);
    all_ratePC = zeros(Ntrial, Nrun);

    for run = 1:Nrun
        fileIO = [rslt_folder sprintf('spk_IO%d_r%d.dat',tc,run)];
        SR = SpikeReader(fileIO);
        spkIO = SR.readSpikes(trial_duration);
        rateIO = mean(spkIO,2);
        all_rateIO(:, run) = rateIO;

        filePC = [rslt_folder sprintf('spk_PC%d_r%d.dat',tc,run)];
        SR = SpikeReader(filePC);
        spkPC = SR.readSpikes(trial_duration);    
        ratePC = mean(spkPC,2);
        all_ratePC(:, run) = ratePC;

        fileCN = [rslt_folder sprintf('spk_CN%d_r%d.dat',tc,run)];
        SR = SpikeReader(fileCN);
        spkCN = SR.readSpikes(trial_duration);    
        rateCN = mean(spkCN,2);
        all_rateCN(:, run) = rateCN;
    end

    all_rateIO = all_rateIO * (1000/trial_duration);
    all_ratePC = all_ratePC * (1000/trial_duration);
    all_rateCN = all_rateCN * (1000/trial_duration);

    all_rateIO = mean(all_rateIO,2);
    all_ratePC = mean(all_ratePC,2);
    all_rateCN = mean(all_rateCN,2);
    
    %plot_metric(all_ratePC,trial);
    %plot_metric(all_rateCN,trial);
    plot_metric(all_rateIO,trial);
    
%     figure
%     subplot(1,3,1)
%     plot(go,all_ratePC(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
%     hold on
%     plot(nogo,all_ratePC(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})
%     ylim([0 90])
%         
%     subplot(1,3,2)
%     plot(go,all_rateCN(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
%     hold on
%     plot(nogo,all_rateCN(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})
%     ylim([15 40])
%     
%     subplot(1,3,3)
%     plot(go,all_rateIO(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
%     hold on
%     plot(nogo,all_rateIO(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})
%     ylim([0 6])
end

return

%% pf-PC weight
all_weight11 = zeros(Ntrial, Nrun);
all_weight22 = zeros(Ntrial, Nrun);

for run = 1:Nrun    
    filewt = [rslt_folder sprintf('wt_PfPC1_r%d.dat',run)];
    CR = ConnectionReader(filewt);
    [~, allWeights] = CR.readWeights();
    weight = mean(allWeights,2);
    all_weight11(:,run) = weight(2:end);
    
    filewt = [rslt_folder sprintf('wt_PfPC2_r%d.dat',run)];
    CR = ConnectionReader(filewt);
    [~, allWeights] = CR.readWeights();
    weight = mean(allWeights,2);
    all_weight22(:,run) = weight(2:end);
end

figure
plot(mean(all_weight11,2),'c-','LineWidth',2)
hold on
plot(mean(all_weight22,2),'m-','LineWidth',2)