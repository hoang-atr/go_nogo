clear 
close all

% add Offline Analysis Toolbox to path
initOAT;

trial = load ('../Release/trial.txt');
go=find(trial==1);
nogo=find(trial==2);

rslt_folder = '../Release/results/';

Nrun = 10;

trial_duration = 500;
Ntrial = 500;
Nneuron = 100;

bin_size = 10;

cc = {'b','r'};

for tc = 1:2 
    run = 2;
    fileIO = [rslt_folder sprintf('spk_IO%d_r%d.dat',tc,run)];
    SR = SpikeReader(fileIO);
	spkIO = SR.readSpikes(-1);
    
    subplot(2,1,tc)
    hold on
    for i = 1:Nneuron
        ix = spkIO(2,:)==(i-1);
        spk = spkIO(1,ix);
        spk_trial = parse_event_timing(spk,Ntrial,trial_duration);
        
        off = 0;
        for n = 1:10
            if trial(n)==1
                plot(spk_trial{n}+off,zeros(size(spk_trial{n}))+i,'b.')
            else
                plot(spk_trial{n}+off,zeros(size(spk_trial{n}))+i,'r.')
            end
            off = off + trial_duration;
        end
    end        
end
return

%% synchrony traces
for tc = 1:2 
    synchrony = zeros(Nrun, Ntrial);
    for run = 1:Nrun
        fileIO = [rslt_folder sprintf('spk_IO%d_r%d.dat',tc,run)];
        SR = SpikeReader(fileIO);
        spkIO = SR.readSpikes(-1);
            
        spk_neuron = cell(Nneuron,Ntrial);
        for i = 1:Nneuron
            ix = spkIO(2,:)==(i-1);
            spk = spkIO(1,ix);
            spk_neuron(i,:) = parse_event_timing(spk,Ntrial,trial_duration);
        end
        
        synchrony(run,:) = extract_synchrony(spk_neuron, bin_size, trial_duration);        
    end
        
    msyn = mean(synchrony,1);
    
    subplot(1,2,tc)   
    plot(go,msyn(go),'o','MarkerSize',8,'MarkerFaceColor',cc{1},'MarkerEdgeColor',cc{1})
    hold on
    plot(nogo,msyn(nogo),'o','MarkerSize',8,'MarkerFaceColor',cc{2},'MarkerEdgeColor',cc{2})
    ylim([0 0.08])
end