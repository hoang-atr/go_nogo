function synchrony = extract_synchrony(spk_neuron, bin_size, trial_duration)
[Nneuron, Ntrial] = size(spk_neuron);

edge=0:bin_size:trial_duration;
NT = numel(edge)-1;

synchrony = zeros(1,Ntrial);
for n = 1:Ntrial    
    count = zeros(1,NT);
    for i = 1:Nneuron
        spk = spk_neuron{i,n};
        count = count + histcounts(spk, edge);
    end
    synchrony(n) = sum(count(count>1));
end
synchrony = synchrony ./ (Nneuron*(Nneuron-1)/2);