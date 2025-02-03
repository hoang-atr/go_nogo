addpath('RW_model/')

load behavior_full_data.mat
models_of_interest = {'base', 'fixed0','fixed1', 'fixed2', 'fixed3'};

for ani = 1:17
    fprintf('Animal = %d\n',ani)
    data = Animal(ani).data;
    RWModel(ani).models = fit_RW_models(models_of_interest,...
        data(:,1), data(:,2), data(:,3));
end
save('behavior_full_model.mat','RWModel','models_of_interest')    
