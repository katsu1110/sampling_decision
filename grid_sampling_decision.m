function [ralf] = grid_sampling_decision

mypath = '/gpfs01/nienborg/group/Katsuhisa/data/sampling_decision/output/';
savepath = '/gpfs01/nienborg/group/Katsuhisa/pupil_project/Figures/Figure6_GridSearch/raw_figs/';
myfiles = {'out_kk12','out_kk13','out_kk14','out_kk15',...
    'out_kk_co9','out_kk_co10'};

for i = 1:length(myfiles)
    E = load([mypath myfiles{i} '.mat']);
    [ralf.colevel(i)] = time_pka(E.(myfiles{i}));
    savefig([savepath 'ralf_co' num2str(ralf.colevel(i).contrast)])
end
save('/gpfs01/nienborg/group/Katsuhisa/pupil_project/Figures/Figure6_GridSearch/data/ralf.mat','ralf')
    
    