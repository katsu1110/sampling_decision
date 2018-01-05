function just_get_more_tr

if mean(ismember('gpfs0',cd))==1
    savepath = '/gpfs01/nienborg/group/Katsuhisa/data/sampling_decision/output/';
else
    savepath = 'Z:\Katsuhisa\data\sampling_decision\output\';
end
% out_kk_co10 = S_Experiment(S_Exp_Para('PK-amplitude10'));
% save([savepath 'out_kk_co10.mat'],'out_kk_co10','-v7.3')
% out_kk_co9 = S_Experiment(S_Exp_Para('PK-amplitude9'));
% out_kk_co9 = rmFieldforMemory(out_kk_co9);
% save([savepath 'out_kk_co9.mat'],'out_kk_co9','-v7.3')
% out_kk_co8 = S_Experiment(S_Exp_Para('PK-amplitude8'));
% out_kk_co8 = rmFieldforMemory(out_kk_co8);
% save([savepath 'out_kk_co8.mat'],'out_kk_co8','-v7.3')
out_kk_co7 = S_Experiment(S_Exp_Para('PK-amplitude7'));
out_kk_co7 = rmFieldforMemory(out_kk_co7);
save([savepath 'out_kk_co7.mat'],'out_kk_co7','-v7.3')
% out_kk_co6 = S_Experiment(S_Exp_Para('PK-amplitude6'));
% save([savepath 'out_kk_co6.mat'],'out_kk_co6','-v7.3')
% out_kk_co5 = S_Experiment(S_Exp_Para('PK-amplitude5'));
% save([savepath 'out_kk_co5.mat'],'out_kk_co5','-v7.3')




