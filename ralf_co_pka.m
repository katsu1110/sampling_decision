function [mat] = ralf_co_pka
%%
% generate a heatmap where x is time, y is co, and z is delta PKA for low
% and high confidence trials
%
% +++++++++++++++++++++++++++++++++++++++++++++++

% pathes
if ismember('gpfs01', cd)==1
    mypath = '/gpfs01/nienborg/group/Katsuhisa/';
else
    mypath = 'Z:/Katsuhisa/';
end
datapath = 'data/sampling_decision/output/';

% load data
k = 5:10;
mat = [];
% tpka_co = cell(1,length(k));
for i = k
    E = load([mypath datapath 'out_kk_co' num2str(i) '.mat']);
%     [tpka_co{i-4}] = time_pka(E.(['out_kk_co' num2str(i)]), [], 0, 0, 1);
%     clc
    disp(['co: ' num2str(E.(['out_kk_co' num2str(i)]).Projection.stimulus_contrast(1))])
    disp(['ntr: ' num2str(size(E.(['out_kk_co' num2str(i)]).Signal,1))])
end

% % save
% save([mypath 'pupil_project/Figures/Figure6_GridSearch/data/tpka_co.mat'],'tpka_co','-v7.3')
% clc
% disp('data saved!')
