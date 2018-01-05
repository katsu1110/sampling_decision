function mat = plot_ralf_co_pka

% pathes
if ismember('gpfs01', cd)==1
    mypath = '/gpfs01/nienborg/group/Katsuhisa/';
else
    mypath = 'Z:/Katsuhisa/';
end
load([mypath 'pupil_project/Figures/Figure6_GridSearch/data/tpka_co.mat'])

% visualization
k = 5:10;
mat = [];
for i = k
    v = tpka_co{i-4}.mat(3,:);
    mat = [mat; v(~isnan(v))];
end
close all;
h = figure;
imagesc(12:size(mat,2), k, mat(:,12:size(mat,2)))
colormap(jet)
colorbar('eastoutside')
xlabel('the number of evidence')
ylabel('contrast')
set(gca,'box','off'); set(gca,'TickDir','out')
set(h, 'Name', 'ralf_contrast_pka', 'NumberTitle', 'off')

% save
savefig([mypath 'pupil_project/Figures/Figure6_GridSearch/raw_figs/ralf_co_pka.fig'])
