function plot_grid_ralf

load('/gpfs01/nienborg/group/Katsuhisa/pupil_project/Figures/Figure6_GridSearch/data/ralf.mat')

co = arrayfun(@(x) x.contrast, ralf.colevel);
[sco, sidx] = sort(co);
mat = nan(length(ralf.colevel), size(ralf.colevel(1).mat,2));
for i = 1:length(ralf.colevel)
    mat(i,:) = ralf.colevel(sidx(i)).mat(2,:);
end
close all;
h = figure;
imagesc(11:size(mat,2), sco, mat(:,11:end))
xlabel('time')
ylabel('signal strength')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
colormap(jet)
colorbar
savefig('/gpfs01/nienborg/group/Katsuhisa/pupil_project/Figures/Figure6_GridSearch/raw_figs/ralf_grid.fig')
 