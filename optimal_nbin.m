function tpka = optimal_nbin(E)
% visualize how the number of time bins change PKA

nbins = [3,4,5,6];
for n = 1:length(nbins)
    [tpka.iter(n)] = time_pka(E, nbins(n), 0, 0, 0);
    clc
end
% close all;
h = figure;
cmap = parula(length(nbins));
for n = 1:length(nbins)
    scatter(tpka.iter(n).mat(2,:), tpka.iter(n).mat(3,:), ...
        80, 'filled', 'markerfacecolor', cmap(n,:), 'markerfacealpha', 0.8)
    text(tpka.iter(n).mat(2,1), tpka.iter(n).mat(3,1)-0.005, [num2str(nbins(n)) ' bins'],...
        'color', cmap(n,:))
    hold on;
end
xx = get(gca, 'XLim');
plot(xx, [0 0], ':k')
xlabel('time (neuronal resplution)')
ylabel('\Delta PKA at the last time bin')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(h, 'Name' ,['SamplingDecision_co' num2str(E.Projection.stimulus_contrast(1))],...
    'NumberTitle', 'off')