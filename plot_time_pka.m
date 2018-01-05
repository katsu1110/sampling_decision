function plot_time_pka(tpka)

savepath = 'Z:\Katsuhisa\data\sampling_decision\time_pka\';

% yellow and green
col = zeros(2,3);
col(2,:) = [0.9576    0.7285    0.2285];
col(1,:) = [0.1059    0.4706    0.2157];

for i = 1:length(tpka.evidence)
    close all;
    h = figure;
    for c = 1:2
        subplot(1,3,1)
        plot(1:length(tpka.evidence(i).pka0), tpka.evidence(i).pkas(c,:), '-', 'color', col(c,:))
        hold on;
        subplot(1,3,2)
        plot(1:length(tpka.evidence(i).pka0), tpka.evidence(i).pkas(c,:)...
            /max(tpka.evidence(i).pka0), '-', 'color', col(c,:))
        hold on;
        subplot(1,3,3)
        plot(1:length(tpka.evidence(i).pka0), tpka.evidence(i).pkas(c,:)...
            /mean(tpka.evidence(i).pkas(:)), '-', 'color', col(c,:))
        hold on;
    end
    % difference in the last bin
    subplot(1,3,1)
    title(tpka.evidence(i).pkas(2,end) - tpka.evidence(i).pkas(1,end))
    xlim([0.5 length(tpka.evidence(i).pka0)+0.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    subplot(1,3,2)
    title((tpka.evidence(i).pkas(2,end) - tpka.evidence(i).pkas(1,end))...
        /max(tpka.evidence(i).pka0))
    xlim([0.5 length(tpka.evidence(i).pka0)+0.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    subplot(1,3,3)
    title((tpka.evidence(i).pkas(2,end) - tpka.evidence(i).pkas(1,end))...
        /mean(tpka.evidence(i).pkas(:)))
    xlim([0.5 length(tpka.evidence(i).pka0)+0.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    % name of the figure
    set(h, 'Name', ['co5_EvidenceNumber' num2str(i)], 'NumberTitle', 'off')
    % save
    savefig([savepath '\raw\' ['co5_EvidenceNumber' num2str(i)] '.fig'])
    print(gcf,'-dpdf',[savepath '\pdf\' ['co5_EvidenceNumber' num2str(i)] '.pdf'], sprintf('-r%d',200))
end
    