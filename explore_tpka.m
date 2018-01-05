function explore_tpka(tpka)

close all;
h = figure;
for i = 1:length(tpka.evidence)
    for j = 1:2
        switch j
            case 1
                delta = tpka.evidence(i).pkas(2,end) - tpka.evidence(i).pkas(1,end);
            case 2
                try
                    delta = mean(tpka.evidence(i).pkas(2,end-2:end)) - mean(tpka.evidence(i).pkas(1,end-2:end));
                catch
                    delta = tpka.evidence(i).pkas(2,end) - tpka.evidence(i).pkas(1,end);
                end
        end
        % difference in the last bin
        subplot(2,3,1+(j-1)*3)        
        scatter(i, delta, 15, 'filled', 'markerfacecolor', 'k',...
            'markeredgecolor', 'k', 'markerfacealpha',0.4, 'markeredgealpha',0.8)
        hold on;
        subplot(2,3,2+(j-1)*3)
        scatter(i, delta/max(tpka.evidence(i).pka0), ...
            15, 'filled', 'markerfacecolor', 'k','markeredgecolor', 'k', 'markerfacealpha',0.4, 'markeredgealpha',0.8)
        hold on;
        subplot(2,3,3+(j-1)*3)
        scatter(i, delta/mean(tpka.evidence(i).pkas(:)), ...
            15, 'filled', 'markerfacecolor', 'k','markeredgecolor', 'k', 'markerfacealpha',0.4, 'markeredgealpha',0.8)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')    
        hold on;
    end
end
tlab = {'raw', '/ max(pka0)', '/ mean(pka1,2)'};
for i = 1:6
    subplot(2,3,i)
    if i < 4
        title(tlab{i})                
    end
    if i==1
        ylabel('last value')
    elseif i==4
        ylabel('mean last 3')
    end
    xlim([5 length(tpka.evidence)+0.5])    
    plot([5 length(tpka.evidence)+0.5],[0 0],':k')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
% name of the figure
set(h, 'Name', 'co5_time_PKAdiff', 'NumberTitle', 'off')
% save
savefig('Z:\Katsuhisa\pupil_project\Figures\Figure6_GridSearch\raw_figs\co5_time_PKAdiff.fig')
print(gcf,'-dpdf',['Z:\Katsuhisa\pupil_project\Figures\Figure6_GridSearch\raw_figs\co5_time_PKAdiff.pdf'], sprintf('-r%d',200))