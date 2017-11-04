function PM_Compute(E)
%%  compute psychometric function by logistic regression

% choice
ch = E.O(:,1,end) - 1;

% start of sampling
n0S=E.InputImage.n_zero_signal; 

% stimulus
s = length(size(E.Signal));
switch s
    case 2
        stm = E.Signal*size(E.Projection.G.G,1)./E.Projection.G.G;
%         stm = mean(E.Signal(:,1), 2);
        x = linspace(-0.1, 0.1, 100);
        xlab = 'signal';
    case 3
        stmmat = zeros(size(E.Signal, 1), size(E.Signal,3));
        for i = 1:size(E.Signal, 1)
            stmmat(i,:) = E.Signal(i,1,:)*size(E.Projection.G.G,1)./E.Projection.G.G;
%             stmmat(i,:) = squeeze(E.Signal(i,ixp,:));
        end
        stmmat = stmmat(:, n0S+1:end);
        stm = sum(stmmat, 2);
        x = linspace(min(stm),max(stm), 100);
%         xlab = 'ratio of vertical stimulus per trial';
        xlab = 'signal';
        
        %  debug stimulus ---------------
        close all;
        figure;
        subplot(1,3,1)
        imagesc(stmmat)
        colorbar('northoutside')
        xlabel('time')
        ylabel('trials')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        subplot(1,3,2)
        histogram(stm)
        xlabel(xlab)
        ylabel('trials')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        subplot(1,3,3)
        plot(stmmat(randi([1 size(E.Signal,1)], 20, 1),:)')
        xlabel('time')
        ylabel('stimulus')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        %  ------------------------------------------

end

% logistic regression
b = glmfit(stm, ch, 'binomial', 'link', 'logit', 'constant', 'on');

% visualize PM
figure;
plot(x, glmval(b, x, 'logit'), '-k')
yy = get(gca, 'YLim');
hold on;
plot([0 0], yy, ':k')
hold on;
plot([min(x) max(x)], [0.5 0.5], ':k')
xlim([min(x) max(x)])
ylim(yy)
xlabel(xlab)
ylabel('P(ch = 2)')
title('psychometric function')
set(gca, 'box','off'); set(gca, 'TickDir', 'out')
