function PM_Compute(E)
%%  compute psychometric function by logistic regression

% choice
ch = E.O(:,1,end) - 1;

% stimulus
s = length(size(E.Signal));
ixp=1;
switch s
    case 2
        stm = mean(E.Signal(:,ixp), 2);
        x = linspace(-0.1, 0.1, 100);
        xlab = 'signal';
    case 3
        stm = zeros(size(E.Signal, 1), 1);
        for i = 1:size(E.Signal, 1)
            stm(i) = sum(E.Signal(i,ixp,:) > 0)/length(E.Signal(i,ixp,:));
        end
        x = linspace(0, 1, 100);
        xlab = 'ratio of vertical stimulus per trial';
end

% %  debug stimulus ---------------
% figure;
% histogram(stm)
% %  ------------------------------------------

% logistic regression
b = glmfit(stm, ch, 'binomial', 'link', 'logit', 'constant', 'on');

% visualize PM
close all;
figure;
plot(x, glmval(b, x, 'logit'), '-k')
yy = get(gca, 'YLim');
hold on;
plot([0 0], yy, ':k')
hold on;
plot([min(x) max(x)], [0.5 0.5], ':k')
ylim(yy)
xlabel(xlab)
ylabel('P(ch = 2)')
title('psychometric function')
set(gca, 'box','off'); set(gca, 'TickDir', 'out')
