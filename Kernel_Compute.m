function Kernel_Compute(E, nsplit, varargin)

if nargin==1
    nsplit = 2;
end

close all;

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

% start of sampling
n0S=E.InputImage.n_zero_signal; 

% % log-odds
% conf_proxy = squeeze(diff(log(E.O(1:size(E.O,1),2:3,:)),[],2));

% posterior
pos = squeeze(E.O(1:size(E.O,1),2,:));
conf_proxy = abs(pos - 0.5) + 0.5;

%%
% debug posterior
figure;
subplot(1,3,1)
imagesc(pos)
colorbar
xlabel('time')
ylabel('trials')
title('posterior')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(1,3,2)
rng(1220);
plot(pos(randi([1, size(E.Signal, 1)], 20, 1),:)')
xlabel('time')
ylabel('trials')
title('posterior')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(1,3,3)
histogram(pos)
xlabel('posterior at the end of a trial')
ylabel('trials')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

%%
% confidence
conf = abs(conf_proxy(:,end));
% med = median(conf);

% % putative reaction time
% rt = zeros(1, s(1));
% for i = 1:s(1)
%     pmax=max(E.O(i,2:3,n0S+1:s(4)),[],2);
%     aux=1+find(pmax<thresh,1,'last');
%     if isempty(aux)
%         rt(i)=n0S;
%     else
%         rt(i)=aux+n0S;
%     end
% end

%%
% debug PK
figure;
pk = getPK(E, n0S);
plot(pk/mean(pk), '-r')
hold on;
pk = getPKbyLogReg(E, n0S);
plot(pk/mean(pk),'--r')
xlabel('time')
ylabel('PK')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

%%
% split PK by confidence level
pivot = min(conf)*ones(1, nsplit);
col = jet(nsplit);
PK_normal = nan(nsplit, length(pk));
PK_logreg = nan(nsplit, length(pk));
for n = 2:nsplit+1
    pivot(n) = percentile(conf, (n-1)*round(100/nsplit));
    E_temp  = E;
    E_temp.Signal = E_temp.Signal(conf >= pivot(n-1) & conf < pivot(n), :, :);
    E_temp.O = E_temp.O(conf >= pivot(n-1) & conf < pivot(n), :, :);
    PK_normal(n-1,:) = getPK(E_temp, n0S);
    PK_logreg(n-1,:) = getPKbyLogReg(E_temp, n0S);
end
figure;
for n = 1:nsplit
    plot(PK_normal(n,:)/mean(PK_normal(:)), '-','color',col(n,:),'linewidth',2)
    hold on;
    plot(PK_logreg(n,:)/mean(PK_logreg(:)), '--','color',col(n,:),'linewidth',2)
    hold on;
end
xlabel('time')
ylabel('PK')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')


%%
function [pk] = getPK(E, n0S)
% the number of V1 neurons
nX =  size(E.X, 2);
% O_pref=1;
ixp=1; ixa=1+nX/2;
idx_pref=(E.O(:,2,end)>0.5);
idx_anti=(E.O(:,3,end)>0.5);
prefpref=mean(E.Signal(idx_pref,ixp,:));
prefanti=mean(E.Signal(idx_pref,ixa,:));
antipref=mean(E.Signal(idx_anti,ixp,:));
antianti=mean(E.Signal(idx_anti,ixa,:));
pk=prefpref-prefanti-antipref+antianti;
pk = squeeze(pk);
pk = pk(n0S+1:end);

function [pk] = getPKbyLogReg(E, n0S)
% the number of V1 neurons
nX =  size(E.X, 2);
% choice
ch = E.O(:,1,end) - 1;
% O_pref=1;
ixp=1; ixa=1+nX/2;
stmmat1 = zeros(size(E.Signal,1),size(E.Signal,3));
stmmat2 = zeros(size(E.Signal,1),size(E.Signal,3));
for n = 1:size(E.Signal,1)
    stmmat1(n,:) = squeeze(E.Signal(n,ixp,:));
    stmmat2(n,:) = squeeze(E.Signal(n,ixa,:));
end
% logistic regression
pk = zeros(1, size(E.Signal, 3));
for n = 1:size(E.Signal,3)
    b1 = glmfit(stmmat1(:,n),ch,'binomial','link','logit','constant','on');
    b2 = glmfit(stmmat2(:,n),ch,'binomial','link','logit','constant','on');
    pk(n) = (abs(b1(2))+abs(b2(2)))/2;
end
pk = pk(n0S+1:end);
