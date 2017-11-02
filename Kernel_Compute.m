function Kernel_Compute(E)

close all;

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

% basic info
s = size(E.G);
n0S = E.InputImage.n_zero_signal;

% log-odds
odds = squeeze(diff(log(E.O(1:size(E.O,1),2:3,:)),[],2));
figure;
imagesc(odds)
colorbar
xlabel('time')
ylabel('trials')
title('log odds')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

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

% % confidence
% conf = abs(odds(:,end));
% med = median(conf);

% % split
% El = E;
% El.Signal = El.Signal(conf < med, :, :);
% El.O = El.O(conf < med, :, :);
% Eh = E;
% Eh.Signal = Eh.Signal(conf > med, :, :);
% Eh.O = Eh.O(conf > med, :, :);

% PK
[pk] = getPK(E);
% [pkh] = getPK(Eh);
% [pkl] = getPK(El);
figure;
plot(pk, '-r')
xlabel('time')
ylabel('PK')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% hold on;
% plot(pkh, '-','color',y)
% hold on;
% plot(pkl, '-','color',g)

% % PK amp
% [pk] = getPKamp(E);
% [pkh] = getPKamp(Eh);
% [pkl] = getPKamp(El);
% figure;
% plot(pk, '-r')
% hold on;
% plot(pkh, '-','color',y)
% hold on;
% plot(pkl, '-','color',g)
% 
% % binned
% figure;
% nbin = 4;
% [pkamp] = pkbin(pk,nbin);
% plot(pkamp, '-r')
% hold on;
% [pkamph] = pkbin(pkh,nbin);
% plot(pkamph, '-','color',y)
% hold on;
% [pkampl] = pkbin(pkl,nbin);
% plot(pkampl, '-','color',g)

%%
function [pk] = getPK(E)
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

% function [pk] = getPKamp(E)
% % the number of V1 neurons
% nX = size(E.X, 2);
% % O_pref=1;
% ixp=1; ixa=1+nX/2;
% ch = E.O(:,1,end) - 1;
% post1 = E.O(:,2,end);
% post2 = E.O(:,3,end);
% idx_pref=find(post1>0.5);
% idx_anti=find(post2>0.5);
% X1p = nan(length(idx_pref), size(E.X, 3));
% X2n = nan(length(idx_pref), size(E.X, 3));
% for i = 1:length(idx_pref)
%     X1p(i,:) = squeeze(E.Signal(idx_pref(i),ixp,:))';
%     X2n(i,:) = squeeze(E.Signal(idx_pref(i),ixa,:))';
% end
% X1n = nan(length(idx_anti), size(E.X, 3));
% X2p = nan(length(idx_anti), size(E.X, 3));
% for i = 1:length(idx_anti)
%     X1n(i,:) = squeeze(E.Signal(idx_anti(i),ixp,:))';
%     X2p(i,:) = squeeze(E.Signal(idx_anti(i),ixa,:))';
% end
% X = [X1p; X1n];
% y = [post1(idx_pref); post1(idx_anti)];
% % X = [X1p; X2n; X1n; X2p];
% % y = [post1(idx_pref); post2(idx_pref); ...
% %     post1(idx_anti); post2(idx_anti)];
% pk = zeros(1, size(E.X, 3));
% for i = 1:size(E.X, 3)
%     b = glmfit(X(:,i), y, 'normal', 'link', 'identity', 'constant', 'on');
%     pk(i) = b(2);
% end

% function binned = pkbin(pkamp, nbin)
% lenv = length(pkamp);
% frameperbin = floor(lenv/nbin);
% binned = zeros(1, nbin);
% begin = 1;
% for n = 1:nbin
%     binned(n) = mean(pkamp(begin:begin+frameperbin-1));
%     begin = begin + frameperbin;
% end
