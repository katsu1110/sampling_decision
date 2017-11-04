function Kernel_Compute(E)

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

% confidence
conf = abs(conf_proxy(:,end));
med = median(conf);

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

% split
El = E;
El.Signal = El.Signal(conf < med, :, :);
El.O = El.O(conf < med, :, :);
Eh = E;
Eh.Signal = Eh.Signal(conf > med, :, :);
Eh.O = Eh.O(conf > med, :, :);

% PK
[pk] = getPK(E, n0S);
[pkh] = getPK(Eh, n0S);
[pkl] = getPK(El, n0S);
figure;
plot(pk, '-r')
xlabel('time')
ylabel('PK')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

figure;
plot(pkh, '-','color',y)
hold on;
plot(pkl, '-','color',g)
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

