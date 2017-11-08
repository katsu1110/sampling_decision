function Kernel_Compute(E, varargin)

nsplit = 2;
cuttime = size(E.O, 3);
resampling_flag = 0;
nbin = cuttime -2;
discretize_flag = 0;
save_flag = 0;

close all;
j = 1;              
while  j <= length(varargin)
    switch varargin{j}
        case 'nsplit'
            nsplit = varargin{j+1};
            j = j + 2;
        case 'cuttime' 
            cuttime = varargin{j+1};
            j = j + 2;
        case 'resample'
            resampling_flag = 1;
            j = j + 1;
        case 'nbin'
            nbin = varargin{j+1};
            j = j + 2;
        case 'discretize'
            discretize_flag = 1;
            j = j + 1;
        case 'save'
            save_flag = 1;
            j = j + 1;
    end
end

if discretize_flag==1
    E = discretize_signal(E);
end

% start of sampling
n0S=E.InputImage.n_zero_signal; 

% % log-odds
% conf_proxy = squeeze(diff(log(E.O(1:size(E.O,1),2:3,:)),[],2));

% posterior
pos = squeeze(E.O(1:size(E.O,1),2,:));
conf = abs(pos(:,cuttime) - 0.5) + 0.5;

% % add noise
% conf = conf + normrnd(median(conf), 0.1*median(conf), size(conf));

%%
% confidence
% conf = abs(conf_proxy);
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
% split PK by confidence level
pivot = min(conf)*ones(1, nsplit);
if nsplit==2
    col = zeros(2,3);
    % yellow and green
    col(2,:) = [0.9576    0.7285    0.2285];
    col(1,:) = [0.1059    0.4706    0.2157];
else
    col = jet(nsplit);
end

% PK
PK_normal = nan(nsplit, E.Projection.n_frames - n0S -1);
PK_logreg = nan(nsplit, E.Projection.n_frames - n0S -1);
pkt = cell(1, nsplit);
serr_pk = zeros(nsplit, E.Projection.n_frames - n0S -1);
serr_logreg = zeros(nsplit, E.Projection.n_frames - n0S -1);
for n = 2:nsplit+1
    pivot(n) = prctile(conf, (n-1)*round(100/nsplit));
    E_temp  = E;
    E_temp.Signal = E_temp.Signal(conf >= pivot(n-1) & conf < pivot(n), :, :);
    E_temp.O = E_temp.O(conf >= pivot(n-1) & conf < pivot(n), :, :);
    PK_normal(n-1,:) = getPK(E_temp, n0S);
    PK_logreg(n-1,:) = getPKbyLogReg(E_temp, n0S);
    if discretize_flag == 1
        [pkt{n-1}] = getDKernel(E_temp);
    end
    if resampling_flag==1
        [serr_pk(n-1,:), serr_logreg(n-1,:)] = resamplePK(E_temp, 1000);
    end
end

%%
% debug posterior
h = figure;
subplot(2,4,1)
imagesc(pos)
colorbar
xlabel('time')
ylabel('trials')
title('posterior')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,4,2)
rng(1220);
plot(pos(randi([1, size(E.Signal, 1)], 100, 1),:)')
hold on;
plot(cuttime*[1 1], [0 1], '--k')
xlabel('time')
ylabel('trials')
title('posterior')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,4,3)
histogram(pos(:,cuttime))
xlabel(['posterior at time: ' num2str(cuttime)])
ylabel('trials')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,4,4)
histogram(conf)
hold on;
yy = get(gca, 'YLim');
plot(median(conf)*[1 1], yy, '-r')
ylim(yy)
xlabel('confidence')
ylabel('trials')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

%%
% debug PK
subplot(2,4,[5 6])
pk0 = getPK(E, n0S);
pk1 = getPKbyLogReg(E, n0S);
stime = 1:length(pk0);
plot([0.5 length(pk0)+0.5],[0 0], ':k')
hold on;
if resampling_flag==1
    hold on
     [err_pk] = resamplePK(E, 100);
     fill_between(stime,(pk0' - err_pk)/mean(pk0), (pk0' + err_pk)/mean(pk0), [1 0 0]);
     hold on;
%      fill_between(stime,(pk1 - err_logreg)/mean(pk1), (pk1 + err_logreg)/mean(pk1), [1 0 0]);
%      hold on;
end
plot(stime, pk0/mean(pk0), '-r')
hold on;
plot(stime, pk1/mean(pk1),'--r')
xlim([0.5 length(pk1)+0.5])
xlabel('time')
ylabel('PK')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

subplot(2,4,[7 8])
plot([0.5 length(pk0)+0.5],[0 0], ':k')
hold on;
for n = 1:nsplit
    if resampling_flag==1
        fill_between(stime, (PK_normal(n,:) - serr_pk(n,:))/mean(PK_normal(:)), ...
            (PK_normal(n,:) + serr_pk(n,:))/mean(PK_normal(:)), col(n,:))
        hold on;
%         fill_between(stime, (PK_logreg(n,:) - serr_logreg(n,:))/mean(PK_logreg(:)), ...
%             (PK_logreg(n,:) + serr_logreg(n,:))/mean(PK_logreg(:)), col(n,:))
%         hold on;
    end
    plot(stime, PK_normal(n,:)/mean(PK_normal(:)), '-','color',col(n,:),'linewidth',1)
    hold on;
    plot(stime, PK_logreg(n,:)/mean(PK_logreg(:)), '--','color',col(n,:),'linewidth',1)
    hold on;
end
yy = get(gca, 'YLim');
plot((E.Sampling.access(cuttime)-n0S-1)*[1 1], yy, '--k')
xlim([0.5 length(pk0)+0.5])
ylim(yy)
xlabel('time')
ylabel('PK')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% figure name
figname =  ['co' num2str(E.Projection.stimulus_contrast(1)) '_alpha' ...
    num2str(E.Projection.alpha) '_kappa' num2str(E.Projection.kappa_O(1)) ...
    '_cut_time' num2str(cuttime) '_nsplit' num2str(nsplit)];
set(h, 'Name', figname, 'NumberTitle','off')
if save_flag==1
    savedir = 'Z:\Katsuhisa\data\sampling_decision\figures\';
    print(h,'-dpdf',[savedir figname],sprintf('-r%d',180))
end

if discretize_flag==1
    figure;
    clim = nan(nsplit, 2);
    binval = [-pi/2, -3*pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2, 7*pi/8, pi];
    for n = 1:nsplit
        subplot(1,nsplit+1,n)
        imagesc(1:size(pkt{nsplit-n+1},2),binval, pkt{nsplit-n+1})
        clim(n,:) = caxis;
    end
    crange = [min(clim(:)) max(clim(:))]
    for n = 1:nsplit
        subplot(1,nsplit+1,n)
        caxis(crange);
        set(gca, 'YTick',[-0.4 2],'YTickLabel',{'0','\pi/2'})
    end

    subplot(1,nsplit+1,nsplit+1)
    plot([0.5 length(pk0)+0.5],[0 0], ':k')
    hold on;
    for n = 1:nsplit
        if resampling_flag==1
            fill_between(stime, (PK_normal(n,:) - serr_pk(n,:))/mean(PK_normal(:)), ...
                (PK_normal(n,:) + serr_pk(n,:))/mean(PK_normal(:)), col(n,:))
            hold on;
        end
        plot(stime, PK_normal(n,:)/mean(PK_normal(:)), '-','color',col(n,:),'linewidth',1)
        hold on;
    end
    yy = get(gca, 'YLim');
    xlim([0.5 length(pk0)+0.5])
    ylim(yy)
    xlabel('time')
    ylabel('PK')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

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
pk = pk(n0S+2:end);

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
b1 = glmfit(stmmat1,ch,'binomial','link','logit','constant','on');
b2 = glmfit(stmmat2,ch,'binomial','link','logit','constant','on');
pk = b1(2:end) - b2(2:end);
pk = pk(n0S+2:end)';

function [err_pk, err_logreg] = resamplePK(E, repeat)
n0S = E.InputImage.n_zero_signal; 
pkrep = nan(repeat, E.Projection.n_frames-n0S-1);
logregrep = nan(repeat, E.Projection.n_frames-n0S-1);
if size(E.Signal,1) > 10000
    sub = 10000;
else
    sub = size(E.Signal,1);
end
for r = 1:repeat
    tr = randi([1, sub], sub, 1);
    E_temp  = E;
    E_temp.Signal = E_temp.Signal(tr, :, :);
    E_temp.O = E_temp.O(tr, :, :);
    pkrep(r,:) = getPK(E_temp, n0S);
    logregrep(r,:) = getPKbyLogReg(E_temp, n0S);
end
err_pk = std(pkrep,[],1);
err_logreg = std(logregrep, [], 1);

% function [binned] = binnize(v, nbin)
% lenv = length(v);
% frameperbin = floor(lenv/nbin);
% binned = nan(1, nbin);
% begin = 1;
% for i = 1:nbin
%     binned(i) = mean(v(begin:begin+frameperbin-1));
%     begin = begin + frameperbin;
% end


function [pkt] = getDKernel(E)
n0S = E.InputImage.n_zero_signal; 
nX =  size(E.X, 2);
% O_pref=1;
ixp=1; ixa=1+nX/2;
% choice
ch = E.O(:,1,end) - 1;
% hdxmat
hdxmat1 = squeeze(E.Signal(:,ixp,n0S+2:end));
hdxmat2 = squeeze(E.Signal(:,ixa,n0S+2:end));
disval1 = unique(hdxmat1(:));
disval2 = unique(hdxmat2(:));
len_d = length(disval1);
svmat1 = zeros(size(hdxmat1,1), len_d);
svmat2 = zeros(size(hdxmat2,1), len_d);
pkt = nan(len_d, E.Projection.n_frames-n0S-1);
for t = 1:size(hdxmat1,2)
    for r = 1:size(hdxmat1,1)
        for d = 1:len_d
            svmat1(r,d) = sum(hdxmat1(r,t)==disval1(d));
            svmat2(r,d) = sum(hdxmat2(r,t)==disval2(d));
        end
    end
    pkt(:,t) = (mean(svmat1(ch==0,:),1) - mean(svmat1(ch==1,:),1))...
        - (mean(svmat2(ch==0,:),1) - mean(svmat2(ch==1,:),1));
end
