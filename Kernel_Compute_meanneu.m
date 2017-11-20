function [pk] = Kernel_Compute_meanneu(out)

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

% basic info
s = size(out.Signal);
ntr = s(1);
% nneuron = s(2);
nneuron = 1;
n_frames = s(3);
disp([num2str(ntr) ' trials, ' num2str(nneuron) ' neurons, ' num2str(n_frames) ' n_frames'])

% trial parameters
ch = out.O(:,1,end);
ch = ch - 1;
posterior = zeros(ntr, n_frames);
for f = 1:n_frames
    posterior(:,f) = abs(out.O(:,2,f) - 0.5) + 0.5;
end
disp(['far-ch tr: ' num2str(sum(ch==1)) ', near-ch tr: ' num2str(sum(ch==0))])

% initialization
nbin = 4;
% v = [1 0.5 0.25 0.125 0 -0.125 -0.25 -0.5 -1];
% v = [-0.3 -0.2 -0.1 0 0.1 0.2 0.3];
v = [-0.5 -0.25 -0.125 -0.0625 0 0.0625 0.125 0.25 0.5];
lenv = length(v);
tkernel = zeros(lenv, nbin);
tkernel_h = zeros(lenv, nbin);
tkernel_l = zeros(lenv, nbin);
for n = 1:nneuron
    hdxmat = zeros(ntr, n_frames);
    for f = 1:n_frames
        hdxmat(:,f) = out.Signal(:,n,f);
    end

    conf = posterior(:,end);
    med = median(conf);    
% med = 0.975;

    % discritize stimuli
    q = 0.001*[-2, -1.25, -0.65, -0.15, 0.15, 0.65, 1.25, 2];     
        
    hdxnew = zeros(size(hdxmat));
    for i = 1:lenv
        if i==1
            range = hdxmat < q(1);
        elseif i==lenv
            range = hdxmat >= q(end);
        else
            range = hdxmat >= q(i-1) & hdxmat < q(i);
        end

        hdxnew(range) = v(i);
    end
    hdxmat = hdxnew;

    % compute PK via difference in choice-triggered stimuli
    tkernel_each = zeros(lenv, nbin);
    tkernel_h_each = zeros(lenv, nbin);
    tkernel_l_each= zeros(lenv, nbin);
    frameperbin = floor(n_frames/nbin);
    begin = 1;
    for k = 1:nbin
        tkernel_each(:,k) = getKernel(hdxmat(:, begin:begin+frameperbin-1), ch);
        tkernel_h_each(:,k) = getKernel(hdxmat(conf > med, begin:begin+frameperbin-1), ch(conf > med));
        tkernel_l_each(:,k) = getKernel(hdxmat(conf < med, begin:begin+frameperbin-1), ch(conf < med));
        begin = begin + frameperbin;
    end
    
    % sum up
    tkernel = tkernel + tkernel_each;
    tkernel_h = tkernel_h + tkernel_h_each;
    tkernel_l = tkernel_l + tkernel_l_each;
end

% average across neurons
tkernel = tkernel/nneuron;
tkernel_h = tkernel_h/nneuron;
tkernel_l = tkernel_l/nneuron;

% PK amplitude
pk0 = mean(tkernel, 2);
pk0_h = mean(tkernel_h, 2);
pk0_l = mean(tkernel_l, 2);
pk = zeros(1, nbin);
pk_h = zeros(1, nbin);
pk_l = zeros(1, nbin);
for n = 1:nbin
    pk(n) = dot(tkernel(:,n), pk0);
    pk_h(n) = dot(tkernel_h(:,n), pk0);
    pk_l(n) = dot(tkernel_l(:,n), pk0);
end

% visualization
close all;
figure;
subplot(1,2,1)
plot(v, pk0, '-', 'color', 'r')
hold on;
plot(v, pk0_h, '-','color',y)
hold on;
plot(v, pk0_l, '-','color',g)
legend('overall','high conf.','low conf.')
title('time-averaged PK')

subplot(1,2,2)
plot(1:nbin, pk, '-','color','r')
hold on;
plot(1:nbin, pk_h, '-','color',y)
hold on;
plot(1:nbin, pk_l, '-','color',g)
title('PK amplitude')

function [pk0] = getKernel(hdxmat, ch)

% trial averaged stimulus distributions
disval = unique(hdxmat);
len_d = length(disval);
svmat = zeros(size(hdxmat,1), len_d);
for r = 1:size(hdxmat,1)
    for d = 1:len_d
        svmat(r,d) = sum(hdxmat(r,:)==disval(d));
    end
end

% compute PK for 0% stimulus
pk0 = mean(svmat(ch==1,:),1) - mean(svmat(ch==0,:),1);
