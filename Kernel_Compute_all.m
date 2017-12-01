function [pk] = Kernel_Compute(out)

s = size(out.Signal);
ntr = s(1);
nneuron = s(2);
% nneuron = 1;
n_frames = s(3);
disp([num2str(ntr) ' trials, ' num2str(nneuron) ' neurons, ' num2str(n_frames) ' n_frames'])

hdxmat = zeros(ntr*nneuron, n_frames);
ch = ones(ntr*nneuron, n_frames);
c = 1;
for i = 1:ntr
    for n = 1:nneuron
        for f = 1:n_frames
            hdxmat(c,f) = out.Signal(i,n,f);
            ch(c,f) = out.O(i,1,f);     
        end
        c = c + 1;
    end
end

ch = ch(:,end);
ch = ch - 1;
disp(['far-ch tr: ' num2str(sum(ch==1)/nneuron) ', near-ch tr: ' num2str(sum(ch==0)/nneuron)])

% discritize stimuli
v = [-2 -1 0 1 2];
q = nan(1,5);
for i = 1:5
    q(i) = quantile(hdxmat(:)', 0.2*i);
end
hdxnew = zeros(size(hdxmat));
for r = 1:size(hdxmat)
    for c = 1:n_frames
        for t = 1:5
            if hdxmat(r,c) >= q(t)-0.2 && hdxmat(r,c) < q(t) 
                hdxnew(r,c) = v(t);
                break
            end
        end
    end
end
hdxmat = hdxnew;
pk0 = getKernel(hdxmat, ch);

% compute PK via difference in choice-triggered stimuli
nbin = 16;
pk = nan(1, nbin);
frameperbin = floor(n_frames/nbin);
begin = 1;
for n = 1:nbin
    pk_temp = getKernel(hdxmat(:, begin:begin+frameperbin-1), ch);
    pk(n) = dot(pk_temp, pk0);
    begin = begin + frameperbin;
end

close all;
figure;
subplot(1,2,1)
plot(pk0)
subplot(1,2,2)
plot(pk)

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
pk0 = mean(svmat(ch==0,:),1) - mean(svmat(ch==1,:),1);

