function explore_SamplingDecision(out)

% how consistent the dynamic stimuli are across neurons?
s = size(out.Signal);
ntr = s(1);
nneuron = s(2);
% nneuron = 1;
n_frames = s(3);
disp([num2str(ntr) ' trials, ' num2str(nneuron) ' neurons, ' num2str(n_frames) ' n_frames'])

stm = nan(nneuron, n_frames);
for n = 1:nneuron
    for f = 1:n_frames
        stm(n,f) = out.Signal(1, n, f);
    end
end

close all;
imagesc(stm)
xlabel('frames')
ylabel('neurons')
colorbar

% stimulus distribution
figure;
histogram(stm(:))
hold on;
% discritize stimuli
% q = 0.001*[-1.3, -0.8, -0.4, -0.1, 0.1, 0.4, 0.8, 1.3]; 
q = 0.001*[-2, -1.25, -0.65, -0.15, 0.15, 0.65, 1.25, 2]; 
yy = get(gca, 'YLim');
for i = 1:length(q)
    plot(q(i)*[1,1], yy, '--r')
    hold on;
end

% posterior (belief): example trials
figure;
N = datasample(1:100, 7);
v = nan(1, 100);
for n = 1:7
    for i = 1:100
        v(i) = abs(out.O(N(n),2,i) - 0.5) + 0.5;
    end
    plot(v)
    hold on;
end
xlabel('frames')
ylabel('posterior')

% V1 activity: example trials
figure;
k = 2;
N = datasample(1:100, 7);
v = nan(1, 100);
for n = 1:7
    for i = 1:100
        v(i) = out.X(N(n),k,i);
    end
    plot(cumsum(v,2))
    hold on;
end
xlabel('frames')
ylabel('V1 res.')

% V1 activity: example neurons
figure;
k = 2;
N = datasample(1:100, 7);
v = nan(1, 100);
for n = 1:7
    for i = 1:100
        v(i) = out.X(k,N(n),i);
    end
    plot(cumsum(v,2))
    hold on;
end
xlabel('frames')
ylabel('V1 res.')

% V1 activity & choice
figure;
ch = out.O(:,1,end);
k = 2;
r1 = zeros(1, 100);
r2 = zeros(1, 100);
for i = 1:100
    r1(i) = mean(out.X(ch==1,k,i));
    r2(i) = mean(out.X(ch==2,k,i));
end
plot(cumsum(r1,2), '-r')
hold on;
plot(cumsum(r2,2), '-b')
xlabel('frames')
ylabel('cumulative V1 res.')
legend('choice 1', 'choice 2')
legend('boxoff')