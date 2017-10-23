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

% out = out_pk;
% 
% % what is 'O'? Orientation, or the stimulus, right?
% v = nan(1, 16);
% for i = 1:16
%     v(i) = out.O(i,1,1);
% end
% 
% figure;
% plot(v)
% 
% % So the third dimension is "choice"
% 
% % is the first dimention stimulus?
% v = nan(1, 16);
% for i = 1:16
%     v(i) = out.O(i,1,1);
% end
% 
% figure;
% plot(v)
% 
% % Signal?
% v = nan(1,76);
% for i = 1:76
%     v(i) = out.Signal(1,1,i);
% end
% 
% plot(v)
% 
% % Signal?
% v = nan(1,16);
% for i = 1:16
%     v(i) = out.Signal(i,1,1);
% end
% 
% plot(v)
% 
% % is the first dimention stimulus?
% v = nan(1, 16);
% vv = nan(1, 16);
% for i = 1:16
%     v(i) = out.O(i,1,1);
%     vv(i) = out.O(i,1,end);
% end
% 
% figure;
% plot(v); hold on; plot(vv)
% 
% % G (higher-order visual areas for decision)?
% for n = 1:256
%     v = zeros(1,100);
%     for i = 2:100
%         v(i) = v(i-1) + out.G(1, 1, n, i);
%         hold on;
%     end
%     disp(v)
%     plot(v); hold on;
% end
% 
% % Style???
% for i = 1:3
%     plot(out.S(i,:))
%     hold on;
% end
% 
% % Pt(D) ++++++++++++++++
% out = out_2afc_corr;
% v = zeros(1,100);
% for n = 1:100
%     a = out.O(n,:,:);
%     for i = 1:100
%         v(i) = a(1,2,i);
%     end
%     plot(v)
%     hold on;
end