function E = merge_output(Es)
%% merge output structures in given cell array
% (assuming trial structures (e.g.'stimulus_contrast', 'n_frames', etc)
% are the same across outputs)
%
% written by Katsuhisa (06.11.17)
% +++++++++++++++++++++++++++++++++++++++++++++++++

% just for PK computation, the following should be enough
% rng(1220);
% rtr = datasample(1:size(Es{1}.Signal,1), round(size(Es{1}.Signal,1)/length(Es)));
E = Es{1};
% E.Signal = E.Signal(rtr,:,:);
% E.O = E.O(rtr,:,:);
for n = 2:length(Es)
%     rtr = datasample(1:size(Es{1}.Signal,1), round(size(Es{1}.Signal,1)/length(Es)));
    E.Signal = cat(1, E.Signal, Es{n}.Signal);
    E.O = cat(1, E.O, Es{n}.O);
end