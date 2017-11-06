function E = merge_output(Es)
%% merge output structures in given cell array
% (assuming trial structures (e.g.'stimulus_contrast', 'n_frames', etc)
% are the same across outputs)
%
% written by Katsuhisa (06.11.17)
% +++++++++++++++++++++++++++++++++++++++++++++++++

% just for PK computation, the following should be enough
E = Es{1};
for n = 2:length(Es)
    E.Signal = cat(3, E.Signal, Es{n}.Signal);
    E.O = cat(3, E.Signal, Es{n}.O);
end
