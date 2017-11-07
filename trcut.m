function E = trcut(E, cuttime)
%% cut trials
%
% written by Katsuhisa (07.11.17)
% +++++++++++++++++++++++++++++++++

% for PK analysis the following should be fine
E.Signal = E.Signal(:, :, 1:E.Sampling.access(cuttime));
E.O = E.O(:, :, 1:cuttime);
E.Projection.n_frames = E.Sampling.access(cuttime);
