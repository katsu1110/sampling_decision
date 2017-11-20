function E = trcut(E, cuttime)
%% cut trials
%
% written by Katsuhisa (07.11.17)
% +++++++++++++++++++++++++++++++++

% for PK analysis the following should be fine
try
    access = E.Sampling.access;
catch
    access = E.access;
end
E.Signal = E.Signal(:, :, 1:access(cuttime));
E.O = E.O(:, :, 1:cuttime);
E.Projection.n_frames = access(cuttime);

