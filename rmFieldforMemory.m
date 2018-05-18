function out = rmFieldforMemory(out)
% remove some irrelevant fields from the output structure 
% for more simulations

% get info
out.ntr = size(out.Signal, 1);
out.nv1 = size(out.X, 2);
out.nvx = size(out.G, 2);
out.access = out.Sampling.access;
out.n0S = out.InputImage.n_zero_signal;

% remove other fields
out = rmfield(out, 'X');
out = rmfield(out, 'G');
out = rmfield(out, 'L');
out = rmfield(out, 'T');
out = rmfield(out, 'S');
out = rmfield(out, 'InputImage');
out = rmfield(out, 'Sampling');