function [E] = concat_all

storepath = 'Z:\Katsuhisa\data\sampling_decision\co5samp100\';
listings = dir(storepath);
listings(1:2) = [];
load([storepath listings(1).name])
E = out;
for i = 2:length(listings)
    try
        load([storepath listings(i).name])
        E.Signal = cat(1, E.Signal, out.Signal);
        E.O = cat(1, E.O, out.O);
        disp([listings(i).name ' was concatenated!'])
    catch
        disp(['failed to concatenate ' listings(i).name '. Next!'])
        continue
    end
end

saveflag = 0;
while saveflag==0
    try
        save([storepath 'E.mat'], 'E', '-v7.3')
        saveflag = 1;
    catch
        E.Signal(1:10000,:,:) = [];
        E.O(1:10000,:,:) = [];
    end
end
    