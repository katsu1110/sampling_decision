function getMoreTr(lab, repeat)

storepath = 'Z:\Katsuhisa\data\sampling_decision\co5samp100\';

for r = 1:repeat
    out = S_Experiment(S_Exp_Para('PK-amplitude'));
    out = rmFieldforMemory(out);
    save([storepath 'out_' num2str(lab) '.mat'], 'out', '-v7.3')
    lab = lab + 1;
end