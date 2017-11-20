function Kernel_explore(E)

for i = 1:10
    Kernel_Compute(E,'nsplit',4,'cuttime',i*10,'save');
end
close all;