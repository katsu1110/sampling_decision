function [tpka] = time_pka(E, nbin, plot_flag, repeat, type)

if nargin<2
    nbin = 4;
elseif nargin<3
    plot_flag = 1;
elseif nargin<4
    repeat = 0;
elseif nargin<5
    type = 0;
end
mat = nan(4, size(E.Signal, 3));
c = 1;
for i = 1:size(E.Signal, 3)
%     ct = find(E.Projection.access==i);
    ct = find(E.Projection.access==i, 1, 'first');
    disp(['time: ' num2str(ct)])
    try
%         pk0_temp = [];
%         pk0_temp_c1 = [];
%         pk0_temp_c2 = [];
%         for k = 1:length(ct)
%             Ecut = trcut(E, ct(k));
%             [temp0, temp1] = Kernel_Compute(Ecut, 'nbin', 4);
%             pk0_temp = [pk0_temp; temp1'];
%             pk0_temp_c1 = [pk0_temp_c1; temp0(1,:)];
%             pk0_temp_c2 = [pk0_temp_c2; temp0(2,:)];
%         end4)4
%         tpka.evidence(c).pka0 = mean(pk0_temp,1);
%         tpka.evidence(c).pkas = [mean(pk0_temp_c1,1); mean(pk0_temp_c2,1)];
        Ecut = trcut(E,ct);
        tpka.evidence(c).time = ct;
        delta = nan;
        sem = nan;
        if type==0
            if mod(size(Ecut.Signal,3)-2, nbin)==0
                [tpka.evidence(c).pkas, tpka.evidence(c).pka0] = Kernel_Compute(Ecut,'nbin',nbin,'resample',repeat);
                delta = (tpka.evidence(c).pkas(2,end) - tpka.evidence(c).pkas(1,end))...
                    /max(tpka.evidence(c).pka0(1,:));
                if repeat > 0
                    sem = mean([tpka.evidence(c).pkas(3,end), tpka.evidence(c).pkas(4,end)])...
                        /max(tpka.evidence(c).pka0(1,:));
                end
                c = c + 1;
            end
        else
            [tpka.evidence(c).pkas, tpka.evidence(c).pka0] = Kernel_Compute(Ecut,'resample',repeat);
            delta = (tpka.evidence(c).pkas(2,:) - tpka.evidence(c).pkas(1,:))...
                /max(tpka.evidence(c).pka0(1,:));
%             delta = mean(delta(end-2:end));
%             delta = delta(end-1);
                delta = mean(delta(end-3:end-1));
            if repeat > 0
                sem = mean([tpka.evidence(c).pkas(3,end-2:end), tpka.evidence(c).pkas(4,end-2:end)])...
                    /max(tpka.evidence(c).pka0(1,:));
            end
            c = c + 1;
        end
          
        mat(1,i) = ct;
        mat(2,i) = i - 2;
        mat(3,i) = delta;
        mat(4,i) = sem;
      
    catch
        disp(['time ' num2str(ct) ' skipped'])
        continue
    end    
end
tpka.mat = mat;
tpka.contrast = E.Projection.stimulus_contrast(1);

if plot_flag==1
    % close all;
    h = figure;
    errorbar(mat(2,:), mat(3,:), mat(4,:), 'linestyle', 'none', 'CapSize', 0, 'color', 'r')
    hold on;
    scatter(mat(2,:), mat(3,:), 30, 'filled', 'markerfacecolor', 'r', 'markerfacealpha', 0.4)
    hold on;
    xx = get(gca, 'XLim');
    plot(xx, [0 0], ':k')
    % hold on;
    % plot(xx, -0.2470*[1 1], '-b')
    xlabel('time (neuronal resplution)')
    ylabel('\Delta PKA at the last time bin')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    set(h, 'Name' ,['SamplingDecision_co' num2str(E.Projection.stimulus_contrast(1))],...
        'NumberTitle', 'off')
end