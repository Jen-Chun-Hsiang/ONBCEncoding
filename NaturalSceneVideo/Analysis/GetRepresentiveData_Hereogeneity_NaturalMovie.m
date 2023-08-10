close all; clear; clc;
%%
CellStr = {'122718_2', '10119_4', '11919_1', '11019_3', '10519_3', '12119_3', '11819_1',...
    '12019_4', '53018_2', '10119_2', '61018_1', '53018_1', '92418_2', '12019_7', '10919_3'};
nCell = length(CellStr);
qualthr = 0.2;
%%
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\BCHeterogeneityLengthConstant';
Topic = 'Ai148_AAV-Grm6Cre';
clear gDdp gDe gDi gDp
for c = 1:nCell
    FileNames = dir(fullfile(FolderPath, '*.mat'));
    FileIds = find(contains({FileNames.name},[Topic '_SponStimLengthConstant_exactestimation_' CellStr{c}]));
    nFile = length(FileIds);
    clear De Di Dap Dai Dbp Dbi Dm
    for f = 1:nFile
        FileName = FileNames(FileIds(f)).name;
        F = load([FolderPath '\' FileName]);
        if median(F.Dp(:), 'omitnan')<qualthr
            continue
        end
        if size(F.Dai, 1) < 2
            continue
        end
        nLen = length(F.LenConsts);
        if f == 1 && c == 1
            Dcp = nan(nCell, nLen);
            Dci = nan(nCell, nLen);
            mDci = nan(nCell, nLen);
            Dcs = nan(nCell, nLen);
            Dch = nan(nCell, nLen);
            Dsp = nan(nCell, nLen);
            Dsi = nan(nCell, nLen);
            Dtab = nan(nCell, 3);
        end
%         rmids = median(F.Dp, 2, 'omitnan')<qualthr;
%         F.Dai(rmids, :, :) = [];
%         F.Dap(rmids, :, :) = [];
%         F.Dbi(rmids, :, :) = [];
%         F.Dbp(rmids, :, :) = [];
%         F.De(rmids, :) = [];
%         F.Di(rmids, :) = [];
%         F.Dp(rmids, :) = [];
%         F.Ddp(rmids) = [];
        ai = permute(F.Dai, [1 3 2]);
        ai = squeeze(median(ai, 2));
        ap = permute(F.Dap, [1 3 2]);
        ap = squeeze(median(ap, 2));
        as = permute(F.Dais, [1 3 2]);
        as = squeeze(median(as, 2));
        ah = permute(F.Daih, [1 3 2]);
        ah = squeeze(median(ah, 2));
        bi = permute(F.Dbi, [1 3 2]);
        bi = squeeze(median(bi, 2));
        bp = permute(F.Dbp, [1 3 2]);
        bp = squeeze(median(bp, 2));
        if f == 1
            De = median(F.De, 2);
            Dm = median(F.Dm, 2);
            Di = median(F.Di, 2);
            Dap = ap;
            Dai = ai;
            Das = as;
            Dah = ah;
            Dbp = bp;
            Dbi = bi;
            
            LenConsts = F.LenConsts;
            Dtab(c, :) = [c F.bctype median(F.Dp(:), 'omitnan')];
            
            gDdp{c} = F.Ddp;
            gDe{c} = F.De;
            gDi{c} = F.Di;
            gDp{c} = F.Dp;
        else
            De = [De; median(F.De, 2)];
            Dm = [Dm; median(F.Dm, 2)];
            Di = [Di; median(F.Di, 2)];
            Dap = [Dap; ap];
            Dai = [Dai; ai];
            Das = [Das; as];
            Dah = [Dah; ah];
            Dbp = [Dbp; bp];
            Dbi = [Dbi; bi];
        end
        clear F
    end
    if ~exist('Dai', 'var')
        continue
    end
    for i = 1:nLen
        kids = ~isnan(Dap(:, i)); 
        Dcp(c, i) = corr(1-Dap(kids, i), 1-De(kids));
        kids = ~isnan(Dai(:, i)); 
        Dci(c, i) = corr(1-Dai(kids, i), 1-De(kids));
        kids = all(~isnan([Dai(:, i) De]), 2); 
        mDci(c, i) = corr(1-Dai(kids, i), 1-Dm(kids));
        kids = ~isnan(Das(:, i)); 
        Dcs(c, i) = corr(1-Das(kids, i), 1-De(kids));
        kids = ~isnan(Dah(:, i)); 
        Dch(c, i) = corr(1-Dah(kids, i), 1-De(kids));
        Dsp(c, i) = corr(1-Dbp(:, i), 1-De(:));
        Dsi(c, i) = corr(1-Dbi(:, i), 1-De(:));
%         Dsp(c, i) = corr(1-reshape(Dbp(:, i, :), [], 1), 1-Di(:));
%         Dsi(c, i) = corr(1-reshape(Dbi(:, i, :), [], 1), 1-Di(:));
    end
    clear Dai Dap De
end
%%
filnam = 'Heterogenity_NatMov_ReconstructedBC.mat';
save(['\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\'...
    filnam], 'Dtab', 'Dcp', 'Dci', 'Dcs', 'Dch', 'mDci');
%%
keyboard;
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCelen = [19.4 18.4 20.8 36.6 14.4 19.7 43.4];
%%
nclip = 11;
Colors = lines(7);
for i = 1:7
    switch i
        case 1
            xtic = 0:10:40;
            xticlab = {'0', '', '20', '', '40'};
        case 2
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
        case 3
            xtic = 0:10:40;
            xticlab = {'0', '', '20', '', '40'};
        case 4
            xtic = 0:10:30;
            xticlab = {'0', '10', '20', '30'};
        case 5
            xtic = 0:12.5:50;
            xticlab = {'0',  '', '25', '', '50'};
        case 6
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
        case 7
            xtic = 0:12.5:50;
            xticlab = {'0',  '', '25', '', '50'};
    end
    close all
    figure;
    bids = find(Dtab(:, 2) == BCTypes(i));
    % 1
    subplot(2, 3, [1 4]); hold on
    x = gDp{bids(1)}(:);
    z1 = repmat(1:11, size(gDp{bids(1)}, 1), 1);
    z2 = repmat((1:size(gDp{bids(1)}, 1))', 1, 11);
    z1 = z1(:);
    z2 = z2(:);
    y1 = gDi{bids(1)}(:);
    y2 = gDe{bids(1)}(:);
    scatter(x, y1, 5, [247 179 32]/255, 'filled');
    scatter(x, y2, 5, [166 30 230]/255, 'filled');
    plot([0 1], [0 1], '--k');
    xticks(0:0.25:1);
    xticklabels({'0', '', '0.5', '', '1'});
    yticks(0:0.25:1);
    yticklabels({'0', '', '0.5', '', '1'});
    xlim([0 1]);
    ylim([0 1]);
    xlabel('Repeat reliability');
    ylabel('Correlation coefficient');
    title('Simulation estimation');
    box off
    keyboard;
    % 2
    subplot(2, 3, [2 5]); hold on
    x = reshape(repmat(gDdp{bids(1)}(:), 1, nclip), [], 1);
    y1 = gDi{bids(1)}(:);
    y2 = gDe{bids(1)}(:);
    scatter(reshape(repmat(gDdp{bids(1)}(:), 1, nclip), [], 1), gDi{bids(1)}(:), 5, [247 179 32]/255, 'filled');
    scatter(reshape(repmat(gDdp{bids(1)}(:), 1, nclip), [], 1), gDe{bids(1)}(:), 5, [166 30 230]/255, 'filled');
    xticks(xtic);
    xticklabels(xticlab);
    yticks(0:0.25:1);
    yticklabels({'0', '', '0.5', '', '1'});
    ylim([0 1]);
%     xlim([0 max(gDdp{bids(1)}(:))])
    xlabel('Path Distance');
    ylabel('Correlation coefficient');
    title('Simulation estimation');
    box off
    keyboard;
    nb = length(bids);
    Colors = [166 30 230]/255;
    Colors = [Colors; zeros(nb-1, 3)];
    for j = 1:nb
        % 3
        subplot(2, 3, 3); hold on
        plot(LenConsts, Dcp(bids(j), :), 'Color', Colors(j ,:));
        if j == 1
            plot(BCelen(i)*ones(1, 2), [0 1], '--', 'Color', 0.4*ones(1, 3));
        end
        xlim([min(LenConsts) max(LenConsts)]);
        yticks(0:0.25:1);
        yticklabels({'0', '', '0.5', '', '1'});
        set(gca, 'XScale', 'log');
        xlabel('Length constant');
        ylabel('Correlation coefficient');
        title('Path distance');
        ylim([-0.05, 1]);
        box off
        % 4
        subplot(2, 3, 6); hold on
        plot(LenConsts, Dci(bids(j), :), 'Color', Colors(j ,:));
        if j == 1
            plot(BCelen(i)*ones(1, 2), [0 1], '--', 'Color', 0.4*ones(1, 3));
        end
        xlim([min(LenConsts) max(LenConsts)]);
        yticks(0:0.25:1);
        yticklabels({'0', '', '0.5', '', '1'});
        xlabel('Length constant');
        ylabel('Correlation coefficient');
        title('Image distance');
        ylim([0, 1]);
        set(gca, 'XScale', 'log');
        box off
    end
    sgtitle(sprintf('%s',  BCTypeLabels{i}));
    %
    keyboard
end
 %%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure6_Heterogeneity_NatMov_DistanceCorr_%s', SaveFolder, BCTypeLabels{i});
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);