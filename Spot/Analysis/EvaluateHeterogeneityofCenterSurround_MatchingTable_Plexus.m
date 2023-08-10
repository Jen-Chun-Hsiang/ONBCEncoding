% run CircleAnalysis_BipolarCellTerminal.m
options.IsPCA = 1;
options.IsDebug = 0;
options.EstimationType = 'DoG';
%[Params, ~, Costv] = estimateCenterSurround(ROIFeatures, options);


options.EstimationType = 'Amp';
[SimpleCSRio, Curves] = estimateCenterSurround(ROIFeatures, options);

options.EstimationType = 'temporal';
[SimpleTemporalRio, TemporalCurves] = estimateCenterSurround(ROIFeatures, options);

options.IsPCA =0;
options.EstimationType = 'Cluster';
[~, ClusterCurves] = estimateCenterSurround(ROIFeatures, options);
figure; imagesc(ClusterCurves'); colorbar;
% save('./Results/CenterSurroundRatio_070822.mat', 'Params', 'Costv', 'SimpleCSRio', 'Curves');
% keyboard;
%clear ROIFeatures
%%
yLabels = {'CenterSurround ratio','Repeat reliability', 'Response function'};
RescaleV = [ones(1, 2)*1000, 1, 1];
SummaryD.Table = [];
SummaryD.RespFunc = [];
SummaryD.TprRespFunc = [];
SummaryD.ClusterFeatures = [];
clear NotProcessing
for i = 1:size(IdTable, 1)
    cIds = find(ROITable(:, 1) == IdTable(i, 1) & ROITable(:, 2) == IdTable(i, 2) &...
        ROITable(:, 4) == IdTable(i, 3));
    cROITable = ROITable(cIds, :);
    UniExp = unique(cROITable(:, [1 2 3 4]), 'rows');
    
    nUniExp = size(UniExp, 1);
    Colors = lines(nUniExp);
    Dots = SimpleCSRio(cIds);
    Ribs = ROIReliab(cIds, :);
    Tprs = SimpleTemporalRio(cIds);
    Raws = ROIRaw(cIds, :);
    %     Amps = ROIDmtAmp(cIds, :);
    Amps = Curves(cIds, :);
    TprAmps = TemporalCurves(cIds, :);
    CluFet = ClusterCurves(cIds, :);
    cAmp = [];
    cTrpAmps = [];
    cClusterFet = [];
    nROIsets = [];
    close all
    figure('visible', 'off');
    %         figure;
    set(gcf, 'Position', get(0, 'Screensize').*[0.5 0.5 0.3 0.7]);
    x = 0;
    for j = 1:nUniExp
        cTable = [];
        eIds = find(cROITable(:, 3) == UniExp(j, 3));
        neIds = length(eIds);
        x = (x(end)+1):(x(end)+neIds);
        plotv = Dots(eIds);
        cTable = [cTable; cROITable(eIds,:)];
        subplot(3, 1, 1)
        scatter(x, plotv, 25, Colors(j, :), 'filled');hold on
        cTable = [cTable plotv];
        subplot(3, 1, 2)
        plotv = mean(Ribs(eIds, 2:4), 2);
        scatter(x, plotv, 25, Colors(j, :), 'filled');hold on
        cTable = [cTable plotv];
        cAmp = [cAmp; Amps(eIds, :)];
        plotv = Tprs(eIds);
        cTable = [cTable plotv];
        cTrpAmps = [cTrpAmps; TprAmps(eIds, :)];
        cClusterFet = [cClusterFet; CluFet(eIds, :)];
        plotv = Raws(eIds, :);
        cTable = [cTable plotv];
        SummaryD.Table = [SummaryD.Table; cTable];
    end
    SummaryD.RespFunc = [SummaryD.RespFunc; cAmp];
    SummaryD.TprRespFunc = [SummaryD.TprRespFunc; cTrpAmps];
    SummaryD.ClusterFeatures = [SummaryD.ClusterFeatures; cClusterFet];
    clear xlabels
    for j = 1:length(Dmt)
        xlabels{j} = sprintf('%d', Dmt(j));
    end
    %
    subplot(3, 1, 1)
    title(sprintf('IPL column: %d-%d ', cROITable(1, 1), cROITable(1, 2)));
    xlim([0 max(x)]);
    box off
    xlabel('#ROI');
    ylabel(yLabels{1});
    %
    subplot(3, 1, 2);
    xlim([0 max(x)]);
    box off
    xlabel('#ROI');
    ylabel(yLabels{2});
    %
    subplot(3, 1, 3)
    imagesc(cAmp); colorbar
    box off
    xticks(1:length(Dmt));
    xticklabels(xlabels);
    ylabel(yLabels{3});
    xlabel('#Spot diameter');
    %         keyboard;
    %%
%     FleNam = sprintf('./Figures/%s/%s_Plexus_SpotCS__reliab_%d_%d', ImgSaveFolder, Topic,...
%         cROITable(1, 1), cROITable(1, 2));
    %     print('-depsc','-painters','-loose', '-r300',FleNam)
    %     saveas(gcf,[FleNam '.fig']);
%     saveas(gcf,[FleNam '.png']);
    %     keyboard;
end
keyboard;
%%
CelTabLabels = {'Day', 'Region', 'Exp', 'Type', 'ROI', 'CelCol', 'CS', 'QL',...
    'TS', 'Rawavg', 'Rawmed'};
assert(length(CelTabLabels) == size(SummaryD.Table, 2));
Collab2num = @(in) find(cellfun(@(x) strcmpi(x, in), CelTabLabels));
SaveFileName = './Results/CenterSurround_Reliability_ROI_Plexus_100722.mat';
% save(SaveFileName, 'SummaryD', 'Dmt', 'CelTabLabels', 'Collab2num');

%%
cIds = SummaryD.Table(:, 7) > 0;
x = abs(SummaryD.Table(cIds, 8));
y = log10(SummaryD.Table(cIds, 7));
figure; scatter(x, y, 15, 'k', 'filled');hold on
ProbThr = 0.521;
CvarThr = 2;
gIds = cIds & SummaryD.Table(:, 8) > ProbThr; % &  log10(SummaryD.Table(:, 7)) > 0;
plot(ProbThr*ones(1,2), [min(y), max(y)], '--b');
title(sprintf('Kept ROIs %d/%d',sum(gIds), length(gIds)));
xlabel('Reliability (R value)');
ylabel('STD Center-surround ratio (log10)');
ylim([-3 5]);

%% Save in mat for batch effect removal
uniD = unique(SummaryD.Table(gIds, 1));
clear BCdat
for i = 1:length(uniD)
    cids = SummaryD.Table(:, 1) == uniD(i);
    val = SummaryD.ClusterFeatures(cids & gIds, :);
    BCdat{i} = (val-min(val, [], 2))./range(val, 2);
end
SaveFileName = './Results/batcheffect_data2python_Plexus_090722.mat';
save(SaveFileName, 'BCdat', 'uniD');
