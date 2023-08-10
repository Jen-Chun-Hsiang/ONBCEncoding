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
% save('./Results/CenterSurroundRatio_070822.mat', 'Params', 'Costv', 'SimpleCSRio', 'Curves');
% keyboard;
clear ROIFeatures
%%
% BCTypes = [5 57 58 6 7 89];
% BCTypelabels = {'BC5', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCTypes = [5 50 51 57 58 6 7 89];
BCTypelabels = {'BC5', 'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
yLabels = {'CenterSurround ratio','Repeat reliability', 'Response function'};
RescaleV = [ones(1, 2)*1000, 1, 1];
SummaryD.Table = [];
SummaryD.RespFunc = [];
SummaryD.TprRespFunc = [];
SummaryD.ClusterFeatures =[];
clear NotProcessing
for i = 1:size(IdTable, 1)
    cIds = find(ROITable(:, 1) == IdTable(i, 1) & ROITable(:, 2) == IdTable(i, 2) &...
        ROITable(:, 4) == IdTable(i, 3));
    cROITable = ROITable(cIds, :);
    UniExp = unique(cROITable(:, 3));
    
    if IdTable(i, 1) < 100000
        SheetName = sprintf('0%d-%d', IdTable(i, 1), IdTable(i, 2));
    else
        SheetName = sprintf('%d-%d', IdTable(i, 1), IdTable(i, 2));
    end
    try
        MatchingT = readtable('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\BCTerminal_MatchingTable_10192022.xlsx',...
            'Sheet', SheetName, 'Range','',...
            'ReadVariableNames',false);
        clc
        fprintf('Load %s (%d/%d) \n', SheetName, i, size(IdTable, 1));
    catch
        NotProcessing{i} = SheetName;
        continue
    end
    
    ExpIds = table2array(MatchingT(1,:));
    getselecids = ismember(ExpIds, UniExp);
    ExpIds = ExpIds(getselecids);
    nUniExp = numel(ExpIds);
    if length(UniExp)<nUniExp
        error('Experiment ids are not matched');
    end
    Colors = lines(nUniExp);
    LindedIds = table2array(MatchingT(2:end,getselecids));
    for j = 1:nUniExp
        cRmvIds = find(RmvIdSet{i, ExpIds(j)});
        while length(cRmvIds) ~= 0
            LindedIds(LindedIds(:, j) == cRmvIds(1), j) = nan;
            discountids = LindedIds(:, j) > cRmvIds(1);
            LindedIds(discountids, j) = LindedIds(discountids, j) - 1;
            cRmvIds = cRmvIds - 1;
            cRmvIds(1) = [];
        end
    end
    if ~(size(LindedIds, 2) == 1)
        CountRepeat = sum(~isnan(LindedIds), 2);
        LindedIds(CountRepeat<2, :) = [];
    end
%     if IdTable(i, 1) == 10119 && IdTable(i, 2) == 2
%         keyboard;
%     end
    nPairs = size(LindedIds, 1);
    %     Colors = [ones(1, 3); Colors];
    %     Colors = [Colors, ones(nUniExp, 1)*0.3];
    %Dots = [Params(cIds, [1 3])'; (Params(cIds, 2)./Params(cIds, 4))'; SimpleCSRio(cIds)'];
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
    %     figure;
    set(gcf, 'Position', get(0, 'Screensize').*[0.5 0.5 0.3 0.7]);
    x = 1:nPairs;
    for j = 1:nUniExp
        cTable = [];
        eIds = find(cROITable(:, 3) == ExpIds(j));
        plotv = Dots(eIds);
        sIds = ~isnan(LindedIds(:, j));
        nROIsets(j) = sum(sIds);
        
        cTable = [cTable; repmat(cROITable(eIds(1), 1:4), sum(sIds), 1),...
            x(sIds)' LindedIds(sIds, j)];
        subplot(3, 1, 1)
        scatter(x(sIds), plotv(LindedIds(sIds, j)), 25, Colors(j, :), 'filled');hold on
        cTable = [cTable plotv(LindedIds(sIds, j))];
        subplot(3, 1, 2)
        plotv = mean(Ribs(eIds, 2:4), 2);
        scatter(x(sIds), plotv(LindedIds(sIds, j)), 25, Colors(j, :), 'filled');hold on
        cTable = [cTable plotv(LindedIds(sIds, j))];
        plotv = Tprs(eIds);
        cTable = [cTable plotv(LindedIds(sIds, j))];
        cAmp = [cAmp; Amps(LindedIds(sIds, j), :)];
        cTrpAmps = [cTrpAmps; TprAmps(LindedIds(sIds, j), :)];
        cClusterFet = [cClusterFet; CluFet(LindedIds(sIds, j), :)];
        plotv = Raws(eIds, :);
        cTable = [cTable plotv(LindedIds(sIds, j), :)];
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
    TypeLabel = BCTypelabels{find(BCTypes == cROITable(1, 4))};
    title(sprintf('Cell: %d-%d (Type: %s) ', cROITable(1, 1), cROITable(1, 2),...
        TypeLabel));
    xlim([0 max(nROIsets)]);
    box off
    xlabel('#ROI');
    ylabel(yLabels{1});
    %
    subplot(3, 1, 2);
    xlim([0 max(nROIsets)]);
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
    %     keyboard;
    %%
%     FleNam = sprintf('./Figures/%s/%s_%s_SpotCS_reliab_%d_%d', ImgSaveFolder, Topic,...
%         TypeLabel, cROITable(1, 1), cROITable(1, 2));
%     print('-depsc','-painters','-loose', '-r300',FleNam)
%     saveas(gcf,[FleNam '.fig']);
%     saveas(gcf,[FleNam '.png']);
%     keyboard;
end
TabLabels = {'Day', 'Region', 'Exp', 'Type', 'ROIpair', 'ROIId', 'CS', 'QL', 'TS', 'RawIntensityMean', 'RawIntensityMedian'};
assert(length(TabLabels) == size(SummaryD.Table, 2));
Name2ColNum = @(in) find(cellfun(@(x) strcmpi(x, in), TabLabels));
keyboard;
%%
SaveFileName = './Results/CenterSurround_Reliability_ROI_BC_122922.mat';
save(SaveFileName, 'SummaryD', 'Dmt');
%%
% figure; scatter(abs(SummaryD.Table(:, 7)), log10(SummaryD.Table(:, 6)), 15, 'k', 'filled');
%%
% figure; gscatter(abs(SummaryD.Table(:, 7)), log10(SummaryD.Table(:, 6)), SummaryD.Table(:, 4), lines(6));
%%
close all
% BCTypes = [5 57 58 6 7 89];
% BCTypelabels = {'BC5', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCTypes = [5 50 51 57 58 6 7 89];
BCTypelabels = {'BC5', 'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
figure;
for i = 1:length(BCTypes)
    subplot(2,length(BCTypes)/2, i);
    histogram(log10(abs(SummaryD.Table(SummaryD.Table(:, 4) == BCTypes(i), 6))), linspace(-3, 2, 30));
    title(sprintf('%s', BCTypelabels{i}));
end
%% Reduce to individual cell and estimate the reliabilities with CS variation across trials
UniCelTab = unique(SummaryD.Table(:, 1:2), 'rows');
nUniCelTab = size(UniCelTab, 1);
SummaryD.CelTable = [];
SummaryD.CelRespf = [];
SummaryD.CelTprRespf = [];
SummaryD.CelClusterFeature = [];
RespFuncL = size(SummaryD.RespFunc, 2);
TprRespFuncL = size(SummaryD.TprRespFunc, 2);
CluFeatL = size(SummaryD.ClusterFeatures, 2);
for i = 1:nUniCelTab
    cIds = find(SummaryD.Table(:, 1) == UniCelTab(i, 1) & SummaryD.Table(:, 2) == UniCelTab(i, 2));
    cSummaryDTab = SummaryD.Table(cIds, :);
    cSummaryDResp = SummaryD.RespFunc(cIds, :);
    cSummaryDTprResp = SummaryD.TprRespFunc(cIds, :);
    cSummaryDCluFeat = SummaryD.ClusterFeatures(cIds, :);
    UniROI = unique(cSummaryDTab(:, Name2ColNum('ROIpair')));
    nUniROI = length(UniROI);
    CelTable= NaN(nUniROI, 15);
    RespFunc = NaN(nUniROI, RespFuncL);
    TprRespFunc = NaN(nUniROI, TprRespFuncL);
    CluFeatFunc = NaN(nUniROI, CluFeatL);
    for j = 1:nUniROI
        eIds = find(cSummaryDTab(:, Name2ColNum('ROIpair')) == UniROI(j));
        CelTable(j, :) = [cSummaryDTab(1, [1 2 4]), UniROI(j), mean(cSummaryDTab(eIds, Name2ColNum('CS')), 1),...
            median(cSummaryDTab(eIds, Name2ColNum('CS')), 1), std(cSummaryDTab(eIds, Name2ColNum('CS')), [], 1),...
            mean(cSummaryDTab(eIds, Name2ColNum('QL')), 1), median(cSummaryDTab(eIds, Name2ColNum('QL')), 1),...
            std(cSummaryDTab(eIds, Name2ColNum('QL')), [], 1), mean(cSummaryDTab(eIds, Name2ColNum('TS')), 1),...
            median(cSummaryDTab(eIds, Name2ColNum('TS')), 1), std(cSummaryDTab(eIds, Name2ColNum('TS')), [], 1),...
            mean(cSummaryDTab(eIds, Name2ColNum('RawIntensityMean')), 1),...
            mean(cSummaryDTab(eIds, Name2ColNum('RawIntensityMedian')), 1)];
        RespFunc(j, :) = mean(cSummaryDResp(eIds, :), 1);
        TprRespFunc(j, :) = mean(cSummaryDTprResp(eIds, :), 1);
        CluFeatFunc(j, :) = mean(cSummaryDCluFeat(eIds, :), 1);
    end
    SummaryD.CelTable = [SummaryD.CelTable; CelTable];
    SummaryD.CelRespf = [SummaryD.CelRespf; RespFunc];
    SummaryD.CelTprRespf = [SummaryD.CelTprRespf; TprRespFunc];
    SummaryD.CelClusterFeature = [SummaryD.CelClusterFeature; CluFeatFunc];
end
CelTabLabels = {'Day', 'Region', 'Type', 'ROI', 'CSavg', 'CSmed', 'CSstd', 'QLavg', 'QLmed', 'QLstd',...
    'TSavg', 'TSmed', 'TSstd', 'RawAmpMeanavg', 'RawAmpMedianavg'};
assert(length(CelTabLabels) == size(SummaryD.CelTable, 2));
Collab2num = @(in) find(cellfun(@(x) strcmpi(x, in), CelTabLabels));

%%
SaveFileName = './Results/CenterSurround_Reliability_individualBC_122922.mat';
save(SaveFileName, 'SummaryD', 'Dmt', 'CelTabLabels', 'Collab2num', 'TabLabels', 'Name2ColNum');
