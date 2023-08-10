%% Correlation based difference
nIter = 2000;
fea = Dfull(1:nBC, :);
typ = dtypes;
TypeDiff_val = nan(nType, nType);
TypeDiff_p = nan(nType, nType);
TypeDiff_map_p = nan(1, size(fea, 2), nType*nType);
DiffMap_ind = nan(nType*nType, 2);
BCMap = nan(nType, size(fea, 2));
BCMap_orig = nan(nType, size(fea, 2));
rng('shuffle');
Count = 1;
for i = 1:nType
    for j = 1:nType
        if j>=i
            cidis = find(typ==i);
            cidjs = find(typ==j);
            ncidi = length(cidis);
            ncidj = length(cidjs);
            if i == j
                BCMap(i, :) = mean(fea(cidis, :), 1);
                BCMap_orig(i, :) = mean(fea(Dtype==BCTypes(i), :), 1);
            end
            %             diff_true = meanR(corrcoef(mean(fea(cidis, :), 1)', mean(fea(cidjs, :), 1)').^2);
            diff_true = pdist([mean(fea(cidis, :), 1); mean(fea(cidjs, :), 1)], 'correlation');
            BCMap_diff_true = mean(fea(cidis, :), 1)-mean(fea(cidjs, :), 1);
            BCMap_diff_sample = nan(nIter, size(BCMap_diff_true, 2));
            TypeDiff_val(i, j) = diff_true;
            diff_sample = nan(nIter, 1);
            for k = 1:nIter
                sid_i = randsample([cidis; cidjs], ncidi, true);
                sid_j = randsample([cidis; cidjs], ncidj, true);
                %                 diff_sample(k) = meanR(corrcoef(mean(fea(sid_i, :), 1)', mean(fea(sid_j, :), 1)').^2);
                diff_sample(k) = pdist([mean(fea(sid_i, :), 1); mean(fea(sid_j, :), 1)], 'correlation');
                BCMap_diff_sample(k, :) = mean(fea(sid_i, :), 1)-mean(fea(sid_j, :), 1);
            end
            
            TypeDiff_p(i, j) = mean(diff_sample > diff_true);
            TypeDiff_map_p(:, :, Count) = 1-mean(BCMap_diff_sample.^2 < BCMap_diff_true.^2, 1);
            %             if j ~= i
            %                 keyboard;
            %             end
            DiffMap_ind(Count, :) = [i, j];
            Count = Count + 1;
        else
            TypeDiff_val(i, j) = TypeDiff_val(j, i);
            TypeDiff_p(i, j) = TypeDiff_p(j, i);
            TypeDiff_map_p(:, :, Count) = TypeDiff_map_p(:, :, DiffMap_ind(:, 1)==j & DiffMap_ind(:, 2)==i);
            DiffMap_ind(Count, :) = [i, j];
            Count = Count + 1;
        end
        clc
        fprintf('progress... %d/%d, %d/%d (%d, %d) \n', i, nType, j, nType, ncidi, ncidj);
    end
end
%%
% close all
% pcorr = TypeDiff_p;
% pcorr(eye(size(TypeDiff_p, 1))==1) = nan;
[~, ~, ~, pcorr] = fdr_bh(TypeDiff_p(:));
pcorr = reshape(pcorr, size(TypeDiff_p));
Targetp = pcorr;
figure;
subplot(1, 2, 1);
imagesc(TypeDiff_val); colorbar
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title('Difference');
ytic = 1:7;

%
enlargefac = 10;
c = imresize(log10(Targetp), enlargefac, 'nearest');
subplot(1, 2, 2);
imagesc(c'); colorbar; hold on
visboundaries(c' < log10(0.05));
% box off
xticks((ytic-0.5)*enlargefac);
xticklabels(BCTypelabels);
yticks((ytic-0.5)*enlargefac);
yticklabels(BCTypelabels);
title('P value (log10)');
% subplot(1, 3, 3);
% imagesc(Targetp<0.05); colorbar
% xticks(1:nType);
% xticklabels(BCTypelabels);
% yticks(1:nType);
% yticklabels(BCTypelabels);
% title('Signficant (p <0.05)');

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure3_MorphologyIdentifiedBCDifference', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
% A = ones(nBC);
% A(eye(nBC) == 1) = 0;
A = ones(nType);
A(eye(nType) == 1) = 0;
B = graph(A);
for i = 1:length(B.Edges.EndNodes)
    B.Edges.Weight(i) = TypeDiff_val(B.Edges.EndNodes(i, 1), B.Edges.EndNodes(i, 2));
end
pos = [1 1.2; 2, 1.1; 3, 1.05; 4, 1; 1 0; 2 0.2; 3 -0.5];

figure;
ReverseWeight = (1.1- B.Edges.Weight./max(B.Edges.Weight));
% ReverseWeight = B.Edges.Weight./max(B.Edges.Weight);
p = plot(B,'Layout','layered','Sources',1:4 ,'Sinks',5:7);
% p = plot(B,'Layout','layered','Sources',1:4 ,'Sinks',5:7, 'XData', pos(:, 1), 'YData', pos(:,2));
% p.LineWidth = ReverseWeight*7./max(ReverseWeight);
p.LineWidth = ReverseWeight*6./max(ReverseWeight);
% p.LineWidth = 3;
p.EdgeColor = 1-(ReverseWeight./max(ReverseWeight))*ones(1, 3);
p.NodeLabel = BCTypelabels;
p.NodeColor = lines(nType);
p.MarkerSize = 12;
p.NodeFontSize = 12;
box off
axis off

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure3_MorphologyIdentifiedBCDifference_graph', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
ReduceFac = 1; % 1: no reduction
IPLDepthProfileType = 2;% 1:our confocal 2: EM
IsNormalized = 1;%
PathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BC diversity morphology data\';
switch IPLDepthProfileType
    case 1 % from our confocal data
        FileName = 'BCTypeIPLProfile.mat';
        IPL = load([PathName FileName], 'IPLdepth', 'IPLprofile', 'DispBC', 'TypeLabels');
        assert(sum(IPL.DispBC-BCTypes) == 0);
        IPLd = nan(size(IPL.IPLprofile));
        for i = 1:nType
            IPLd(i, :) = IPL.IPLprofile(IPL.DispBC == BCTypes(i), :);
        end
        IPLx = IPL.IPLdepth;
    case 2 % from other EM reconstruction data
        for i = 1:nType
            FileName = sprintf('EMBC_IPLProfile_%d.mat', BCTypes(i));
            IPL = load([PathName FileName], 'IPLprofile', 'DepthVec');
            if i == 1
                IPLd = nan(nType, length(IPL.DepthVec));
            end
            IPLd(i, :) = IPL.IPLprofile;
        end
        IPLx = IPL.DepthVec/100;
end

%%
IPLd = smoothdata(IPLd, 2, 'movmean', round(0.1*100));
% Adjust probability to set boundary
IPLd = IPLd - quantile(IPLd, 0.25, 2);
IPLd(IPLd<0) = 0;
v = corr(IPLd');
%%
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(7);
figure;
ax1 = subplot(1, 2, 1);
imagesc(v); colorbar;
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title(sprintf('max:%0.3G, min:%0.3G', max(v(:)), min(v(:))));
box off
ax2 = subplot(1, 2, 2);
imagesc((1:7)); colorbar;
colormap(ax1, 'gray');
colormap(ax2, Colors);
box off
axis off

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure1_EMBCType_DepthCorr', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
IPLd = smoothdata(IPLd, 2, 'movmean', round(0.1*100));
% Adjust probability to set boundary
IPLd = IPLd - quantile(IPLd, 0.25, 2);
IPLd(IPLd<0) = 0;
if ~IsNormalized
    cIPLd = sum(IPLd, 1);
    IPLd = IPLd./cIPLd;
else
    IPLd = IPLd ./ max(IPLd, [], 2);
    IPLd = IPLd./sum(IPLd, 1);
end
IPLd(:, IPLx< 0.47, :) = 0;
IPLd(:, IPLx> 0.95, :) = 0;
IPLd = IPLd*ReduceFac+(1-ReduceFac);
figure; plot(repmat(IPLx(:), 1, nType), IPLd');
legend(BCTypelabels);
box off
xlabel('Depth');
xlim([0.3 1]);
ylabel('Weight');
title(sprintf('Depth%d-Profile%d-Normal%d', IsConsideringDepth, IPLDepthProfileType, IsNormalized));

%%
SliceResp_type = Dtype(Dtype ~= 81);
nBC = size(SliceResp_type, 1);
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
X = Dfull([(1:nBC)'; remids(:)+nBC], :);
nP = size(X, 1);
DistanceType = 'correlation';
Pairdist = squareform(pdist(X, DistanceType));
Pairdist(eye(nP)== 1) = nan;
nP = size(X, 1);
pd = fitdist(min(Pairdist, [], 2),'gamma');
distprob = 1-cdf(pd, Pairdist);
depthprob = IPL2depthprob(IPLd);
% depthprob = depthprob./sum(depthprob, 2, 'omitnan');
PlexusTypeGround = nan(nBC, 2);
awd = nan(nBC, nBC);
for i = 1:nBC
    depthweight =  depthprob(dtypes(i),  dtypes);
    [PlexusTypeGround(i, 1), PlexusTypeGround(i, 2)] = max(distprob(i, 1:nBC),...
        [], 'omitnan');
    [PlexusTypeGround(i, 3), PlexusTypeGround(i, 4)] = max(distprob(i, 1:nBC).*depthweight,...
        [], 'omitnan');
    awd(i, :) = depthweight;
end
PlexusTypeGround(:, 5) = dtypes(PlexusTypeGround(:, 2));
PlexusTypeGround(:, 6) = dtypes(PlexusTypeGround(:, 4));
PlexusTypeGround(:, 7) = dtypes;
LeaveOneOutTest = nan(nType, 3);
classerrors = nan(nType, nType);
classerrors_depth = nan(nType, nType);
for i = 1:nType
    cids = PlexusTypeGround(:, 7) == i;
    LeaveOneOutTest(i, :) = [mean(PlexusTypeGround(cids, 5)==PlexusTypeGround(cids, 7)),...
        mean(PlexusTypeGround(cids, 6)==PlexusTypeGround(cids, 7)) (1/nType)];
    c = PlexusTypeGround(cids, 5)==(1:7);
    classerrors(i, :) = mean(c, 1, 'omitnan');
    c = PlexusTypeGround(cids, 6)==(1:7);
    classerrors_depth(i, :) = mean(c, 1, 'omitnan');
end
%%
close all
figure;
subplot(1, 3, 1); hold on
plot(LeaveOneOutTest(:, 1), 'Color', [74 118 253]/255);
plot(LeaveOneOutTest(:, 2), 'k');
plot(LeaveOneOutTest(:, 3), '--', 'Color', 0.3*ones(1, 3));
xticks(1:7);
xticklabels({'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'});
box off
legend({'No depth', 'Include depth'});
ylabel('Leave one out accuracy');
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylim([0 1]);
% title(sprintf('%s (k=%d)', DistanceType, nk));
subplot(1, 3, 2);
imagesc(classerrors, [0 1]);colorbar;
colormap(gray);
xticks(1:7);
xticklabels({'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'});
yticks(1:7);
yticklabels({'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'});
xlabel('Classified');
yticks(1:7);
title('Classification error');
subplot(1, 3, 3);
imagesc(classerrors_depth, [0 1]);colorbar;
xticks(1:7);
xticklabels({'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'});
yticks(1:7);
yticklabels({'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'});
xlabel('Classified');
yticks(1:7);
title('Classification error with depth');
%%
keyboard;
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure3_MorphologyIdentifiedBCDifference_classifyerror', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
Fz = 9.47;
StmDur = 3; % in second
BasDur = 0.5;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm;
t = round(linspace(-BasDurFrm/Fz, StmDurFrm/Fz, length(1:1:SpnDurFrm))*100)/100;
t = t(4:33);
%%
figure; 
nDmt = 7;
enlargefac = 10;
Dmt = [10 25 50 100 200 300 400];
nDmt = length(Dmt);
ytic = 1:7;
yticlab = num2cell(Dmt);
xtic = interp1(t, 1:30, [0 1 2 3], 'linear', 'extrap');
xticlab = {'0', '1', '2', '3'};
for i = 1:nType
    subplot(1, nType, i);
    c = BCMap(i, :);
    c = reshape(c, [], nDmt);
    c = imresize(c, enlargefac, 'nearest');
    imagesc(flipud(c'), [-0.5 1]);
    box off
    title(BCTypelabels{i});
    
    yticks((ytic-0.5)*enlargefac);
    yticklabels(fliplr(yticlab));
    ylabel('Radius (um)');
    if i == nType
        xticks((xtic-0.5)*enlargefac);
        xticklabels(xticlab);
        xlabel('Time (s)');
    else
        
        h = gca;
        h.XAxis.Visible = 'off';
    end
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure3_MorphologyIdentifiedBCDifference_HeatmapResponse', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);

%%
figure; 
c = BCMap(1, :);
c = reshape(c, [], nDmt);
t =  linspace(-0.5, 3, 33);
t = t(4:33);
s = fea(typ==1, :);
s = std(s, [], 1, 'omitnan')/sqrt(sum(typ==1));
s = reshape(s, [], nDmt);
Colors = lines(7);
for i = 1:nDmt
    subplot(1, nDmt, i); hold on
%     plot(t, c(:, i));
    shadePlot(t, c(:, i), s(:, i), Colors(1, :), 0.3, 2);
    plot([0 0 ], [-0.2, 1], '--', 'Color', 0.4*ones(1, 3));
    rectangle('Position', [-0.2, 1, 0.2, 0.05], 'FaceColor', 0.5*ones(1, 3), 'EdgeColor', 0.5*ones(1, 3));
    rectangle('Position', [0, 1, 1.5, 0.05], 'FaceColor', [253 237 74]/255, 'EdgeColor', [253 237 74]/255);
    rectangle('Position', [1.5, 1, 1.5, 0.05], 'FaceColor', zeros(1, 3), 'EdgeColor', zeros(1, 3));
    ylim([-0.2 1.05]);
    box off
    h = gca;
    h.YAxis.Visible = 'off';
    h.XAxis.Visible = 'off';
    title(sprintf('%d (um)', Dmt(i)));
    if i == 1
        plot([0.5 1.5], [0.7 0.7], 'k');
        plot([0.5 0.5], [0.7 0.9], 'k');
    end
    xlim([min(t(:)) max(t(:))]);
end

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure3_MorphologyIdentifiedBCDifference_SpotResponseTrace_BC5o', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);