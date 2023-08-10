%%
options.EstimationType = 'surround';
[SimpleCSRio, Curves] = estimateCenterSurround(SliceResp_full, options);
options.EstimationType = 'surroundreduction';
[CSreduction, Curves] = estimateCenterSurround(SliceResp_full, options);
options.EstimationType = 'linearity';
[SimpleONOFFRio, Curves] = estimateCenterSurround(SliceResp_full, options);
% SimpleCSRio(SimpleCSRio<-0.4) = -0.4;
options.EstimationType = 'transience';
[SimpleTemporalRio, TemporalCurves] = estimateCenterSurround(SliceResp_full, options);
%% Transience
x = SimpleTemporalRio;
[~, g] = max(SliceResp_type == BCTypes, [], 2);
[p,tbl,stats] = kruskalwallis(x, g);
c = multcompare(stats);
%% compare two block separately
[p,tbl,stats] = kruskalwallis(x(ismember(g, 1:4)), g(ismember(g, 1:4)));
c = multcompare(stats);

%%
[p,tbl,stats] = kruskalwallis(x(ismember(g, 5:7)), g(ismember(g, 5:7)));
c = multcompare(stats);

%% Surround strength
x = SimpleCSRio;
[~, g] = max(SliceResp_type == BCTypes, [], 2);
[p,tbl,stats] = kruskalwallis(x, g);
c = multcompare(stats);

%% compare two block separately
[p,tbl,stats] = kruskalwallis(x(ismember(g, 1:4)), g(ismember(g, 1:4)));
c = multcompare(stats);

%%
[p,tbl,stats] = kruskalwallis(x(ismember(g, 5:7)), g(ismember(g, 5:7)));
c = multcompare(stats);

%% test on two groups (above or below chatband
clc
x = SimpleTemporalRio;
[~, g] = max(SliceResp_type == BCTypes, [], 2);
% x(g==3) = [];
% g(g==3) = [];
g = ismember(g, [5 6 7]);
[p,tbl,stats] = ranksum(x(g==0), x(g==1), 'tail', 'right')
% c = multcompare(stats);
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};

dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(7);
mksize = 25;
%%
y = SimpleCSRio;
[~, g] = max(SliceResp_type == BCTypes, [], 2);
figure; subplot(2, 1, 1); hold on
swarmchart(g, y, mksize, Colors(g, :), 'filled', 'XJitter', 'density');
for i = 1:7
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y(g == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
% xticks(1:nType);
% xticklabels(BCTypelabels);
set(gca,'xtick',[])
ylabel('Spatial contrast');

 subplot(2, 1, 2); hold on
 y = SimpleTemporalRio;
swarmchart(g, y, mksize, Colors(g, :), 'filled', 'XJitter', 'density');
for i = 1:7
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y(g == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
xticks(1:nType);
xticklabels(BCTypelabels);
ylabel('Temporal contrast');
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure3_SpotContrastDifference_SwarmChart', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%% Correlation based difference
nIter = 2000;
fea = SliceResp_full;
[~, typ] = max(SliceResp_type == BCTypes, [], 2);
nType = length(BCTypes);
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
                BCMap_orig(i, :) = mean(fea(cidis, :), 1);
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
cTypeDiff_p = TypeDiff_p;
cTypeDiff_val = TypeDiff_val;
cTypeDiff_p(eye(size(cTypeDiff_p, 1))==1) = nan;
[~, ~, ~, pcorr] = fdr_bh(TypeDiff_p(:));
pcorr = reshape(pcorr, size(TypeDiff_p));
Targetp = pcorr;
figure;
subplot(1, 2, 1); hold on

imagesc(flipud(cTypeDiff_val)); colorbar
pv = nan(size(Targetp));
pv(Targetp<0.05) = 1;
pv(Targetp<0.01) = 2;
pv(Targetp<0.001) = 3;
for i = 1:size(Targetp, 1)
    for j = 1:size(Targetp, 2)
        if ~isnan(pv(i, j)) & j<i
            switch pv(i, j)
                case 1
                    signi = '*';
                case 2
                    signi = '**';
                case 3
                    signi = '***';
            end
            text(i, 8-j-0.1, signi, 'Color', 'r', 'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'middle', 'FontSize', 15);
        end
    end
end
xlim([0.5 7.5]);
ylim([0.5 7.5]);
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(fliplr(BCTypelabels));
colormap(gray);
title('Difference');



subplot(1, 2, 2);
ytic = 1:7;
enlargefac = 10;
c = imresize(log10(Targetp), enlargefac, 'nearest');
imagesc(c'); colorbar; hold on
visboundaries(c' < log10(0.05));
% box off
xticks((ytic-0.5)*enlargefac);
xticklabels(BCTypelabels);
yticks((ytic-0.5)*enlargefac);
yticklabels(BCTypelabels);
title('P value (log10)');

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure3_SpotResponseBCTypeDifference_Heatmap2', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
close all
figure;
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(7);
Fz = 9.47;
StmDur = 3; % in second
BasDur = 0.5;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm;
t = round(linspace(-BasDurFrm/Fz, StmDurFrm/Fz, length(1:1:SpnDurFrm))*100)/100;
t = t(4:33);
for j = 1:7
    c = BCMap(j, :);
    c = reshape(c, [], nDmt);
    s = fea(typ==j, :);
    s = std(s, [], 1, 'omitnan')/sqrt(sum(typ==j));
    s = reshape(s, [], nDmt);
    for i = 1:nDmt
        subplot(1, nDmt, i); hold on
        %     plot(t, c(:, i));
        shadePlot(t, c(:, i)-j+1, s(:, i), Colors(j, :), 0.3, 2);
        if i == 1
            plot([min(t(:)) max(t(:))], [0 0]-j+1, '--', 'Color', 0.4*ones(1, 3));
        end
        if j == 1
%             plot([0 0 ], [-6.2, 1], '--', 'Color', 0.4*ones(1, 3));
            rectangle('Position', [-0.2, 1, 0.2, 0.05], 'FaceColor', 0.5*ones(1, 3), 'EdgeColor', 0.5*ones(1, 3));
            rectangle('Position', [0, 1, 1.5, 0.05], 'FaceColor', [253 237 74]/255, 'EdgeColor', [253 237 74]/255);
            rectangle('Position', [1.5, 1, 1.5, 0.05], 'FaceColor', zeros(1, 3), 'EdgeColor', zeros(1, 3));
        end
        ylim([-6.3 1.05]);
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
end
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_SpotResponseTrace_All', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);