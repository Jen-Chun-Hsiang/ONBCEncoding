
%%
X = Dfull([(1:nBC)'; remids(:)+nBC], :);
options.EstimationType = 'surround';
[SimpleCSRio, Curves] = estimateCenterSurround(X, options);
options.EstimationType = 'linearity';
[SimpleONOFFRio, Curves] = estimateCenterSurround(X, options);
% SimpleCSRio(SimpleCSRio<-0.4) = -0.4;
options.EstimationType = 'transience';
[SimpleTemporalRio, TemporalCurves] = estimateCenterSurround(X, options);
%%
tabd = Dtab(remids(:)+nBC, :);
depd = Ddepth(remids(:));
typ = [PlexusTypeGround(:, 4); PlexusTypeD(remids, 3)];
Pairdist = squareform(pdist(X, DistanceType));
mdstype = 2;

switch mdstype
    case 1
        [Y,stress] = mdscale(Pairdist, 20);
        r_recont = nan(1, size(Y, 2));
        for i = 1:size(Y, 2)
            Pairdist_reconstruct = squareform(pdist(Y(:,1:i), 'euclidean'));
            r_recont(i) = corr(Pairdist(triu(ones(size(Pairdist)), 1)==1),...
                Pairdist_reconstruct(triu(ones(size(Pairdist_reconstruct)), 1)==1)).^2;
        end
    case 2
        opt = statset('MaxIter', 2000);
        weightmask = reshape(~isoutlier(Pairdist(:)), size(Pairdist, 1), []); % ; Pairdist<1
        [Y,stress, disparity] = mdscale(Pairdist,16, 'criterion', 'sstress', 'weights', weightmask, 'start', 'random', 'Options', opt);
        r_recont = nan(1, size(Y, 2));
        for i = 1:size(Y, 2)
            Pairdist_reconstruct = squareform(pdist(Y(:,1:i), 'euclidean'));%
            r_recont(i) = corr(Pairdist(triu(ones(size(Pairdist)), 1)==1 & weightmask),...
                Pairdist_reconstruct(triu(ones(size(Pairdist_reconstruct)), 1)==1& weightmask)).^2;
        end
end
%% 
violationLevel = testTriangleInequality(Pairdist);

%%
figure; hold on
x = 1:size(Y, 2);
y = [r_recont(1), diff(r_recont)];
h = bar( x(1:16), y(1:16));
h.EdgeColor = 'w';
h.FaceColor = 0.3*ones(1, 3);
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylabel('Variance explained');

% yyaxis right
x = 1:size(Y, 2);
y = r_recont;
plot(x(1:16), y(1:16), 'b');
% yticks(0.4:0.2:1);
% yticklabels({'0.4', '', '0.8', '', '1'});
box off
xlabel('MDS components');
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
%%
keyboard;
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure4_SpotProjection_ExplainedVariance', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%% Show type distribution in the space
close all
figure;
MarkSize = 20;
subplot(2, 4, 1); hold on
ndot = size(Y, 1);
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(nType);
for i = 1:ndot
    plot3(Y(i, 1), Y(i, 2), Y(i, 3), '.', 'color', Colors(typ(i), :), 'markersize', MarkSize);
end
for i = 1:nType
    text(i*0.27-1.3,0.85,BCTypelabels{i},'Color',Colors(i, :),'FontSize',11);
end
ylabel('MDS component 2');
xlabel('MDS component 1');

for j = 1:nType
    subplot(2, 4, 1+j); hold on
    tids = find(typ == j);
    scatter(Y(:, 1), Y(:, 2), MarkSize, 0.7*ones(1, 3), 'filled');
    scatter(Y(tids, 1), Y(tids, 2), MarkSize, Colors(j, :), 'filled');
    text(0.5,0.75,BCTypelabels{j},'Color',Colors(j, :),'FontSize',15);
    ylabel('MDS component 2');
    xlabel('MDS component 1');
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_CellTypeDistribution', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
% close all
figure; 
MarkSize = 20;
ndot = size(Y, 1);
tits = {'Surround strength', 'Transience', 'Linearity'};
mY = median(Y, 1);
% scalefac = 0.5./sqrt(sum(sumv(:, 1:2).^2, 2));
for j = 1:3
    subplot(1, 3, j); hold on
    switch j
        case 1
            [OutColor, OutIds] = GroupScatterColor(SimpleCSRio, parula(256));
        case 2
            [OutColor, OutIds] = GroupScatterColor(SimpleTemporalRio, parula(256));
        case 3
            [OutColor, OutIds] = GroupScatterColor(SimpleONOFFRio, parula(256));
    end
    nids = length(OutIds);
    for i = 1:nids
        scatter(Y(OutIds(i), 1), Y(OutIds(i), 2), MarkSize, OutColor(i, :), 'filled');
%         plot3(Y(OutIds(i), 1), Y(OutIds(i), 2), Y(OutIds(i), 3), '.', 'color', OutColor(i, :), 'markersize', MarkSize);
    end
    ylabel('MDS component 2');
    xlabel('MDS component 1');
    title(tits{j});
%     quiver3(mY(1), mY(2), mY(3), scalefac(j)*sumv(j, 1), scalefac(j)*sumv(j, 2), scalefac(j)*sumv(j, 3),...
%         'Color', [1 1 0], 'LineWidth', 2);
%     quiver3(mY(1), mY(2), mY(3), -scalefac(j)*sumv(j, 1), -scalefac(j)*sumv(j, 2), -scalefac(j)*sumv(j, 3),...
%         'Color', [0 0 1], 'LineWidth', 2);
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_FeatureProjection_vector', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
