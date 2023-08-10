clear; close all; clc;
%% Load Excel Experimental Document
% FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTypeDendriteAxonData.xlsx';
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\BCTypeDendriteAxonData.xlsx'; 
[DocNum, DocStr] = xlsread(FilNam,'Sheet1', 'A:F');
ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));
DocStr = DocStr(1, :);

%%
BCTypes = [50 51 57 58 6 7 89];
nType = length(BCTypes);

%%
x = [];
y = [];
for i = 1:nType
    cids = DocNum(:, ColName2Ind('Type')) == BCTypes(i);
    v = DocNum(cids, ColName2Ind('Dendrite'));
    v = v(~isnan(v));
    y = [y; v(:)];
    x = [x; i*ones(length(v), 1)];
end

%%
x2 = [];
y2 = [];
for i = 1:nType
    cids = DocNum(:, ColName2Ind('Type')) == BCTypes(i);
    v = DocNum(cids, ColName2Ind('Axon'));
    v = v(~isnan(v));
    y2 = [y2; v(:)];
    x2 = [x2; i*ones(length(v), 1)];
end
%%
x3 = [];
y3 = [];
for i = 1:nType
    cids = DocNum(:, ColName2Ind('Type')) == BCTypes(i);
    v = DocNum(cids, ColName2Ind('cone-contacted'));
    v = v(~isnan(v));
    y3 = [y3; v(:)];
    x3 = [x3; i*ones(length(v), 1)];
end
%%
mksize = 30;
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
close all
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(nType);
figure; subplot(1, 3, 1); hold on
swarmchart(x, y, mksize, Colors(x, :), 'filled', 'XJitter', 'density');
for i = 1:nType
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y(x == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(200:200:600);
ylim([0 700]);
ylabel('Dendritic territory');
subplot(1, 3, 2); hold on
swarmchart(x2, y2, mksize, Colors(x2, :), 'filled');
for i = 1:nType
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y2(x2 == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(200:200:600);
ylim([0 700]);
ylabel('Axonal territory');
xlim([0.5 nType+1]);

subplot(1, 3, 3); hold on
rng('shuffle');
cids = randperm(length(x3));
swarmchart(x3(cids), y3(cids), mksize, Colors(x3(cids), :), 'filled', 'XJitter', 'density', 'XJitterWidth', 0.3);
for i = [1 2]
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y3(x3 == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
xticks(1:nType);
xticklabels(BCTypelabels(1:2));
yticks(0:4:8);
ylim([0 8]);
ylabel('# cone contacted');
xlim([0.5 2+0.5]);
%%
keyboard;
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure1_MorphologyQuantification_Axon&Dend_Upper', SaveFolder);
FleNam = sprintf('%sFigure1_MorphologyQuantification_Axon&Dend_Lower', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
Colors = dazure(nType);
figure; subplot(1, 2, 1); hold on
swarmchart(x, y, mksize, Colors(x, :), 'filled', 'XJitter', 'density');
for i = 1:nType
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y(x == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
xticks(1:nType);
xticklabels(BCTypelabels);
yticks([1000 2000]);
yticklabels({'1000', '2000'});
ylim([0 3000]);
ylabel('Dendritic territory');
subplot(1, 2, 2); hold on
swarmchart(x2, y2, mksize, Colors(x2, :), 'filled');
for i = 1:nType
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*mean(y2(x2 == i)), 'Color', Colors(i, :), 'LineWidth', 2);
end
xticks(1:nType);
xticklabels(BCTypelabels);
ylabel('Axonal territory');
yticks([1000 2000]);
yticklabels({'1000', '2000'});
ylim([0 3000]);
xlim([0.5 nType+1]);
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure1_MorphologyQuantification_Axon&Dend_Top', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);