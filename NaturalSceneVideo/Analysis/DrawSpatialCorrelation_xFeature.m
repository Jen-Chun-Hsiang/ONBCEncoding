clear clc close all
%%
Folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\SpatialCorrection_110123';
LM = load([Folder '\SpatialCorrelation_luminance_Intercept.mat']);
SR = load([Folder '\SpatialCorrelation_Surround_Intercept.mat']);
TS = load([Folder '\SpatialCorrelation_Transience_Intercept.mat']);
MT = load([Folder '\SpatialCorrelation_Motion_Intercept.mat']);
IC = load([Folder '\SpatialCorrelation_In-centercontrast_Intercept.mat']);

%%
QT95ints = [LM.QT95int SR.QT95int TS.QT95int MT.QT95int IC.QT95int];
mksize = 15;
figure; hold on
x = repmat(1:5, 11, 1);
x = x(:);
y = QT95ints(:);
swarmchart(x(:), y, mksize, 'k', 'filled', 'XJitter', 'density');
for i = 1:5
    plot(0.7*[-0.5 0.5]+i, ones(1, 2)*median(y(x == i)), 'Color', 'k', 'LineWidth', 2);
end
xticks(1:5);
xticklabels({'Luminance', 'Spatial', 'Temporal', 'Motion', 'In-center'});
ylim([0.6 1]);
yticks(0.6:0.1:1);
yticklabels({'0.6', '0.7', '0.8', '0.9', '1'});
ylabel('Spatial corr. at 95%');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sSupFigure6_FeatureSpatialCorrelation_Summary', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);