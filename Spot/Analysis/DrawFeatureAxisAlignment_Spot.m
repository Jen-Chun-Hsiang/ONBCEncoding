
ndim = size(Y, 2);
% ndim = 3;
%%
x = SimpleCSRio;
ev_s = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    ev_s(i) = corr(x, y).^2;
end

x = SimpleTemporalRio;
ev_t = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    ev_t(i) = corr(x, y).^2;
end

%%
% close all
figure;
subplot(1, 2, 1); hold on
b = bar(ev_s);
b.EdgeColor = 'w';
b.FaceColor = 0.6*ones(1, 3);
ylim([0 0.6]);
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
xlim([0.5 16.5]);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
box off
title('Spatial contrast');
xlabel('MDS dimension');
ylabel('Explained variance');


subplot(1, 2, 2);
b = bar(ev_t);
b.EdgeColor = 'w';
b.FaceColor = 0.6*ones(1, 3);
ylim([0 0.6]);
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
xlim([0.5 16.5]);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
box off
title('Temporal contrast');
xlabel('MDS dimension');
ylabel('Explained variance');

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure4_SpotProjection_FeatureMDScomponents', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);