% show the angle distribution
np = size(Y, 1);
ndim = 3;
rng('shuffle');
%% Surround 
nIter = 1000;
gSize = 90;
targetv = SimpleCSRio;
pairids = nchoosek(1:np,2);
npair = size(pairids, 1);
pairv_s = nan(npair, ndim);
for i = 1:npair
    v1 = Y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_s(i, :) = v1*w1;
end
suma_s = nan(nIter, 1);
bsv_s = median(pairv_s, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_s(i) = angle_between_vectors(median(pairv_s(ids, :), 1, 'omitnan'), bsv_s);
end
%%
% Transience
targetv = SimpleTemporalRio;
pairv_t = nan(npair, ndim);
for i = 1:npair
    v1 = Y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_t(i, :) = v1*w1;
end
suma_t = nan(nIter, 1);
bsv_t = median(pairv_t, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_t(i) = angle_between_vectors(median(pairv_t(ids, :), 1, 'omitnan'), bsv_t);
end
%% linearity
targetv = SimpleONOFFRio;
pairv_n = nan(npair, ndim);
for i = 1:npair
    v1 = Y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_n(i, :) = v1*w1;
end
suma_n = nan(nIter, 1);
bsv_n = median(pairv_n, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_n(i) = angle_between_vectors(median(pairv_n(ids, :), 1, 'omitnan'), bsv_n);
end

%% Between
suma_st = nan(nIter, 1);
suma_sn = nan(nIter, 1);
suma_tn = nan(nIter, 1);
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_st(i) = angle_between_vectors(median(pairv_s(ids, :), 1, 'omitnan'),...
        median(pairv_t(ids, :),1, 'omitnan'));
    suma_sn(i) = angle_between_vectors(median(pairv_s(ids, :), 1, 'omitnan'),...
        median(pairv_n(ids, :),1));
    suma_tn(i) = angle_between_vectors(median(pairv_t(ids, :), 1, 'omitnan'),...
        median(pairv_n(ids, :),1, 'omitnan'));
end

%%
close all
figure; 
subplot(1, 3, 1);hold on
[f,x,flo,fup] = ecdf(suma_s, 'Function','survivor','Alpha',0.01,'Bounds','on');
shadePlot(x, f, abs(fup-f), 0.3*ones(1, 3));
% plot(x, f, 'k');
% plot(x, flo, ':b');
% plot(x, fup, ':b');
box off
xlim([0 180]);
xticks(0:90:180);
xticklabels({'0', '90', '180'});
ylim([0 1]);
% yticks(0:0.03:0.06);
% yticklabels({'0', '0.03', '0.06'});

subplot(1, 3, 2)
[f,x,flo,fup] = ecdf(suma_t, 'Function','survivor','Alpha',0.01,'Bounds','on');
shadePlot(x, f, abs(fup-f), 0.3*ones(1, 3));
xlim([0 180]);
xticks(0:90:180);
xticklabels({'0', '90', '180'});
ylim([0 1]);
% yticks(0:0.03:0.06);
% yticklabels({'0', '0.03', '0.06'});

subplot(1, 3, 3);hold on
h = histogram(suma_st, 0:6:180);
h.Normalization = 'Probability';
box off
h.EdgeColor = 'w';
h.FaceColor = 0.3*ones(1, 3);
xlim([-0 180]);
xticks(0:90:180);
xticklabels({'0', '90', '180'});
suma_st = suma_st
plot(median(suma_st, 'omitnan')*ones(1, 2), [0 0.005], '--k');
% yticks(0:0.03:0.06);
% yticklabels({'0', '0.03', '0.06'});
%%
hrange = 0:5:160;
hlims = [min(hrange) max(hrange)];
close all
figure; 
subplot(2, 3, 1)
h = histogram(suma_s, hrange);
h.Normalization = 'Probability';
box off
h.EdgeColor = 'w';
h.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.54]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('Surround strength');

subplot(2, 3, 2)
h2 = histogram(suma_t, hrange);
h2.Normalization = 'Probability';
box off
h2.EdgeColor = 'w';
h2.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.54]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('Transience');

subplot(2, 3, 3)
h2 = histogram(suma_n, hrange);
h2.Normalization = 'Probability';
box off
h2.EdgeColor = 'w';
h2.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.54]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('Nonlinearity');

% COMPARISON
hrange = linspace(0, 120, 25);
subplot(2, 3, 4);hold on
h12 = histogram(suma_st, hrange);
h12.Normalization = 'Probability';
box off
h12.EdgeColor = 'w';
h12.FaceColor = 0.3*ones(1, 3);
xlim([0 120]);
xticks(0:30:120);
xticklabels({'0', '', '60', '', '120'});
mval = median(suma_st, 'omitnan');
plot(mval*ones(1, 2), [0 0.14], '--k');
yticks(0:0.07:0.14);
yticklabels({'0', '0.07', '0.14'});
ylim([0 0.14]);
xlabel('Angle');
ylabel('Prob.');
title('Surround vs Transience');
title(sprintf('Surround vs Transience %0.3G', mval));

subplot(2, 3, 5);hold on
h12 = histogram(suma_sn, hrange);
h12.Normalization = 'Probability';
box off
h12.EdgeColor = 'w';
h12.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
plot(median(suma_sn, 'omitnan')*ones(1, 2), [0 0.2], '--k');
% yticks(0:0.1:0.2);
% yticklabels({'0', '0.1', '0.2'});
xlabel('Angle');
ylabel('Prob.');
title('Surround vs nonlinearity');

subplot(2, 3, 6);hold on
h12 = histogram(suma_tn, hrange);
h12.Normalization = 'Probability';
box off
h12.EdgeColor = 'w';
h12.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
plot(median(suma_tn, 'omitnan')*ones(1, 2), [0 0.2], '--k');
% yticks(0:0.1:0.2);
% yticklabels({'0', '0.1', '0.2'});
xlabel('Angle');
ylabel('Prob.');
title('Transience vs nonlinearity');

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure4_SpotProjection_FeatureAlignment_Angles_ndim%d', SaveFolder, ndim);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
