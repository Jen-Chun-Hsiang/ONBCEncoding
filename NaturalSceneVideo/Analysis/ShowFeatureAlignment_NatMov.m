% show the angle distribution
np = size(Y, 1);
ndim = 3;
rng('shuffle');
ChooseClip = 1:11;
%% for axes alignment
cids = Qual < qualthr;
%% Surround 
nIter = 1000;
gSize = 90;
% 
a = CSctr;
a(cids) = nan;
b = BSctr;
b(cids) = nan;
targetv = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
targetv = targetv(:);
pairids = nchoosek(1:np,2);
npair = size(pairids, 1);
pairv_s = nan(npair, ndim);
for i = 1:npair
    v1 = Y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_s(i, :) = v1*w1;
end
suma_s = nan(nIter, 1);
bsv_s = mean(pairv_s, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_s(i) = angle_between_vectors(mean(pairv_s(ids, :), 1, 'omitnan'),...
        bsv_s);
end
%%
% Transience
a = TSHctr;
a(cids) = nan;
b = TSLctr;
b(cids) = nan;
targetv = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
targetv = targetv(:);
pairv_t = nan(npair, ndim);
for i = 1:npair
    v1 = Y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_t(i, :) = v1*w1;
end
suma_t = nan(nIter, 1);
bsv_t = mean(pairv_t, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_t(i) = angle_between_vectors(mean(pairv_t(ids, :), 1, 'omitnan'), bsv_t);
end
%% Coherent motion
a = MSAdctr;
a(cids) = nan;
targetv = log(median(1./a(:, ChooseClip), 2, 'omitnan'));
targetv = targetv(:);
pairv_m = nan(npair, ndim);
for i = 1:npair
    v1 = Y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_m(i, :) = v1*w1;
end
suma_m = nan(nIter, 1);
bsv_m = mean(pairv_m, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_m(i) = angle_between_vectors(mean(pairv_m(ids, :), 1, 'omitnan'),...
        bsv_m);
end

%% Nonlinearity
a = squeeze(fNLmvctr(:, 1, :));
a(squeeze(fNLmvctr(:, 2, :)).^2 < qualthr) = nan;
targetv = median(log10(a(:, ChooseClip)), 2, 'omitnan');
targetv = targetv(:);
rmids = isnan(targetv);
y = Y(~rmids, :);
pairids = nchoosek(1:size(y, 1),2);
npair = size(pairids, 1);
pairv_n = nan(npair, ndim);

for i = 1:npair
    v1 = y(pairids(i, 1), 1:ndim)-Y(pairids(i, 2), 1:ndim);
    w1 = targetv(pairids(i, 1))-targetv(pairids(i, 2));
    pairv_n(i, :) = v1*w1;
end
suma_n = nan(nIter, 1);
bsv_n = mean(pairv_n, 1, 'omitnan');
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_n(i) = angle_between_vectors(mean(pairv_n(ids, :), 1, 'omitnan'), bsv_n);
end
%% Between
suma_st = nan(nIter, 1);
suma_tm = nan(nIter, 1);
suma_sm = nan(nIter, 1);
for i = 1:nIter
    ids = randsample(1:npair, gSize);
    suma_st(i) = angle_between_vectors(mean(pairv_s(ids, :), 1, 'omitnan'),...
        mean(pairv_t(ids, :),1, 'omitnan'));
    suma_tm(i) = angle_between_vectors(mean(pairv_t(ids, :), 1, 'omitnan'),...
        mean(pairv_m(ids, :),1, 'omitnan'));
    suma_sm(i) = angle_between_vectors(mean(pairv_s(ids, :), 1, 'omitnan'),...
        mean(pairv_m(ids, :),1, 'omitnan'));
end

%%
hrange = 0:5:160;
hlims = [min(hrange) max(hrange)];
close all
figure; 
subplot(2, 4, 1)
h = histogram(suma_s, hrange);
h.Normalization = 'Probability';
box off
h.EdgeColor = 'w';
h.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.45]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('Surround strength');

subplot(2, 4, 2)
h2 = histogram(suma_t, hrange);
h2.Normalization = 'Probability';
box off
h2.EdgeColor = 'w';
h2.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.45]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('Transience');

subplot(2, 4, 3);hold on
h3 = histogram(suma_m, hrange);
h3.Normalization = 'Probability';
box off
h3.EdgeColor = 'w';
h3.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.45]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('Coherent motion');

subplot(2, 4, 4);hold on
h3 = histogram(suma_n, hrange);
h3.Normalization = 'Probability';
box off
h3.EdgeColor = 'w';
h3.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:40:160);
xticklabels({'0', '', '80', '', '160'});
ylim([0 0.45]);
yticks(0:0.2:0.4);
yticklabels({'0', '0.2', '0.4'});
xlabel('Angle');
ylabel('Prob.');
title('In-center contrast');

% COMPARISON

hrange = linspace(0, 120, 25);
hlims = [min(hrange) max(hrange)];
subplot(2, 4, 5);hold on
h3 = histogram(suma_st, hrange);
h3.Normalization = 'Probability';
box off
h3.EdgeColor = 'w';
h3.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:30:120);
xticklabels({'0', '', '60', '', '120'});
mval = mean(suma_st, 'omitnan');

switch ndim
    case 2
        plot(mval*ones(1, 2), [0 0.16], '--k');
        yticks(0:0.08:0.16);
        yticklabels({'0', '0.08', '0.16'});
        ylim([0 0.16]);
    case 3
        plot(mval*ones(1, 2), [0 0.2], '--k');
        yticks(0:0.1:0.2);
        yticklabels({'0', '0.1', '0.2'});
        ylim([0 0.2]);
end
xlabel('Angle');
ylabel('Prob.');
title(sprintf('Surround vs Transience %0.3G', mval));


subplot(2, 4, 6);hold on
h3 = histogram(suma_sm, hrange);
h3.Normalization = 'Probability';
box off
h3.EdgeColor = 'w';
h3.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:30:120);
xticklabels({'0', '', '60', '', '120'});
mval = mean(suma_sm, 'omitnan');

switch ndim
    case 2
        plot(mval*ones(1, 2), [0 0.16], '--k');
        yticks(0:0.08:0.16);
        yticklabels({'0', '0.08', '0.16'});
        ylim([0 0.16]);
    case 3
        plot(mval*ones(1, 2), [0 0.21], '--k');
        yticks(0:0.1:0.2);
        yticklabels({'0', '0.1', '0.2'});
        ylim([0 0.21]);
end
xlabel('Angle');
ylabel('Prob.');
title(sprintf('Surround vs coherent motion %0.3G', mval));

subplot(2, 4, 7);hold on
h3 = histogram(suma_tm, hrange);
h3.Normalization = 'Probability';
box off
h3.EdgeColor = 'w';
h3.FaceColor = 0.3*ones(1, 3);
xlim(hlims);
xticks(0:30:120);
xticklabels({'0', '', '60', '', '120'});
mval = mean(suma_tm, 'omitnan');

switch ndim
    case 2
        plot(mval*ones(1, 2), [0 0.3], '--k');
        yticks(0:0.1:0.3);
        yticklabels({'0', '0.1', '0.2', '0.3'});
        ylim([0 0.3]);
    case 3
        plot(mval*ones(1, 2), [0 0.24], '--k');
        yticks(0:0.1:0.2);
        yticklabels({'0', '0.1', '0.2'});
        ylim([0 0.24]);
end
xlabel('Angle');
ylabel('Prob.');
title(sprintf('Transience vs coherent motion %0.3G', mval));

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure6_NatMovieProjection_FeatureAlignment_Angles_ndim%d', SaveFolder, ndim);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);