%% for axes alignment
cids = Qual < qualthr;
ndim = size(Y, 2);
% 
a = CSctr;
a(cids) = nan;
b = BSctr;
b(cids) = nan;
x = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
x = x-x';
x(eye(size(x, 1))==1) = nan;
x = x(:);
x(isnan(x)) = [];

ev_s = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    y = y-y';
    y(eye(size(y, 1))==1) = nan;
    y = y(:);
    y(isnan(y)) = [];
    ev_s(i) = corr(x, y).^2;
end
%
a = TSHctr;
a(cids) = nan;
b = TSLctr;
b(cids) = nan;
x = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
x = x-x';
x(eye(size(x, 1))==1) = nan;
x = x(:);
x(isnan(x)) = [];
ev_t = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    y = y-y';
    y(eye(size(y, 1))==1) = nan;
    y = y(:);
    y(isnan(y)) = [];
    ev_t(i) = corr(x, y).^2;
end
%
a = MSAdctr;
a(cids) = nan;
x = log(median(1./a(:, ChooseClip), 2, 'omitnan'));
x = x-x';
x(eye(size(x, 1))==1) = nan;
x = x(:);
x(isnan(x)) = [];
ev_c = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    y = y-y';
    y(eye(size(y, 1))==1) = nan;
    y = y(:);
    y(isnan(y)) = [];
    ev_c(i) = corr(x, y).^2;
end
%
a = squeeze(fNLmvctr(:, 1, :));
a(squeeze(fNLmvctr(:, 2, :)).^2 < qualthr) = nan;

x = median(log10(a(:, ChooseClip)), 2, 'omitnan');
rmids = isnan(x);
x(rmids) = [];
x = x-x';
x(eye(size(x, 1))==1) = nan;
x = x(:);
x(isnan(x)) = [];
ev_n = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    y(rmids) = [];
    y = y-y';
    y(eye(size(y, 1))==1) = nan;
    y = y(:);
    y(isnan(y)) = [];
    ev_n(i) = corr(x, y).^2;
end
%%
% close all
figure;
subplot(1, 4, 1); hold on

b = bar(ev_s);
b.EdgeColor = 'w';
b.FaceColor = 0.6*ones(1, 3);
ylim([0 0.7]);
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
xlim([0.5 16.5]);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
box off
title('Spatial contrast');
xlabel('MDS dimension');
ylabel('Explained variance');

subplot(1, 4, 2);
b = bar(ev_t);
b.EdgeColor = 'w';
b.FaceColor = 0.6*ones(1, 3);
ylim([0 0.7]);
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
xlim([0.5 16.5]);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
box off
title('Temporal contrast');
xlabel('MDS dimension');
ylabel('Explained variance');

subplot(1, 4, 3);
b = bar(ev_c);
b.EdgeColor = 'w';
b.FaceColor = 0.6*ones(1, 3);
ylim([0 0.7]);
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
xlim([0.5 16.5]);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
box off
title('Coherent motion');
xlabel('MDS dimension');
ylabel('Explained variance');

subplot(1, 4, 4);
b = bar(ev_n);
b.EdgeColor = 'w';
b.FaceColor = 0.6*ones(1, 3);
ylim([0 0.7]);
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
xlim([0.5 16.5]);
yticks(0:0.3:0.6);
yticklabels({'0', '0.3', '0.6'});
box off
title('Nonlinearity');
xlabel('MDS dimension');
ylabel('Explained variance');

%% for axes alignment (without minus) and run the above section to get the plot
cids = Qual < qualthr;
ndim = size(Y, 2);
% 
a = CSctr;
a(cids) = nan;
b = BSctr;
b(cids) = nan;
x = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
x = x(:);

ev_s = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    ev_s(i) = corr(x, y).^2;
end
%
a = TSHctr;
a(cids) = nan;
b = TSLctr;
b(cids) = nan;
x = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
x = x(:);
ev_t = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    ev_t(i) = corr(x, y).^2;
end
%
a = MSAdctr;
a(cids) = nan;
x = log(median(1./a(:, ChooseClip), 2, 'omitnan'));
ev_c = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    ev_c(i) = corr(x, y).^2;
end
%
a = squeeze(fNLmvctr(:, 1, :));
a(squeeze(fNLmvctr(:, 2, :)).^2 < qualthr) = nan;
x = median(log10(a(:, ChooseClip)), 2, 'omitnan');
rmids = isnan(x);
x(rmids) = [];
ev_n = nan(ndim, 1);
for i = 1:ndim
    y = Y(:, i);
    y(rmids) = [];
    ev_n(i) = corr(x, y).^2;
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_NatMovEncodingSpace_noise', SaveFolder);
FleNam = sprintf('%sFigure5_NatMovEncodingSpace_featureaxisalignment', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);