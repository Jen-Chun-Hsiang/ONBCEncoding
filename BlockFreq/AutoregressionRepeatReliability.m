%% Recalculate the repeat reliability and the shuffle value
Diameter = 150;
dids = find(ismember(ROITable(:, 8), Diameter));
ndids = length(dids);
nFz = size(DMat, 2);
Quals = nan(ndids, nFz);
rQuals = nan(ndids, nFz);
Qualtyp = nan(ndids, 1);
mPhis = nan(ndids, nFz);
stda = std(DMat(:));
rng('shuffle');
for i = 1:ndids
    clc
    fprintf('progress...%d/%d \n', i, ndids);
    Qualtyp(i) = ROITable(dids(i), 4);
    for j = 1:nFz
        a = squeeze(DMat(:, j, :, dids(i)));
        Quals(i, j) = estimateRepeatReliability(a);
        [sa, mPhis(i, j)] = autoregressionSim(a);
        rQuals(i, j) = estimateRepeatReliability(sa);
    end
end

%% Cell report
Celltypeids = [50 51 57 58 6 7 89];
a = unique(ROITable(ismember(ROITable(:, 8), Diameter) & ismember(ROITable(:, 4), Celltypeids), 1:2), 'rows')
%% Phi - the composition of the noise
Fzs = [Freqs];
Color = [240 45 108]/255;
close all
Passed = nan(7, 2);
figure;
for i = 1:6
    subplot(1, 6, i); hold on
    a = mPhis(:, i);
    [f,x,flo,fup] = ecdf(a,'Alpha',0.01,'Bounds','on');
    shadePlot(x, f, abs(flo-f), 'k', 0.5);
    
    ylim([0 1]);
    title(sprintf('%3G Hz', Fzs(i)));
    ylabel('Probability');
    xlabel('Noise composition');
end
%%
Fzs = [Freqs];
Color = [240 45 108]/255;
Celltypeids = [50 51 57 58 6 7 89];
% thr = quantile(sQuals(:, end), 0.95);
close all
Passed = nan(7, 2);
PassedTyp = nan(6, 7);
figure;
for i = 1:6
    thr = quantile(rQuals(:, i), 0.95);
    subplot(1, 6, i); hold on
    bins = linspace(0, 1, 25);
    a = Quals(:, i);
    [f,x,flo,fup] = ecdf(a,'Alpha',0.01,'Bounds','on');
    shadePlot(x, f, abs(flo-f), 'k', 0.5);
    [f,x,flo,fup] = ecdf(rQuals(:, i),'Alpha',0.01,'Bounds','on');
    shadePlot(x, f, abs(flo-f), Color, 0.5);
    Passed(i, :) = [mean(a>thr) mean(rQuals(:, i)>thr)];
    for j = 1:7
        a = Quals(Qualtyp == Celltypeids(j), i);
        PassedTyp(i,j) = mean(a>thr);
    end
    ylim([0 1]);
    xlim([0 1]);
    title(sprintf('%3G Hz', Fzs(i)));
    ylabel('Probability');
    xlabel('Repeat reliability (R2)');
    yticks(0:0.5:1);
    yticklabels({'0', '0.5', '1'});
end
%%
c = PassedTyp./max(PassedTyp, [], 1);
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_FrequencyResponseRepeatReliability_ECDF', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
figure; hold on 
plot(Passed(1:6, 1), '-ok', 'MarkerFaceColor',zeros(1, 3));
plot(Passed(1:6, 2), '-o', 'Color', Color, 'MarkerFaceColor',Color);
xticks(1:6);
xticklabels(num2cell(Fzs));
ylabel('Above chance probability');
xlabel('Frequency (Hz)');
legend({'ROI traces', 'Simulated ROI traces'})
ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0', '0.5', '1'});
xlim([1 6]);
xticks(1:6);
xticklabels({'0.5', '1', '2', '4', '8', '16'});
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure2_FrequencyResponseRepeatReliability_ChanceCurve', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
figure;
for i = 1:7
    subplot(1, 7, i);
    bins = linspace(-0.2, 1, 25);
    h = histogram(ROIQuality(:, i), bins);
    h.Normalization = 'Probability';
end

%% Test autoregression
time_series = squeeze(DMat(:, 7, :, 10));
a = time_series(1, :)';

n = length(a);
X = [ones(n-1, 1) a(1:end-1)];
y = a(2:end);
beta_hat = X \ y;
phi = beta_hat(2);

epsilon = std(a)*randn(n, 1);
b = zeros(size(a));
for t = 2:n
    b(t) = phi * b(t-1) + epsilon(t);
end

figure; hold on
plot(a, 'k');
plot(b, 'r');