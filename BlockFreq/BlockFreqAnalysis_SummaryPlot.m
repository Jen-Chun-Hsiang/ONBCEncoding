addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');

BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCTypes = [50 51 57 58 6 7 89];
nType = length(BCTypes);
SpotSize = {[150 300],  800};
NumSpotSize = length(SpotSize);
T_detail = linspace(min(T(:)), max(T(:)), (length(T)-1)*5);
SinW = sin(2*pi*[Freqs 32]'*T_detail(:)');
SinW(:, T_detail(:)<0) = 0;
SinW(:, T_detail(:)>2) = 0;
AngS = [];
AngL = [];
typenums = nan(nType, NumDmt-1, NumSpotSize);
for k = 1:nType
    for i = 1:NumDmt-1
        for j = 1:NumSpotSize
            ExpIds = uniCel(:, 4) == BCTypes(k) & ismember(uniCel(:, 3),SpotSize{j});
            uniExp = unique(uniCel(ExpIds, 1:2), 'rows');
            NumExp = size(uniExp, 1);
            typenums(k, i, j) = NumExp;
            if i == 1 && j == 1 && k == 1
                f0mean = nan(20, NumDmt, NumSpotSize, nType);
                f1pow = nan(20, NumDmt, NumSpotSize, nType);
                f1pha = nan(20, NumDmt, NumSpotSize, nType);
                freqTrace_mean = nan(nType, size(CelMeanNrm, 2), NumDmt, NumSpotSize);
                freqTrace_sem = nan(nType, size(CelMeanNrm, 2), NumDmt, NumSpotSize);
            end
            cResp = [];
            cPow = [];
            cPha = [];
            for p = 1:NumExp
                ExpIds = uniCel(:, 4) == BCTypes(k) & ismember(uniCel(:, 3),SpotSize{j}) &...
                    uniCel(:, 1) == uniExp(p, 1) & uniCel(:, 2) == uniExp(p, 2);
                cResp = [cResp; squeeze(mean(CelMeanNrm(i, :, ExpIds), 3))];
                cPow = [cPow; squeeze(mean(mean(CellFreqPw(i, fsets{i}, ExpIds), 2, 'omitnan'), 3, 'omitnan'))];
                cPha = [cPha; meanAngle(CellFreqPh(i, fsets{i}, ExpIds))];
                if j == 1
                    AngS = [AngS; reshape(CellFreqPh(i, fsets{i}, ExpIds), [], 1)];
                else
                    AngL = [AngL; reshape(CellFreqPh(i, fsets{i}, ExpIds), [], 1)];
                end
            end
            freqTrace_mean(k, :, i, j) = mean(cResp, 1);
            freqTrace_sem(k, :, i, j) = std(cResp, [], 1)/sqrt(NumExp);
            f0mean(1:NumExp, i, j, k) = mean(cResp(:, T>0), 2);
            f1pow(1:NumExp, i, j, k) = mean(cPow, 2);
            f1pha(1:NumExp, i, j, k) = cPha;
        end
    end
end
%%
close all
figure; 
hold on
v = [reshape(f1pha(:, freqid, 1, :), [], 1) reshape(f1pha(:, freqid, 2, :), [], 1)];
v(any(isnan(v), 2), :) = [];
[h, p] = kstest2(v(:, 1), v(:, 2));

[f,x_values, flo, fup] = ecdf(v(:, 1));
shadePlot(x_values,f, (fup-f)/norminv(0.975), [144 39 142]/255);
% plot(x_values, flo, '--b');
% plot(x_values, fup, '--b');
[f,x_values, flo, fup] = ecdf(v(:, 2));
shadePlot(x_values,f, (fup-f)/norminv(0.975), [246 146 30]/255);
% K = plot(x_values,f,'r');
% plot(x_values, flo, '--r');
% plot(x_values, fup, '--r');
% set(J,'LineWidth',2);
% set(K,'LineWidth',2);
% legend([J K],'Center','Surround');
xlim([-180 180])
title('Paired difference');
xlabel('Phase (angle)');
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
xticks(-180:90:180);
xticklabels({'-180', '', '0', '', '180'});
ylabel('Cumulative Probability'); 
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_FrequencyPhaseDifference', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
[h, p, k] = kstest2(v(:, 1), v(:, 2));
%%
IsCorrGC6f = 1;
Sim = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Simulation\GCaMP6f\simf1powere.mat', 'samplefreqs', 'f1pow');
% Sim.f1pow = Sim.f1pow./sum(Sim.f1pow);
% freqPwf0 = nan(NumDmt-1, :, k);
freqPwf1_mean = nan(NumDmt-1, NumSpotSize, nType);
freqPwf1_std = nan(NumDmt-1, NumSpotSize, nType);
freqPwf1_diff_mean = nan(NumDmt-1, nType);
freqPwf1_diff_std = nan(NumDmt-1, nType);
freqPhf1_diff_mean = nan(NumDmt-1, nType);
freqPhf1_diff_std = nan(NumDmt-1, nType);
y = [];
yc = [];
ya = [];
gtyp = [];
gfrq = [];
gtypc = [];
gfrqc = [];
gsiz = [];
for k = 1:nType
    target = f1pow(:, 1:6, :, k);
    target = target./max(target, [], 2);
%     target = target./sum(target,  2);
    if IsCorrGC6f
        target = target(:, 1:6, :)./Sim.f1pow';
    end
    NumExp = sum(~isnan(target(:, 1, 1)));
    targetmean = squeeze(mean(log10(target(:, 1:6, :)), 1, 'omitnan'));
    targetstd = squeeze(std(log10(target(:, 1:6, :)), [], 1, 'omitnan'))/sqrt(NumExp);
    freqPwf1_mean(:, :, k) = targetmean;
    freqPwf1_std(:, :, k) = targetstd;
    targetmean = squeeze(mean(log10(target(:, 1:6, 2))-log10(target(:, 1:6, 1)), 1, 'omitnan'));
    targetstd = squeeze(std(abs(log10(target(:, 1:6, 2))-log10(target(:, 1:6, 1))), [], 1, 'omitnan'))/sqrt(NumExp);
    freqPwf1_diff_mean(:, k) = targetmean;
    freqPwf1_diff_std(:, k) = targetstd;
    for i = 1:5
        yc = [yc; log10(target(:, i, 2))-log10(target(:, i, 1))];
        gtypc = [gtypc; k*ones(size(target, 1), 1)];
        gfrqc = [gfrqc; i*ones(size(target, 1), 1)];
        for j = 1:2
            %             ya = [ya; reshape(f1pha(:, i, j, k), [], 1)];
            y = [y; log10(target(:, i, 1))];
            gtyp = [gtyp; k*ones(size(target, 1), 1)];
            gfrq = [gfrq; i*ones(size(target, 1), 1)];
            gsiz = [gsiz; j*ones(size(target, 1), 1)];
        end
    end
    %%
    target = f1pha(:, 1:6, :, k);
    v = normAngle(target(:, 1:6, 2)-target(:, 1:6, 1));
    
    targetmean = squeeze(mean(v, 1, 'omitnan'));
    targetstd = squeeze(std(v, [], 1, 'omitnan'))/sqrt(NumExp);
    freqPhf1_diff_mean(:, k) = targetmean;
    freqPhf1_diff_std(:, k) = targetstd;
end

%% Test for XBC against the rest
testtypeid = 4; 
y1 = y(gsiz == 1 & ismember(gfrq, [3 4 5]));
g1 = gtyp(gsiz == 1 & ismember(gfrq, [3 4 5]));
g2 = gfrq(gsiz == 1 & ismember(gfrq, [3 4 5]));
keepids = all(~isnan([y1, g1, g2]), 2);
y1 = y1(keepids);
g1 = g1(keepids);
g2 = g2(keepids);
yx = mean([y1(g1 == testtypeid & g2==3) y1(g1 == testtypeid & g2==4) y1(g1 == testtypeid & g2==5)], 2);
yy = mean([y1(g1~= testtypeid & g2==3) y1(g1 ~= testtypeid & g2==4) y1(g1 ~= testtypeid & g2==5)], 2);
[p, ~, stats] = ranksum(yx, yy, 'tai', 'right')
%% Test for BC5t against the rest
testtypeid = 3; 
y1 = y(gsiz == 1 & ismember(gfrq, [3 4 5]));
g1 = gtyp(gsiz == 1 & ismember(gfrq, [3 4 5]));
g2 = gfrq(gsiz == 1 & ismember(gfrq, [3 4 5]));
keepids = all(~isnan([y1, g1, g2]), 2);
y1 = y1(keepids);
g1 = g1(keepids);
g2 = g2(keepids);
yx = mean([y1(g1 == testtypeid & g2==3) y1(g1 == testtypeid & g2==4) y1(g1 == testtypeid & g2==5)], 2);
yy = mean([y1(g1~= testtypeid & g2==3) y1(g1 ~= testtypeid & g2==4) y1(g1 ~= testtypeid & g2==5)], 2);
[p, ~, stats] = ranksum(yx, yy, 'tail', 'left')
%% Test for BC5t Surround enhancement against the rest
testtypeid = 7;
y1 = yc(ismember(gfrqc, [4 5]));
g1 = gtypc(ismember(gfrqc, [4 5]));
g2 = gfrqc(ismember(gfrqc, [4 5]));
keepids = all(~isnan([y1, g1, g2]), 2);
y1 = y1(keepids);
g1 = g1(keepids);
g2 = g2(keepids);
yx = mean([y1(g1 == testtypeid & g2==4) y1(g1 == testtypeid & g2==5)], 2);
yy = mean([y1(g1 ~= testtypeid & g2==4) y1(g1 ~= testtypeid & g2==5)], 2);
[p, ~, stats] = ranksum(yx, yy, 'tail', 'right')

%% Test for Surround enhancement against itself in lower frequency, for BC5t, BC7, and BC8/9
y1 = yc(ismember(gfrqc, [1 2 3 4 5]));
g1 = gtypc(ismember(gfrqc, [1 2 3 4 5]));
g2 = gfrqc(ismember(gfrqc, [1 2 3 4 5]));
keepids = all(~isnan([y1, g1, g2]), 2);
y1 = y1(keepids);
g1 = g1(keepids);
g2 = g2(keepids);
testtypeid = 3;
yx = mean([y1(g1 == testtypeid & g2==1) y1(g1 == testtypeid & g2==2) y1(g1 == testtypeid & g2==3)], 2);
yy = mean([y1(g1 == testtypeid & g2==4) y1(g1 == testtypeid & g2==5)], 2);
[p, ~, stats] = signrank(yx, yy, 'tail', 'left')
%
testtypeid = 6;
yx = mean([y1(g1 == testtypeid & g2==1) y1(g1 == testtypeid & g2==2) y1(g1 == testtypeid & g2==3)], 2);
yy = mean([y1(g1 == testtypeid & g2==4) y1(g1 == testtypeid & g2==5)], 2);
[p, ~, stats] = signrank(yx, yy, 'tail', 'left')
%
testtypeid = 7;
yx = mean([y1(g1 == testtypeid & g2==1) y1(g1 == testtypeid & g2==2) y1(g1 == testtypeid & g2==3)], 2);
yy = mean([y1(g1 == testtypeid & g2==4) y1(g1 == testtypeid & g2==5)], 2);
[p, ~, stats] = signrank(yx, yy, 'tail', 'left')
%
yx = mean([y1(g2==1) y1(g2==2) y1(g2==3)], 2);
yy = mean([y1(g2==4) y1(g2==5)], 2);
[p, ~, stats] = signrank(yx, yy, 'tail', 'left')

%%
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
% close all
Colors = dazure(7);
sizetype =1;
figure; 
for i = 1:nType
    ctrace = squeeze(freqTrace_mean(i, :, :, sizetype));
    ctrace = ctrace-min(ctrace(:));
    scaling = 0.8/max(ctrace(:));
    ctrace = ctrace*scaling;
    for j = 1:NumDmt-1
        subplot(1, NumDmt-1, j);hold on
        shadePlot(T,  ctrace(:, j)-i+1,  squeeze(freqTrace_sem(i, :, j, sizetype))*scaling,...
                Colors(i, :), 0.3, 2); hold on;
        ylim([-6 1.5]);
        box off
        axis off
        if i == 1
            plot(T_detail, 0.5*0.5*(SinW(j, :)+1)+1, 'k');
            if j == 1
               plot([0.2 0.2], 0.8*[0 0.2]+0.2, 'k'); 
               plot([0.2 0.7], [0 0]+0.2, 'k'); 
            end
        end
        if i == nType
            plot([0 0], [-6 1.5], '--k');
        end
        
    end
    
end
subplot(1, NumDmt-1, 1); hold on
switch sizetype
    case 1
        tit = 'Center';
    case 2
        tit = 'Surround';
end
sgtitle(tit);
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_FrequencyResponseTrace_%s', SaveFolder, tit);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);


%%
cosmic = @(n) [linspace(51, 255, n); linspace(51, 0, n); linspace(153, 204, n)]'/255;
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
close all
Colors = dazure(7);
figure;
subplot(1, 4, 1); hold on
textpos = linspace(log(0.5), log(8), 7);
for i = 1:nType
    errorbar([0.5 1 2 4 8], freqPwf1_mean(1:5, 1, i), freqPwf1_std(1:5, 1, i), 'CapSize', 0,...
        'Color', Colors(i, :));
    text(exp(textpos(i)), 1, BCTypelabels{i}, 'Color',  Colors(i, :))
end
set(gca, 'XScale','log')
box off;
xlabel('Frequency (Hz)');
ylim([-1.5, 1.3]);
xlim([0.4, 10]);
yticks(-1:1:1);
yticklabels({'-1', '0', '1'});
ylabel('F1 power (log10 norm.)');
namelabel = 'f1';

subplot(1, 4, 2); hold on
for i = 1:nType
%     errorbar([0.5 1 2 4 8 16], freqPwf1_mean(:, 2, i), freqPwf1_std(:, 2, i), 'CapSize', 0,...
%         'Color', Colors(i, :));
    errorbar([0.5 1 2 4 8], freqPwf1_diff_mean(1:5, i), freqPwf1_diff_std(1:5, i), 'CapSize', 0,...
        'Color', Colors(i, :));
end
set(gca, 'XScale','log')
box off;
xlabel('Frequency (Hz)');
ylim([-1.5, 1.3]);
xlim([0.4, 10]);
yticks(-1:1:1);
yticklabels({'-1', '0', '1'});
ylabel({'Surround-Center F1 power';' (log10 norm.)'});
namelabel = 'f1';


subplot(1, 4, 3); hold on
for i = 1:nType
%     errorbar([0.5 1 2 4 8 16], freqPwf1_mean(:, 2, i), freqPwf1_std(:, 2, i), 'CapSize', 0,...
%         'Color', Colors(i, :));
    errorbar([0.5 1 2 4 8], freqPhf1_diff_mean(1:5, i), freqPhf1_diff_std(1:5, i), 'CapSize', 0,...
        'Color', Colors(i, :));
end
set(gca, 'XScale','log')
box off;
xlabel('Frequency (Hz)');
% ylim([-1.5, 1.3]);
% xlim([0.4, 10]);
% yticks(-1:1:1);
% yticklabels({'-1', '0', '1'});
ylabel({'Phase difference (S-C)';' (Angle)'});
namelabel = 'f1';

subplot(1, 4, 4); hold on
for i = 1:nType
    errorbar([0.5 1 2 4 8], freqPwf1_mean(1:5, 2, i), freqPwf1_std(1:5, 2, i), 'CapSize', 0,...
        'Color', Colors(i, :));
    text(exp(textpos(i)), 1, BCTypelabels{i}, 'Color',  Colors(i, :))
end
set(gca, 'XScale','log')
box off;
xlabel('Frequency (Hz)');
ylim([-1.5, 1.3]);
% xlim([0.4, 10]);
yticks(-1:1:1);
yticklabels({'-1', '0', '1'});
ylabel({'F1 power (log10 norm.)'});
namelabel = 'f1';
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_FrequencyResponseSummary', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);