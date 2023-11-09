addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');

BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCTypes = [50 51 57 58 6 7 89];
nType = length(BCTypes);
SpotSize = [1 2];% small (1) and large (2)
SpotSizeLabel = {'Small', 'Large'};
NumSpotSize = length(SpotSize);
T_detail = linspace(min(T(:)), max(T(:)), (length(T)-1)*5);
SinW = sin(2*pi*[Freqs 32]'*T_detail(:)');
SinW(:, T_detail(:)<0) = 0;
SinW(:, T_detail(:)>2) = 0;
AngS = [];
AngL = [];
typenums = nan(nType, NumDmt-1, NumSpotSize);
Celtab = [];
for k = 1:nType
    for i = 1:NumDmt-1
        for j = 1:NumSpotSize
            ExpIds = uniCel(:, 4) == BCTypes(k) & uniCel(:, 3) == SpotSize(j);
            uniExp = unique(uniCel(ExpIds, 1:2), 'rows');
            NumExp = size(uniExp, 1);
            typenums(k, i, j) = NumExp;
            if i == 1 && j == 1 && k == 1
                f0mean = nan(20, NumDmt, NumSpotSize, nType);
                f1pow = nan(20, NumDmt, NumSpotSize, nType);
                Celtab = nan(20, NumDmt, NumSpotSize, nType);
                freqTrace_mean = nan(nType, size(CelMeanNrm, 2), NumDmt, NumSpotSize);
                freqTrace_sem = nan(nType, size(CelMeanNrm, 2), NumDmt, NumSpotSize);
                freqrepqual = nan(nType, NumDmt, NumSpotSize);
            end
            cResp = [];
            cPow = [];
            cPha = [];
            cQual = [];
            for p = 1:NumExp
                ExpIds = uniCel(:, 4) == BCTypes(k) &  uniCel(:, 3) == SpotSize(j) &...
                    uniCel(:, 1) == uniExp(p, 1) & uniCel(:, 2) == uniExp(p, 2);
                rNumExp = sum(~isnan(CelMeanNrm(i, 1, ExpIds)));
                cQual = [cQual; ExpQual(ExpIds, i)];
                cResp = [cResp; squeeze(mean(CelMeanNrm(i, :, ExpIds), 3, 'omitnan'))];
                cPow = [cPow; squeeze(mean(mean(CellFreqPw(i, fsets{i}, ExpIds), 2, 'omitnan'), 3, 'omitnan'))];
            end
            freqTrace_mean(k, :, i, j) = mean(cResp, 1, 'omitnan');
            freqTrace_sem(k, :, i, j) = std(cResp, [], 1, 'omitnan')/sqrt(rNumExp);
            f0mean(1:NumExp, i, j, k) = mean(cResp(:, T>0), 2, 'omitnan');
            f1pow(1:NumExp, i, j, k) = mean(cPow, 2, 'omitnan');
            freqrepqual(k, i, j) = mean(cQual, 'omitnan');
            Celtab(1:NumExp, i, j, k) = k;
        end
    end
end


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
            y = [y; log10(target(:, i, j))];
            gtyp = [gtyp; k*ones(size(target, 1), 1)];
            gfrq = [gfrq; i*ones(size(target, 1), 1)];
            gsiz = [gsiz; j*ones(size(target, 1), 1)];
        end
    end
    %% for phase
%     target = f1pha(:, 1:6, :, k);
%     v = normAngle(target(:, 1:6, 2)-target(:, 1:6, 1));
%     
%     targetmean = squeeze(mean(v, 1, 'omitnan'));
%     targetstd = squeeze(std(v, [], 1, 'omitnan'))/sqrt(NumExp);
%     freqPhf1_diff_mean(:, k) = targetmean;
%     freqPhf1_diff_std(:, k) = targetstd;
end
%%
c = tabulate(gtyp(gfrq==3 & gsiz ==1 & ~isnan(y)));
%%
keepids = all(~isnan([y(:) gtyp(:) gsiz(:) gfrq(:)]), 2);
gtyp = gtyp(keepids);
gsiz = gsiz(keepids);
gfrq = gfrq(keepids);
y = y(keepids);

%%
%% Statistical test
cid = 7;
keepids = all(~isnan([y, gtyp, gfrq, gsiz]), 2);
y1 = y(keepids & gtyp == cid);
g1 = gfrq(keepids & gtyp == cid);
g2 = gsiz(keepids & gtyp == cid);
[~,~,stats] = anovan(y1,{g1,g2},"Model","interaction", ...
    "Varnames",["freq","size"]);
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);


%% Trace comparison
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
close all
Colors = dazure(7);
figure; 
for i = 1:nType
    for j = 1:NumDmt-2
        subplot(1, NumDmt-1, j);hold on
        ctrace = squeeze(freqTrace_mean(i, :, :, 1));
        ctrace = ctrace-min(ctrace(:));
        scaling = 0.8/max(ctrace(:));
        ctrace = ctrace*scaling;
        shadePlot(T,  ctrace(:, j)-(i-1)*2,  squeeze(freqTrace_sem(i, :, j, 1))*scaling,...
            Colors(i, :), 0.3, 2);
        ylim([-13 1.5]);
        if j == 1
            text(-0.2, -(i-1)*2+0.5, 'S');
        end
        ctrace = squeeze(freqTrace_mean(i, :, :, 2));
        ctrace = ctrace-min(ctrace(:));
        scaling = 0.8/max(ctrace(:));
        ctrace = ctrace*scaling;
        shadePlot(T,  ctrace(:, j)-i*2+1,  squeeze(freqTrace_sem(i, :, j, 2))*scaling,...
            Colors(i, :)*0.7, 0.3, 2);
        if j == 1
            text(-0.2, -i*2+1.5, 'L');
        end
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
            plot([0 0], [-13 1.5], '--k');
        end
        
    end
    
end
% subplot(1, NumDmt-1, 1); hold on
sgtitle('Comparison ROI diameter size');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sSupFig2X_RecordingSplitFrequencyResponseTrace_%s', SaveFolder);
FleNam = sprintf('%sSupFig2X_RecordingSplitFrequencyResponseTrace_%s', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);


%%
cosmic = @(n) [linspace(51, 255, n); linspace(51, 0, n); linspace(153, 204, n)]'/255;
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
% close all
Colors = dazure(7);
figure;
for i = 1:nType
    subplot(2, 4, i); hold on
    errorbar([0.5 1 2 4 8], freqPwf1_mean(1:5, 1, i), freqPwf1_std(1:5, 1, i), 'CapSize', 0,...
        'Color', Colors(i, :));
    errorbar([0.5 1 2 4 8], freqPwf1_mean(1:5, 2, i), freqPwf1_std(1:5, 2, i), 'CapSize', 0,...
        'Color', Colors(i, :)*0.7);
    text(1, 1, BCTypelabels{i}, 'Color',  Colors(i, :))
    set(gca, 'XScale','log')
    box off;
    xlabel('Frequency (Hz)');
    ylim([-1.5, 1.3]);
    xlim([0.4, 10]);
    yticks(-1:1:1);
    yticklabels({'-1', '0', '1'});
    ylabel('F1 power (log10 norm.)');
    namelabel = 'f1';
    legend({'small', 'large'});
end

%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sSupFig2_RecordingSplitFrequencyResponseSummary', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
keyboard;
%% Test across all cell type
%%
IsCorrGC6f = 1;
Sim = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Simulation\GCaMP6f\simf1powere.mat', 'samplefreqs', 'f1pow');
% Sim.f1pow = Sim.f1pow./sum(Sim.f1pow);
% freqPwf0 = nan(NumDmt-1, :, k);
afreqPwf1_mean = nan(NumDmt-1, NumSpotSize);
afreqPwf1_std = nan(NumDmt-1, NumSpotSize);
y = [];
gfrq = [];
gsiz = [];
gtyp = [];
for j = 1:NumSpotSize
    target = squeeze(f1pow(:, 1:6, j, :));
    target = permute(target, [1 3 2]);
    target = reshape(target, [], size(target, 3));
    rmids = all(isnan(target(:, 1:3)), 2);
    ctab = squeeze(Celtab(:, 1:6, j, :));
    ctab = permute(ctab, [1 3 2]);
    ctab = reshape(ctab, [], size(ctab, 3));
    ctab(rmids, :) = [];
    target(rmids, :) = [];
    target = target./max(target, [], 2);
    %     target = target./sum(target,  2);
    if IsCorrGC6f
        target = target(:, 1:6)./Sim.f1pow';
    end
    NumExp = sum(~isnan(target(:, 1)));
    keyboard;
    targetmean = squeeze(mean(log10(target(:, 1:6)), 1, 'omitnan'));
    targetstd = squeeze(std(log10(target(:, 1:6)), [], 1, 'omitnan'))/sqrt(NumExp);
    afreqPwf1_mean(:, j) = targetmean;
    afreqPwf1_std(:, j) = targetstd;
    for i = 1:5
        y = [y; log10(target(:, i))];
        gfrq = [gfrq; i*ones(size(target, 1), 1)];
        gsiz = [gsiz; j*ones(size(target, 1), 1)];
        gtyp = [gtyp; ctab];
    end
    fprintf('%d n=%d \n', j, NumExp);
end
%%
figure; hold on
for i = 1:nType
    errorbar([0.5 1 2 4 8], afreqPwf1_mean(1:5, 1), afreqPwf1_std(1:5, 1), 'CapSize', 0,...
        'Color', [196 0 255]/255);
    errorbar([0.5 1 2 4 8], afreqPwf1_mean(1:5, 2), afreqPwf1_std(1:5, 2), 'CapSize', 0,...
        'Color', [255 192 0]/255);
    set(gca, 'XScale','log')
    box off;
    xlabel('Frequency (Hz)');
    ylim([-1.5, 1.3]);
    xlim([0.4, 10]);
    yticks(-1:1:1);
    yticklabels({'-1', '0', '1'});
    ylabel('F1 power (log10 norm.)');
    legend({'small', 'large'});
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure2_SmallvsLargeROIFreqResponse_all', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
keepids = all(~isnan([y, gsiz, gfrq]), 2);
[~,tbl,stats] = anovan(y(keepids),{gsiz(keepids),gfrq(keepids)},"Model","interaction", ...
    "Varnames",["size","freq"]);
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

%% Test if Size difference in BC6 is significant
testtypeid = 4; 
frqwin = 3:5;
y1 = y(ismember(gfrq, frqwin));
g1 = gfrq(ismember(gfrq, frqwin));
g2 = gsiz(ismember(gfrq, frqwin));
% y1 = y(ismember(gfrq, frqwin));
% g1 = gfrq(ismember(gfrq, frqwin));
% g2 = gsiz(ismember(gfrq, frqwin));
keepids = all(~isnan([y1, g1, g2]), 2);
y1 = y1(keepids);
g1 = g1(keepids);
g2 = g2(keepids);
yx = mean([y1(g1 == 3 & g2==1) y1(g1 == 4 & g2==1) y1(g1 == 5 & g2==1)], 2, 'omitnan');
yy = mean([y1(g1 == 3 & g2==2) y1(g1 == 4 & g2==2) y1(g1 == 5 & g2==2)], 2, 'omitnan');
[p, ~, stats] = ranksum(yx, yy, 'tai', 'left')