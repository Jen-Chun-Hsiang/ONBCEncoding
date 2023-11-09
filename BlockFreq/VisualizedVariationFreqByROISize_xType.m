a = sqrt(ROIdia(ismember(ROITable(:, 4), celltype) & ROISizepix >= 5 & max(ROIQuality.^2, [], 2) >0.1));
ma = median(a);
figure; 
h = histogram(a, linspace(0, 5, 25));
h.Normalization = 'Probability';


%% Test if larger ROI size induce difference in frequency tuning
% split the data of each experiments by ROIs size (large and small)
Twin = size(DMat, 3);
ExpMeanNrmL = nan(7, Twin, nExpTab);
ExpMeanNrmS = nan(7, Twin, nExpTab);
ExpSizeTab = nan(nExpTab, 6);
ROISizeS = [];
ROISizeL = [];
ExpQual = nan(nExpTab, 7, 2);
for i = 1:nExpTab
    aids= find(ROITable(:, 1) == ExpTab(i, 1) & ROITable(:, 2) == ExpTab(i, 2) &...
        ROITable(:, 3) == ExpTab(i, 3) & ROISizepix >= 5 & max(ROIQuality.^2, [], 2) >0.1);
    
    Sig = squeeze(DMat(:, :, :, aids));
%     cqual = [];
%     for k = 1:length(aids)
%         cqual = [cqual; estimateRepeatReliability(Sig(:, :, k))];
%     end
%     aids = aids(cqual > 0.1);
    if isempty(aids)
        continue
    end
    c = sqrt(ROIdia(aids));
    sidsS = c <= ma;
    sidsL = c > ma;
    ExpSizeTab(i, :) = [ExpTab(i, :) ROITable(aids(1), 4) sum(sidsS) sum(sidsL)];
    ROISizeS = [ROISizeS(:); c(sidsS)];
    ROISizeL = [ROISizeL(:); c(sidsL)];
    for j = 1:7
        if ~isempty(sidsS)
            SigS = mean(squeeze(Sig(:, j, :, sidsS)), 3, 'omitnan');
            ExpQual(i, j, 1) = estimateRepeatReliability(SigS);
            b = squeeze(DMat(:, j, :, aids));
            ExpMeanNrmS(j, :, i) = mean(b(:, :, sidsS), [1 3], 'omitnan');
        end
        if ~isempty(sidsL)
            SigL = mean(squeeze(Sig(:, j, :, sidsL)), 3, 'omitnan');
            ExpQual(i, j, 2) = estimateRepeatReliability(SigL);
            b = squeeze(DMat(:, j, :, aids));
            ExpMeanNrmL(j, :, i) = mean(b(:, :, sidsL), [1 3], 'omitnan');
        end
    end
end
%%
sum(ExpSizeTab(ismember(ExpSizeTab(:, 4), celltype), 5:6))
%%
figure; hold on
h1 = histogram(ROISizeS, linspace(0, 5, 25));
h1.EdgeColor = 'k';
h1.FaceColor = [196 0 255]/255;
h1.Normalization = 'Probability';
h2 = histogram(ROISizeL, linspace(0, 5, 25));
h2.EdgeColor = 'k';
h2.FaceColor = [255 192 0]/255;
h2.Normalization = 'Probability';
ylim([0 0.35]);
yticks(0:0.1:0.3);
yticklabels({'0', '0.1', '0.2', '0.3'});
xlim([0 5]);
xticks(0:2:4);
xticklabels({'0', '2', '4'});
ax = gca;
ax.YColor = 'w';
ax.XColor = 'w';
xlabel('ROI slice diameter (um)');
ylabel('Probability');
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_ROISizeFixedSplit_histogram', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%% Parameter should have saved
StmDur = 2.5; % in second
BasDur = 0.5;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm;
T = round(linspace(-BasDurFrm/Fz, StmDurFrm/Fz, length(1:1:SpnDurFrm))*100)/100;

NumDmt = length(Dmt);
%% Normalize for each recording
ExpMeanNrmS = ExpMeanNrmS-min(ExpMeanNrmS, [], 1:2);
ExpMeanNrmS = ExpMeanNrmS./range(ExpMeanNrmS, 1:2);
ExpMeanNrmS = ExpMeanNrmS-mean(ExpMeanNrmS(:, 1:BasDurFrm, :), 2);
ExpMeanNrmL = ExpMeanNrmL-min(ExpMeanNrmL, [], 1:2);
ExpMeanNrmL = ExpMeanNrmL./range(ExpMeanNrmL, 1:2);
ExpMeanNrmL = ExpMeanNrmL-mean(ExpMeanNrmL(:, 1:BasDurFrm, :), 2);
%%
UpFac = 5;
TI = linspace(T(1), T(end), length(T)*UpFac);
n = length(TI);  
fs = Fz*UpFac;
f = (0:n-1)*(fs/n); 
fl = f(f<Fz/2);
nf = length(fl);
nuniExp = size(ExpMeanNrmS, 3);
ExpFreqPwS = nan(NumDmt, nf, nuniExp);
ExpFreqPwL = nan(NumDmt, nf, nuniExp);
for i = 1:nuniExp
    for j = 1:NumDmt
        if all(~isnan(ExpMeanNrmS(j, :, i)))
            YI = interp1(T, ExpMeanNrmS(j, :, i), TI, 'pchip');
            YI = abs(fft(YI)).^2/n;
            ExpFreqPwS(j, :, i) = YI(1:nf);
        end
        if all(~isnan(ExpMeanNrmL(j, :, i)))
            YI = interp1(T, ExpMeanNrmL(j, :, i), TI, 'pchip');
            YI = abs(fft(YI)).^2/n;
            ExpFreqPwL(j, :, i) = YI(1:nf);
        end
    end
end
Ti = T;
%% Combine
uniCel = [ExpSizeTab(:, 1:2) ones(size(ExpSizeTab, 1), 1) ExpSizeTab(:, 4);
    ExpSizeTab(:, 1:2) 2*ones(size(ExpSizeTab, 1), 1) ExpSizeTab(:, 4)];
CelMeanNrm = cat(3, ExpMeanNrmS, ExpMeanNrmL);
CellFreqPw = cat(3, ExpFreqPwS, ExpFreqPwL);
ExpQual = cat(1, ExpQual(:, :, 1), ExpQual(:, :, 2));