ExpTab = unique(ROITable(ismember(ROITable(:, 8), 150), 1:3) , 'rows');
nExpTab = size(ExpTab, 1);
%% Get ROI size
if exist('ROIdiameterMeasurement.mat', 'file')
    load('ROIdiameterMeasurement.mat');
else
    nonexistfile = [];
    iex = 1;
    ROISize = nan(size(ROITable, 1), 1);
    ROISizepix = nan(size(ROITable, 1), 1);
    ROIdia = nan(size(ROITable, 1), 1);
    ScanningSize = [...
        32,40, 1, 1;
        64,40, 1, 0.5;
        64,80, 1, 1;
        128,40, 1, 0.25;
        128,80, 1, 0.5;
        256,40, 1, 0.125;
        256,80, 1, 0.25;
        512,512, 1, 1];
    for i = 1:nExpTab
        Day = num2str(ExpTab(i, 1));
        if numel(Day) == 5
            Day = ['0' Day];
        end
        Folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BlockFreq\Analysis\ROIs';
        FileName = sprintf('Ai148_AAV-Grm6Cre_ROIs_MorphSeg_%s_%d%03d.mat', Day, ExpTab(i, 2), ExpTab(i, 3));
        if ~exist([Folder '\' FileName], 'file')
            nonexistfile{iex} = FileName;
            iex = iex + 1;
            continue
        else
            I = load([Folder '\' FileName]);
            RmIds = any(isnan(I.wROISig), 2);
            I.SlcROI(:, :, RmIds) = [];
        end
        Folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BlockFreq\Analysis\Translation';
        FileName = sprintf('Ai148_AAV-Grm6Cre_Translation_%s_%d%03d.mat', Day, ExpTab(i, 2), ExpTab(i, 3));
        Q = load([Folder '\' FileName]);
        %% Process spatial
        Folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BlockFreq\Analysis\ROIs\CrossSection';
        FileName = sprintf('Ai148_AAV-Grm6Cre_CrossSection_%s_%d%03d.mat', Day, ExpTab(i, 2), ExpTab(i, 3));
        if ~exist([Folder '\' FileName], 'file')
            close all
            Sections = nan(2, 2, size(I.SlcROI, 3));
            figure;
            subplot(1, 2, 1);
            imagesc(std(Q.I, [], 3));
            subplot(1, 2, 2);
            for k = 1:size(I.SlcROI, 3)
                imagesc(I.SlcROI(:, :, k));
                [x,y] = ginput(2);
                Sections(:, :, k) = [x y];
            end
            %%
            figure;
            subplot(1, 2, 1);
            imagesc(std(Q.I, [], 3));
            subplot(1, 2, 2);
            imagesc(sum(I.SlcROI, 3)); hold on
            for k = 1:size(I.SlcROI, 3)
                plot(Sections(:, 1, k), Sections(:, 2, k), '-r');
                text(mean(Sections(:, 1, k)), mean(Sections(:, 2, k)), sprintf('%d', k));
            end
            save([Folder '\' FileName], 'Sections');
            keyboard;
        else
            load([Folder '\' FileName]);
        end
        %%
        
        aids= ROITable(:, 1) == ExpTab(i, 1) & ROITable(:, 2) == ExpTab(i, 2) &...
            ROITable(:, 3) == ExpTab(i, 3);
        a = ROITable(aids, :);
        assert(size(a, 1) == size(I.SlcROI, 3));
        
        %
        CurDocNum = DocNum(DocNum(:,ColName2Ind('Day')) == ExpTab(i, 1) &...
            DocNum(:, ColName2Ind('Cel')) == ExpTab(i, 2) &...
            DocNum(:,ColName2Ind('Exp')) == ExpTab(i, 3), :);
        cX = CurDocNum(ColName2Ind('SizeX'));
        cY = CurDocNum(ColName2Ind('SizeY'));
        cRio = ScanningSize(ScanningSize(:, 1) == cX & ScanningSize(:, 2) == cY, 3:4);
        wd = cRio(1)*300/(cX*CurDocNum(ColName2Ind('Zoom')));
        hd = cRio(2)*300/(cY*CurDocNum(ColName2Ind('Zoom')));
        ROISize(aids) = squeeze(sum(I.SlcROI, 1:2))*wd*hd;
        sSection = Sections.*[wd hd];
        ROISizepix(aids) = squeeze(sum(I.SlcROI, 1:2));
        ROIdia(aids) = squeeze(sum((sSection(1, :, :)-sSection(2, :, :)).^2, 2));
        clc
        fprintf('Progress...%d/%d \n', i, nExpTab);
    end
    save('ROIdiameterMeasurement.mat', 'ROISize', 'ROISizepix', 'ROIdia');
end
%%
celltype = [50 51 57 58 6 7 89];
celltypelabel = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC8/9'};
nType = numel(celltype);
nL = 10;
fsets = {2:3, 4, 6:7, 12:13, 23:24, 44:46};
%%
c = ROITable(~isnan(ROITable(:, 1)) & ROITable(:, 4) ~=5, :);
c = unique(c(:, 1:2), 'rows');
%%
keyboard;
%%
mRad = nan(7, 1);
mSize = nan(7, 1);
clear Vdata
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(nType);
keepids = zeros(size(ROIdia, 1), 1);
figure; 
for i = 1:7
    subplot(2, 4, i); hold on
    cids = find(ROITable(:, 4) == celltype(i) & ROISizepix >= 5 &...
        max(ROIQuality.^2, [], 2) >0.1);
    a = sqrt(ROIdia(cids));
    keepids(cids(~isoutlier(a))) = 1;
    Vdata{i} = a(~isoutlier(a));
    mRad(i) = median(Vdata{i} );
    mSize(i) = length(Vdata{i} );
    h = histogram(Vdata{i} , linspace(0, 5, 25));
    h.Normalization = 'Probability';
    h.EdgeColor = 'k';
    h.FaceColor = Colors(i, :);
    plot(mRad(i)*ones(1, 2), [0 0.24], '--w');
    title(sprintf('%s \n', celltypelabel{i}));
    box off
    xlabel('ROI slice diameter (um)');
    ylabel('Probability');
    xlim([0 5]);
    xticks(0:2:4);
    xticklabels({'0', '2', '4'});
    yticks(0:0.1:0.2);
    yticklabels({'0', '0.1', '0.2'});
    ylim([0 0.24]);
    ax = gca;
    ax.YColor = 'w';
    ax.XColor = 'w';
    
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure2_ROISizexBCType', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
close all
figure; 
violin(Vdata, 'facecolor', Colors, 'mc', 'w', 'medc', '--w');
box off
xticks(1:7);
xticklabels(celltypelabel);
yticks(0:2:4);
yticklabels({'0', '2', '4'});
ylabel('ROI slice diameter (um)');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure2_ROISizexBCType_violinplot', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource');
pVals = nan(7, 7);
for i = 1:7
    for j = 1:7
        if j >= i
            a = sqrt(ROIdia(ROITable(:, 4) == celltype(i) & ROISizepix >= 5 &...
                max(ROIQuality.^2, [], 2) >0.1));
            b = sqrt(ROIdia(ROITable(:, 4) == celltype(j) & ROISizepix >= 5 &...
                max(ROIQuality.^2, [], 2) >0.1));
            [~, pVals(i, j)] = kstest2(a, b, 'Tail', 'smaller');
        else
            a = sqrt(ROIdia(ROITable(:, 4) == celltype(i) & ROISizepix >= 5 &...
                max(ROIQuality.^2, [], 2) >0.1));
            b = sqrt(ROIdia(ROITable(:, 4) == celltype(j) & ROISizepix >= 5 &...
                max(ROIQuality.^2, [], 2) >0.1));
            [~, pVals(i, j)] = kstest2(a, b, 'Tail', 'smaller');
        end
    end
end
adj_pVals = pVals;
adj_pVals(eye(7)==1) = nan;
[~, ~, ~, adj_pVals]=fdr_bh(adj_pVals,0.05,'pdep','yes');
figure; 
imagesc(log10(adj_pVals)<log10(0.05)); colorbar;
% imagesc(log10(adj_pVals)); colorbar;
%% 
y1 = [];
g1 = [];
for i = 1:7
    ids = ROITable(:, 4) == celltype(i) & ROISizepix >= 5 &...
                max(ROIQuality.^2, [], 2) >0.1;
    y1 = [y1; sqrt(ROIdia(ids))];
    g1 = [g1; i*ones(sum(ids), 1)];
end

[p, tbl, stats] = kruskalwallis(y1, g1);
c = multcompare(stats);
%% Test if larger ROI size induce difference in frequency tuning
% split the data of each experiments by ROIs size (large and small)
Twin = size(DMat, 3);
ExpMeanNrmL = nan(7, Twin, nExpTab);
ExpMeanNrmS = nan(7, Twin, nExpTab);
ExpSizeTab = nan(nExpTab, 4);
ROISizeS = [];
ROITypeS = [];
ROISizeL = [];
ROITypeL = [];
ExpQual = nan(nExpTab, 7, 2);
for i = 1:nExpTab
    aids= find(ROITable(:, 1) == ExpTab(i, 1) & ROITable(:, 2) == ExpTab(i, 2) &...
        ROITable(:, 3) == ExpTab(i, 3) & ROISizepix >= 5 & max(ROIQuality.^2, [], 2) >0.1);
    if isempty(aids)
        continue
    end
    ExpSizeTab(i, :) = [ExpTab(i, :) ROITable(aids(1), 4)];
    
%     a = ROITable(aids, :);
    Sig = squeeze(DMat(:, :, :, aids));
%     cqual = [];
%     for k = 1:length(aids)
%         cqual = [cqual; estimateRepeatReliability(Sig(:, :, k))];
%     end
%     aids = aids(cqual > 0.1);
    nids = length(aids);
    nidsS = round(nids/2);
    c = sqrt(ROIdia(aids));
    [~, sids] = sort(c);
    sidsS = sids(1:nidsS);
    sidsL = sids((nidsS+1):end);
    ROISizeS = [ROISizeS(:); c(sidsS)];
    ROITypeS = [ROITypeS(:); ExpSizeTab(i, 4)*ones(length(sidsS), 1)];
    ROISizeL = [ROISizeL(:); c(sidsL)];
    ROITypeL = [ROITypeL(:); ExpSizeTab(i, 4)*ones(length(sidsL), 1)];
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
keyboard;
%%
close all
ids = find(ExpSizeTab(:, 4) ==7);
ids = ids(randperm(length(ids)));
figure;
for j = 1:5
    subplot(1, 5, j); hold on
    for i = 1:min([6, length(ids)])
        plot(T, ExpMeanNrmS(j, :, i)-(i-1)*2+1, 'k');
        plot(T, ExpMeanNrmL(j, :, i)-(i-1)*2, 'r');
    end
    ylim([-11 2]);
end
%%
figure; 
for i = 1:7
    subplot(2, 4, i); hold on
    h1 = histogram(ROISizeS(ROITypeS == celltype(i)), linspace(0, 5, 25));
    h1.EdgeColor = 'w';
    h1.FaceColor = 'b';
    h1.Normalization = 'Probability';
    h2 = histogram(ROISizeL(ROITypeL == celltype(i)), linspace(0, 5, 25));
    h2.EdgeColor = 'w';
    h2.FaceColor = 'r';
    h2.Normalization = 'Probability';
    xlabel('ROI slice diameter (um)');
    ylabel('Probability');
end
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
ylim([0 0.25]);
yticks(0:0.1:0.2);
yticklabels({'0', '0.1', '0.2'});
xlim([0 5]);
xticks(0:2:4);
xticklabels({'0', '2', '4'});
ax = gca;
ax.YColor = 'w';
ax.XColor = 'w';
xlabel('ROI slice diameter (um)');
ylabel('Probability');
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure2_ROISizeRecordingSplit_histogram', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%% Parameter should have saved
StmDur = 2.5; % in second
BasDur = 0.5;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm;
T = round(linspace(-BasDurFrm/Fz, StmDurFrm/Fz, length(1:1:SpnDurFrm))*100)/100;

NumDmt = length(Dmt);
%% Get FFT Power
[uniCel, cids] = unique(ExpSizeTab(:, [1:2]), 'rows');
uniCel = [uniCel ExpSizeTab(cids, 4)];
nuniCel = size(uniCel, 1);
CelMeanNrmS = nan(size(ExpMeanNrmS, 1), size(ExpMeanNrmS, 2), nuniCel);
CelMeanNrmL = nan(size(ExpMeanNrmL, 1), size(ExpMeanNrmL, 2), nuniCel);
for i = 1:nuniCel
    cids = ExpSizeTab(:, 1) == uniCel(i, 1) & ExpSizeTab(:, 2) == uniCel(i, 2);
    CelMeanNrmS(:, :, i) = mean(ExpMeanNrmS(:, :, cids), 3, 'omitnan');
    CelMeanNrmL(:, :, i) = mean(ExpMeanNrmL(:, :, cids), 3, 'omitnan');
end
CelMeanNrmS = CelMeanNrmS-min(CelMeanNrmS, [], 1:2);
CelMeanNrmS = CelMeanNrmS./range(CelMeanNrmS, 1:2);
CelMeanNrmS = CelMeanNrmS-mean(CelMeanNrmS(:, 1:BasDurFrm, :), 2);
CelMeanNrmL = CelMeanNrmL-min(CelMeanNrmL, [], 1:2);
CelMeanNrmL = CelMeanNrmL./range(CelMeanNrmL, 1:2);
CelMeanNrmL = CelMeanNrmL-mean(CelMeanNrmL(:, 1:BasDurFrm, :), 2);

%%
UpFac = 5;
TI = linspace(T(1), T(end), length(T)*UpFac);
n = length(TI);  
fs = Fz*UpFac;
f = (0:n-1)*(fs/n); 
fl = f(f<Fz/2);
nf = length(fl);
CellFreqPwS = nan(NumDmt, nf, nuniCel);
CellFreqPwL = nan(NumDmt, nf, nuniCel);
for i = 1:nuniCel
    for j = 1:NumDmt
        if all(~isnan(CelMeanNrmS(j, :, i)))
            YI = interp1(T, CelMeanNrmS(j, :, i), TI, 'pchip');
            YI = abs(fft(YI)).^2/n;
            CellFreqPwS(j, :, i) = YI(1:nf);
        end
        if all(~isnan(CelMeanNrmL(j, :, i)))
            YI = interp1(T, CelMeanNrmL(j, :, i), TI, 'pchip');
            YI = abs(fft(YI)).^2/n;
            CellFreqPwL(j, :, i) = YI(1:nf);
        end
    end
end
Ti = T;
%% Combine
uniCel = [uniCel(:, 1:2) ones(size(uniCel, 1), 1) uniCel(:, 3);
    uniCel(:, 1:2) 2*ones(size(uniCel, 1), 1) uniCel(:, 3)];
CelMeanNrm = cat(3, CelMeanNrmS, CelMeanNrmL);
CellFreqPw = cat(3, CellFreqPwS, CellFreqPwL);
ExpQual = cat(1, ExpQual(:, :, 1), ExpQual(:, :, 2));