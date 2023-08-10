clear; close all; clc;
%%
% from RFMapping_RF_DiffEncoding_batch
% Goal: summary plot of Different encoding properties of individual BC to
%   to exam if spatial information of noise play an role in heterogenous
%   responses in the subcellular units
% FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\SpotPathDistance\';
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\StimSponBCCorrelation_Image';
FileNames = dir(fullfile(FolderPath, '*.mat'));
nFile = length(FileNames);
%% Load Excel Experimental Document
% FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
% [DocNum, DocStr] = xlsread(FilNam,'Sheet1', 'A:G');
% ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));
% DocStr = DocStr(1, :);
%% get data from image without path reconstruction
Topic = 'Ai148_AAV-Grm6Cre';
Dtab = [];
for f = 1:nFile
    F = load([FolderPath '\' FileNames(f).name], 'bctype', 'Day', 'Cel','Pla','Dp', 'nROI');
    Dtab = [Dtab; F.Day, F.Cel, F.Pla, F.bctype, median(F.Dp(:), 'omitnan'), F.nROI];
end
%%
figure; histogram(Dtab(:, 5));
qualthr = 0.2;
FileNames(Dtab(:,5)<qualthr)=[];
Dtab(Dtab(:,5)<qualthr, :) = [];
uniPla = unique(Dtab(:, [1 2 3]), 'rows');
nuniPla = size(uniPla, 1);
%%
clear Dai Dbi De Di gDai gDbi gDtab gmDai
rng('shuffle')
for i = 1:nuniPla
    fids = find(Dtab(:, 1) == uniPla(i, 1) & Dtab(:, 2) == uniPla(i, 2) & Dtab(:, 3) == uniPla(i, 3));
    nfid = length(fids);
    for f = 1:nfid
        F = load([FolderPath '\' FileNames(fids(f)).name], 'bctype', 'Dai', 'Dbi','De', 'Di', 'Dais', 'Dbis', 'LenConsts',...
            'Daih', 'Dbih', 'Dp', 'Dm');
        if size(F.Dai, 1) < 2
            continue
        end
        a = permute(F.Dai, [1 3 2]);
        a = squeeze(median(a, 2));
        b = permute(F.Dbi, [1 3 2]);
        b = squeeze(median(b, 2));
        as = permute(F.Dais, [1 3 2]);
        as = squeeze(median(as, 2));
        bs = permute(F.Dbis, [1 3 2]);
        bs = squeeze(median(bs, 2));
        ah = permute(F.Daih, [1 3 2]);
        ah = squeeze(median(ah, 2));
        bh = permute(F.Dbih, [1 3 2]);
        bh = squeeze(median(bh, 2));
        if f == 1
            LenConsts = F.LenConsts;
            nLen = length(LenConsts);
            bctype = F.bctype;
            Dai = a;
            Dbi = b;
            Dais = as;
            Dbis = bs;
            Daih = ah;
            Dbih = bh;
            De = median(F.De, 2);
            Di = median(F.Di, 2);
            Dp = median(F.Dp, 2);
            Dm = median(F.Dm, 2);
        else
            Dai = [Dai; a];
            Dbi = [Dbi; b];
            Dais = [Dais; as];
            Dbis = [Dbis; bs];
            Daih = [Daih; ah];
            Dbih = [Dbih; bh];
            De = [De; median(F.De, 2)];
            Di = [Di; median(F.Di, 2)];
            Dp = [Dp; median(F.Dp, 2)];
            Dm = [Dm; median(F.Dm, 2)];
        end
        clear F
    end
    gDtab(i, :) = [i, bctype];
    if ~exist('Dai', 'var')
        continue
    end
    for q = 1:nLen
        kids = all(~isnan([Dai(:, q) Dm(:) De(:)]), 2);
        gDai(i, q) = corr(1-Dai(:, q), 1-De(:));
        gmDai(i, q) = corr(1-Dai(kids, q), 1-Dm(kids));
        gDbi(i, q) = corr(1-Dbi(:, q), 1-Di(:));
        gDais(i, q) = corr(1-Dais(:, q), 1-De(:));
        gDbis(i, q) = corr(1-Dbis(:, q), 1-De(:));
        kids = ~isnan(Daih(:, q));
        gDaih(i, q) = corr(1-Daih(kids, q), 1-De(kids));
        kids = ~isnan(Dbih(:, q));
        gDbih(i, q) = corr(1-Dbih(kids, q), 1-De(kids));
    end
    clear Dai Dbi Dais Dbis Daih Dbih De Di
end
%% get data from path reconstruction
% first run @GetRepresentiveData_Hereogeneity_NaturalMovie to get the mat
filnam = 'Heterogenity_NatMov_ReconstructedBC.mat';
load(['\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\'...
    filnam], 'Dtab', 'Dcp', 'Dci', 'Dcs', 'Dch', 'mDci');
rmids = isnan(Dtab(:, 1));
Dtab(rmids, :) = [];
Dcp(rmids, :) = [];
Dci(rmids, :) = [];
Dcs(rmids, :) = [];
Dch(rmids, :) = [];
mDci(rmids, :) = [];
%%
keyboard;
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCelen = [19.4 18.4 20.8 36.6 14.4 19.7 43.4];
%%
% a = [Dci; gDai];
a = gDai;
lng_m = mean(a, 1);
lng_se = std(a, [], 1)/sqrt(size(a, 1));
% a = [Dcs; gDais];
a = gDais;
upb_m = mean(a, 1);
upb_se = std(a, [], 1)/sqrt(size(a, 1));
% a = [Dch; gDaih];
% lwb_m = mean(a, 1);
% lwb_se = std(a, [], 1)/sqrt(size(a, 1));
% a = [mDci; gmDai];
a = gmDai;
lwb_m = mean(a, 1, 'omitnan');
lwb_se = std(a, [], 1, 'omitnan')/sqrt(size(a, 1));
R = corr(gDai', gmDai');
%%
thr = 0.5;
close all
figure; 
subplot(1, 3, 1); hold on
shadePlot(LenConsts, lng_m, lng_se, 0*ones(1, 3));
shadePlot(LenConsts, upb_m, upb_se, [52 177 235]/255);
shadePlot(LenConsts, lwb_m, lwb_se, [26 23 209]/255);
set(gca, 'XScale', 'log');
xlim([min(LenConsts) max(LenConsts)]);
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
xlabel('Length constant');
ylabel('Correlation coefficient (r)');
text(5, 0.95, 'experiment', 'Color', zeros(1, 3)/255);
text(50, 0.95, 'location shuffle', 'Color', [52 177 235]/255);
text(500, 0.95, 'optimal local', 'Color', [26 23 209]/255);

a = gDai;
b = gDais;
c = gmDai; 
rmids = any([all(gDai==0, 2) all(gDais==0, 2) all(gmDai==0, 2)], 2);
a(rmids, :) = [];
b(rmids, :) = [];
c(rmids, :) = [];

ab = corr(a', b');
ab = ab(eye(size(ab, 1))==1);
ac = corr(a', c');
ac = ac(eye(size(ac, 1))==1);

nm = size(a, 2);
a = a-min(a, [], 2);
b = b-min(b, [], 2);
c = c-min(c, [], 2);
a = cumsum(a, 2);
b = cumsum(b, 2);
c = cumsum(c, 2);
a = a./max(a, [], 2)+linspace(0, 1, nm)*10e-6;
b = b./max(b, [], 2)+linspace(0, 1, nm)*10e-6;
c = c./max(c, [], 2)+linspace(0, 1, nm)*10e-6;
nd = size(a, 1);
al = nan(nd, 1);
bl = nan(nd, 1);
cl = nan(nd, 1);
for i = 1:nd
    al(i) = interp1(a(i, :), log10(LenConsts), thr, 'pchip');
    bl(i) = interp1(b(i, :), log10(LenConsts), thr, 'pchip');
    cl(i) = interp1(c(i, :), log10(LenConsts), thr, 'pchip');
end
subplot(1, 3, 2); hold on
swarmchart([ones(nd, 1) 2*ones(nd, 1) 3*ones(nd, 1)],...
    10.^[al(:), bl(:), cl(:)], 20, [0 0 0; 52 177 235; 26 23 209]/255, 'filled');
title(sprintf('50 percent cummulative threshold'));
xticks(1:3);
xticklabels({'experiment', 'location shuffle', 'optimal local'});
set(gca, 'YScale', 'log');
ylabel('Length constant (um)');


subplot(1, 3, 3); hold on
swarmchart([ones(length(ab), 1) 2*ones(length(ac), 1)],...
    [ab(:), ac(:)], 20, [52 177 235; 26 23 209]/255, 'filled');
xticks(1:2);
xticklabels({'location shuffle', 'optimal local'});
yticks(-1:0.5:1);
yticklabels({'-1', '', '0', '', '1'});
ylabel('Correlatio coefficient (r)');


%%
keyboard;
%%
%[247 179 32]/255 [166 30 230]/255
close all
figure;
for i = 1:7
    subplot(2, 4, i); hold on
    switch i
        case 1
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
           
        case 2
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
          
        case 3
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
          
        case 4
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
     
        case 5
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
          
        case 6
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
      
        case 7
            xtic = 0:12.5:50;
            xticlab = {'0', '', '25', '', '50'};
        
    end    
    bids = find(gDtab(:, 2) == BCTypes(i));

    nbid = length(bids);
    y = [gDai(bids, :); Dci(Dtab(:, 2) == BCTypes(i), :)];
    
    ym = mean(y, 1);
    ny = size(y, 1);
    ye = std(y, [], 1)/ny;
    shadePlot(LenConsts, ym, ye, 0*ones(1, 3));
%     for j = 1:nbid
%         
%         plot(LenConsts,y , 'Color', 0*ones(1, 3));
%         
%     end
    bids = find(Dtab(:, 2) == BCTypes(i));
    
%     for j = 1:nbid
%         y = Dci(bids(j), :);
%         ym = mean(y, 1);
%         ny = size(y, 1);
%         ye = std(y, [], 1)/ny;
%         shadePlot(LenConsts, ym, ye, 0*ones(1, 3));
% %         if j == 1
% %             plot(LenConsts,y , 'Color', [166 30 230]/255);
% %         else
% %             plot(LenConsts,y , 'Color', 0*ones(1, 3));
% %         end
%     end
    plot(BCelen(i)*ones(1, 2), [0 1], '--', 'Color', 0.4*ones(1, 3));
    ylim([0, 1]);
    set(gca, 'XScale', 'log');
    xlabel('Length constant');
    ylabel('Correlation coefficient');
    title('Image distance');
    box off
    yticks(0:0.25:1);
    yticklabels({'0', '', '0.5', '', '1'});
    xlim([min(LenConsts) max(LenConsts)]);
    title(BCTypeLabels{i});
    keyboard;
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure6_Heterogeneity_NatMov_DistanceCorr_Image_Shade', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam);
saveas(gcf,[FleNam '.png']);

%%
%[247 179 32]/255 [166 30 230]/255
close all
figure;
for i = 1:7
    subplot(2, 4, i); hold on
    switch i
        case 1
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
           
        case 2
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
          
        case 3
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
          
        case 4
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
     
        case 5
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
          
        case 6
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
      
        case 7
            xtic = 0:12.5:50;
            xticlab = {'0', '', '25', '', '50'};
        
    end    
    bids = find(Dtab(:, 2) == BCTypes(i));

    nbid = length(bids);
    y = Dcp(bids, :);
    ym = mean(y, 1);
    ny = size(y, 1);
    ye = std(y, [], 1)/ny;
    shadePlot(LenConsts, ym, ye, 0*ones(1, 3));
%     for j = 1:nbid
%         
%         plot(LenConsts,y , 'Color', 0*ones(1, 3));
%         
%     end
    bids = find(Dtab(:, 2) == BCTypes(i));
    
%     for j = 1:nbid
%         y = Dci(bids(j), :);
%         ym = mean(y, 1);
%         ny = size(y, 1);
%         ye = std(y, [], 1)/ny;
%         shadePlot(LenConsts, ym, ye, 0*ones(1, 3));
% %         if j == 1
% %             plot(LenConsts,y , 'Color', [166 30 230]/255);
% %         else
% %             plot(LenConsts,y , 'Color', 0*ones(1, 3));
% %         end
%     end
    plot(BCelen(i)*ones(1, 2), [0 1], '--', 'Color', 0.4*ones(1, 3));
    ylim([0, 1]);
    set(gca, 'XScale', 'log');
    xlabel('Length constant');
    ylabel('Correlation coefficient');
    title('Image distance');
    box off
    yticks(0:0.25:1);
    yticklabels({'0', '', '0.5', '', '1'});
    xlim([min(LenConsts) max(LenConsts)]);
    title(BCTypeLabels{i});
end
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure6_Heterogeneity_NatMov_DistanceCorr_Path_Shade', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam);
saveas(gcf,[FleNam '.png']);