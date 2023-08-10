clear; close all; clc;
%%
CellStr = {'53018_2', '10119_4', '61018_1', '53018_1', '12019_7', ....
           '10919_3', '11819_1', '12019_4', '122718_2', '10119_2', ....
           '92418_2', '11019_3', '10519_3', '12119_3', '11919_1'};
CelFils = [ 1 9 1 3 3,...
           16 3 4 4 4,...
            3 3 9 4 9];
nCell = length(CellStr);
%%
% from RFMapping_RF_DiffEncoding_batch
% Goal: summary plot of Different encoding properties of individual BC to
%   to exam if spatial information of noise play an role in heterogenous
%   responses in the subcellular units
% FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\SpotPathDistance\';
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\SpotPathDistance_woQualtiy\';

FilContain = 'SpotStimLengthConstant';
FileNames = dir(fullfile(FolderPath, '*.mat'));
%% Load Excel Experimental Document
% FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
% [DocNum, DocStr] = xlsread(FilNam,'Sheet1', 'A:G');
% ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));
% DocStr = DocStr(1, :);
%%
clear De Dap Dai Dcp Dci Dtab Dabsp Dabsi Dabtp Dabti
Topic = 'Ai148_AAV-Grm6Cre';
clear gDdp gDdi gDbs gDbt gDq gDe
for c = 1:nCell
    FileIds = find(contains({FileNames.name},[Topic '_' FilContain, '_', CellStr{c} '_' num2str(CelFils(c))]));
    nFile = length(FileIds);
    for f = 1:nFile
        FileName = FileNames(FileIds(f)).name;
        F = load([FolderPath '\' FileName], 'LenConsts', 'Da', 'De', 'bctype', 'Db', 'Dq',...
            'Ds', 'Dac', 'Ddp', 'Ddi', 'bctype');
        LenConsts = F.LenConsts;
        nLen = length(LenConsts);
        if f == 1
            De = F.De;
            Dap = squeeze(F.Da(:, :, 1));
            Dai = squeeze(F.Dac(:, :, 1));
            Ddp = F.Ddp;
            Ddi = F.Ddi;
            Dbs = F.Db(:, 1);
            Dbt = F.Db(:, 2);
            Dq = F.Dq;
            bctype = F.bctype;
        else
            De = [De; F.De];
            Dap = [Dap; squeeze(F.Da(:, :, 1))];
            Dai = [Dai; squeeze(F.Dac(:, :, 1))];
            Ddp = [Ddp; F.Ddp];
            Ddi = [Ddi; F.Ddi];
            Dbs = [Dbs; F.Db(:, 1)];
            Dbt = [Dbt; F.Db(:, 2)];
            Dq = [Dq; F.Dq];
        end
        clear F
    end
    gDdp{c} = Ddp;
    gDdi{c} = Ddi;
    gDbs{c} = Dbs;
    gDbt{c} = Dbt;
    gDq{c} = Dq;
    gDe{c} = De;
    gDap{c} = Dap;
    gDai{c} = Dai;
    Dtab(c, :) = [c bctype];
    
    for i = 1:nLen
        Dcp(c, i) = corr(Dap(:, i), De(:));
        Dci(c, i) = corr(Dai(:, i), De(:));
        Dabsp(c, i) = corr(1-Dap(:, i), Dbs(:));
        Dabtp(c, i) = corr(1-Dap(:, i), Dbt(:));
        Dabsi(c, i) = corr(1-Dai(:, i), Dbs(:));
        Dabti(c, i) = corr(1-Dai(:, i), Dbt(:));
    end
end
%%
SelectSpotRecordingFiles;
FilContain = 'SpotStimLengthConstantImageDist';
nCell = size(SelRec4Het, 1);
iDci = nan(nCell, nLen);
iDtab = nan(nCell, 1);
for c = 1:nCell
    Daystr = num2str(SelRec4Het(c, 1));
    FileIds = find(contains({FileNames.name},[Topic '_' FilContain, '_', Daystr '_' num2str(SelRec4Het(c, 2)),...
        '_' num2str(SelRec4Het(c, 3)) '.']));
    assert(length(FileIds)==1);
    FileName = FileNames(FileIds).name;
    F = load([FolderPath '\' FileName], 'LenConsts', 'Dac', 'De', 'Dq', 'bctype');
    iDtab(c) = F.bctype;
    LenConsts = F.LenConsts;
    nLen = length(LenConsts);
    for i = 1:nLen
        iDci(c, i) = corr(squeeze(F.Dac(:, i, 1)), F.De);
    end
    %     De = F.De;
    %     Dai = squeeze(F.Dac(:, :, 1));
    %     igDe{c} = De;
    %     igDai{c} = Dai;
end
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCelen = [19.4 18.4 20.8 36.6 14.4 19.7 43.4];
%%
%[247 179 32]/255 [166 30 230]/255
datacounts = nan(7, 2);
for i = 1:7
    switch i
        case 1
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
            ytic1 =0:0.1:0.3;
            yticlab1 = ({'0', '0.1', '0.2', '0.3'});
            ytic2 =0:0.1:0.3;
            yticlab2 = ({'0', '0.1', '0.2', '0.3'});
            ytic3 =0:0.05:0.1;
            yticlab3 = ({'0', '0.05', '0.1'});
        case 2
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
            ytic1 =0:0.1:0.3;
            yticlab1 = ({'0', '0.1', '0.2', '0.3'});
            ytic2 =0:0.1:0.3;
            yticlab2 = ({'0', '0.1', '0.2', '0.3'});
            ytic3 =0:0.1:0.2;
            yticlab3 = ({'0', '0.1', '0.2'});
        case 3
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
            ytic1 =0:0.1:0.2;
            yticlab1 = ({'0', '0.1', '0.2'});
            ytic2 =0:0.5:1;
            yticlab2 = ({'0', '0.5', '1'});
            ytic3 =0:0.1:0.2;
            yticlab3 = ({'0', '0.1', '0.2'});
        case 4
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
            ytic1 =0:0.2:0.4;
            yticlab1 = ({'0', '0.2', '0.4'});
            ytic2 =0:0.1:0.2;
            yticlab2 = ({'0', '0.1', '0.2'});
            ytic3 =0:0.1:0.2;
            yticlab3 = ({'0', '0.1', '0.2'});
            ylms2 = [0 0.2];
            ylms3 = [0 0.2];
        case 5
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
            ytic1 =0:0.05:0.15;
            yticlab1 = ({'0', '0.05', '0.1', '0.15'});
            ytic2 =0:0.2:0.4;
            yticlab2 = ({'0', '0.2', '0.4'});
            ytic3 =0:0.05:0.1;
            yticlab3 = ({'0', '0.05', '0.1'});
        case 6
            xtic = 0:15:60;
            xticlab = {'0', '', '30', '', '60'};
            ytic1 =0:0.1:0.2;
            yticlab1 = ({'0', '0.1', '0.2'});
            ytic2 =0:0.4:0.8;
            yticlab2 = ({'0', '0.4', '0.8'});
            ytic3 =0:0.01:0.03;
            yticlab3 = ({'0', '0.01', '0.02', '0.03'});
        case 7
            xtic = 0:12.5:50;
            xticlab = {'0', '', '25', '', '50'};
            ytic1 =0:0.2:0.4;
            yticlab1 = ({'0', '0.2', '0.4'});
            ytic2 =0:0.6:1.2;
            yticlab2 = ({'0', '0.6', '1.2'});
            ytic3 =0:0.1:0.2;
            yticlab3 = ({'0', '0.1', '0.2'});
    end
    close all
    figure; hold on
    bids = find(Dtab(:, 2) == BCTypes(i));
%     bids = bids(1);
    x1 = []; x2 = []; y = []; q = [];
    for j = 1:length(bids)
        x1 = [x1; gDdp{bids(j)}(:)];
        x2 = [x2; gDdi{bids(j)}(:)];
        y = [y; gDbs{bids(j)}(:)];
        q = [q; gDq{bids(j)}(:)];
    end
    kids = q>0.3 & ~isoutlier(y); %
    %%
    subplot(2, 4, [1 5]); hold on
    xa = x1(kids);
    ya = y(kids);
    scatter(xa, ya, 5, [247 179 32]/255, 'filled');
    datacounts(i, 1) = sum(kids);
    R1 = corr(xa, ya);
    X = [ones(sum(kids),1) xa];
    b1 = X\y(kids);
    yCalc2 = X*b1;
    plot(xa,yCalc2,'Color', [247 179 32]/255);
    xb = x2(kids);
    scatter(xb, ya, 5, [166 30 230]/255, 'filled');
    R2 = corr(xb, ya);
    X = [ones(sum(kids),1) xb];
    b1 = X\y(kids);
    yCalc2 = X*b1;
    plot(xb,yCalc2,'Color', [166 30 230]/255);
    title(sprintf('%0.3G (o) %0.3G (p)', R2^2, R1^2));
    xlabel('Distance (um)');
    ylabel('Diff. Surround(arb.)');
    xticks(xtic);
    xticklabels(xticlab);
    yticks(ytic1);
    yticklabels(yticlab1);
        
    
    %%
    keyboard;
    
    %%
    y = []; 
    for j = 1:length(bids)
        y = [y; gDbt{bids(j)}(:)];
    end
    kids = q>0.3 & ~isoutlier(y);
    xa = x1(kids);
    ya = y(kids);
    
    subplot(2, 4, [2 6]); hold on
    scatter(xa, ya, 5, [247 179 32]/255, 'filled');
    datacounts(i, 2) = sum(kids);
    R1 = corr(xa, ya);
    X = [ones(sum(kids),1) xa];
    b1 = X\y(kids);
    yCalc2 = X*b1;
    plot(xa,yCalc2,'Color', [247 179 32]/255);
    
    xb = x2(kids);
    scatter(xb, ya, 5, [166 30 230]/255, 'filled');
    R2 = corr(xb, ya);
    X = [ones(sum(kids),1) xb];
    b1 = X\y(kids);
    yCalc2 = X*b1;
    plot(xb,yCalc2,'Color', [166 30 230]/255);
    title(sprintf('%0.3G (o) %0.3G (p)', R2^2, R1^2));
    xlabel('Distance (um)');
    ylabel('Diff. Transience (arb.)');
    xticks(xtic);
    xticklabels(xticlab);
    yticks(ytic2);
    yticklabels(yticlab2);
%     ylim(ylms3);   
    
    y = [];
    yde = [];
    xap = [];
    xai = [];
    for j = 1:length(bids)
        y = [y; 1-gDe{bids(j)}(:);];
        yde = [yde; gDe{bids(j)}(:)];
        xap = [xap; gDap{bids(j)}];
        xai = [xai; gDai{bids(j)}];
    end
    kids = q>0.3 & ~isoutlier(y);
    
    %%
    keyboard;
    %%
    subplot(2, 4, [3 7]); hold on
    scatter(x1(kids), y(kids), 5, [247 179 32]/255, 'filled');
    R1 = corr(x1(kids), y(kids));
    X = [ones(sum(kids),1) x1(kids)];
    b1 = X\y(kids);
    yCalc2 = X*b1;
    plot(x1(kids),yCalc2,'Color', [247 179 32]/255);
    scatter(x2(kids), y(kids), 5, [166 30 230]/255, 'filled');
    R2 = corr(x2(kids), y(kids));
    X = [ones(sum(kids),1) x2(kids)];
    b1 = X\y(kids);
    yCalc2 = X*b1;
    plot(x2(kids),yCalc2,'Color', [166 30 230]/255);
    title(sprintf('%0.3G (o) %0.3G (p)', R2^2, R1^2));
    nb = length(bids);
    Colors = [166 30 230]/255;
    Colors = [Colors; zeros(nb-1, 3)];
    xlabel('Distance (um)');
    ylabel('Diff. Corr. (arb.)');
%     ylim(ylms3);
    xticks(xtic);
    xticklabels(xticlab);
    yticks(ytic3);
    yticklabels(yticlab3);
    
%     subplot(2, 4, 4); hold on
%     dcp = nan(nLen, 1);
%     dci = nan(nLen, 1);
%     kids = q>0.3 & ~isoutlier(y);%
%     for j = 1:nLen
%         dcp(j) = corr(xap(kids, j), yde(kids));
%         dci(j) = corr(xai(kids, j), yde(kids));
%     end
%     plot(LenConsts, dcp, 'Color', 'k');
%     set(gca, 'XScale', 'log');
%     subplot(2, 4, 8); hold on
%     plot(LenConsts, dci, 'Color', 'k');
%     set(gca, 'XScale', 'log');
    %%
    for j = 1:nb
        % 3
        subplot(2, 4, 4); hold on
        plot(LenConsts, Dcp(bids(j), :), 'Color', 'k');
        if j == 1
            plot(BCelen(i)*ones(1, 2), [-0.18 1], '--', 'Color', 0.4*ones(1, 3));
        end
       
        xlim([min(LenConsts) max(LenConsts)]);
        ylim([-0.18, 1]);
        set(gca, 'XScale', 'log');
        xlabel('Length constant');
        ylabel('Correlation coefficient');
        title('Path distance');
        
        box off
        
        subplot(2, 4, 8); hold on
        y = Dci(bids(j), :);
        plot(LenConsts, y, 'Color', 'k');
        if j ==1
            plot(BCelen(i)*ones(1, 2), [-0.18 1], '--', 'Color', 0.4*ones(1, 3));
        end
        
    end
    
    subplot(2, 4, 8); hold on
    bids = find(iDtab == BCTypes(i));
    nbid = length(bids);
    for j = 1:nbid
        y = iDci(bids(j), :);
        plot(LenConsts,y , 'Color', 0.3*ones(1, 3));
    end
    xlim([min(LenConsts) max(LenConsts)]);
    ylim([-0.18, 1]);
    set(gca, 'XScale', 'log');
    xlabel('Length constant');
    ylabel('Correlation coefficient');
    title('Image distance');
    box off
    
    sgtitle(sprintf('%s',  BCTypeLabels{i}));
    keyboard;
end

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure7_Heterogeneity_SpotParametersPathDistance_%s', SaveFolder, BCTypeLabels{i});
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
keyboard;
