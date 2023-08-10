close all; clear; clc;
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis');
MainPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\';
load([MainPath 'Results\SplitROILevelSummary\SplitData_sampled_02072023.mat']);
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
%%
MatchROIPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\ROIMatchingTable';
Topic = 'Ai148_AAV-Grm6Cre';
Subject = 'ROIMatchingTable';
FileNames = dir(fullfile(MatchROIPath, '*.mat'));
FileIds = find(contains({FileNames.name},[Topic '_ROIMatchingTable']));
% FileIds = find(contains({FileNames.name},[Topic '_ROIMatchingTable_10119_4_1']));
nFile = length(FileIds);
IsDisplay = 0;
%%
Fz = 100;
cliplines = cumsum(ClipLens-90)/Fz;
clipvec = nan(1, cliplines(end)*Fz);
nClip = length(cliplines);
for i = 1:nClip
    if i == 1
        clipvec(1:round(cliplines(i)*Fz)) = i;
    else
        clipvec(round(cliplines(i-1)*Fz+1):round(cliplines(i)*Fz)) = i;
    end
end
%%
skipids = [];
for f = 1:nFile
    FileName = FileNames(FileIds(f)).name;
    load([MatchROIPath '\' FileName]);
    
    Words = strsplit(FileName, ['_' Subject '_']);
    Words = strsplit(Words{2}, '_');
    Daystr = Words{1};
    Day = str2double(Words{1});
    Cel = str2double(Words{2});
    Words = strsplit(Words{3}, '.');
    Pla = str2double(Words{1});
    % look for natural movie files
    if anckor(3) == 3
        expids = anckor(2);
    else
        expids = [];
    end
    a = unique(ROIMatchingTab(:, [3 4]), 'rows');
    expids = [expids; a(a(:, 2) == 3, 1)];
    nexpid = length(expids);
    if nexpid < 2
        skipids = [skipids f];
        continue
    end
    % Generate paired matching ROI table with the fixed
    pairids = nchoosek(expids, 2);
    npairid = size(pairids, 1);
    for p = 1:npairid
        clc
        fprintf('progress...%d/%d, %d/%d \n', f, nFile, p, npairid);
        expids = pairids(p, :);
        if any(expids == anckor(2))
            cexp = expids(expids~=anckor(2));
            cTab = ROIMatchingTab(ROIMatchingTab(:, 3) == cexp(1), 5:6);
            expids = [anckor(2) cexp(1)];
        else
            cTab1 = ROIMatchingTab(ROIMatchingTab(:, 3) == expids(1), 5:6);
            cTab2 = ROIMatchingTab(ROIMatchingTab(:, 3) == expids(2), 5:6);
            cTab = crossMatchingROI(cTab1, cTab2);
        end
        cTab = unique(cTab, 'rows');
        %
        cids = find(DTab(:, 1) == Day & DTab(:, 2) == Cel & DTab(:, 5) == expids(1));
        targetD = 0.5*(DTra1(cids, :)+DTra2(cids, :));
        ampD1 = targetD(cTab(:, 1), :);
        targetD = targetD - min(targetD, [], 2);
        targetD1 = targetD./max(targetD, [], 2);
        targetD1 = targetD1(cTab(:, 1), :);
        targetD1 = targetD1 + 0.01*std(targetD1(:))*rand(size(targetD1));
        %
        cids = find(DTab(:, 1) == Day & DTab(:, 2) == Cel & DTab(:, 5) == expids(2));
        targetD = 0.5*(DTra1(cids, :)+DTra2(cids, :));
        ampD2 = targetD(cTab(:, 2), :);
        targetD = targetD - min(targetD, [], 2);
        targetD2 = targetD./max(targetD, [], 2);
        targetD2 = targetD2(cTab(:, 2), :);
        targetD2 = targetD2 + 0.01*std(targetD2(:))*rand(size(targetD2));
        bctype = BCTypes(DTab(cids(1), 4));
%         if bctype ~= 50
%             continue
%         end
        %
        nROI = size(cTab, 1);
        RepRli = nan(nROI, nClip);
        RespSTD = nan(nROI, nClip, 2);
        for i = 1:nROI
            for j = 1:nClip
                RepRli(i, j) = sqrt(corr(targetD1(i, clipvec==j)', targetD2(i, clipvec==j)').^2);
                if i == 1
                    RespSTD(:, j, 1) = std(ampD1(:, clipvec==j), [], 2, 'omitnan');
                    RespSTD(:, j, 2) = std(ampD2(:, clipvec==j), [], 2, 'omitnan');
                end
            end
        end
        PirRepRli = nan(nROI, nROI, nClip);
        PirSimlar1 = nan(nROI, nROI, nClip);
        PirSimlar2 = nan(nROI, nROI, nClip);
        PirSimlarx = nan(nROI, nROI, nClip);
        PirSimlari = nan(nROI, nROI, nClip);
        for i = 1:nROI
            for j = 1:nROI
                if i < j
                    for k = 1:nClip
                        x1 = targetD1(i, clipvec==k)';
                        x2 = targetD2(i, clipvec==k)';
                        y1 = targetD1(j, clipvec==k)';
                        y2 = targetD2(j, clipvec==k)';
                        [PirSimlarx(i, j, k), PirRepRli(i, j, k), PirSimlari(i, j, k)] = estimatepairedcorr(x1, x2, y1, y2);
                        
                    end
                elseif i == j
                end
            end
        end
        
        %%
        if IsDisplay
            figure('visible','off'); hold on
            scatter(PirRepRli(:), PirSimlarx(:), 5, 'b', 'filled');
            scatter(PirRepRli(:), PirSimlari(:), 5, 'k', 'filled');
            plot([0 1], [0 1], '--', 'Color', 0.3*ones(1, 3));
        end
        incIds = ~isnan(PirRepRli) & ~isnan(PirSimlarx);
        y = PirSimlarx(incIds);
        x = PirRepRli(incIds);
        pcor(1) = corr(x, y);
        %         [beta, residue] = leastsquarereg(PirRepRli(incIds), PirSimlarx(incIds));
        %         plot([0 1], [residue beta+residue], '--', 'Color', [49 142 247]/255);
        X = [ones(length(x),1) x];
        b1 = X\y;
        yCalc2 = X*b1;
        if IsDisplay, plot(x,yCalc2,'--b'); end
        
        incIds = ~isnan(PirRepRli) & ~isnan(PirSimlari);
        y = PirSimlari(incIds);
        x = PirRepRli(incIds);
        pcor(2) = corr(x, y);
        X = [ones(length(x),1) x];
        b2 = X\y;
        yCalc2 = X*b2;
        tit = sprintf('%d-%d-%d_%d_%d (R=%.03G, R=%.03G)', Day, Cel, Pla,...
            expids(1), expids(2), pcor(1) , pcor(2));
        if IsDisplay
            plot(x,yCalc2,'--k');
            xlabel('Repeat reliability');
            ylabel('Correlation coeff. (paired)');
            xlim([0 1]);
            ylim([0 1]);
            title(tit);
        end
        
        %%
%         keyboard;
        
        FleNam = [MainPath '\Results\StimSponBCCorrelation\' sprintf('%d-%d-%d_%d_%d.mat',...
            Day, Cel, Pla, expids(1), expids(2))];
        save(FleNam, 'pcor', 'b1', 'b2', 'Day', 'Cel', 'Pla', 'expids', 'bctype',...
            'PirRepRli', 'RepRli', 'PirSimlari', 'PirSimlarx', 'RespSTD');
        %%
        if IsDisplay
            FolderName = BCTypeLabels{DTab(cids(1), 4)};
            FleNam = [MainPath '\Figures\ROIHeterogeneity\' FolderName '\' tit];
            saveas(gcf,[FleNam '.png']);
            close all
        end
    end
end