% identified the path distance
% estimate spontaneous and stimulus compositions
% Variables:
% (1) length contants for spontaneous and stimulus
% (2) ratio between sponatneous and stimulus
%
% additional version
% (1) fixed length contants from previous spot response analysis
% Problem (1) ratio between spon. and stim. responses may vary across movie
% clips and worse if spon. has different origins (from different terminals)
close all; clear; clc;
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis');
MainPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\';
load([MainPath 'Results\SplitROILevelSummary\SplitData_sampled_02072023.mat']);
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
%% Load Excel Experimental Document
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
[DocNum, DocStr] = xlsread(FilNam,'ScanLocMatch', 'A:H');
ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));
%%
Topic = 'Ai148_AAV-Grm6Cre';
LenConsts = exp(1:0.35:8); % in um
nLen = length(LenConsts);

nIter = 100;
NatMovPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis';
BCTerminalPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis';
savefilefolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\BCHeterogeneityLengthConstant';
%%
[Opt, Met]             = imregconfig('Multimodal');
Opt.Epsilon            = 1.5e-6;
Opt.MaximumIterations  = 300;
Opt.InitialRadius      = 6.25e-4;
%%
FilName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\Results\BCHeterogeneityLengthConstant';
load([FilName, '\EqualNumPointSimulation_052023.mat'], 'a', 'E');
E = median(E, 2:3);
getnp = @(x) round(interp1(E, a, x, 'linear', 'extrap'));
%%
q2a = @(xq) 1./(1+sqrt(xq.^-2-1));
%%
for ii = 13
    switch ii
        case 1
            Day = 12019; Cel = 4; Exp = 19; % BC5o
        case 2
            Day = 53018; Cel = 2; Exp = 27; % BC5o
        case 15
            Day = 122718; Cel = 2; Exp = 19; % BC5o
        case 3
            Day = 10119; Cel = 2; Exp = 20; % BC5i
        case 4
            Day = 10119; Cel = 4; Exp = 22; % BC5i
        case 5
            Day = 11919; Cel = 1; Exp = 18; % BC5t
        case 6
            Day = 61018; Cel = 1; Exp = 28; % BC5t
        case 7
            Day = 53018; Cel = 1; Exp = 23; % XBC
        case 16
            Day = 11019; Cel = 3; Exp = 16; %  XBC
        case 8
            Day = 92418; Cel = 2; Exp = 13; % XBC
        case 9
            Day = 10519; Cel = 3; Exp = 18;% BC6
        case 10
            Day = 12019; Cel = 7; Exp = 14; % BC6
        case 11
            Day = 10919; Cel = 3; Exp = 21; % BC7
        case 12
            Day = 12119; Cel = 3; Exp = 16; % BC7
        case 13
            Day = 11819; Cel = 1; Exp = 26; % BC89
        case 14
            Day = 12219; Cel = 1; Exp = 15; % BC89
    end
    run(sprintf('%s/SkeletonGraph/Ai148_AAVGrm6Cre_SkeletonGraph_%d_%d_%d', BCTerminalPath, Day, Cel, Exp));
    if ~exist('Nodes', 'var')
        Nodes = [Nodezxy(:, [3 2]) Nodezxy(:,1)];
    end
    %%
    nedge = size(Edges, 1);
    %     PlaneIds = 126:134;
    npage = length(PlaneIds);
    Daystr = num2str(Day);
    if length(Daystr) == 5
        Daystr = ['0' Daystr];
    end
    %% Align the scaned imaging to the 3d zstack
    % register to the specific set of plane in zstack
    % project the ROI in two 2d to the 3D structure (nodes and edges) by the
    % minimum distance
    PathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\CalciumZstack';
    SkeletonName = sprintf('%s/%s_%s_%d%03d_1_processed.tif', PathName, Topic, Daystr, Cel, Exp);
    
    %     Pla = DocNum(DocNum(:, ColName2Ind('Day')) == Day & DocNum(:, ColName2Ind('Cel')) == Cel &...
    %         DocNum(:, ColName2Ind('Exp')) == exps(f), ColName2Ind('Pla'));
    
    
    
    %% get the natural movie recordings
    MatchROIPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\ROIMatchingTable';
    Topic = 'Ai148_AAV-Grm6Cre';
    Subject = 'ROIMatchingTable';
    FileNames = dir(fullfile(MatchROIPath, '*.mat'));
    FileIds = find(contains({FileNames.name},[Topic '_ROIMatchingTable_' num2str(Day) '_' num2str(Cel)]));
    nFile = length(FileIds);
    IsDisplay = 0;
    skipids = [];
    for f = 1:nFile
        FileName = FileNames(FileIds(f)).name;
        load([MatchROIPath '\' FileName]);
        
        Words = strsplit(FileName, ['_' Subject '_']);
        Words = strsplit(Words{2}, '_');
        %     Daystr = Words{1};
        %     Day = str2double(Words{1});
        Cel = str2double(Words{2});
        Words = strsplit(Words{3}, '.');
        Pla = str2double(Words{1});
        %
        rids = DocNum(:, ColName2Ind('Day')) == Day & DocNum(:, ColName2Ind('Cel')) == Cel &...
            DocNum(:, ColName2Ind('Pla')) == Pla;
        pids = DocNum(rids & ~isnan(DocNum(:, ColName2Ind('PlaneStart'))),...
            [ColName2Ind('PlaneStart'), ColName2Ind('PlaneEnd')]);
        if ~isempty(pids)
            PlaneIds = pids(1, 1):pids(1, 2);
        end
        Info = imfinfo(SkeletonName);
        TifLink = Tiff(SkeletonName, 'r');
        anchorImg = nan(Info(1).Height, Info(1).Width, npage);
        for i = 1:npage
            TifLink.setDirectory(PlaneIds(i));
            anchorImg(:, :, i) = TifLink.read();
        end
        TifLink.close();
        anchorImg = mean(anchorImg, 3);
        anchorImg = anchorImg./quantile(anchorImg(:), 0.995);
        Fxd = anchorImg;
        
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
            close all
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
            %
            cids = find(DTab(:, 1) == Day & DTab(:, 2) == Cel & DTab(:, 5) == expids(2));
            targetD = 0.5*(DTra1(cids, :)+DTra2(cids, :));
            ampD2 = targetD(cTab(:, 2), :);
            targetD = targetD - min(targetD, [], 2);
            targetD2 = targetD./max(targetD, [], 2);
            targetD2 = targetD2(cTab(:, 2), :);
            bctype = BCTypes(DTab(cids(1), 4));
            
            % get paired path distance
            GetPairedPathDistance
            %
            cPairD1 = PairD1(cTab(:, 1), cTab(:, 1));
            cPairD2 = PairD2(cTab(:, 2), cTab(:, 2));
            PairD = 0.5*(cPairD1+cPairD2);
            dtest = 0.5*(targetD1+targetD2);
            X = dtest+1e-6*randn(size(dtest));
            Pairdist = squareform(pdist(X, 'correlation'));
            %%
            figure;
            scatter(PairD(:), Pairdist(:), 5, 'k', 'filled');
            %%
            keyboard;
        end
    end
end