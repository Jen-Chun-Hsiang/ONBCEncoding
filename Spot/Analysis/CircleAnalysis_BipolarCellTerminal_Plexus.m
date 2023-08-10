% run CircleAnalysis_ResponseScreening072122_MorphSeg.m or
% CircleAnalysis_ResponseScreening071222_MorphSeg.m first
clear all; clc;
GroupData.Days = {...
%   1         2        3         4         5         6          7        8         9         10
'092519', '093019', '100319', '100719', '103118', '110318', '110818',...
};
GroupData.Cels{1} = 1:3;
GroupData.Cels{2} = [1 3]; %2 doesn't have good response to spot
GroupData.Cels{3} = 1:3;
GroupData.Cels{4} = 2;
GroupData.Cels{5} = 3:4;
GroupData.Cels{6} = 1:2;
GroupData.Cels{7} = 1:3;
GroupData.Exps{1, 1} = [4 7 10];
GroupData.Exps{1, 2} = [1 5 9];
GroupData.Exps{1, 3} = 5;
GroupData.Exps{2, 1} = [5 9];
GroupData.Exps{2, 3} = [5 9];
GroupData.Exps{3, 1} = [5 9];
GroupData.Exps{3, 2} = 10;
GroupData.Exps{3, 3} = [5 9];
GroupData.Exps{4, 2} = 9;
GroupData.Exps{5, 3} = [1 9];
GroupData.Exps{5, 4} = [1 2 3 8 13 17];
GroupData.Exps{6, 1} = [4 12 19];%
GroupData.Exps{6, 2} = [2 3 8 9 13 14 18];
GroupData.Exps{7, 1} = [1:3 8 12 16];
GroupData.Exps{7, 2} = [2 3 8 10 16];
GroupData.Exps{7, 3} = [2 3 8 12 16];%
Topic = 'Ai148_Grm6Cre';
ImgSaveFolder = 'BCTerminal_070822';
ResultType = 'AutoCluster_';

%% Load Excel Experimental Document
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
[DocNum, DocStr] = xlsread(FilNam,'Sheet1', 'A:G');
ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));

%%
NumDays = length(GroupData.Days);
IdTable = [];
MeanNrm = [];
MeanRawNrm = [];
ExpTable = [];
ExpMeanNrm = [];
ExpMeanRawNrm = [];
ROITable = [];
ROIFeatures = [];
ROIRaw = [];
ROIReliab = [];
ROIDmtAmp = [];
clear RmvIdSet 
ROINrmResizeFac = 1;
RespWindow = 4:33;
CellIds = 1;
for dl = 1:NumDays
    if strcmp(GroupData.Days{dl}, '030918') | strcmp(GroupData.Days{dl}, '022218')
        continue;
    end
    Cels = GroupData.Cels{dl};
    NumCel = length(Cels);
    for cl = 1:NumCel
        Exps = GroupData.Exps{dl, Cels(cl)};
        NumExp = length(Exps);
        
        Cell.MeanNrm = [];
        Cell.MeanRawNrm = [];        
        for e = 1:NumExp
            load(sprintf('./Results/%s_Result_%s%s_%d%03d.mat',Topic, ResultType, GroupData.Days{dl}, Cels(cl), Exps(e)), 'Data', 'Fz', 'Dmt', 'RmIds');
%             load(sprintf('./ProcessedData/%s_Result_%s%s_%d%03d.mat',Topic, ResultType, GroupData.Days{dl}, Cels(cl), Exps(e)));
            CurDocNum = DocNum(DocNum(:,ColName2Ind('Day')) == str2num(GroupData.Days{dl}) & DocNum(:, ColName2Ind('Cel')) == Cels(cl) & DocNum(:,ColName2Ind('Exp')) == Exps(e), :);
            Cell.MeanNrm = cat(3, Cell.MeanNrm, Data.MeanNrm);
            Cell.MeanRawNrm = cat(3, Cell.MeanRawNrm, Data.MeanRawNrm);
            NumROI = size(Data.ROIsNrm{1, 1}, 1);
            ROITable = [ROITable; repmat(CurDocNum([1:3 7]), NumROI, 1) (1:NumROI)' ones(NumROI, 1)*CellIds];
            NumDmt = length(Dmt);
            Features = [];
            Amp = NaN(NumROI, NumDmt);
            Reliabilities = NaN(NumROI, NumDmt);
            Raw = [];
            for i = 1:NumDmt
                ROIsNrm = Data.ROIsNrm{i, 1};
                ROIsReliab = max([Data.ROIsReliab{i, 1} Data.ROIsReliab{i, 1}], [], 2);
                DmtFeatures = zeros(NumROI, length(RespWindow)/ROINrmResizeFac);
                for r = 1:NumROI
                    DmtFeatures(r, :) = mean(reshape(ROIsNrm(r, RespWindow), ROINrmResizeFac, []), 1);
                end
                Features = [Features DmtFeatures];
                Reliabilities(:, i) = ROIsReliab;
                Amp(:, i) = mean(ROIsNrm(:, 4:17), 2);
                Raw = [Raw Data.ROIsRaw{i, 1}];
            end     
            ROIRaw = [ROIRaw; mean(Raw, 2), median(Raw, 2)];
            ROIFeatures = [ROIFeatures; Features];
            ROIReliab = [ROIReliab; Reliabilities];
            ROIDmtAmp = [ROIDmtAmp; Amp];
            RmvIdSet{CellIds, Exps(e)} = RmIds;
        end
        MeanNrm = cat(3, MeanNrm, mean(Cell.MeanNrm, 3));
        MeanRawNrm = cat(3, MeanRawNrm, mean(Cell.MeanRawNrm, 3));
        ExpMeanNrm = cat(3, ExpMeanNrm, Cell.MeanNrm);
        ExpMeanRawNrm = cat(3, ExpMeanRawNrm, Cell.MeanRawNrm);
        IdTable = [IdTable; CurDocNum([1:2 7])];        
        ExpTable = [ExpTable; repmat(CurDocNum([1:2 7]), NumExp,1),  Exps'];
        CellIds = CellIds+1;
    end
end

% Check
if size(IdTable, 1) ~= size(MeanNrm, 3)
    error('Wrong!');
end
if size(ROITable, 1) ~= size(ROIFeatures, 1)
    error('Wrong!');
end
keyboard;
%%
EvaluateHeterogeneityofCenterSurround_MatchingTable_Plexus



