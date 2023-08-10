% run CircleAnalysis_ResponseScreening072122_MorphSeg.m or
% CircleAnalysis_ResponseScreening071222_MorphSeg.m first
clear all; clc;
GroupData.Days = {...
%   1         2         3         4         5        6          7         8         9        10
'021818', '021918', '022218', '022718', '030618', '030918', '031818', '031918', '041618', '042218',...
'053018', '061018', '061118', '061918', '062518', '062618', '073018', '080218', '081518', '091918',...
'092418', '103118', '111218', '121418', '122518', '122718', '122918', '123118', '010119', '010219',...
'010319', '010519', '010919', '011019', '011319', '011819', '011919', '012019', '012119', '012219',...
'081618', '091818'...
};
GroupData.Cels{1} = 1;
GroupData.Cels{2} = 1;
GroupData.Cels{3} = 2; % 022218-1 is noisy
GroupData.Cels{4} = 1:2;
GroupData.Cels{5} = [3 5];
GroupData.Cels{6} = 2;
GroupData.Cels{7} = 3;
GroupData.Cels{8} = 2;
GroupData.Cels{9} = [1 31 32]; 
GroupData.Cels{10} = [2];
GroupData.Cels{11} = 1:2;
GroupData.Cels{12} = 1:2;
GroupData.Cels{13} = [3 11 12];
GroupData.Cels{14} = 1;
GroupData.Cels{15} = 1;
GroupData.Cels{16} = 1;
GroupData.Cels{17} = 1;
GroupData.Cels{18} = [1 3];
GroupData.Cels{19} = 3;
GroupData.Cels{20} = 1:2;
GroupData.Cels{21} = 2;
GroupData.Cels{22} = 2;
GroupData.Cels{23} = 1;
GroupData.Cels{24} = 1:2;
GroupData.Cels{25} = 4;
GroupData.Cels{26} = 2;
GroupData.Cels{27} = [2 11 12];
GroupData.Cels{28} = 1:2;
GroupData.Cels{29} = [2 4];
GroupData.Cels{30} = 1:2;
GroupData.Cels{31} = 3;
GroupData.Cels{32} = 3;
GroupData.Cels{33} = 3;
GroupData.Cels{34} = 3;
GroupData.Cels{35} = 4;
GroupData.Cels{36} = 1;
GroupData.Cels{37} = [1 5];
GroupData.Cels{38} = [4 7];
GroupData.Cels{39} = 2:3;
GroupData.Cels{40} = [1 5];
GroupData.Cels{41} = 5;
GroupData.Cels{42} = 4;
GroupData.Exps{1, 1} = [3 8 9 12]; 
GroupData.Exps{2, 1} = [5 10]; 
GroupData.Exps{3, 1} = [1:3]; 
GroupData.Exps{3, 2} = [1 6];
GroupData.Exps{4, 1} = [2]; % 1 10 bad
GroupData.Exps{4, 2} = [2];
GroupData.Exps{5, 3} = [1 11];
GroupData.Exps{5, 5} = [1];
GroupData.Exps{6, 2} = [1 2 11 13];
GroupData.Exps{7, 3} = [1 2 4 13 15 22 23];
GroupData.Exps{8, 2} = [2 3 4 15 16];
GroupData.Exps{9, 1} = [9 15 21];
GroupData.Exps{9, 31} = [21];
GroupData.Exps{9, 32} = [21];
GroupData.Exps{10, 2} = [12 13 21];
GroupData.Exps{10, 11} = [4 18];
GroupData.Exps{10, 12} = [4 18];
GroupData.Exps{11, 1} = [2 3 16 22];
GroupData.Exps{11, 2} = [1 2 3 4 17 18];
GroupData.Exps{12, 1} = [1 2 3 26 27];
GroupData.Exps{12, 2} = [2 3];
GroupData.Exps{13, 3} = [3 6];
GroupData.Exps{13, 11} = [1 2 3];
GroupData.Exps{13, 12} = [1 2 3];
GroupData.Exps{14, 1} = [2 3 14 15];
GroupData.Exps{15, 1} = [2 3 15 20 21];
GroupData.Exps{16, 1} = [2 3 4];
GroupData.Exps{17, 1} = [10 13 19];
GroupData.Exps{18, 1} = [2 3 4 5 16 17 22];
GroupData.Exps{18, 3} = [2 4 13 14];
GroupData.Exps{19, 3} = [14 19]; %5 6 21
GroupData.Exps{20, 1} = [13 19];
GroupData.Exps{20, 2} = [2 3 4];% 12 16
GroupData.Exps{21, 2} = [3 4 12];
GroupData.Exps{22, 2} = [5 6 13 17 21];
GroupData.Exps{23, 1} = [8 12 16 20];
GroupData.Exps{24, 1} = [10 14 18 19];
GroupData.Exps{24, 2} = [4 9 13];
GroupData.Exps{25, 4} = [3 8 13];
GroupData.Exps{26, 2} = [4 12 16];
GroupData.Exps{27, 2} = [8 13];
GroupData.Exps{27, 11} = [3 8 12 16];
GroupData.Exps{27, 12} = [3 8 12 16];
GroupData.Exps{28, 1} = [8 12];
GroupData.Exps{28, 2} = [8 12 16];
GroupData.Exps{29, 2} = [4 9 13 15];
GroupData.Exps{29, 4} = [9 11 15 19];
GroupData.Exps{30, 1} = [9 14 17];
GroupData.Exps{30, 2} = [6 13 19 22];
GroupData.Exps{31, 3} = [8 12 14];
GroupData.Exps{32, 3} = [9 13 16];
GroupData.Exps{33, 3} = [16 17];
GroupData.Exps{34, 3} = [3 4 10];% 9(stimulus not aligned)
GroupData.Exps{35, 4} = [9 14];
GroupData.Exps{36, 1} = [3 8 13 14 15 20];% 21 24 are the trunk of axon
GroupData.Exps{37, 1} = [9 14 16];
GroupData.Exps{37, 5} = [2 3 8 13 14];
GroupData.Exps{38, 4} = [4 9 14 15];
GroupData.Exps{38, 7} = [3 4 9];
GroupData.Exps{39, 2} = [10 15];
GroupData.Exps{39, 3} = [4 9 14]; % Set Align2Stimuli = 1
GroupData.Exps{40, 1} = [3 4 9 14];
GroupData.Exps{40, 5} = [12 17];
GroupData.Exps{41, 5} = [5 13 18 20];
GroupData.Exps{42, 4} = [16];

Topic = 'Ai148_AAV-Grm6Cre';
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
RespWindow = 4:33;% 4:33
CellIds = 1;
for dl = 1:NumDays
    if strcmp(GroupData.Days{dl}, '030918') % | strcmp(GroupData.Days{dl}, '022218')
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
%             if strcmpi(GroupData.Days{dl}, '010119') && Cels(cl) == 2 && Exps(e) == 4
%                 keyboard;
%             end
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
% keyboard;
%%
EvaluateHeterogeneityofCenterSurround_MatchingTable
keyboard;









