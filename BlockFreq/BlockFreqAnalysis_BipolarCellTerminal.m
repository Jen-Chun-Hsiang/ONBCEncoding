clear all; clc;
GroupData.Days = {...
    %   1         2         3         4        5         6         7         8         9        10
    '053018', '061118', '061918', '062518', '062618', '073018', '080218', '081518', '081618', '091918',...
    '092418', '103118', '111218', '121418', '122518', '122718', '122918', '123118', '010119', '010219',...
    '010319', '010519', '010919', '011019', '011319', '011819', '011919', '012019', '012119', '012219',...
    '031918', '042218', '061018', '091818'};
GroupData.Cels{1} = 1:2;
GroupData.Cels{2} = [3 11 12];
GroupData.Cels{3} = 1;
GroupData.Cels{4} = 1;
GroupData.Cels{5} = [1];
GroupData.Cels{6} = 1;
GroupData.Cels{7} = [1 3];
GroupData.Cels{8} = 3;
GroupData.Cels{9} = 5;
GroupData.Cels{10} = 1:2;
GroupData.Cels{11} = 2; % traces look bad, better start at 7
GroupData.Cels{12} = 2;
GroupData.Cels{13} = 1;
GroupData.Cels{14} = 1:2;
GroupData.Cels{15} = 4;
GroupData.Cels{16} = 2;
GroupData.Cels{17} = [12 2]; %1 has XBC & BC5o
GroupData.Cels{18} = 1:2;
GroupData.Cels{19} = [2 4];
GroupData.Cels{20} = 1:2;
GroupData.Cels{21} = 3;
GroupData.Cels{22} = 3;
GroupData.Cels{23} = 3;
GroupData.Cels{24} = 3;
GroupData.Cels{25} = 4;
GroupData.Cels{26} = 1;
GroupData.Cels{27} = [1 5];
GroupData.Cels{28} = [4 7];
GroupData.Cels{29} = [2 3];
GroupData.Cels{30} = [1 5];
GroupData.Cels{31} = [2];
GroupData.Cels{32} = [2];
GroupData.Cels{33} = [1];
GroupData.Cels{34} = [4];
GroupData.Exps{1, 1} = [5:11];
GroupData.Exps{1, 2} = [5:8];
GroupData.Exps{2, 11} = [4 5 6 7 16 17 18 19];
GroupData.Exps{2, 12} = [6 7 18 19];
GroupData.Exps{2, 3} = [4 5];
GroupData.Exps{3, 1} = [4:7];
GroupData.Exps{4, 1} = 4:10;
GroupData.Exps{5, 1} = 5:9;
GroupData.Exps{6, 1} = [4:6 11 12];
GroupData.Exps{7, 1} = [6:9];
GroupData.Exps{7, 3} = [5:9 18:20];
GroupData.Exps{8, 3} = 7:10;
GroupData.Exps{9, 5} = 6:9;
GroupData.Exps{10, 1} = 5:8;
GroupData.Exps{10, 2} = 5:8;
GroupData.Exps{11, 2} = 5:8;
GroupData.Exps{12, 2} = 7:12;
GroupData.Exps{13, 1} = 4:7;
GroupData.Exps{14, 1} = 6:9;
GroupData.Exps{14, 2} = 5:8; %5-7 are bad
GroupData.Exps{15, 4} = 4:7; % 4 doesn't look good
GroupData.Exps{16, 2} = 5:8;
GroupData.Exps{17, 12} = 4:7;
GroupData.Exps{17, 2} = 4:7;
GroupData.Exps{18, 1} = 4:7;
GroupData.Exps{18, 2} = 4:7;
GroupData.Exps{19, 2} = 5:8;
GroupData.Exps{19, 4} = 5:8;
GroupData.Exps{20, 1} = 5:8;
GroupData.Exps{20, 2} = 7:12;
GroupData.Exps{21, 3} = 4:7;
GroupData.Exps{22, 3} = 4:7;
GroupData.Exps{23, 3} = 7:10; %bad in 7, 8
GroupData.Exps{24, 3} = 5:8;
GroupData.Exps{25, 4} = 5:8;
GroupData.Exps{26, 1} = 4:7; % bad in 6
GroupData.Exps{27, 1} = 5:8;
GroupData.Exps{27, 5} = 4:7;
GroupData.Exps{28, 4} = 5:8;
GroupData.Exps{28, 7} = 5:8;
GroupData.Exps{29, 2} = 5:9;
GroupData.Exps{29, 3} = 5:8;
GroupData.Exps{30, 1} = 5:8;
GroupData.Exps{30, 5} = 7:11;
GroupData.Exps{31, 2} = [8 9];
GroupData.Exps{32, 2} = [5 6];
GroupData.Exps{33, 1} = [15 16 18 19];
GroupData.Exps{34, 4} = [5:8];

Topic = 'Ai148_AAV-Grm6Cre';
ImgSaveFolder = 'BCTerminal_031219';
ResultType = 'AutoCluster_';

%% Load Excel Experimental Document
FilNam = '//storage1.ris.wustl.edu/kerschensteinerd/Active/Emily/Recording_Document/BCTerminal_DataMatrix.xlsx';
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
ExpFreqResp = [];
ExpFreqRespNrm = [];
ExpFreqRespVar = [];
ExpFreqRespNrmVar = [];
ExpFreqPw = [];
ROIFreqPw = [];
ROITable = [];
ROIQuality = [];
DMat = [];
ROIFeatures = [];
ROINrmResizeFac = 1;
ROIFreqResp = [];
ROISpotSize = [];
CellIds = 1;
ExpQuality = [];
ExpSelfC = [];
for dl = 1:NumDays
    fprintf('%d/%d \n', dl, NumDays);
%     if strcmp(GroupData.Days{dl}, '062618')
%         continue;
%     end
    Cels = GroupData.Cels{dl};
    NumCel = length(Cels);
    for cl = 1:NumCel
        Exps = GroupData.Exps{dl, Cels(cl)};
        NumExp = length(Exps);
        
        Cell.MeanNrm = [];
        Cell.MeanRawNrm = [];
        ExpSpotSize = [];
        
        for e = 1:NumExp
            load(sprintf('./Results/%s_Result_%s%s_%d%03d.mat',Topic, ResultType, GroupData.Days{dl}, Cels(cl), Exps(e)), 'Data', 'Fz', 'Dmt', 'ROIReliability');
            %             load(sprintf('./ProcessedData/%s_Result_%s%s_%d%03d.mat',Topic, ResultType, GroupData.Days{dl}, Cels(cl), Exps(e)));
            CurDocNum = DocNum(DocNum(:,ColName2Ind('Day')) == str2num(GroupData.Days{dl}) & DocNum(:, ColName2Ind('Cel')) == Cels(cl) & DocNum(:,ColName2Ind('Exp')) == Exps(e), :);
            Cell.MeanNrm = cat(3, Cell.MeanNrm, Data.MeanNrm);
            Cell.MeanRawNrm = cat(3, Cell.MeanRawNrm, Data.MeanRawNrm);
            NumROI = size(Data.ROIsNrm{1, 1}, 1);
            QdatFileName = sprintf('Quality/%s_SelfCorrQuality_%s_%d%03d.mat', Topic, GroupData.Days{dl}, Cels(cl), Exps(e));
            load(QdatFileName, 'CorrV');
            assert(size(CorrV, 1) == NumROI);
            ROITable = [ROITable; repmat(CurDocNum([1:3 7]), NumROI, 1) (1:NumROI)' ones(NumROI, 1)*CellIds,...
                CorrV];
            NumDmt = length(Dmt);
            cQuality = mean(reshape(ROIReliability, NumROI, NumDmt, []), 3);
            ROIQuality = [ROIQuality; cQuality];
            Features = [];
            cDMat = [];
            for i = 1:NumDmt
%                 ROIsNrm = Data.ROIsRaw{i, 1};
                ROIsNrm = Data.ROIsNrm{i, 1};
                ROIAll = Data.ROIsMatAll{i, 1};
                nRepeat = size(ROIAll, 2);
                cDMat = cat(2, cDMat, reshape(permute(ROIAll, [2 3 1]), nRepeat, 1, [], NumROI));
                Features = cat(1, Features, reshape(ROIsNrm', 1, [], NumROI));
            end
            DMat = cat(4, DMat, cDMat); %[repeat, fz, trace, roi]
            ROIFeatures = cat(3, ROIFeatures, Features);
            clear Features ROIsNrm
            load(sprintf('./Results/%s_Result_FreqAnalysis_%s_%d%03d.mat',Topic, GroupData.Days{dl}, Cels(cl), Exps(e)),...
                'MainFreqResp', 'NrmMainFreqResp', 'FreqPw', 'fl', 'SpotSize', 'Freqs');
            ExpFreqResp = [ExpFreqResp; mean(MainFreqResp, 1)];
%             NrmMainFreqResp = NrmMainFreqResp./NrmMainFreqResp(:, 1);
            NrmMainFreqResp = (NrmMainFreqResp-NrmMainFreqResp(:, 1))./range(NrmMainFreqResp, 2);
            ROIFreqResp = [ROIFreqResp; NrmMainFreqResp];
            ROISpotSize = [ROISpotSize; SpotSize*ones(NumROI, 1)];
            ExpFreqRespNrm = [ExpFreqRespNrm; mean(NrmMainFreqResp, 1)];
            ExpFreqRespVar = [ExpFreqRespVar; mad(MainFreqResp, [], 1)];
            ExpFreqRespNrmVar = [ExpFreqRespNrmVar; mad(NrmMainFreqResp, [], 1)];
            ROIFreqPw = cat(3, ROIFreqPw, FreqPw);
            ExpFreqPw = cat(3, ExpFreqPw, mean(FreqPw, 3, 'omitnan'));
            ExpSpotSize = [ExpSpotSize SpotSize];
            ExpQuality = [ExpQuality; median(cQuality, 1)];
            ExpSelfC = [ExpSelfC; median(CorrV)];
        end
        MeanNrm = cat(3, MeanNrm, mean(Cell.MeanNrm, 3));
        MeanRawNrm = cat(3, MeanRawNrm, mean(Cell.MeanRawNrm, 3));
        ExpMeanNrm = cat(3, ExpMeanNrm, Cell.MeanNrm);
        ExpMeanRawNrm = cat(3, ExpMeanRawNrm, Cell.MeanRawNrm);
        IdTable = [IdTable; CurDocNum([1:2 7])];
        ExpTable = [ExpTable; repmat(CurDocNum([1:2 7]), NumExp,1),  Exps', ExpSpotSize'];
        CellIds = CellIds+1;
    end
end
ROITable = [ROITable, ROISpotSize];
clear ROISpotSize

% Check
if size(IdTable, 1) ~= size(MeanNrm, 3)
    error('Wrong!');
end
% if size(ROITable, 1) ~= size(ROIFeatures, 1)
%     error('Wrong!');
% end
%% Correct stimulus doc
ExpTable(ExpTable(:, 1) == 31918 & ExpTable(:, 2) == 2 & ExpTable(:, 4) == 8, 5) = 800;
ExpTable(ExpTable(:, 1) == 31918 & ExpTable(:, 2) == 2 & ExpTable(:, 4) == 9, 5) = 300;
ExpTable(ExpTable(:, 1) == 42218 & ExpTable(:, 2) == 2 & ExpTable(:, 4) == 5, 5) = 150;
ExpTable(ExpTable(:, 1) == 42218 & ExpTable(:, 2) == 2 & ExpTable(:, 4) == 6, 5) = 800;

ROITable(ROITable(:, 1) == 31918 & ROITable(:, 2) == 2 & ROITable(:, 3) == 8, 8) = 800;
ROITable(ROITable(:, 1) == 31918 & ROITable(:, 2) == 2 & ROITable(:, 3) == 9, 8) = 300;
ROITable(ROITable(:, 1) == 42218 & ROITable(:, 2) == 2 & ROITable(:, 3) == 5, 8) = 150;
ROITable(ROITable(:, 1) == 42218 & ROITable(:, 2) == 2 & ROITable(:, 3) == 6, 8) = 800;

%% Parameter should have saved
StmDur = 2.5; % in second
BasDur = 0.5;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm;
T = round(linspace(-BasDurFrm/Fz, StmDurFrm/Fz, length(1:1:SpnDurFrm))*100)/100;

NumDmt = length(Dmt);
%%
% keyboard;
%% Get FFT Power
[uniCel, cids] = unique(ROITable(:, [1:3 5]), 'rows');
uniCel = [uniCel ROITable(cids, [4 8]) ROIQuality(cids, :)];
nuniCel = size(uniCel, 1);
CelMeanNrm = nan(size(ROIFeatures, 1), size(ROIFeatures, 2), nuniCel);
CelDMat = nan(size(DMat, 1), size(DMat, 2), size(DMat, 3), nuniCel);
for i = 1:nuniCel
    cids = ROITable(:, 1) == uniCel(i, 1) & ROITable(:, 2) == uniCel(i, 2) &...
        ROITable(:, 3) == uniCel(i, 3) & ROITable(:, 5) == uniCel(i, 4);
    CelMeanNrm(:, :, i) = mean(ROIFeatures(:, :, cids), 3);
    CelDMat(:, :, :, i) = mean(DMat(:, :, :, cids), 4);
end
CelMeanNrm = CelMeanNrm-min(CelMeanNrm, [], 1:2);
CelMeanNrm = CelMeanNrm./range(CelMeanNrm, 1:2);
CelMeanNrm = CelMeanNrm-mean(CelMeanNrm(:, 1:BasDurFrm, :), 2);
cids = any(isnan(CelMeanNrm), 1:2);
uniCel(cids, :) = [];
CelMeanNrm(:, :, cids) = [];
CelDMat(:, :, :, cids) = [];
%%
cids = uniCel(:, 5) == 57;
figure; imagesc(mean(CelMeanNrm(:, :, cids), 3)); colorbar;
%%
UpFac = 5;
TI = linspace(T(1), T(end), length(T)*UpFac);
n = length(TI);  
fs = Fz*UpFac;
f = (0:n-1)*(fs/n); 
fl = f(f<Fz/2);
nf = length(fl);
CellFreqPw = nan(NumDmt, nf, nuniCel);
for i = 1:nuniCel
    for j = 1:NumDmt
        YI = interp1(T, CelMeanNrm(j, :, i), TI, 'pchip');
%         YI = abs(fft(YI))/n;
        YI = abs(fft(YI)).^2/n;
        CellFreqPw(j, :, i) = YI(1:nf);
    end
end
Ti = T;
%%
BCTypes = [5 50 51 57 58 6 7 89];
BCTypelabels = {'BC5', 'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
nType = numel(BCTypes);
nL = 10;
fsets = {2:3, 4, 6:7, 12:13, 23:24, 44:46};
nfset = length(fsets);
TargetMatrix = nan(nuniCel, nfset);
for i = 1:nfset
    TargetMatrix(:, i) = reshape(mean(log10(CellFreqPw(i, fsets{i}, :)), 2), [], nuniCel)';
end
cids = ismember(uniCel(:, 6), [150 300]);
TargetMatrix = TargetMatrix(cids, 1:5);
Dtab =  uniCel(cids, 1:11);
Dtrace = CelMeanNrm(1:5, :, cids);
Datrace = CelDMat(:, 1:5, :, cids);
cids = any(TargetMatrix < -10, 2);
TargetMatrix(cids, :) = [];
Dtab(cids, :) = [];
Dtrace(:, :, cids) = [];
Datrace(:, :, :, cids) = [];
% run @BlockFreqAnalysis_ROIIndividualTraces
%% Get FFT Power
[uniCel, cids] = unique(ExpTable(:, [1:2 5]), 'rows');
uniCel = [uniCel ExpTable(cids, 3)];
nuniCel = size(uniCel, 1);
CelMeanNrm = nan(size(ExpMeanNrm, 1), size(ExpMeanNrm, 2), nuniCel);
for i = 1:nuniCel
    cids = ExpTable(:, 1) == uniCel(i, 1) & ExpTable(:, 2) == uniCel(i, 2) & ExpTable(:, 5) == uniCel(i, 3);
    CelMeanNrm(:, :, i) = mean(ExpMeanNrm(:, :, cids), 3);
end
CelMeanNrm = CelMeanNrm-min(CelMeanNrm, [], 1:2);
CelMeanNrm = CelMeanNrm./range(CelMeanNrm, 1:2);
CelMeanNrm = CelMeanNrm-mean(CelMeanNrm(:, 1:BasDurFrm, :), 2);
% CelMeanNrm = CelMeanNrm./mean(CelMeanNrm(:, 1:BasDurFrm, :), 1:2);
%%
UpFac = 5;
TI = linspace(T(1), T(end), length(T)*UpFac);
n = length(TI);  
fs = Fz*UpFac;
f = (0:n-1)*(fs/n); 
fl = f(f<Fz/2);
nf = length(fl);
CellFreqPw = nan(NumDmt, nf, nuniCel);
CellFreqPh = nan(NumDmt, nf, nuniCel);
for i = 1:nuniCel
    for j = 1:NumDmt
        YI = interp1(T, CelMeanNrm(j, :, i), TI, 'pchip');
        fx = fft(YI);
%         fz = fftshift(fx);
%         fz(abs(fz) < 1e-6) = 0;
        YI = abs(fx).^2/n;
        CellFreqPw(j, :, i) = YI(1:nf);
        CellFreqPh(j, :, i) = angle(fx(1:nf))*180/pi;
    end
end