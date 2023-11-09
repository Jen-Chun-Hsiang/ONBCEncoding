% Goal Visualize if spike events is in an individual trials
% combine @CircleAnalysis_BipolarCellTerminal  &
% @CircAnalysis_ResponseScreening072122_MorphSeg &
% @ReanalyzeVG3Circle_withoutintensitycorrection

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
NumDays = length(GroupData.Days);
Fz = 9.47;

celltype = [50 51 57 58 6 7 89];
celltypelabel = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC8/9'};
%% Load Excel Experimental Document
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
[DocNum, DocStr] = xlsread(FilNam,'Sheet1', 'A:G');
ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));
%%
StmDur = 1.75; % in second
BasDur = 0.25;
StmDurFrm = round(StmDur*Fz);
BasDurFrm = round(BasDur*Fz);
SpnDurFrm = StmDurFrm+BasDurFrm+1;

%%
IsDisplayAlign = 0;
Folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Circle\Analysis\Translation';
if exist('ONBCProj_110623_revise.mat', 'file') & ~IsDisplayAlign %'ONBCProj_102423_revise.mat'
    load('ONBCProj_110623_revise.mat');
    clc
    disp('Loading...');
else
    % Sig = nan(200, 1600);
    Typ = nan(200, 1);
    % DataROIs = nan(19, 7, 2, 200);
    DataROIsDf = nan(20, 5, 7, 2, 200); % [time, repeat, size, contrast, roi]
    ROItab = nan(200, 4);
    % DataROIsDff = nan(19, 7, 2, 200);
    ii = 1;
    iex = 1;
    nonexistfile = [];
    for dl = 1:NumDays
        Cels = GroupData.Cels{dl};
        NumCel = length(Cels);
        for cl = 1:NumCel
            Exps = GroupData.Exps{dl, Cels(cl)};
            if dl == 39 && cl == 2
                Align2Stimuli = 1;
            else
                Align2Stimuli = 0;
            end
            for e = 1:length(Exps)
                FileName = sprintf('%s_Translation_%s_%d%03d.mat', Topic, GroupData.Days{dl}, Cels(cl), Exps(e));
                if ~exist([Folder '/' FileName], 'file')
                    nonexistfile{iex} = FileName;
                    iex = iex +1;
                    continue
                else
                    load([Folder '/' FileName]);
                end
                % Get type
                CurDocNum = DocNum(DocNum(:,ColName2Ind('Day')) == str2num(GroupData.Days{dl}) & DocNum(:, ColName2Ind('Cel')) == Cels(cl) & DocNum(:,ColName2Ind('Exp')) == Exps(e), :);
                CurType = CurDocNum(ColName2Ind('Type'));
                Typ(ii) = CurType;
                ROItab(ii, :) = [str2num(GroupData.Days{dl}) Cels(cl) Exps(e) CurType];
                I = medfilt3(I);
                stdI = std(I, [], 3);
                % determine the threshold for ROI regions
                mask = stdI > 3;
                % normalized to the number of pixels of the ROIs
                Resp = reshape(I, [], size(I, 3))'*mask(:)/sum(mask(:));
                %             Sig(ii, :) = Resp;
                % Get aligned signals
                Pixel_x = size(I, 1);
                Pixel_y = size(I, 2);
                NumFrm  = size(I, 3);
                % Load Stimulus
                load(sprintf('./DataStimulus/%s_%s_%d_%d.mat',Topic, GroupData.Days{dl}, Cels(cl), Exps(e)));
                OUT = circleOut;
                
                % Projector Imaging
                FilNam = sprintf('%s_%s_%d%03d_4',Topic, GroupData.Days{dl}, Cels(cl), Exps(e) );
                if exist(['./Preprocessed/MatFile/' FilNam '.mat'], 'file')
                    disp('Loading preprocessed prejector imaging dataset...');
                    PI = load(['./Preprocessed/MatFile/' FilNam '.mat']);
                    PI = PI.I;
                else
                    keyboard;
                end
                
                % Stimulus - convert to a time line
                ImgS = nan(Pixel_x*Pixel_y*NumFrm,1);
                for i = 1:NumFrm
                    Img = PI(:,:,i);
                    Img = Img';
                    ImgS((i-1)*Pixel_x*Pixel_y+1 : i*Pixel_x*Pixel_y) = Img(:);
                end
                
                %% Time lock to the stimulus
                FilNam = sprintf('./Preprocessed/TriggerWindow/%s_TriggerWindow_%s_%d%03d.mat',Topic, GroupData.Days{dl}, Cels(cl), Exps(e));
                if exist(FilNam, 'file')
                    load(FilNam);
                else
                    keyboard
                end
                
                Ind = find(ImgS > Thr);
                if Align2Stimuli
                    IndStr = Ind(Ind > Jag);
                    StrPnt = IndStr(1);
                    IndEnd = Ind(Ind > Jag2);
                    EndPnt = IndEnd(1);
                else
                    if exist('Jag2','var')
                        Ind = Ind(Ind > Jag & Ind <Jag2);
                        disp('Jag2');
                    else
                        Ind = Ind(Ind > Jag);
                    end
                    StrPnt = Ind(1);
                    EndPnt = Ind(end);
                end
                if isfield(OUT, 'Proj_End')
                    TrgType = 2; % 1: 12 pulses
                elseif isfield(OUT, 'stop')
                    TrgType = 1; % 1: 12 pulses
                end
                if Align2Stimuli
                    TrgType = 3;
                end
                    
                NumTri = numel(OUT.stim);
                switch TrgType
                    case 1
                        TimSpn = OUT.stop - OUT.start;
                    case 2
                        TimSpn = OUT.Proj_End - OUT.Proj_Start;
                    case 3
                        MaxDia = OUT.DataTable(OUT.DataTable(:, 2)>= max(OUT.DataTable(:, 2)), :);
                        TimSpn = MaxDia(end-1, 3) - MaxDia(1, 3);
                end
                StrFrm = round(StrPnt/(Pixel_x*Pixel_y)); % The frame number that links to trigger "Start"
                EndFrm = round(EndPnt/(Pixel_x*Pixel_y));
                
                switch TrgType
                    case 1
                        NumSecPfr = TimSpn/(EndFrm-StrFrm+1);
                        StrTim = OUT.start - NumSecPfr*(mod(StrPnt, (Pixel_x*Pixel_y))/(Pixel_x*Pixel_y));
                    case 2 %for turning on the "StimOn" & "StimOff"
                        NumSecPfr = TimSpn/(EndFrm-StrFrm+1);
                        StrTim = OUT.Proj_Start - NumSecPfr*(mod(StrPnt, (Pixel_x*Pixel_y))/(Pixel_x*Pixel_y));
                    case 3
                        NumSecPfr = TimSpn/((EndPnt-StrPnt)/(Pixel_x*Pixel_y));
                        TimeErrGap = mean([MaxDia(1, 3) - StrFrm*NumSecPfr,...
                            MaxDia(end-1, 3) - EndFrm*NumSecPfr]);
                        StrTim = OUT.Proj_Start-TimeErrGap;
                        StrFrm = round(StrTim/NumSecPfr);
                        EndFrm = round((OUT.Proj_End-TimeErrGap)/NumSecPfr);
                        MaxEndFrm = floor(length(ImgS)/(Pixel_x*Pixel_y));
                        if EndFrm > MaxEndFrm
                            EndFrm = MaxEndFrm;
                        end
                    otherwise
                        error('No such Trigger Type!');
                end
                
                % Combine On and Off Trials
                Dmt = [OUT.diameter; OUT.diameter];
                Cts = [ones(size(OUT.diameter)); -ones(size(OUT.diameter))];
                TriTim = [OUT.stim(:,:,1); OUT.stim(:,:,2)];
                TriTim = TriTim(:);
                Dmt = Dmt(:);
                Cts = Cts(:);
                [TriTim, IdT] = sort(TriTim);
                Dmt = Dmt(IdT);
                Cts = Cts(IdT);
                
                NumSigTrn = length(StrFrm:EndFrm);
                StmTim = TriTim(:) - StrTim;
                StmFrm = round(StmTim/NumSecPfr);
                DmtMat = Dmt';
                CtsMat = Cts';
                
                Dmt = unique(DmtMat);
                Cts = unique(CtsMat);
                Cts = flip(Cts);
                NumDmt = length(Dmt);
                NumCts = length(Cts);
                
                ROISigTrn = Resp(StrFrm:EndFrm);
                DmtVecFrm = zeros(1, NumSigTrn);
                CtsVecFrm = zeros(1, NumSigTrn);
                DmtVecFrm(StmFrm) = DmtMat;
                CtsVecFrm(StmFrm) = CtsMat;
                
                %
                stim = squeeze(mean(PI, 1:2));
                stim = stim(StrFrm:EndFrm);
                
                %%
                
                if IsDisplayAlign
                    %                 if Typ(ii) ~= 58 && Typ(ii) ~= 5 && ii > 109
                    close all
                    figure;hold on;
                    T = (0:NumSigTrn-1)/Fz;
%                     T = 1:NumSigTrn;
                    sig = ROISigTrn';
                    sig = sig - mean(ROISigTrn(35:266));
                    sig = sig/max(sig);
                    plot(T', sig, 'k');
                    RngRawSig = range(sig);
                    NrmDmtVecFrm = zeros(size(DmtVecFrm));
                    for kl = 1:7
                        ids = find(DmtVecFrm == Dmt(kl));
                        for i = 1:5
                            NrmDmtVecFrm(ids((i-1)*2+1):ids(2*i)) = kl;
                        end
                    end
%                     NrmDmtVecFrm = NrmDmtVecFrm*RngRawSig/range(NrmDmtVecFrm)- 1.5*RngRawSig;
                    NrmDmtVecFrm = NrmDmtVecFrm/range(NrmDmtVecFrm);
                    NrmDmtVecFrm(isnan(NrmDmtVecFrm)) = 0;
                    plot(T', 0.2*NrmDmtVecFrm-0.35, 'b');
%                     NrmCtsVecFrm = CtsVecFrm*RngRawSig/range(CtsVecFrm)- 1.5*RngRawSig;
%                     plot(T', 0.2*NrmCtsVecFrm);
%                     NrmStimFrm = stim*RngRawSig/range(stim)-1*RngRawSig;
%                     plot(T', 0.1*NrmStimFrm);
                    title(sprintf('%s-%d-%d (%s).mat', GroupData.Days{dl}, Cels(cl), Exps(e),...
                        celltypelabel{find(celltype == CurType)}));
                    
%                     pages =  [309 551 608; 524 881 1035]'
%                     pages =  [1092 524 1035; 1220 1306 609]'%
%                     for i = 1:numel(pages)
%                         plot(T(pages(i)), sig(pages(i))+0.05, 'vr');
%                         text(T(pages(i)), sig(pages(i))+0.12, num2str(i), 'HorizontalAlignment', 'center');
%                     end
%                     xlim([T(266) T(1378)]);
                    plot([45 50], 0.5*ones(1, 2), 'k');
                    plot([45 45], [0.5 0.7], 'k');
                    
%                     box off
%                     axis off
                    
%                     yticks([-0.35+0.2*(1:7)/7 0 0.5 1])
%                     yticklabels({'20', '50', '100', '200', '400', '600', '800', '0', '0.5', '1'});
                    %                 end
                    %%
%                     SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
%                     FleNam = sprintf('%sSupFig3_NoSpikeExampleTraceXBC', SaveFolder);
%                     print('-depsc','-painters','-loose', '-r300',FleNam)
%                     saveas(gcf,[FleNam '.png']);
                end
                %%
                
                Colors = [0 0 0; parula(256)];
                if IsDisplayAlign
%                                     pages = [309 551 608; 524 881 1035]'; % XBC 053018-1-3
                                    pages =  [1092 524 1035; 1220 1306 609]' %XBC 053018-1-3
                    %                 pages = [309 425 609; 342 483 908]'; % XBC 062518-1-20
                    %                 pages = [608 1335 1122; 445 553 1309]'; % XBC 080218-3-2
                    %                 pages = [339 481 638; 452 851 312]'; % BC5o 081518-3-19
%                     pages = [315 428 677; 556 855 368]'; % BC7 041618-1-21
                    % pages = []; % BC7 012119 3 4
                    figure;
                    colormap(Colors)
                    for i = 1:numel(pages)
                        subplot(2, 3, i);
                        t = pages(i)+StrFrm-1;
                        imagesc(mean(I(:, :, t+[-1 0 1]), 3), [0 530]);colorbar
                        title(sprintf('%d', i));
                        axis off
                    box off
                    end
                    %%
%                     SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
%                     FleNam = sprintf('%sSupFig3_NoSpikeExampleSpatialMapXBC', SaveFolder);
%                     print('-depsc','-painters','-loose', '-r300',FleNam)
%                     saveas(gcf,[FleNam '.png']);
                end
               %%
               if IsDisplayAlign
                   keyboard;
               end
                %%
                IsDisplayEpoch = 0;
                
                for d = 1:NumDmt
                    CurDmt = Dmt(d);
                    for c = 1:NumCts
                        % Epoch
                        CurCts = Cts(c);
                        CurInd = find(DmtVecFrm==CurDmt & CtsVecFrm==CurCts);
                        %         EpcFrm = bsxfun(@plus, repmat(CurInd(:), 1, SpnDurFrm), (1:SpnDurFrm)-(BasDurFrm+1))';
                        EpcFrm = repmat(CurInd(:), 1, SpnDurFrm) + ((1:SpnDurFrm)-(BasDurFrm+1));
                        CurROIMat = ROISigTrn(EpcFrm');
                        CurROIMat = reshape(CurROIMat, SpnDurFrm, [])';
                        CurROIMat_f  = mean(CurROIMat(:, 1:BasDurFrm), 2);
                        ErrDtc = CurROIMat_f == 0;
                        if any(ErrDtc)
                            ErrInd = find(ErrDtc);
                            CurROIMat_f(ErrInd, :) = [];
                            CurROIMat(ErrInd, :)   = [];
                        end
                        %                     CurROIMat_Df = CurROIMat; %raw
                        CurROIMat_Df = CurROIMat - CurROIMat_f;
                        %                     CurROIMat_Df = CurROIMat./CurROIMat_f;
                        %                     CurROIMat_Dff = CurROIMat_Df ./ CurROIMat_f;
                        
                        DataROIsDf(:, 1:size(CurROIMat_Df, 1), d, c, ii) = CurROIMat_Df';
                        
                        if d == 3 && c == 1
                            if IsDisplayEpoch
                                figure(2);
                                subplot(1, 3, 1);
                                imagesc(CurROIMat);colorbar
                                subplot(1, 3, 2);
                                imagesc(CurROIMat_Df);colorbar
                                subplot(1, 3, 3);
                                imagesc(CurROIMat_Dff);colorbar
                                keyboard;
                            end
                        end
                    end
                end
                %%
                ii = ii + 1;
                clc
                fprintf('progress...%d \n', ii);
            end
        end
    end
    save('ONBCProj_110623_revise.mat', 'DataROIsDf', 'BasDurFrm', 'StmDurFrm', 'Fz', 'ROItab', 'Typ', 'Dmt',...
        'NumCts', 'NumDmt');
end

%%
roi = 1;
figure;
for k = 1:2
    for i = 1:7
        subplot(2, 7, (k-1)*7+i);
        imagesc(squeeze(DataROIsDf(:, :, i, k, roi))');
    end
end
%% Get Repeat Reliability and normalized signals
nROIs = size(DataROIsDf, 5);
Maxes = nan(nROIs, 1);
DataROIsDfnrm = nan(size(DataROIsDf));
for i = 1:nROIs
    Maxes(i) = max(DataROIsDf(:, :, :, :, i), [], 'all');
    DataROIsDfnrm(:, :, :, :, i) = DataROIsDf(:, :, :, :, i)/Maxes(i);
end
RepReli = nan(7, 2, 200);
for k = 1:200
    for i = 1:2
        for j = 1:7
            RepMat = squeeze(DataROIsDfnrm(:, :, j, i, k))';
            if isnan(RepMat(1, 1))
                continue
            end
            RepReli(j, i, k) = estimateRepeatReliability(RepMat);
        end
    end
end

%%
celltype = [50 51 57 58 6 7 89];
celltypelabel = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC8/9'};
lighttype = [1 2];
lighttypelabel = {'ON', 'OFF'};