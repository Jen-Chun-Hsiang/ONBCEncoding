close all; clear; clc;
% Get plexus data and its assigned BCTypes
%%
NatMovPath = '//storage1.ris.wustl.edu/kerschensteinerd/Active/Emily/NaturalSceneVideo/Analysis/';
SpotMovPath = '//storage1.ris.wustl.edu/kerschensteinerd/Active/Emily/Circle/Analysis/';
NoisePath = '//storage1.ris.wustl.edu/kerschensteinerd/Active/Emily/Noise/Analysis/';
MainPath = ' \\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\';
BCROIFolder = 'BipolarCellTerminal\Analysis\ROIMatchingTable\';
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis');
%%
IsBaselineSub = 0;
Txcorr = [5 95];
Fileradius_NL = 50;
Fileradius_fNL = 40;
Fileradius_MS = 100;
Fileradius_CS = 300;
Fileradius_HP = 50;
lockoffset = [80];%80 20
IsCellAlignment = 0;
%%
Topic = 'Ai148_Grm6Cre';
LoadDay = '030723';
% from spot, and functional alignment
LoadFileName = [NatMovPath 'Results/BCtype_augdata_aligned_assigned_' LoadDay '.mat'];
load(LoadFileName);

LoadFileName = [SpotMovPath 'Results/CenterSurround_Reliability_ROI_Plexus_100722.mat'];
PL = load(LoadFileName);
PLcolab = PL.Collab2num;
PL = PL.SummaryD;
%% PATCH create id number for each experiments
uniExp = unique(PLTab(:, 1:3), 'rows');
nuniExp = size(uniExp, 1);
for i = 1:nuniExp
    cids = PLTab(:, 1) == uniExp(i, 1) & PLTab(:, 2) == uniExp(i, 2) & PLTab(:, 3) == uniExp(i, 3);
    PLTab(cids, 9) = 1:sum(cids);
end
%% Get the minimum length for align NaturalMovie Response altogether
LoadDay = '041623';
PathName = [NatMovPath '/ProcessedData/'];
ClipLenFileName = sprintf('%s/ClipLen_%s.mat', PathName, LoadDay);
if ~exist(ClipLenFileName)
    Subject = 'ROISignalAlign';
    FileNames = dir(fullfile(PathName, '*.mat'));
    FileIds = find(contains({FileNames.name},'Ai148_Grm6Cre_Natural_ROISignalAlign'));
    nFile = length(FileIds);
    Dcliplen = nan(11, nFile);
    for f = 1:nFile
        FileName = FileNames(FileIds(f)).name;
        Words = strsplit(FileName, ['_' Subject '_']);
        Words = strsplit(Words{2}, '_');
        Daystr = Words{1};
        Day = str2double(Words{1});
        Words = strsplit(Words{2}, '.');
        Cel = floor(str2double(Words{1})/1000);
        Exp = mod(str2double(Words{1}),1000);
        D = load([PathName FileName]);
        D.TemporalRange = [-1 0];
        SelectedTW = [];
        NumClips = length(getClipIdinOrderFromVec(D.ClipVecFrm));
        ClipIds = nan(NumClips, 3);
        for i = 1:NumClips
            cStartT = find(D.ClipVecFrm==i);
            cEndT = cStartT(end);
            cStartT = cStartT(1)+round(D.UpSamplingFz*abs(D.TemporalRange(1)));
            cEndT = cEndT+round(D.UpSamplingFz*abs(D.TemporalRange(1)));
            ClipIds(i, :) = [cStartT cEndT length(SelectedTW)+1];
            SelectedTW = [SelectedTW cStartT:cEndT];
        end
        ClipIds = [ClipIds ClipIds(:, 2)-ClipIds(:, 1)+1];
        Dcliplen(:, f) = ClipIds(:, 4);
        FilName = sprintf('./ProcessedData/NatMovAlignmentLength/%s_ExpPlexus_%d_%d_%d.mat',...
            Topic, Day, Cel, Exp);
        save(FilName, 'ClipIds');
    end
    save(ClipLenFileName, 'Dcliplen');
else
    load(ClipLenFileName);
end


%%
ClipLens = round(median(Dcliplen, 2));
nClip = length(ClipLens);
%% Get Data from individual BC
cTopic = 'Ai148_AAV-Grm6Cre';
BCTypes = [50 51 57 58 6 7 89];
%Folders = {'Noise', 'Circle', 'NaturalSceneVideo', 'BlockFreq'};
SpotPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Circle\Analysis\';
LoadFileName = [SpotPath 'Results/CenterSurround_Reliability_individualBC_122122.mat']; % for deconvolution -0.5 baseline
BC = load(LoadFileName);
BC = BC.SummaryD.CelTable(:, 1:3);
[~, BC(:, 4)] = max(BC(:, 3) == BCTypes, [], 2);
% load spatial-nonlinear trace of BC
NL = load(sprintf('./ProcessedData/NatMovRFLNestimation_gc6f_radius%d.mat', Fileradius_NL),...
    'SharedVar', 'MovieSize', 'RFctr_L_gc6f', 'RFctr_NL_gc6f','RFctr_dNL_gc6f');
fNL = load(sprintf('./ProcessedData/NatMovRFSpatialDiffPosition_gc6f_radius%d.mat', Fileradius_fNL),...
    'RFctr_L_gc6f', 'RFctr_NL_gc6f', 'RFctr_ClipIds','posi');
MS = load(sprintf('./ProcessedData/NatMovMotionOpticFlow_gc6f_radius%d.mat', Fileradius_MS),...
    'FlowV_gc6f', 'FlowVA_gc6f', 'FlowVAd_gc6f');
CS = load(sprintf('./ProcessedData/NatMovCSestimation_gc6f_radius%d.mat', Fileradius_CS),...
    'RFctr_CS_gc6f');
HP = load(sprintf('./ProcessedData/NatMovHPestimation_gc6f_radius%d.mat', Fileradius_HP),...
    'TSctr_HP_gc6f', 'TSctr_LP_gc6f');
PathName = ['\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Analysis\'...
    'ProcessedData'];
FileNames = dir(fullfile(PathName, '*.mat'));
FileIds = find(contains({FileNames.name}, [cTopic '_Natural_ROISignalAlign']));
nFile = length(FileIds);
DYexp = [];
DYtab = nan(nFile, 4);
Dclipid = [];
% DNL = nan(nFile, 4, nClip);
for f = 1:nFile
    clc;
    fprintf('progress...%d/%d \n', f, nFile);
    FilName = FileNames(FileIds(f)).name;
    D = load([PathName '\' FilName]);
    Words = strsplit(FilName, '_Natural_ROISignalAlign_');
    Words = strsplit(Words{2}, '_');
    Daystr = Words{1};
    Day = str2double(Words{1});
    Words = strsplit(Words{2}, '.');
    Cel = floor(str2double(Words{1})/1000);
    Exp = mod(str2double(Words{1}),1000);
    FilName = sprintf('./ProcessedData/NatMovAlignmentLength/Ai148_AAV-Grm6Cre_spotBioModel_%d_%d_%d.mat',...
        Day, Cel, Exp);
    load(FilName, 'ClipIds');
    cexp = [];
    YExp = mean(D.ROISig, 1);
    for c = 1:nClip
        x = 1:ClipIds(c, 4);
        v = YExp(ClipIds(c, 1):(ClipIds(c, 2)));
        xq = linspace(1, ClipIds(c, 4), ClipLens(c)+1);
        Dsig = interp1(x, v, xq, 'pchip');
        if f == 1
            x = 1:NL.MovieSize(c, 3);
            xq = linspace(1, NL.MovieSize(c, 3), ClipLens(c)+1);
            % refine linear and nonlinear
            for j = 1:size(fNL.RFctr_L_gc6f, 1)
                v = fNL.RFctr_L_gc6f(j, fNL.RFctr_ClipIds == c);
                vq = interp1(x, v, xq, 'pchip');
                if j == 1 
                    fNL.RF_L{c} = nan(size(fNL.RFctr_L_gc6f, 1), ClipLens(c)+1);
                end
                fNL.RF_L{c}(j, :) = vq;
                
                v = fNL.RFctr_NL_gc6f(j, fNL.RFctr_ClipIds == c);
                vq = interp1(x, v, xq, 'pchip');
                if j == 1
                    fNL.RF_NL{c} = nan(size(fNL.RFctr_L_gc6f, 1), ClipLens(c)+1);
                end
                fNL.RF_NL{c}(j, :) = vq;
            end
            % linear
            v = NL.RFctr_L_gc6f{c};
            NL.RF_L{c} = interp1(x, v, xq, 'pchip');
            % nonlinear
            v = NL.RFctr_NL_gc6f{c};
            NL.RF_NL{c} = interp1(x, v, xq, 'pchip');
            % difference nonlinear
            v = NL.RFctr_dNL_gc6f{c};
            NL.RF_dNL{c} = interp1(x, v, xq, 'pchip');
            % Motion sensitivity (simple average)
            v = smoothdata(MS.FlowV_gc6f{c}, 'gaussian', [0 12]);
            MS.intMSV{c} = interp1(x, v, xq, 'pchip');
            % Motion sensitivity (weighted average)
            v = smoothdata(MS.FlowVA_gc6f{c}, 'gaussian', [0 12]);
            MS.intMSVA{c} = interp1(x, v, xq, 'pchip');
            %
            v = smoothdata(MS.FlowVAd_gc6f{c}, 'gaussian', [0 12]);
            MS.intMSVAd{c} = interp1(x, v, xq, 'pchip');
            % linear
            v = CS.RFctr_CS_gc6f{c};
            NL.RF_CS{c} = interp1(x, v, xq, 'pchip');
            % Transience - high pass
            v = HP.TSctr_HP_gc6f{c};
            HP.intHP{c} = interp1(x, v, xq, 'pchip');
            % Transience - low pass
            v = HP.TSctr_LP_gc6f{c};
            HP.intLP{c} = interp1(x, v, xq, 'pchip');
        end
        if IsBaselineSub
            Dsig = Dsig - mean(Dsig(65:90));
        end
        if f == nFile
            %             Dclipid = [Dclipid; c*ones(ClipLens(c)-90, 1)];
            Dclipid = [Dclipid; c*ones(ClipLens(c)+1, 1)];
        end
        %         cexp = [cexp Dsig(91:end)];
        cexp = [cexp Dsig(1:end)];
    end
    DYexp = [DYexp; cexp];
    cids = find(BC(:, 1) == Day & BC(:, 2) == Cel);
    DYtab(f, :) = BC(cids(1), :);
end
%% non calculating here (outdated)
% figure;
% for i = 1:11
%     subplot(2, 6, i); hold on
%     scatter(squeeze(DNL(:,3, :)), squeeze(DNL(:,4, :)), 10, 'k', 'filled');
% end
%% Get mean responses
mResp = median(DYexp, 1, 'omitnan');
%%

uniCel = unique(DYtab(:, 1:2), 'rows');
nuniCel = size(uniCel, 1);
MovBCexp = nan(nuniCel, size(DYexp, 2));
MovBCtab = nan(nuniCel, 3);
MovBCrpt = nan(nuniCel, 1);
MovBCrep = nan(nuniCel, nClip);
MovBCspNL = nan(nuniCel, 6, nClip);
MovBCfNLnoi = nan(nuniCel, 4);
MovBCfNLmov = nan(nuniCel, 4,  nClip);
MovBCspms = nan(nuniCel, 6, nClip);
MovBCmrp = nan(nuniCel, nClip);
MovBCcs = nan(nuniCel, 2, nClip);
MovBClag = nan(nuniCel, 2);
MovBCts = nan(nuniCel, 2);
NLRF = [];
NLRFNL = [];
NLRFdNL = [];
NLRFCS = [];
MSV = [];
MSVA = [];
MSVAd = [];
TSHP = [];
TSLP = [];
for i = 1:nClip
    NLRF = [NLRF NL.RF_L{i}(:)'];
    NLRFNL = [NLRFNL NL.RF_NL{i}(:)'];
    NLRFdNL = [NLRFdNL NL.RF_dNL{i}(:)'];
    NLRFCS = [NLRFCS NL.RF_CS{i}(:)'];
    MSV = [MSV MS.intMSV{i}(:)'];
    MSVA = [MSVA MS.intMSVA{i}(:)'];
    MSVAd = [MSVAd MS.intMSVAd{i}(:)'];
    TSHP = [TSHP HP.intHP{i}(:)'];
    TSLP = [TSLP HP.intLP{i}(:)'];
end
% For getting noise RF location 
FileType = 'NoiseSpatialNonlinear';
PathName = [NoisePath 'ProcessedData'];
FileNames = dir(fullfile([PathName '/RFSpatialNL_04242023/'], '*.mat'));


%
for i = 1:nuniCel
    cids = find(DYtab(:, 1) == uniCel(i, 1) & DYtab(:, 2) == uniCel(i, 2));
    a = DYexp(cids, :);
    [MovBCexp(i, :), MovBCrep(i, :)] = AssembleClipbyQuality(a, Dclipid,'outmedian');
    MovBCtab(i, :) = DYtab(cids(1), [1:2 4]);
    MovBCrpt(i) = length(cids);
    % Get the nonlinear weights of nonlinear filter
    Daystr = num2str(DYtab(cids(1), 1));
    if length(Daystr)<6
        Daystr = ['0' Daystr];
    end
    FileIds = find(contains({FileNames.name},[cTopic '_' FileType]) &...
        contains({FileNames.name},[Daystr '_' num2str(DYtab(cids(1), 2))]));
    nFile = length(FileIds);
    if nFile > 0
        cMovBCfNLnoi = nan(nFile, 4);
        clipr = nan(nFile, 4, nClip);
        aclipr = nan(nFile, 4, nClip, 441);
        for k = 1:nFile
            FileName = FileNames(FileIds(k)).name;
            load([PathName '/RFSpatialNL_04242023/' FileName], 'RF_locs', 'Optm_m');
            
            cMovBCfNLnoi(k, 1:2) = Optm_m(1:2);
            % Get weighted position
            w = RF_locs(:, 4, 1);
            rmids = ~isnan(w) & w >0.2;
            if sum(rmids) > 1
                w = w(rmids);
                w = w.^2;
                w = w/sum(w);
            else
                rmids = ~isnan(w);
                w = w(rmids);
                w = (w-min(w))/range(w);
                w = w/sum(w);
            end
            mposi = [w'*RF_locs(rmids, 5, 1) w'*RF_locs(rmids, 6, 1)];
            mposi = min([40*ones(size(mposi)); mposi], [], 1);
            mposi = max([-40*ones(size(mposi)); mposi], [], 1);
            mposiId = find(fNL.posi(:, 1) == round(mposi(1)/4)*4 &...
            fNL.posi(:, 2) == round(mposi(2)/4)*4);
            cMovBCfNLnoi(k, 3:4) = mposi;
            % calculate nonlinearity in natural movies
            if ~IsCellAlignment
                if ~isempty(lockoffset)
                    for j = 1:nClip
                        b = MovBCexp(i, Dclipid'==j);
                        b = b(1:length(b)-lockoffset);
                        bL = fNL.RF_L{j}(mposiId, :);
                        bNL = fNL.RF_NL{j}(mposiId, :);
                        bL = bL((lockoffset+1):end);
                        bNL = bNL((lockoffset+1):end);
                        [clipr(k, 1, j), clipr(k, 2, j)] = LinearNLCombine(bL, bNL, b);
                        clipr(k, 3, j) = corr(bL(:), b(:));
                        clipr(k, 4, j) = corr(bNL(:), b(:));
%                         for w = 1:441
%                             bL = fNL.RF_L{j}(w, :);
%                             bNL = fNL.RF_NL{j}(w, :);
%                             bL = bL((lockoffset+1):end);
%                             bNL = bNL((lockoffset+1):end);
%                             [aclipr(k, 1, j, w), aclipr(k, 2, j, w)] = LinearNLCombine(bL, bNL, b);
%                             aclipr(k, 3, j, w) = corr(bL(:), b(:));
%                             aclipr(k, 4, j, w) = corr(bNL(:), b(:));
%                         end
                    end
                end
            end
            % 
            
        end
        % decide which recording wins
        if nFile > 1
            [~, maxid] = max(mean(clipr(:, 2, :).^2, 3));
            MovBCfNLmov(i, :, :) = clipr(maxid, :, :);
            [~, maxid] = max(cMovBCfNLnoi(:, 2));
            MovBCfNLnoi(i, :) = cMovBCfNLnoi(maxid, :);
        else
            MovBCfNLnoi(i, :) = squeeze(cMovBCfNLnoi(1, :));
            MovBCfNLmov(i, :, :) = squeeze(clipr(1, :, :));
        end
    end
    % correlate with visual features
    if IsCellAlignment
        a = MovBCexp(i, :);
        %         for j = 1:nClip
        %             b = a(Dclipid'==j);
        %             a(Dclipid'==j) = b-mean(b);
        %         end
        [r, lag] = xcorr(NLRF, a, D.UpSamplingFz, 'coeff');
        r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
        [MovBClag(i, 1), mid] = max(r);
        MovBClag(i, 2) = 1-lag(mid)/D.UpSamplingFz;
        MovBClag(i, 1) = MovBClag(i, 1)^2;
        lockoffset = mid-100;
        for j = 1:nClip
            b = a(Dclipid'==j);
            b = b(1:length(b)-lockoffset);
            for k = 1:10
                switch k
                    case 1
                        c = mResp(:, Dclipid'==j);
                    case 2
                        c = NLRF(:, Dclipid'==j);
                    case 3
                        c = NLRFNL(:, Dclipid'==j);
                    case 4
                        c = NLRFCS(:, Dclipid'==j);
                    case 5
                        c = MSV(:, Dclipid'==j);
                    case 6
                        c = MSVA(:, Dclipid'==j);
                    case 7
                        c = TSHP(:, Dclipid'==j);
                    case 8
                        c = TSLP(:, Dclipid'==j);
                    case 9
                        c = NLRFdNL(:, Dclipid'==j);
                    case 10
                        c = MSVAd(:, Dclipid'==j);
                end
                c = c((lockoffset+1):end);
                r = corr(c(:), b(:));
                switch k
                    case 1
                        MovBCmrp(i, j) = r;
                    case 2
                        MovBCspNL(i,1, j) = r;
                        MovBCspNL(i,2, j) = MovBClag(i, 2);
                    case 3
                        MovBCspNL(i,3, j) = r;
                        MovBCspNL(i,4, j) = MovBClag(i, 2);
                    case 4
                        MovBCcs(i,1, j) = r;
                        MovBCcs(i,2, j) = MovBClag(i, 2);
                    case 5
                        MovBCspms(i,1, j) = r;
                        MovBCspms(i,2, j) = MovBClag(i, 2);
                    case 6
                        MovBCspms(i,3, j) = r;
                        MovBCspms(i,4, j) = MovBClag(i, 2);
                    case 7
                        MovBCts(i,1, j) = r;
                        MovBCts(i,3, j) = MovBClag(i, 2);
                    case 8
                        MovBCts(i,2, j) = r;
                        MovBCts(i,4, j) = MovBClag(i, 2);
                    case 9
                        MovBCspNL(i,5, j) = r;
                        MovBCspNL(i,6, j) = MovBClag(i, 2);
                    case 10
                        MovBCspms(i,5, j) = r;
                        MovBCspms(i,6, j) = MovBClag(i, 2);
                end
            end
        end
    end
    for j = 1:nClip
        if ~isempty(lockoffset)
            b = MovBCexp(i, Dclipid'==j);
            b = b(1:length(b)-lockoffset);
            for k = 1:10
                switch k
                    case 1
                        a = mResp(Dclipid'==j);
                    case 2
                        a = NL.RF_L{j};
                    case 3
                        a = NL.RF_NL{j};
                    case 4
                        a = NL.RF_CS{j};
                    case 5
                        a = MS.intMSV{j};
                    case 6
                        a = MS.intMSVA{j};
                    case 7
                        a = HP.intHP{j};
                    case 8
                        a = HP.intLP{j};
                    case 9
                        a = NL.RF_dNL{j};
                    case 10
                        a = MS.intMSVAd{j};
                end
                a = a((lockoffset+1):end);
                r = corr(a(:), b(:));
                switch k
                    case 1
                        MovBCmrp(i, j) = r;
                    case 2
                        MovBCspNL(i,1, j) = r;
                        MovBCspNL(i,2, j) = 1-lockoffset/D.UpSamplingFz;
                    case 3
                        MovBCspNL(i,3, j) = r;
                        MovBCspNL(i,4, j) = 1-lockoffset/D.UpSamplingFz;
                    case 4
                        MovBCcs(i,1, j) = r;
                        MovBCcs(i,2, j) = 1-lockoffset/D.UpSamplingFz;
                    case 5
                        MovBCspms(i,1, j) = r;
                        MovBCspms(i,2, j) = 1-lockoffset/D.UpSamplingFz;
                    case 6
                        MovBCspms(i,3, j) = r;
                        MovBCspms(i,4, j) = 1-lockoffset/D.UpSamplingFz;
                    case 7
                        MovBCts(i,1, j) = r;
                        MovBCts(i,3, j) = 1-lockoffset/D.UpSamplingFz;
                    case 8
                        MovBCts(i,2, j) = r;
                        MovBCts(i,4, j) = 1-lockoffset/D.UpSamplingFz;
                    case 9
                        MovBCspNL(i,5, j) = r;
                        MovBCspNL(i,6, j) = 1-lockoffset/D.UpSamplingFz;
                    case 10
                        MovBCspms(i,5, j) = r;
                        MovBCspms(i,6, j) = 1-lockoffset/D.UpSamplingFz;
                end
            end
        else
            % Control average response
            [r, lag] = xcorr(mResp(Dclipid'==j), MovBCexp(i, Dclipid'==j), D.UpSamplingFz, 'coeff');
            [MovBCmrp(i, j), mid] = max(r);
            % Linear intensity
            [r, lag] = xcorr(NL.RF_L{j}, MovBCexp(i, Dclipid'==j), D.UpSamplingFz, 'coeff');
            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
            [MovBCspNL(i,1, j), mid] = max(r);
            MovBCspNL(i,3, j) = 1-lag(mid)/D.UpSamplingFz;
            % Nonlinear intensity
            [r, lag] = xcorr(NL.RF_NL{j}, MovBCexp(i, Dclipid'==j), D.UpSamplingFz, 'coeff');
            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
            [MovBCspNL(i,2, j), mid] = max(r);
            MovBCspNL(i,4, j) = 1-lag(mid)/D.UpSamplingFz;
            % Motion sensitivity
            [r, lag] = xcorr(MS.intMSV{j}, MovBCexp(i, Dclipid'==j), D.UpSamplingFz, 'coeff');
            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
            [MovBCspms(i,1, j), mid] = max(r);
            MovBCspms(i,3, j) = 1-lag(mid)/D.UpSamplingFz;
            % Motion sensitivity weighted
            [r, lag] = xcorr(MS.intMSVW{j}, MovBCexp(i, Dclipid'==j), D.UpSamplingFz, 'coeff');
            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
            [MovBCspms(i,2, j), mid] = max(r);
            MovBCspms(i,4, j) = 1-lag(mid)/D.UpSamplingFz;
            % Center surround
            [r, lag] = xcorr(NL.RF_CS{j}, MovBCexp(i, Dclipid'==j), D.UpSamplingFz, 'coeff');
            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
            [MovBCcs(i,1, j), mid] = max(abs(r));
            MovBCcs(i,2, j) = 1-lag(mid)/D.UpSamplingFz;
        end
    end
end

%%
% figure;
% for i = 1:nClip
%     subplot(2, 6, i); hold on
%     scatter(squeeze(MovBCspNL(:,1, i)).^2, squeeze(MovBCspNL(:,2, i)).^2, 10, 'k', 'filled');
%     plot([0 1], [0 1], '--k');
% end
% figure;
% for i = 1:nClip
%     subplot(2, 6, i); hold on
%     scatter(squeeze(MovBCspms(:,2, i)).^2, squeeze(MovBCspNL(:,1, i)).^2, 10, 'k', 'filled');
%     plot([0 1], [0 1], '--k');
% end
% figure;
% for i = 1:nClip
%     subplot(2, 6, i); hold on
%     scatter(MS.intMSVA{i}, NL.RF_NL{i}, 10, 'k', 'filled');
%     plot([0 1], [0 1], '--k');
% end
% MovBCnl = (squeeze(MovBCspNL(:,2, :)).^2-squeeze(MovBCspNL(:,1, :)).^2)./...
%     (squeeze(MovBCspNL(:,2, :)).^2+squeeze(MovBCspNL(:,1, :)).^2);
MovBCbs = squeeze(MovBCspNL(:,1, :)).^2;
MovBCnl = squeeze(MovBCspNL(:,3, :)).^2;
MovBCdnl = squeeze(MovBCspNL(:,5, :)).^2;
MovBCms = squeeze(MovBCspms(:,1, :)).^2;
MovBCmsa = squeeze(MovBCspms(:,3, :)).^2;
MovBCmsad = squeeze(MovBCspms(:,5, :)).^2;
MovBCmrp = MovBCmrp.^2;
MovBCcs = squeeze(MovBCcs(:, 1, :).^2);
MovBCtshp = squeeze(MovBCts(:, 1, :).^2);
MovBCtslp = squeeze(MovBCts(:, 2, :).^2);
% MovBCnl = squeeze(MovBCspNL(:,1, :)).^2;
% MovBCms = MovBCrep.^2;
% MovBCnl = squeeze(MovBCspNL(:,1, :)).^2 ./(squeeze(MovBCspNL(:,1, :)).^2+squeeze(MovBCspNL(:,2, :)).^2-NL.SharedVar'.^2);
%% check number of unique types
% c = tabulate(PLTab(:, 8));
%% Get Plexus data based on plane alignment
typesids = 1:7;
uniExp = unique(PLTab(ismember(PLTab(:, 8), typesids), [1:2 5]), 'rows');
nuniExp = size(uniExp, 1);
MovResp = [];
MovRpt = [];
MovTab = [];
MovQual = [];
MovspNL = [];
MovspMS = [];
MovspDV = [];
MovspCS = [];
SpotTab = [];
NoiTab = [];
MovfNL = [];
% {'Noise', 'Circle', 'NaturalSceneVideo', 'BlockFreq'};
%     1         2               3               4
clear NoAccessFiles
NoAccessFiles{1} = 'no access';
for i = 1:nuniExp
    clc
    fprintf('progress... %d/%d \n', i, nuniExp);
    cids = find(PLTab(:, 1) == uniExp(i, 1) & PLTab(:, 2) == uniExp(i, 2) & PLTab(:, 5) == uniExp(i, 3) &...
        ismember(PLTab(:, 8), typesids));
    Day = num2str(uniExp(i, 1));
    if length(Day) ~= 6
        Day = ['0' Day];
    end
    sids = PLTab(cids, 9);
    % get the plane data as well as related recordings (Natural Movie)
    FileName = sprintf('%s_%s_%d_%d_%d.mat', Topic, 'ROIMatchingTable',...
        uniExp(i, 1),  uniExp(i, 2), uniExp(i, 3));
    try
        load([MainPath BCROIFolder FileName], 'ROIMatchingTab','anckor');
    catch
        NoAccessFiles{end+1} = FileName;
        continue
    end
    % find circle to align, if no, skip
    if anckor(3) == 2 || any(ROIMatchingTab(:, 4)==2)
        % find Natural Movie to align, if no, skip
        if anckor(3) == 3 || any(ROIMatchingTab(:, 4)==3)
            spotid = [];
            natmovid = [];
            if anckor(3) == 2
                spotid = [spotid; anckor(2)];
            end
            if anckor(3) == 3
                natmovid = [natmovid; anckor(2)];
            end
            if any(ROIMatchingTab(:, 4)==2)
                spotid = [spotid; unique(ROIMatchingTab(ROIMatchingTab(:, 4)==2, 3))];
            end
            if any(ROIMatchingTab(:, 4)==3)
                natmovid = [natmovid; unique(ROIMatchingTab(ROIMatchingTab(:, 4)==3, 3))];
            end
        else
            continue
        end
    else
        continue
    end
    if anckor(3) == 1
        noivids = anckor(2);
    else
        noivids = [];
    end
    if any(ROIMatchingTab(:, 4)==1)
        noivids = [noivids; unique(ROIMatchingTab(ROIMatchingTab(:, 4)==1, 3))];
    end
    
    % use circle to locate corresponding ROI in natural movie
    %     [Nids, Sids, Fids] = NatMovROIMap(ROIMatchingTab, natmovid, spotid);
    addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\');
    [Nids, Sids, Oids, Fids] = MovSpotNoiROIMap(ROIMatchingTab, natmovid, spotid, noivids);
    
    % Get Spot info
    cSpotTab = nan(size(Sids, 1), 2, size(Sids, 2));
    for j = 1:length(spotid)
        cids = find(PLTab(:, 1) == uniExp(i, 1) & PLTab(:, 2) == uniExp(i, 2) & PLTab(:, 3) == spotid(j));
        if ~isempty(cids)
            cSpotTab(~isnan(Sids(:, j)), :, j) = PLTab(cids(Sids(~isnan(Sids(:, j)), j)), [7 8]);
        end
    end
    spotL = 0;
    if ndims(cSpotTab) == 2
        SpotTab = [SpotTab; cSpotTab];
        spotL = spotL + 1;
    else
        for j = 1:length(spotid)
            if ~all(isnan(cSpotTab(:, 1, j)), 1)
                SpotTab = [SpotTab; cSpotTab(:, :, j)];
                spotL = spotL + 1;
            end
        end
    end
    if spotL == 0
        keyboard;
    end
    % Get Noise info
    cNoiTab = nan(size(Nids, 1), 6, length(noivids));
    if ~isempty(noivids)
        clipr = nan(length(noivids), 4, nClip);
        for j = 1:length(noivids)
            % (1) get the nonlinear weights and RF locations
            FileName = sprintf('%sProcessedData/RFSpatialNL_04242023/%s_NoiseSpatialNonlinear_radius%d_%s_%d_%d.mat',...
                NoisePath, Topic, 40, Day, uniExp(i, 2), noivids(j));
            if exist(FileName, 'file')
                D = load(FileName, 'RF_locs', 'Rd_posi');
                rmvids = find(all(isnan(D.RF_locs(:, :, 2)), 2));
                D.RF_locs(rmvids, :, :) = [];
                D.Rd_posi(:, :, rmvids) = [];
                nkpids = size(D.Rd_posi, 3);
                MaxV = nan(nkpids, 4);
                for c = 1:nkpids
                    [~, mid] = max(D.Rd_posi(:, 2, c));
                    MaxV(c, :) = D.Rd_posi(mid, :, c);
                end
                % Optimal weight, max L+NL r, max L r, max NL r, fine x loc,
                % fine y loc
                cNoiTab(~isnan(Oids(:, j)), :, j) = [MaxV(Oids(~isnan(Oids(:, j)), j), :)...
                    D.RF_locs(Oids(~isnan(Oids(:, j)), j), 5:6, 1)];
                % (2) for convolution, requires to use RF locations in the same plane
                % as movies
            end
        end
        cNoiTab = winnertakesall(cNoiTab, 2);
    else
        cNoiTab = nan(size(Nids, 1), 6, 1);
    end
    NoiTab = [NoiTab; repmat(cNoiTab, spotL, 1)];
    
    % Get Natural Movie Response and Info
    cTab = nan(size(Nids));
    clipids = [];
    expresp = nan(size(Nids, 1), sum(ClipLens)+nClip, length(natmovid));
    clipr = nan(length(natmovid), size(Nids, 1), 4, nClip);
    for j = 1:length(natmovid)
        % get the ROI type to the natural movie ROIs
        FileName = sprintf('%sProcessedData/%s_Natural_ROISignalAlign_%s_%d%03d.mat', NatMovPath, Topic,...
            Day, uniExp(i, 2), natmovid(j));
        D = load(FileName);
        FilName = sprintf('%sProcessedData/NatMovAlignmentLength/%s_ExpPlexus_%d_%d_%d.mat',...
            NatMovPath, Topic, uniExp(i, 1), uniExp(i, 2), natmovid(j));
        load(FilName, 'ClipIds');
        nROI = size(D.ROISig, 1);
        FilName = sprintf('./ROIs/%s_ROIs_MorphSeg_%s_%d%03d.mat',Topic, Day, uniExp(i, 2), natmovid(j));
        load(FilName, 'wROISig');
        RmIds = any(isnan(wROISig), 2);
        wROISig(RmIds, :) = [];
        assert(nROI == size(wROISig, 1));
        mInt =  std(wROISig(:, :), [], 2);
        cTab(~isnan(Nids(:, j)), j) = mInt(Nids(~isnan(Nids(:, j)), j));
        cexp = [];
        for c = 1:nClip
            Dsig = nan(nROI, ClipLens(c)+1);
            for h = 1:nROI
                x = 1:ClipIds(c, 4);
                v = D.ROISig(h, ClipIds(c, 1):(ClipIds(c, 2)));
                xq = linspace(1, ClipIds(c, 4), ClipLens(c)+1);
                Dsig(h, :) = interp1(x, v, xq, 'pchip');
            end
            seq = Dsig(:, 1:end);
            if j == 1
                clipids = [clipids c*ones(1, size(seq, 2))];
            end
            cexp = [cexp seq];
        end
        expresp(~isnan(Nids(:, j)), :, j) = cexp(Nids(~isnan(Nids(:, j)),j), :);
        % Extend to nonlinear 
        cNoiTab(:, 5:6) = round(cNoiTab(:, 5:6)/4)*4;
        cNoiTab(:, 5:6) = min(cat(3, 40*ones(size(cNoiTab(:, 5:6))), cNoiTab(:, 5:6)), [], 3);
        cNoiTab(:, 5:6) = max(cat(3, -40*ones(size(cNoiTab(:, 5:6))), cNoiTab(:, 5:6)), [], 3);
        cNoiTab(isnan(cNoiTab(:, 1)), 5:6) = nan;
        for c = 1:nClip
            if ~isempty(lockoffset)
                for h = 1:size(Fids, 1)
                    if ~isnan(cNoiTab(h, 5)) && ~isnan(Nids(h, j))                        
                        b = expresp(h, clipids==c);
                        b = b(1:length(b)-lockoffset);
                        if all(~isnan(b))
                            mposiId = find(fNL.posi(:, 1) == round(cNoiTab(h, 5)/4)*4 &...
                                fNL.posi(:, 2) == round(cNoiTab(h, 6)/4)*4);
                            bL = fNL.RF_L{c}(mposiId, :);
                            bNL = fNL.RF_NL{c}(mposiId, :);
                            bL = bL((lockoffset+1):end);
                            bNL = bNL((lockoffset+1):end);
                            [clipr(j, h, 1, c), clipr(j, h, 2, c)] = LinearNLCombine(bL, bNL, b);
                            clipr(j, h, 3, c) = corr(bL(:), b(:));
                            clipr(j, h, 4, c) = corr(bNL(:), b(:));
                        end
                    end
                end
            end
        end
    end
    if length(natmovid) > 1
        cclipr = nan(size(Nids, 1), 4, nClip);
        for h = 1:size(Fids, 1)
            if ~isnan(cNoiTab(h, 5)) && ~isnan(Nids(h, j))
                [~, maxid] = max(mean(squeeze(clipr(:, h, 2, :).^2), 2));
                cclipr(h, :, :) = squeeze(clipr(maxid, h, :, :));
            end
        end
        clipr = cclipr;
        clear cclipr
    else
        clipr = squeeze(clipr);
    end
    MovfNL = cat(1, MovfNL, repmat(clipr, spotL, 1, 1));
    cResp = nan(size(Fids, 1), size(expresp, 2));
    cQual = nan(size(Fids, 1), nClip);
    crpt = nan(size(Fids, 1), 1);
    for j = 1:size(Fids, 1)
        cids = ~isnan(Nids(j, :));
        crpt(j) = sum(cids);
        if sum(cids) > 1
            [cResp(j, :),  cQual(j, :)] = AssembleClipbyQuality(squeeze(expresp(j, :, cids))', clipids, 'outmedian');
            
        else
            cResp(j, :) = squeeze(expresp(j, :, cids))';
        end
    end
    MovResp = [MovResp; repmat(cResp, spotL, 1)];
    MovQual = [MovQual; repmat(cQual, spotL, 1)];
    MovRpt = [MovRpt; repmat(crpt, spotL, 1)];
    cTab = [uniExp(i, 1:3).*ones(size(cTab, 1), 1) mean(cTab, 2, 'omitnan')];
    MovTab = [MovTab; repmat(cTab, spotL, 1)];
    %     cMovspNL = nan(size(Fids, 1), 4, nClip);
    %     cMovspMS = nan(size(Fids, 1), 4, nClip);
    %     cMovmrp = nan(size(Fids, 1), nClip);
    %     cMovcs = nan(size(Fids, 1), nClip);
    %     for k = 1:size(Fids, 1)
    %         for j = 1:nClip
    %             % Control average response
    %             [r, lag] = xcorr(mResp(Dclipid'==j), cResp(k, Dclipid'==j), D.UpSamplingFz, 'coeff');
    %             [cMovmrp(k, j), mid] = max(r);
    %             % loca contrast
    %             [r, lag] = xcorr(NL.RF_L{j}, cResp(k, clipids'==j), D.UpSamplingFz, 'coeff');
    %             r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    %             [cMovspNL(k, 1, j), mid] = max(r);
    %             cMovspNL(k, 3, j) = 1-lag(mid)/D.UpSamplingFz;
    %             [r, lag] = xcorr(NL.RF_NL{j}, cResp(k, clipids'==j), D.UpSamplingFz, 'coeff');
    %             r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    %             [cMovspNL(k, 2, j), mid] = max(r);
    %             cMovspNL(k, 4, j) = 1-lag(mid)/D.UpSamplingFz;
    %             % Motion sensitivity
    %             [r, lag] = xcorr(MS.intMSV{j}, cResp(k, Dclipid'==j), D.UpSamplingFz, 'coeff');
    %             r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    %             [cMovspMS(k,1, j), mid] = max(r);
    %             cMovspMS(k,3, j) = 1-lag(mid)/D.UpSamplingFz;
    %             [r, lag] = xcorr(MS.intMSVW{j}, cResp(k, Dclipid'==j), D.UpSamplingFz, 'coeff');
    %             r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    %             [cMovspMS(k,2, j), mid] = max(r);
    %             cMovspMS(k,4, j) = 1-lag(mid)/D.UpSamplingFz;
    %             % center surround
    %             [r, lag] = xcorr(NL.RF_CS{j}, cResp(k, clipids'==j), D.UpSamplingFz, 'coeff');
    %             r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    %             [cMovcs(k, j), mid] = max(abs(r));
    %         end
    %     end
    %     MovspNL = cat(1, MovspNL, repmat(cMovspNL, spotL, 1, 1));
    %     MovspMS = cat(1, MovspMS, repmat(cMovspMS, spotL, 1, 1));
    %     MovspDV = cat(1, MovspDV, repmat(cMovmrp, spotL, 1));
    %     MovspCS = cat(1, MovspCS, repmat(cMovcs, spotL, 1));
    % use repeat natural movie to measure its quality
    % get the all repeat match ROIs
end

% Movnl = squeeze(MovspNL(:,2, :)).^2;
% Movms = squeeze(MovspMS(:,1, :)).^2;
% Movdv = squeeze(MovspDV(:,:)).^2;
% Movcs = squeeze(MovspCS(:,:)).^2;

clear r lag cMovspNL cMovspMS cMovmrp cMovcs cTab cResp cQual Nids Sids Fids ROIMatchingTab
% clear NS CS MS
%%
loadDay = '40423';
loadFileName = ['\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Circle\Analysis\Results\IdentifiedBCtype_Response' loadDay '.mat'];
load(loadFileName, 'SliceResp_info', 'SliceResp_full', 'SliceResp_type');
ntab = size(MovBCtab, 1);
Mov2Spot = nan(ntab, 1);
for i = 1:ntab
    cid = find(SliceResp_info(:, 1) == MovBCtab(i, 1) & SliceResp_info(:, 2) == MovBCtab(i, 2));
    if ~isempty(cid)
        assert(SliceResp_info(cid, 3) == BCTypes(MovBCtab(i, 3)));
        Mov2Spot(i) = cid;
    end
end
%%
clear wROISig seq expresp cexp a YExp SelectedTW SliceResp_full FileNames Dfull Dsig DYexp D
%%
Features = [];
FeatureNames = {'RF_L', 'RF_NL', 'RF_dNL', 'RF_CS', 'MS_V', 'MS_VA', 'MS_VAd', 'TS_HP', 'TS_LP', 'MResp'};
for i = 1:nClip
    cfea = [NL.RF_L{i}(:)'; NL.RF_NL{i}(:)'; NL.RF_dNL{i}(:)'; NL.RF_CS{i}(:)'; MS.intMSV{i}(:)';...
        MS.intMSVA{i}(:)'; MS.intMSVAd{i}(:)'; HP.intHP{i}(:)'; HP.intLP{i}(:)'];
    Features = [Features cfea];
end
Features = [Features; mResp];
%%
% clear cfea mResp
%%
% SaveDay = '40423';
% SaveFileName = ['./ProcessedData/BCNatMovResponse_AssignedROI' SaveDay '.mat'];
% save(SaveFileName, 'SpotTab', 'Mov2Spot', 'Features', 'FeatureNames', 'clipids',...
%     'MovBCexp','MovBCrep','MovBCtab', 'MovBCnl','MovBCcs', 'MovBCms','MovBCtshp','MovBCtslp', ...
%     'MovResp', 'MovQual', 'MovTab', 'MovBCspNL','MovBCbs','MovBCmrp','MovBCts',...
%     'lockoffset','Txcorr', 'IsCellAlignment', 'MovBClag');%, 'Movdv', 'Movnl', 'Movms', 'Movcs'

% @CheckAvaliableData
keyboard;
