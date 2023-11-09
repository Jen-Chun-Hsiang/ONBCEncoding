% Run after @GenerateExpectedGC6fResponseNoise
% Get the filename for the batch loop
clear; close all; clc;
ChopClipL = 60;
UpSamplingFz = 100;
Txcorr = [-95 -5];
radius = 40; %
Addnoise = 0.002;
IsDisplay = 0;
%% load comparison grids
Topic = 'Ai148_Grm6Cre';
Day = '110318';
Cell = 1;
Exp = 5;
NoisePath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\';
FileName = sprintf('%s%s%s_%s_%d_%d', NoisePath, 'DataStimulus\', Topic, Day, Cell, Exp);
AnkGrids = load(FileName);
AnkGrids = AnkGrids.DG_OUT.Grids;
AnkL = size(AnkGrids, 3);
clear topic Day Cell Exp FolderPath
%%
saveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\ProcessedData\RFDiffEncoding_050923\';
%%
position = -40:4:40;
% radius = 40; % in  um
% pixel2um = 2;
% nonlineartype = 1; % std
[X, Y] = meshgrid(position, position);
XY = [X(:), Y(:)];
%%
% broken files 031818-3-1
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\ProcessedData\';
Topic = 'Ai148_AAV-Grm6Cre';
% Topic = 'Ai148_Grm6Cre';
FilContain = 'ROISignalAlign';
FileNames = dir(fullfile(FolderPath, '*.mat'));
FileIds = find(contains({FileNames.name},[Topic '_' FilContain]));
nFile = length(FileIds);

Positions = nan(nFile, 4);
for f = 1:nFile
    clc
    fprintf('progress... %d/%d \n', f, nFile);
    %%
    FileName = FileNames(FileIds(f)).name;
    Words = strsplit(FileName, [Topic '_' FilContain '_']);
    Words = strsplit(Words{2}, '.');
    Words = strsplit(Words{1}, '_');
    Day = Words{1};
    Cel = floor(str2double(Words{2})/1000);
    Exp = mod(str2double(Words{2}), 1000);
    %%
    saveFileName = sprintf('%s_RFDiffEncodingSingle_%s_%d_%d.mat',Topic, Day, Cel, Exp);
%     if exist([saveFolder saveFileName], 'file')
%         continue
%     end
    %%
    load([FolderPath FileName]);
    assert(UpSamplingFz == 100);
    nL = size(Grids, 3);
    if size(Grids, 1) == 40
        GenType = 1;
    elseif size(Grids, 3) == 600
        GenType = 3;
    elseif sum(abs(Grids(:) - AnkGrids(:))) < 1e-6 && nL == AnkL
        GenType = 2; % Ai148_Grm6Cre since 2019
    elseif sum(abs(Grids(:) - AnkGrids(:))) > 1e-6
        GenType = 3;
    else
        error('No such file available');
    end
    
    clear RFctr_L_gc6f posi
    switch GenType
        case 1
            load(sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d_2019.mat', radius), 'MovieSize',...
                'RFctr_L_gc6f', 'radius', 'pixel2um', 'posi'); %, 'RFctr_NL_gc6f','RFctr_dNL_gc6f'
        case 2
            load(sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d.mat', radius), 'MovieSize',...
                'RFctr_L_gc6f', 'radius', 'pixel2um', 'posi');%, 'RFctr_NL_gc6f','RFctr_dNL_gc6f'
        case 3
            loadFileName = sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d_%s_%d_%d.mat',...
                radius, Day, Cel, Exp);
            if ~exist(loadFileName, 'file')
                FailedFiles = [FailedFiles; f 1];
                continue
            else
                load(loadFileName, 'RFctr_L_gc6f', 'radius', 'pixel2um', 'posi');
            end
        otherwise
            FailedFiles = [FailedFiles; f 2];
            continue
    end
    %%
    RmIds = any(isnan(ROISig), 2) | all(ROISig == 0, 2);
    ROISig(RmIds, :) = [];
    nT = size(ROISig, 2);
    nROI = size(ROISig, 1);
    weights = nan(3, 1);
    %%
    FilNam = sprintf('ROIs/%s_ROIs_MorphSeg_%s_%d%03d.mat',Topic, Day, Cel, Exp);
    M = load([NoisePath FilNam]);
    if isfield(M, 'Cent')
        Cent = M.Cent;
    else
        Cent = nan(size(M.wROISig, 1), 2);
        for i = 1:size(M.wROISig, 1)
            [Cent(i, 2), Cent(i, 1)] = CenterMass(M.wSlcROI(:, :, i));
        end
    end
    RmIds = any(isnan(M.wROISig), 2) | all(M.wROISig == 0, 2);
    Cent(RmIds, :) = [];
    assert(size(Cent, 1) == nROI);
    clear M
    %% upsample stimulus
    nf = size(RFctr_L_gc6f, 2);
    
    post_c_id = find(posi(:, 1) == 0 & posi(:, 2) == 0); % Use center of scanning as the start
    x = [1 StmFrm(:)' nT];
    v = [0.5 RFctr_L_gc6f(post_c_id, :) RFctr_L_gc6f(post_c_id, end)];
    mResp = median(ROISig, 1, 'omitnan');
    xq = linspace(1, nT, nT);
    Sig = interp1(x, v, xq, 'pchip');
    %% Draw figures
    if IsDisplay
        if ishandle(1)
            close(1);
        end
        figure(1);
        t = (0:nT-1)/UpSamplingFz;
        subplot(2, 3, 1); hold on
        plot(t, Sig, 'b');
        plot(t, mResp/max(mResp), 'k');
        sgtitle(sprintf('%s %s-%d-%d',  Topic, Day, Cel, Exp));
    end
    %% Temporal align
    [r, lag] = xcorr(Sig(StmFrm(1):StmFrm(end)), mResp(StmFrm(1):StmFrm(end)), UpSamplingFz, 'coeff');
    r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    [var, mid] = max(r);
    a = mResp(StmFrm(1):StmFrm(end));
    b = Sig(StmFrm(1):StmFrm(end));
    a = a(-lag(mid):end);
    b = b(1:length(a));
    b = b+Addnoise*randn(size(b));
    t = (0:length(a)-1)/UpSamplingFz;
    r = corr(a', b');
    weights(1) = corrchop(a', b');
    
    if IsDisplay
        subplot(2, 3, 2:3); hold on
        plot(t, a/max(a), 'b');
        plot(t, b, 'k');
        title(sprintf('%0.3G %0.3G (%0.3G)', lag(mid)/UpSamplingFz, weights(1), r));
    end
    %% Get the best location
    nposi = size(posi, 1);
    Rd = nan(nROI, 1);
    for i =1:nposi
        v = [0.5 RFctr_L_gc6f(i, :) RFctr_L_gc6f(i, end)];
        Sig = interp1(x, v, xq, 'pchip');
        b = Sig(StmFrm(1):StmFrm(end));
        b = b(1:length(a));
        b = b+0.01*randn(size(b));
        b = b+Addnoise*randn(size(b));
        Rd(i) = corrchop(b(:), a(:));
    end
    %%
    posi_ids = posi/4;
    posi_ids = posi_ids-min(posi_ids, [], 1)+1;
    pxs = unique(posi_ids(:, 1));
    pys = unique(posi_ids(:, 2));
    npxs = length(pxs);
    npys = length(pys);
    r_posi = nan(npxs, npys);
    for j = 1:nposi
        r_posi(posi_ids(j, 1), posi_ids(j, 2)) = Rd(j);
    end
    if IsDisplay
        subplot(2, 3, 4);
        imagesc(r_posi); colorbar
    end
    
    %% Realign temporal again
    [maxv, posi_max] = max(Rd);
    Positions(f, :) = [posi_max, XY(posi_max, :), maxv.^2];
    
end

%%
% save('ExtractRFCenterPosition.mat', 'Positions');
%%
load('ExtractRFCenterPosition.mat', 'Positions');
%%
x = Positions(:, 4); % quality
y = sqrt(sum(Positions(:, 2:3).^2, 2));
y = y(x>0.1);
x = x(x>0.1);
y = y - mean(y);
%%
Bound95 = quantile(y, 0.95);
SD = std(y, [], 'omitnan');
SD2 = 1.96*SD/sqrt(length(x));
