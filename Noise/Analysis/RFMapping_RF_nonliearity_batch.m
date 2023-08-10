% Run after @GenerateExpectedGC6fResponseNoise
% Get the filename for the batch loop
clear; close all; clc;
ChopClipL = 60;
UpSamplingFz = 100;
Txcorr = [-95 -5];
radius = 40; %
Addnoise = 0.002;
%% load comparison grids
Topic = 'Ai148_Grm6Cre';
Day = '110318';
Cell = 1;
Exp = 5;
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\DataStimulus\';
FileName = sprintf('%s%s_%s_%d_%d', FolderPath, Topic, Day, Cell, Exp);
AnkGrids = load(FileName);
AnkGrids = AnkGrids.DG_OUT.Grids;
AnkL = size(AnkGrids, 3);
clear topic Day Cell Exp FolderPath

%% parameters
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\';
Topic = 'Ai148_Grm6Cre';
% Topic = 'Ai148_AAV-Grm6Cre';
FileType = 'ROISignalAlign';

PathName = [FolderPath 'ProcessedData/'];
FileNames = dir(fullfile(PathName, '*.mat'));
FileIds = find(contains({FileNames.name},[Topic '_' FileType]));
nFile = length(FileIds);
%%
FailedFiles = [];
rng('shuffle');
fshuffle = randperm(nFile);
% keyboard;
for f = 37%;1:nFile
% for f = fshuffle
    FileName = FileNames(FileIds(f)).name;
    Words = strsplit(FileName, [Topic '_' FileType '_']);
    Words = strsplit(Words{2}, '.');
    Words = strsplit(Words{1}, '_');
    Day = Words{1};
    Cel = floor(str2double(Words{2})/1000);
    Exp = mod(str2double(Words{2}), 1000);
    
    
    load([PathName FileName]);
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
    saveFileName = sprintf('./ProcessedData/RFSpatialNL_04242023/%s_NoiseSpatialNonlinear_radius%d_%s_%d_%d.mat',...
        Topic, radius, Day, Cel, Exp);
%     if exist(saveFileName, 'file')
%         continue
%     end
    switch GenType
        case 1
            load(sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d_2019.mat', radius), 'MovieSize',...
                'RFctr_L_gc6f', 'RFctr_NL_gc6f','RFctr_dNL_gc6f', 'radius', 'pixel2um', 'posi');
        case 2
            load(sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d.mat', radius), 'MovieSize',...
                'RFctr_L_gc6f', 'RFctr_NL_gc6f','RFctr_dNL_gc6f', 'radius', 'pixel2um', 'posi');
        case 3
            loadFileName = sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d_%s_%d_%d.mat',...
                radius, Day, Cel, Exp);
            if ~exist(loadFileName, 'file')
                FailedFiles = [FailedFiles; f 1];
                continue
            else
                load(loadFileName);
            end
        otherwise
            FailedFiles = [FailedFiles; f 2];
            continue
    end
    %% upsample stimulus
    nf = size(RFctr_L_gc6f, 2);
    nT = size(ROISig, 2);
    nROI = size(ROISig, 1);
    post_c_id = find(posi(:, 1) == 0 & posi(:, 2) == 0); % Use center of scanning as the start
    x = [1 StmFrm(:)' nT];
    v = [0.5 RFctr_L_gc6f(post_c_id, :) RFctr_L_gc6f(post_c_id, end)];
    mResp = median(ROISig, 1, 'omitnan');
    xq = linspace(1, nT, nT);
    Sig = interp1(x, v, xq, 'pchip');
    %% Draw figures
    t = (0:nT-1)/UpSamplingFz;
    close all
%     figure('visible', 'off');
        figure;
    subplot(2, 3, 1); hold on
    plot(t, Sig, 'b');
    plot(t, mResp/max(mResp), 'k');
    sgtitle(sprintf('%s %s-%d-%d',  Topic, Day, Cel, Exp));
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
    subplot(2, 3, 2:3); hold on
    plot(t, a/max(a), 'b');
    plot(t, b, 'k');
    r = corr(a', b');
    R = corrchop(a', b');
    title(sprintf('%0.3G %0.3G (%0.3G)', lag(mid)/UpSamplingFz, R, r));
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
    subplot(2, 3, 4); 
    imagesc(r_posi); colorbar
    %% Realign temporal again
    [~, posi_max] = max(Rd);
    v = [0.5 RFctr_L_gc6f(posi_max, :) RFctr_L_gc6f(posi_max, end)];
    Sig = interp1(x, v, xq, 'pchip');
    [r, lag] = xcorr(Sig(StmFrm(1):StmFrm(end)), mResp(StmFrm(1):StmFrm(end)), UpSamplingFz, 'coeff');
    r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
    [var, mid] = max(r);
    a = mResp(StmFrm(1):StmFrm(end));
    b = Sig(StmFrm(1):StmFrm(end));
    a = a(-lag(mid):end);
    b = b(1:length(a));
    b = b+Addnoise*randn(size(b));
    t = (0:length(a)-1)/UpSamplingFz;
    subplot(2, 3, 5:6); hold on
    plot(t, a/max(a), 'b');
    plot(t, b, 'k');
    r = corr(a', b');
    R = corrchop(a', b');
    title(sprintf('%0.3G %0.3G (%0.3G)', lag(mid)/UpSamplingFz, R, r));
    %%
    a = ROISig(:, StmFrm(1):StmFrm(end))';
    a = a(-lag(mid):end, :);
    Quality = corrchop(a, b(:));
    keyboard;
    %% Run @TestSubCellularNoiseDifferenceEncodingSpace
%     TestSubCellularNoiseDifferenceEncodingSpace
%     keyboard;
%     continue
    %% Get the optimal combination of linear and nonlinear
    vL = [0.5 RFctr_L_gc6f(posi_max, :) RFctr_L_gc6f(posi_max, end)];
    vNL = [0.5 RFctr_NL_gc6f(posi_max, :) RFctr_NL_gc6f(posi_max, end)];
    SigL = interp1(x, vL, xq, 'pchip');
    SigNL = interp1(x, vNL, xq, 'pchip');
    bL = SigL(StmFrm(1):StmFrm(end));
    bNL = SigNL(StmFrm(1):StmFrm(end));
    bL = bL(1:length(a));
    bNL = bNL(1:length(a));
    bL = bL+Addnoise*randn(size(b));
    bNL = bNL+Addnoise*randn(size(b));
    %     Optm = nan(nROI, 4);
    %     for j = 1:nROI
    %         ca = a(:, j);
    %         if ~all(isnan(ca))
    %             [Optm(j, 1), Optm(j, 2)] = LinearNLCombine(bL, bNL, ca);
    %             Optm(j, 3) = corrchop(bL(:), ca);
    %             Optm(j, 4) = corrchop(bNL(:), ca);
    %         end
    %     end
    mResp = median(a, 2, 'omitnan');
    Optm_m = nan(1, 5);
    [Optm_m(1), Optm_m(2)] = LinearNLCombine(bL, bNL, mResp);
    Optm_m(3) = corrchop(bL(:), mResp);
    Optm_m(4) = corrchop(bNL(:), mResp);
    Optm_m(5) = posi_max;
    %%
    Rd_posi = nan(nposi, 4, nROI);
    for i =1:nposi
        vL = [0.5 RFctr_L_gc6f(i, :) RFctr_L_gc6f(i, end)];
        vNL = [0.5 RFctr_NL_gc6f(i, :) RFctr_NL_gc6f(i, end)];
        SigL = interp1(x, vL, xq, 'pchip');
        SigNL = interp1(x, vNL, xq, 'pchip');
        bL = SigL(StmFrm(1):StmFrm(end));
        bNL = SigNL(StmFrm(1):StmFrm(end));
        bL = bL(1:length(a));
        bNL = bNL(1:length(a));
        bL = bL+0.01*randn(size(b));
        bNL = bNL+0.01*randn(size(b));
        for j = 1:nROI
            ca = a(:, j);
            if ~all(isnan(ca)) && ~all(ca==0)
                [Rd_posi(i, 1, j), Rd_posi(i, 2, j)] = LinearNLCombine(bL, bNL, ca);
                Rd_posi(i, 3, j) = corrchop(bL(:), ca);
                Rd_posi(i, 4, j) = corrchop(bNL(:), ca);
            end
        end
        clc
        fprintf('progress...%d/%d', i , nposi');
    end
    
    %%
    IsDisplay = 0;
    %     close all
    RF_locs = nan(nROI, 6, 3);
    exts = -44:4:44;
    exts_fine = -44:44;
    for j = 1:3
        for i = 1:nROI
            v = reshape(squeeze(Rd_posi(:, j+1, i)), sqrt(nposi), []);
            if any(isnan(v(:)))
                continue;
            end
            v = [v(:, 1) v v(:, end)];
            v = [v(1, :); v; v(end, :)];
            [x, y] = meshgrid(exts, exts);
            [xx, yy] =meshgrid(exts_fine, exts_fine);
            Vq = interp2(x, y,...
                v, xx, yy, 'spline');
            [RF_locs(i, 4, j), maxid] = max(Vq(:));
            RF_locs(i, 5, j) = xx(maxid);
            RF_locs(i, 6, j) = yy(maxid);
            [RF_locs(i, 1, j), maxid] = max(v(:));
            RF_locs(i, 2, j) = x(maxid);
            RF_locs(i, 3, j) = y(maxid);
            
            if IsDisplay
                close all
                figure;
                subplot(1, 2, 1); hold on
                imagesc(v); colorbar
                plot(find(exts == RF_locs(i, 2, j)), find(exts == RF_locs(i, 3, j)), 'ok');
                
                subplot(1, 2, 2); hold on
                imagesc(Vq); colorbar;
                plot(find(exts_fine == RF_locs(i, 5, j)), find(exts_fine == RF_locs(i, 6, j)), 'ok');
                title(sprintf('ROI: %d', i));
                keyboard;
            end
        end
    end
    
    %%
    TimeLag = lag(mid)/UpSamplingFz;
    save(saveFileName, 'Rd_posi','Optm_m', 'TimeLag','RF_locs', 'radius', 'pixel2um', 'posi');
    
    %%
    keyboard;
end
