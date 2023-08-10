IsMakeVideo = 0;
if IsMakeVideo
    FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
    nClip = 11;
    MovieSize = nan(nClip, 3);
    SharedVar = nan(nClip, 1);
    clear RFctr_L RFctr_NL
    for i = 1:nClip
        
        FileName = sprintf('%s0131182_UVProj_%d.mat', FolderPath, i);
        load(FileName);
        VidObj = VideoWriter(sprintf('%s0131182_UVProj_%d.avi', FolderPath, i), 'Uncompressed AVI'); %set your file name and video compression
        VidObj.FrameRate = 52; %set your frame rate
        open(VidObj);
        nframe = size(RsizeFilm, 3);
        for f = 1:nframe  %T is your "100x150x75" matrix
            writeVideo(VidObj,RsizeFilm(:, :, f));
            clc
            fprintf('progress...%d/%d, %d/%d \n', i, nClip, f, nframe);
        end
        close(VidObj);
    end
end

%% Once the avi movies were generated
clear;close all;clc;
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
nClip = 11;
opticFlow = opticalFlowHS('Smoothness', 3);
MovieSize = nan(nClip, 3);
radius = 100; % in  um
pixel2um = 2;
clear FlowV FlowVW
for i = 1:nClip
    FileName = sprintf('%s0131182_UVProj_%d.mat', FolderPath, i);
    load(FileName);
    h = size(RsizeFilm, 1);
    w = size(RsizeFilm, 2);
    f = size(RsizeFilm, 3);
    MovieSize(i, :) = [h, w, f];
    [X, Y] = meshgrid(1:h, 1:w);
    radpix = radius/pixel2um;
    Mask = sqrt((X-0.5*(1+h)).^2+(Y-0.5*(1+w)).^2) <  radpix;
    Mask = Mask';
    gaufilt = mvnpdf([X(:) Y(:)],0.5*[1+h 1+w],2*pi*radpix*eye(2));
    gaufilt = reshape(gaufilt, w, h)';
    gaufilt = gaufilt/max(gaufilt(:));
    cMask = Mask;
    Mask = Mask.*gaufilt;
    Mask = Mask(:);
    weight = Mask(cMask(:)==1);
    
    vidReader = VideoReader(sprintf('%s0131182_UVProj_%d.avi', FolderPath, i));
    count = 1;
    cFlowV = nan(f, 4);
    cFlowVA = nan(f, 4);
    while hasFrame(vidReader)
        frameRGB = readFrame(vidReader);
        frameGray = im2gray(frameRGB);
        flow = estimateFlow(opticFlow,frameGray);
        vx = flow.Vx(cMask == 1);
        vy = flow.Vy(cMask == 1);
        cFlowV(count, 1:2) = [mean(vx) mean(vy)];
        cFlowV(count, 3:4) = [sqrt(sum(cFlowV(count, 1:2).^2))...
            atan2d(cFlowV(count, 2), cFlowV(count, 1))];
%         cFlowVW(count, 1:2) = [weight'*vx weight'*vy];
        cFlowVA(count, 1:2) = [mean(abs(vx)) mean(abs(vy))];
        cFlowVA(count, 3:4) = [sqrt(sum(cFlowVA(count, 1:2).^2))...
            atan2d(cFlowVA(count, 2), cFlowVA(count, 1))];
        count = count + 1;
        clc
        fprintf('progress...%d/%d, %d/%d \n', i, nClip, count, f);
    end
    cFlowVAd = cFlowVA(:, 3) - cFlowV(:, 3);
%     figure; hold on
%     plot(smoothdata(cFlowV(:, 3), 'gaussian', [0 10]), 'k');
%     plot(smoothdata(cFlowVA(:, 3), 'gaussian', [0 10]), 'b');
%     plot(smoothdata(cFlowVAd, 'gaussian', [0 10]), 'm');
%         keyboard;
    FlowV{i} = cFlowV;
    FlowVA{i} = cFlowVA;
    FlowVAd{i} = cFlowVAd;
    clear vidReader flow
end

%%
% save(sprintf('./ProcessedData/NatMovMotionOpticFlow_radius%d.mat', radius), 'MovieSize',...
%     'FlowV', 'FlowVW','radius', 'pixel2um');

%%
Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r)) /(Tau_f - Tau_r);
SamplingRate = 1000; % times per second
FrameHz = 52;
T = -1:(1/SamplingRate):3;
DecayTau = 0.142;
RiseTau = 0.045;
figure; plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)));
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<1));
clear FlowV_gc6f FlowVA_gc6f FlowVAd_gc6f
normalize = @(x) (x-min(x(:)))/range(x(:));
for i = 1:nClip
    x = 1:MovieSize(i, 3);
    xq = linspace(1, MovieSize(i, 3), round(SamplingRate*(MovieSize(i, 3)/FrameHz)));
    
    for j = 1:3
        switch j
            case 1
                wave = normalize(FlowV{i}(:, 3));
            case 2
                wave = normalize(FlowVA{i}(:, 3));
            case 3
                wave = normalize(FlowVAd{i});
        end
        wave = interp1(x, wave, xq, 'pchip');
        wave = [wave(1)*ones(1, SamplingRate) wave(:)' zeros(1, SamplingRate)];
        Response = conv(GC6f, wave,'full')/SamplingRate;
        nL = length(wave);
        resp = resample(Response(1:nL), FrameHz, SamplingRate);
        %
        t = (0:nL-1)/SamplingRate;
        rt = resample(t, FrameHz, SamplingRate, 0);
        rt(end) = rt(end-1)+median(diff(rt(1:end-1)));
%         figure; hold on
%         plot(t, wave, 'k');
%         plot(t, Response(1:nL), 'm');
%         plot(rt, resp, 'b');
%         keyboard;
        resp = resp(rt>=1);
        resp = resp(1:MovieSize(i, 3));
        switch j
            case 1
                FlowV_gc6f{i} = resp;
            case 2
                FlowVA_gc6f{i} = resp;
            case 3
                FlowVAd_gc6f{i} = resp;
        end
    end
end

%%
save(sprintf('./ProcessedData/NatMovMotionOpticFlow_gc6f_radius%d.mat', radius), 'MovieSize',...
    'FlowV_gc6f', 'FlowVA_gc6f','FlowVAd_gc6f', 'radius', 'pixel2um');