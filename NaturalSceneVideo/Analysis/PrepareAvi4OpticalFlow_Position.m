
clear;close all;clc;
%%
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

%%
position = -40:4:40;
radius = 40; % in  um
pixel2um = 2;
nonlineartype = 1; % std
[X, Y] = meshgrid(position, position);
gridposi = [X(:), Y(:)]/pixel2um;
ngridposi = size(gridposi, 1);

Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r)) /(Tau_f - Tau_r);
SamplingRate = 1000; % times per second
FrameHz = 52;
T = -1:(1/SamplingRate):3;
DecayTau = 0.142;
RiseTau = 0.045;
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<1));
normalize = @(x) (x-min(x(:)))/range(x(:));
%% Once the avi movies were generated
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
nClip = 11;
opticFlow = opticalFlowHS('Smoothness', 3);
MovieSize = nan(nClip, 3);
radius = 100; % in  um
pixel2um = 2;
clear FlowV FlowVW
if exist('SpatialCorrelationOpticalFlow.mat', 'file')
    load('SpatialCorrelationOpticalFlow.mat');
    disp('Loading...');
else
    FlowV_pos = [];
    FlowVA_pos = [];
    FlowVAd_pos = [];
    RFctr_ClipIds = [];
    for i = 9:nClip
        FileName = sprintf('%s0131182_UVProj_%d.mat', FolderPath, i);
        load(FileName);
        h = size(RsizeFilm, 1);
        w = size(RsizeFilm, 2);
        f = size(RsizeFilm, 3);
        MovieSize(i, :) = [h, w, f];
        [X, Y] = meshgrid(1:h, 1:w);
        radpix = radius/pixel2um;
        cFlowV_pos = nan(ngridposi, f);
        cFlowVA_pos = nan(ngridposi, f);
        cFlowVAd_pos = nan(ngridposi, f);
        for k = 1:ngridposi
            Mask = sqrt((X-0.5*(1+h)-gridposi(k, 1)).^2+(Y-0.5*(1+w)-gridposi(k, 2)).^2) <  radpix;
            Mask = Mask';
            gaufilt = mvnpdf([X(:) Y(:)],0.5*[1+h 1+w]+gridposi(k, :),2*pi*radpix*eye(2));
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
                fprintf('progress...%d/%d, %d/%d, %d/%d \n', i, nClip, k, ngridposi, count, f);
            end
            cFlowVAd_pos(k, :) = cFlowVA(:, 3) - cFlowV(:, 3);
            %     figure; hold on
            %     plot(smoothdata(cFlowV(:, 3), 'gaussian', [0 10]), 'k');
            %     plot(smoothdata(cFlowVA(:, 3), 'gaussian', [0 10]), 'b');
            %     plot(smoothdata(cFlowVAd, 'gaussian', [0 10]), 'm');
            %         keyboard;
            cFlowV_pos(k, :) = cFlowV(:, 3);
            cFlowVA_pos(k, :) = cFlowVA(:, 3);
        end
        FlowV_pos = [FlowV_pos cFlowV_pos];
        FlowVA_pos = [FlowVA_pos cFlowVA_pos];
        FlowVAd_pos = [FlowVAd_pos cFlowVAd_pos];
        RFctr_ClipIds = [RFctr_ClipIds i*ones(1, f)];
        clear vidReader flow
        save(sprintf('SpatialCorrelationOpticalFlow_%d.mat', i), 'FlowV_pos', 'FlowVA_pos', 'FlowVAd_pos',...
            'RFctr_ClipIds', 'MovieSize');
    end
   
end
%%
nClip = numel(unique(RFctr_ClipIds));

%%
Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r)) /(Tau_f - Tau_r);
SamplingRate = 1000; % times per second
FrameHz = 52;
T = -1:(1/SamplingRate):3;
DecayTau = 0.142;
RiseTau = 0.045;
figure; plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)));
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<1));
normalize = @(x) (x-min(x(:)))/range(x(:));
FlowV_pos_gc6f = [];
FlowVA_pos_gc6f = [];
FlowVAd_pos_gc6f = [];
RFctr_ClipIds_gc6f = [];
for i = 1:nClip
    clipids = RFctr_ClipIds == i;
    nclipids = sum(clipids);
    x = 1:nclipids;
    xq = linspace(1, nclipids, round(SamplingRate*nclipids/FrameHz));
    
    cFlowV_pos_gc6f =  nan(ngridposi, nclipids);
    cFlowVA_pos_gc6f = nan(ngridposi, nclipids);
    cFlowVAd_pos_gc6f = nan(ngridposi, nclipids);
    for k = 1:ngridposi
        for j = 1:3
            switch j
                case 1
                    wave = normalize(FlowV_pos(k, clipids));
                case 2
                    wave = normalize(FlowVA_pos(k, clipids));
                case 3
                    wave = normalize(FlowVAd_pos(k, clipids));
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
            resp = resp(1:nclipids);
            switch j
                case 1
                    cFlowV_pos_gc6f(k, :) = resp;
                case 2
                    cFlowVA_pos_gc6f(k, :) = resp;
                case 3
                    cFlowVAd_pos_gc6f(k, :) = resp;
            end
        end
    end
    FlowV_pos_gc6f = [FlowV_pos_gc6f cFlowV_pos_gc6f];
    FlowVA_pos_gc6f = [FlowVA_pos_gc6f cFlowVA_pos_gc6f];
    FlowVAd_pos_gc6f = [FlowVAd_pos_gc6f cFlowVAd_pos_gc6f];
    RFctr_ClipIds_gc6f = [RFctr_ClipIds_gc6f i*ones(1, nclipids)];
end
%%
save('SpatialCorrelationOpticalFlow_gc6f.mat', 'FlowV_pos_gc6f', 'FlowVA_pos_gc6f',...
    'FlowVAd_pos_gc6f', 'RFctr_ClipIds_gc6f');
%% Examine the spatial correlation
clear R_L;
for i = 1:nClip
    R_L(:, :, i) = corr(FlowVAd_pos_gc6f(:, RFctr_ClipIds_gc6f ==i)').^2;
end
%%
[X, Y] = meshgrid(position, position);
Coords = [X(:), Y(:)];
Distans = sqrt((Coords(:, 1)'-Coords(:, 1)).^2 + (Coords(:, 2)'-Coords(:, 2)).^2);
%%
figure;
scatter(Distans(:), reshape(R_L(:, :, 2), [], 1), 5, 'k', 'filled');
%%
xv = linspace(0, 100, 100)';
v = reshape(R_L(:, :, 2), [], 1);
distu = unique(Distans(:));
ndistu = length(distu);
RepSam = nan(ndistu, 7);
st = linspace(0, 1, 7);
for i = 1:ndistu
    cv = v(Distans(:) == distu(i));
    for j = 1:7
        RepSam(i, j) = quantile(cv, st(j));
    end
end
x = repmat(distu(:), 1, 7);
gprMdl1 = fitrgp(x(:), RepSam(:));
[ypred1,~,yint1] = predict(gprMdl1,xv);
%%
ypred = nan(nClip, 55);
yerr = nan(nClip, 55);
for k = 1:nClip
    v = reshape(R_L(:, :, k), [], 1);
    RepSam = nan(50, 3);
    st = [2.5 50 97.5]/100;
    for i = 1:55
        cv = v(Distans(:)>= (i-1)*2 & Distans(:)<i*2);
        for j = 1:3
            RepSam(i, j) = quantile(sqrt(cv), st(j));
        end
    end
    
    ypred(k, :) = RepSam(:, 2);
    yerr(k, :) = 0.5*(RepSam(:, 3)-RepSam(:, 1));
end
xv = 1:2:109;
%%
QT95 = 18.2615;
SD2 = 9.9436*2;
QT95int = nan(11, 1);
SDint =  nan(11, 1);
for i = 1:11
    ids = ~isnan(ypred(i, :));
    QT95int(i) = interp1(xv(ids), ypred(i, ids), QT95);
    SDint(i) = interp1(xv(ids), ypred(i, ids), SD2/2);
end
save('./Results/SpatialCorrection_110123/SpatialCorrelation_Motion_Intercept.mat',...
    'QT95', 'SD2', 'QT95int', 'SDint');
%%
figure;
shadePlot(xv, mean(ypred, 1), std(ypred, [], 1)/sqrt(nClip), 'k', 0.5);
plot(QT95*ones(1, 2), [0 1], '--k');
text(QT95+1, 0.5, '95%');
xlabel('Distance (um)');
ylabel('Correlation coefficient');
title('Spatial correlation of coherent motion');
xlim([1 76]);
xticks(0:25:75);
xticklabels({'0', '25', '50', '75'});
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylim([0 1.01]);
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sSupFigure6_FeatureSpatialCorrelation_CoherentMotion', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);