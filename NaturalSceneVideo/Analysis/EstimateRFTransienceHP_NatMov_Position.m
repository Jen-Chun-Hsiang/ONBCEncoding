clear; close all; clc;
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
%% Loop through all the clip and generate the traces and evaluation (correlation)
radius_c = 50; % in  um
PassHz = 4;
pixel2um = 2;
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
nClip = 11;
MovieSize = nan(nClip, 3);
SharedVar = nan(nClip, 1);

if exist('SpatialCorrelationTransience.mat', 'file')
    disp('Loading...');
    load('SpatialCorrelationTransience.mat');
    
else
    RFctr_HPf = [];
    RFctr_ClipIds = [];
    for i = 1:nClip
        FileName = sprintf('%s0131182_UVProj_%d.mat', FolderPath, i);
        load(FileName);
        h = size(RsizeFilm, 1);
        w = size(RsizeFilm, 2);
        f = size(RsizeFilm, 3);
        MovieSize(i, :) = [h, w, f];
        [X, Y] = meshgrid(1:h, 1:w);
        % center
        radpix = radius_c/pixel2um;
        cRFctr_HP = nan(ngridposi, f);
        for k = 1:ngridposi
            clear RFctr_L RFctr_NL RFctr_dNL
            Mask = sqrt((X-0.5*(1+h)-gridposi(k, 1)).^2+(Y-0.5*(1+w)-gridposi(k, 2)).^2) <  radpix;
            Mask = Mask';
            gaufilt = mvnpdf([X(:) Y(:)],0.5*[1+h 1+w]+gridposi(k, :),2*pi*radpix*eye(2));
            gaufilt = reshape(gaufilt, w, h)';
            gaufilt = gaufilt/max(gaufilt(:));
            cMask = Mask;
            Mask = Mask.*gaufilt;
            movie = double(reshape(RsizeFilm, [], f))/255-0.5;
            contrat_C = squeeze(movie'*Mask(:)/sum(Mask(:)));
            cRFctr_HP(k, :) = contrat_C;
            clc
            fprintf('progress...%d/%d (%d/%d) \n', i, nClip, k, ngridposi);
        end
        RFctr_HPf = [RFctr_HPf cRFctr_HP];
        RFctr_ClipIds = [RFctr_ClipIds i*ones(1, f)];
    end
    %%
    save('SpatialCorrelationTransience.mat', 'RFctr_HPf', 'RFctr_ClipIds', 'MovieSize');
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
Ypf = @(w) (gaussmf(T(T>=0 & T<1), [w(1) w(3)])*w(5)-gaussmf(T(T>=0 & T<1), [w(2) w(4)])*w(6))+w(7);
freqf_t = Ypf([0.05 0.12 0.08 0.12 1 0.5 0]);
freqf_s = Ypf([0.10 0.24 0.24 0.30 1 0.001 0]);
freqf_t = freqf_t./sum(abs(freqf_t));
freqf_s = freqf_s./sum(abs(freqf_s));
%%
figure; hold on
plot(T(T>=0 & T<1), freqf_t); hold on
plot(T(T>=0 & T<1), freqf_s); hold on
%%
normalize = @(x) (x-min(x(:)))/range(x(:));
figure;
TSctr_HP_gc6f =  [];
TSctr_LP_gc6f =  [];
RFctr_ClipIds_gc6f = [];
for i = 1:nClip
    x = 1:MovieSize(i, 3);
    xq = linspace(1, MovieSize(i, 3), round(SamplingRate*(MovieSize(i, 3)/FrameHz)));
    clipids = RFctr_ClipIds == i;
    nclipids = sum(clipids);
    cTSctr_HP_gc6f =  nan(ngridposi, nclipids);
    cTSctr_LP_gc6f =  nan(ngridposi, nclipids);
    for k = 1:ngridposi
        wave = RFctr_HPf(k, clipids);
        wave = interp1(x, wave, xq, 'pchip');
        wave = [wave(1)*ones(1, SamplingRate) wave(:)' zeros(1, SamplingRate)];
        bpwave_t = conv(freqf_t, wave,'full')/SamplingRate;
        bpwave_s = conv(freqf_s, wave,'full')/SamplingRate;
        Response_t = conv(GC6f, bpwave_t,'full')/SamplingRate;
        Response_s = conv(GC6f, bpwave_t-bpwave_s,'full')/SamplingRate;
        
        nL = length(wave);
        resp_t = resample(Response_t(1:nL), FrameHz, SamplingRate);
        resp_s = resample(Response_s(1:nL), FrameHz, SamplingRate);
        % 
        t = (0:nL-1)/SamplingRate;
        rt = resample(t, FrameHz, SamplingRate, 0);
        rt(end) = rt(end-1)+median(diff(rt(1:end-1)));
        resp_t = resp_t(rt>=1);
        resp_s = resp_s(rt>=1);
        resp_t = resp_t(1:nclipids);
        resp_s = resp_s(1:nclipids);
        cTSctr_HP_gc6f(k, :) = resp_t;
        cTSctr_LP_gc6f(k, :) = resp_s;
    end
    TSctr_HP_gc6f = [TSctr_HP_gc6f cTSctr_HP_gc6f];
    TSctr_LP_gc6f = [TSctr_LP_gc6f cTSctr_LP_gc6f];
    RFctr_ClipIds_gc6f = [RFctr_ClipIds_gc6f i*ones(1, nclipids)];
end
%%
% save('SpatialCorrelationTransience_gc6f.mat', 'TSctr_HP_gc6f', 'TSctr_LP_gc6f','RFctr_ClipIds_gc6f');
%% Examine the spatial correlation
clear R_L;
for i = 1:nClip
    R_L(:, :, i) = corr(TSctr_HP_gc6f(:, RFctr_ClipIds_gc6f ==i)').^2;
end
%%
[X, Y] = meshgrid(position, position);
Coords = [X(:), Y(:)];
Distans = sqrt((Coords(:, 1)'-Coords(:, 1)).^2 + (Coords(:, 2)'-Coords(:, 2)).^2);

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
save('./Results/SpatialCorrection_110123/SpatialCorrelation_Transience_Intercept.mat',...
    'QT95', 'SD2', 'QT95int', 'SDint');
%%
figure;
shadePlot(xv, mean(ypred, 1), std(ypred, [], 1)/sqrt(nClip), 'k', 0.5);
plot(QT95*ones(1, 2), [0 1], '--k');
text(QT95+1, 0.5, '95%');
xlabel('Distance (um)');
ylabel('Correlation coefficient');
title('Spatial correlation of Transience');
xlim([1 76]);
xticks(0:25:75);
xticklabels({'0', '25', '50', '75'});
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylim([0 1.01]);
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sSupFigure6_FeatureSpatialCorrelation_TemporalContrast', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);