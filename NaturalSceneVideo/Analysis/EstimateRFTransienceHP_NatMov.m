%% Loop through all the clip and generate the traces and evaluation (correlation)
radius_c = 50; % in  um
PassHz = 4;
pixel2um = 2;
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
nClip = 11;
MovieSize = nan(nClip, 3);
SharedVar = nan(nClip, 1);

clear RFctr_HP
 
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
    Mask = sqrt((X-0.5*(1+h)).^2+(Y-0.5*(1+w)).^2) <  radpix;
    Mask = Mask';
    gaufilt = mvnpdf([X(:) Y(:)],0.5*[1+h 1+w],2*pi*radpix*eye(2));
    gaufilt = reshape(gaufilt, w, h)';
    gaufilt = gaufilt/max(gaufilt(:));
    cMask = Mask;
    Mask = Mask.*gaufilt;
    movie = double(reshape(RsizeFilm, [], f))/255-0.5;
    contrat_C = squeeze(movie'*Mask(:)/sum(Mask(:)));
    RFctr_HP{i} =contrat_C;
end

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
clear TSctr_HP_gc6f TSctr_LP_gc6f
normalize = @(x) (x-min(x(:)))/range(x(:));
figure; 
for i = 1:nClip
    x = 1:MovieSize(i, 3);
    xq = linspace(1, MovieSize(i, 3), round(SamplingRate*(MovieSize(i, 3)/FrameHz)));
    
    wave = RFctr_HP{i};
    wave = interp1(x, wave, xq, 'pchip');
    wave = [wave(1)*ones(1, SamplingRate) wave(:)' zeros(1, SamplingRate)];
    bpwave_t = conv(freqf_t, wave,'full')/SamplingRate;
    bpwave_s = conv(freqf_s, wave,'full')/SamplingRate;
    Response_t = conv(GC6f, bpwave_t,'full')/SamplingRate;
    Response_s = conv(GC6f, bpwave_t-bpwave_s,'full')/SamplingRate;
%     figure; hold on
% %     plot(wave, 'k');
%     plot(bpwave_t, 'b');
%     plot(bpwave_s, 'm');
%     plot(bpwave_t - bpwave_s, 'k');
%     keyboard;
    
    nL = length(wave);
    resp_t = resample(Response_t(1:nL), FrameHz, SamplingRate);
    resp_s = resample(Response_s(1:nL), FrameHz, SamplingRate);
    %
    t = (0:nL-1)/SamplingRate;
    rt = resample(t, FrameHz, SamplingRate, 0);
    rt(end) = rt(end-1)+median(diff(rt(1:end-1)));
    resp_t = resp_t(rt>=1);
    resp_s = resp_s(rt>=1);
    resp_t = resp_t(1:length(RFctr_HP{i}));
    resp_s = resp_s(1:length(RFctr_HP{i}));
    TSctr_HP_gc6f{i} = resp_t;
    TSctr_LP_gc6f{i} = resp_s;
    
    subplot(2, 6, i);hold on
    plot(TSctr_HP_gc6f{i});
    plot(TSctr_LP_gc6f{i});
end

%%
save(sprintf('./ProcessedData/NatMovHPestimation_gc6f_radius%d.mat', radius_c), 'MovieSize',...
    'TSctr_HP_gc6f','TSctr_LP_gc6f','PassHz','radius_c', 'pixel2um');