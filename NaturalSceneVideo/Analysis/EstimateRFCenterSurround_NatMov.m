%% Loop through all the clip and generate the traces and evaluation (correlation)
radius_c = 50; % in  um
radius_s = 300;
pixel2um = 2;
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
nClip = 11;
MovieSize = nan(nClip, 3);
SharedVar = nan(nClip, 1);
clear RFctr_CS

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
    
    % surround
    radpix = radius_s/pixel2um;
    Mask = sqrt((X-0.5*(1+h)).^2+(Y-0.5*(1+w)).^2) <  radpix;
    Mask = Mask';
    gaufilt = mvnpdf([X(:) Y(:)],0.5*[1+h 1+w],2*pi*radpix*eye(2));
    gaufilt = reshape(gaufilt, w, h)';
    gaufilt = gaufilt/max(gaufilt(:));
    cMask = Mask;
    Mask = Mask.*gaufilt;
    movie = double(reshape(RsizeFilm, [], f))/255-0.5;
    contrat_S= squeeze(movie'*Mask(:)/sum(Mask(:)));
    
    
    
%     plot(contrat_C, 'k');
%     plot(contrat_S, 'm');
    RFctr_CS{i} = contrat_S-contrat_C;
%     contrat_C = (contrat_C-min(contrat_C))./range(contrat_C);
%     contrat_S = (contrat_S-min(contrat_S))./range(contrat_S);
%     RFctr_CS{i} = (contrat_C-contrat_S)./(contrat_C+contrat_S+eps);
    
end
%%
% keyboard;
%%
Act = @(MaxAmp, Alpha, Tau_r, Tau_f, t)  MaxAmp * Alpha * (exp(- t/Tau_f) - exp(- t/Tau_r)) /(Tau_f - Tau_r);
SamplingRate = 1000; % times per second
FrameHz = 52;
T = -1:(1/SamplingRate):3;
DecayTau = 0.142;
RiseTau = 0.045;
figure; plot(T(T>=0 & T<2 ), Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<2)));
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<1));
clear RFctr_CS_gc6f
normalize = @(x) (x-min(x(:)))/range(x(:));
% figure; 
for i = 1:nClip
    x = 1:MovieSize(i, 3);
    xq = linspace(1, MovieSize(i, 3), round(SamplingRate*(MovieSize(i, 3)/FrameHz)));
    
    wave = normalize(RFctr_CS{i});
    wave = interp1(x, wave, xq, 'pchip');
    wave = [wave(1)*ones(1, SamplingRate) wave(:)' zeros(1, SamplingRate)];
    Response = conv(GC6f, wave,'full')/SamplingRate;
    nL = length(wave);
    resp = resample(Response(1:nL), FrameHz, SamplingRate);
    %
    t = (0:nL-1)/SamplingRate;
    rt = resample(t, FrameHz, SamplingRate, 0);
    rt(end) = rt(end-1)+median(diff(rt(1:end-1)));
%     figure; hold on
%     plot(t, wave, 'k');
%     plot(t, Response(1:nL), 'm');
%     plot(rt, resp, 'b');
%     keyboard;
    resp = resp(rt>=1);
    resp = resp(1:length(RFctr_CS{i}));
    RFctr_CS_gc6f{i} = resp;
    
%     subplot(2, 6, i);hold on
%     plot(RFctr_CS_gc6f{i});
end
%%
save(sprintf('./ProcessedData/NatMovCSestimation_gc6f_radius%d.mat', radius_s), 'MovieSize',...
    'RFctr_CS_gc6f','radius_s','radius_c', 'pixel2um');