clear; close all; clc; 
% modification from @GenerateExpectedGC6fResponsetoNoise in noise folder
% and from @TestSpatialDecorrelation_position_NatMov in NaturalSceneVideo
% folder
%%
FolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\NaturalSceneVideo\Preprocessed\Blue\';
nClip = 11;
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

%%
RFctr_L_gc6f = [];
RFctr_NL_gc6f = [];
RFctr_dNL_gc6f = [];
RFctr_ClipIds = [];  
for i = 1:nClip
    MovieSize = nan(nClip, 3);
    FileName = sprintf('%s0131182_UVProj_%d.mat', FolderPath, i);
    load(FileName);
    h = size(RsizeFilm, 1);
    w = size(RsizeFilm, 2);
    f = size(RsizeFilm, 3);
    MovieSize(i, :) = [h, w, f];
    [X, Y] = meshgrid(1:h, 1:w);
    radpix = radius/pixel2um;
    cRFctr_L_gc6f = nan(ngridposi, f);
    cRFctr_NL_gc6f = nan(ngridposi, f);
    cRFctr_dNL_gc6f = nan(ngridposi, f);
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
        
        contrat_lin = squeeze(movie'*Mask(:)/sum(Mask(:)));
        switch nonlineartype
            case 1
                mInt = squeeze(movie'*cMask(:)/sum(cMask(:)));
                contrat_non = squeeze(max(cat(3, movie-mInt', zeros(size(movie))), [], 3)'*Mask(:)/sum(Mask(:)));
            case 2
                contrat_non = movie;
                contrat_non =contrat_non(cMask(:)==1, :);
                avg = mean(contrat_non, 1);
                Mask = Mask(:);
                weight = Mask(cMask(:)==1);
                contrat_non = sqrt(weight'*(contrat_non-avg).^2/(sum(weight)-1))';
            case 3
                nDiv = 4;
                Divs = nan(nDiv, 2);
                dMask = nan(numel(X), nDiv);
                for j = 1:nDiv
                    Divs(j, :) = radpix*[cos((j-1)*2*pi/nDiv) sin((j-1)*2*pi/nDiv)];
                    dMask(:, j) = sqrt((X(:)-0.5*(1+h)-Divs(j, 1)).^2+(Y(:)-0.5*(1+w)-Divs(j, 2)).^2);
                end
                [~, dMask] = min(dMask, [], 2);
                dMask = reshape(dMask, w, h)';
                contrat_non = zeros(size(contrat_lin));
                for j = 1:nDiv
                    cdMask = reshape((dMask==j).*Mask, [], 1);
                    contrat_non = contrat_non + max([movie'*cdMask/sum(cdMask)-contrat_lin, zeros(size(contrat_lin))], [], 2);
                end
        end
        RFctr_L = contrat_lin;
        RFctr_NL = contrat_non;
        RFctr_dNL = contrat_non-contrat_lin;
                
        x = 1:MovieSize(i, 3);
        xq = linspace(1, MovieSize(i, 3), round(SamplingRate*(MovieSize(i, 3)/FrameHz)));
        
        for j = 1:3
            switch j
                case 1
                    wave = normalize(RFctr_L);
                case 2
                    wave = normalize(RFctr_NL);
                case 3
                    wave = normalize(RFctr_dNL);
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
            resp = resp(rt>=1);
            resp = resp(1:length(RFctr_L));
            switch j
                case 1
                    cRFctr_L_gc6f(k, :) = resp;
                case 2
                    cRFctr_NL_gc6f(k, :) = resp;
                case 3
                    cRFctr_dNL_gc6f(k, :) = resp;
            end
        end
        clc
        fprintf('progress...%d/%d (%d/%d) \n', i, nClip, k, ngridposi);
    end
    RFctr_L_gc6f = [RFctr_L_gc6f cRFctr_L_gc6f];
    RFctr_NL_gc6f = [RFctr_NL_gc6f cRFctr_NL_gc6f];
    RFctr_dNL_gc6f = [RFctr_dNL_gc6f cRFctr_dNL_gc6f];
    RFctr_ClipIds = [RFctr_ClipIds i*ones(1, f)];
end
%%
posi = gridposi*pixel2um;
saveFileName = sprintf('./ProcessedData/NatMovRFSpatialDiffPosition_gc6f_radius%d.mat', radius);
save(saveFileName, 'MovieSize','RFctr_L_gc6f', 'RFctr_NL_gc6f','RFctr_dNL_gc6f', 'RFctr_ClipIds', 'radius', 'pixel2um', 'posi');
A = load(saveFileName);