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
saveFileName = sprintf('./ProcessedData/NatMovRFSpatialDiffPosition_gc6f_radius%d.mat', radius);
%%
if exist(saveFileName, 'file')
    load(saveFileName);
    disp('Loading...');
else
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
        posi = gridposi*pixel2um;
        save(saveFileName, 'MovieSize','RFctr_L_gc6f', 'RFctr_NL_gc6f','RFctr_dNL_gc6f', 'RFctr_ClipIds', 'radius', 'pixel2um', 'posi');
    end
end

%% Examine the spatial correlation
nClip = numel(unique(RFctr_ClipIds));
clear R_L;
for i = 1:nClip
    R_L(:, :, i) = corr(RFctr_NL_gc6f(:, RFctr_ClipIds ==i)').^2;
end
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
save('./Results/SpatialCorrection_110123/SpatialCorrelation_In-centercontrast_Intercept.mat',...
    'QT95', 'SD2', 'QT95int', 'SDint');
%%
figure;
shadePlot(xv, mean(ypred, 1), std(ypred, [], 1)/sqrt(nClip), 'k', 0.5);
plot(QT95*ones(1, 2), [0 1], '--k');
text(QT95+1, 0.5, '95%');
xlabel('Distance (um)');
ylabel('Correlation coefficient');
title('Spatial correlation of in-center contrast');
xlim([1 76]);
xticks(0:25:75);
xticklabels({'0', '25', '50', '75'});
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylim([0 1.01]);

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sSupFigure6_FeatureSpatialCorrelation_InCenterContrast', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
keyboard;
%% Examine the spatial correlation - luminance
nClip = numel(unique(RFctr_ClipIds));
clear R_L;
for i = 1:nClip
    R_L(:, :, i) = corr(RFctr_L_gc6f(:, RFctr_ClipIds ==i)').^2;
end
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
% save('./Results/SpatialCorrection_110123/SpatialCorrelation_luminance_Intercept.mat',...
%     'QT95', 'SD2', 'QT95int', 'SDint');
%%
figure;
shadePlot(xv, mean(ypred, 1), std(ypred, [], 1)/sqrt(nClip), 'k', 0.5);
plot(QT95*ones(1, 2), [0 1], '--k');
text(QT95+1, 0.5, '95%');
xlabel('Distance (um)');
ylabel('Correlation coefficient');
title('Spatial correlation of Luminance');
xlim([1 76]);
xticks(0:25:75);
xticklabels({'0', '25', '50', '75'});
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylim([0 1.01]);

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sSupFigure6_FeatureSpatialCorrelation_Luminance', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);