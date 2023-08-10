clear; close all; clc;

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
Topic = 'Ai148_AAV-Grm6Cre';
FileType = 'ROISignalAlign';
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
T = -1:(1/SamplingRate):3;
DecayTau = 0.142;
RiseTau = 0.045;
GC6f = Act(1, 1, RiseTau, DecayTau, T(T>=0 & T<1));
normalize = @(x) (x-min(x(:)))/range(x(:));
%%
PathName = [FolderPath 'ProcessedData/'];
FileNames = dir(fullfile(PathName, '*.mat'));
FileIds = find(contains({FileNames.name},[Topic '_' FileType]));
nFile = length(FileIds);
IsSelfGen = 0;
for f = 1:nFile
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
    if nL ~= AnkL
        IsSelfGen = 1;
    elseif size(Grids, 1) == 40
        IsSelfGen = 0; % Ai148_Grm6Cre since 2019
    elseif sum(abs(Grids(:) - AnkGrids(:))) > 1e-6
        IsSelfGen = 1;
    end
    if IsSelfGen
        saveFileName = sprintf('./ProcessedData/NoiseRFSpatialDiffPosition_gc6f_radius%d_%s_%d_%d.mat', radius, Day, Cel, Exp);
        if exist(saveFileName, 'file')
            continue
        end
        StimFolderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Noise\Analysis\DataStimulus\';
        FileName = sprintf('%s%s_%s_%d_%d', StimFolderPath, Topic, Day, Cel, Exp);
        load(FileName);
        h = round(size(DG_OUT.Grids, 2)*DG_IN.lengthEachCell/pixel2um);
        w = round(size(DG_OUT.Grids, 1)*DG_IN.lengthEachCell/pixel2um);
        nf = size(DG_OUT.Grids, 3);
        Movie = nan(h, w, nf);
        for i = 1:nf
            Movie(:, :, i) = imresize(squeeze(DG_OUT.Grids(:, :, i))', [h, w], 'nearest');
        end
        RFctr_L_gc6f = nan(ngridposi, nf);
        RFctr_NL_gc6f = nan(ngridposi, nf);
        RFctr_dNL_gc6f = nan(ngridposi, nf);
        
        for i = 1:ngridposi
            movie = Movie;
            [X, Y] = meshgrid(1:h, 1:w);
            radpix = radius/pixel2um;
            Mask = sqrt((X-0.5*(1+h)-gridposi(i, 1)).^2+(Y-0.5*(1+w)-gridposi(i, 2)).^2) <  radpix;
            Mask = Mask';
            gaufilt = mvnpdf([X(:) Y(:)],0.5*[1+h 1+w]+gridposi(i, :),2*pi*radpix*eye(2));
            gaufilt = reshape(gaufilt, w, h)';
            gaufilt = gaufilt/max(gaufilt(:));
            cMask = Mask;
            Mask = Mask.*gaufilt;
            movie = double(reshape(movie, [], nf))-0.5;
            
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
            clear contrat_lin contrat_non
            
           
            FrameHz = 1/mean(diff(DG_OUT.DataTable));
            [p, q] = rat(FrameHz/SamplingRate, 1e-3);
            
            
            x = 1:nf;
            xq = linspace(1, nf, round(SamplingRate*(nf/FrameHz)));
            
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
                resp = resample(Response(1:nL), p, q);
                %
                t = (0:nL-1)/SamplingRate;
                rt = resample(t, p, q, 0);
                rt(end) = rt(end-1)+median(diff(rt(1:end-1)));
                resp = resp(rt>=1);
                resp = resp(1:length(RFctr_L));
                switch j
                    case 1
                        RFctr_L_gc6f(i, :) = resp;
                    case 2
                        RFctr_NL_gc6f(i, :) = resp;
                    case 3
                        RFctr_dNL_gc6f(i, :) = resp;
                end
            end
            clc
            fprintf('Progress...%d/%d (%d/%d) \n', f, nFile, i, ngridposi);
        end
        %%
        MovieSize = size(DG_OUT.Grids);
        posi = gridposi*pixel2um;
        
        save(saveFileName, 'MovieSize','RFctr_L_gc6f', 'RFctr_NL_gc6f','RFctr_dNL_gc6f', 'radius', 'pixel2um', 'posi');
    end
end