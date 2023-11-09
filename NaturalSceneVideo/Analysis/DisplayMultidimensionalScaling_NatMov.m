IsCombiMultiROIs = 0; % 0 for luminance repeat measurement, 1 for all other encoding space
HasProcessed = 0;
clear PL PLTab
nclip = 11;
nType = 7;
% lockoffset = [80];%76
% IsCellAlignment = 1;
% qualthr = 0.2; % R-squared
addpath(genpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource\ExpVecEDM'));

roitype = 3;%1: 'individual', 2: 'clusters' 3:'cluster rerun feature correlation'
IsRemove5Hz = 1;
UpSamplingFz = 100;
DistanceType = 'correlation';
cids = ismember(SpotTab(:, 1), 1:max(SpotTab(:, 1), [], 1,'omitnan')) &...
    ismember(SpotTab(:, 2), 1:nType);
savetimemark = '2023042702';
%%
if ~HasProcessed
    switch roitype
        case 1
            Pairdist = nan(sum(cids), sum(cids), nclip);
            Resp = MovResp(cids, :);
            Qual = MovQual(cids, :);
            typ = SpotTab(cids, 2);
            oids = SpotTab(cids, 1);
        case 3
            roiids = unique(SpotTab(cids, 1));
            nroiids = length(roiids);
            Resp = nan(nroiids, size(MovResp, 2));
            Qual = nan(nroiids, nclip);
            rept = nan(nroiids, 1);
            NLctr = nan(nroiids, nclip);
            fNLmvctr = nan(nroiids, 2, nclip);
            fNLnsctr = nan(nroiids, 2);
            dNLctr = nan(nroiids, nclip);
            BSctr = nan(nroiids, nclip);
            MSctr = nan(nroiids, nclip);
            MSActr = nan(nroiids, nclip);
            MSAdctr = nan(nroiids, nclip);
            DVctr = nan(nroiids, nclip);
            CSctr = nan(nroiids, nclip);
            TSHctr = nan(nroiids, nclip);
            TSLctr = nan(nroiids, nclip);
            lagD = nan(nroiids, nclip, 5);
            BSlag = nan(nroiids, 2);
            typ = nan(nroiids, 1);
            Plane = nan(nroiids, 2);
            oids = roiids;
            nrepeats = nan(nroiids, 1);
            % find columns
            findcol = @(in) find(cellfun(@(x) strcmpi(x, in), FeatureNames));
            for i = 1:nroiids
                clc
                fprintf('progress...%d/%d \n', i, nroiids);
                ids = find(cids & SpotTab(:, 1) == roiids(i));
                
                if IsCombiMultiROIs
                    qualthr = 0.2; % R-squared
                    [Resp(i, :), Qual(i, :)] = AssembleClipbyQuality(MovResp(ids, :),clipids, 'outmedian');
                    rept(i) = median(MovRpt(ids));
                    %                     if rept(i) > 1
                    %                         Qual(i, :) = median(MovQual(ids, :), 1);
                    %                     end
                else
                    qualthr = 0.1; % R-squared
                    [~, maxid] = max(median(MovQual(ids, :), 2, 'omitnan'));
                    Qual(i, :) = MovQual(ids(maxid), :);
                    Resp(i, :) = MovResp(ids(maxid), :);
                    rept(i) = MovRpt(ids(maxid));
                end
                for j = 1:length(ids)
                    % for fNLmvctr
                    for k = 1:nclip
                        dids = ids(MovfNL(ids, 2, k) > 0.2);
                        if isempty(dids)
                            continue
                        end
                        w = MovfNL(dids, 2, k).^2;
                        w = w./sum(w);
                        fNLmvctr(i, 1, k) = 10^(w(:)'*log10(MovfNL(dids, 1, k)));
                        fNLmvctr(i, 2, k) = sqrt(w(:)'*MovfNL(dids, 2, k).^2);
                    end
                    % for fNLnsctr
                    dids = ids(NoiTab(ids, 2) > 0.2);
                    if isempty(dids)
                        continue
                    end
                    w = NoiTab(dids, 2).^2;
                    w = w./sum(w);
                    fNLnsctr(i, 1) = 10^(w(:)'*log10(NoiTab(dids, 1)));
                    fNLnsctr(i, 2) = sqrt(w(:)'*NoiTab(dids, 2).^2);
                end
                %                 [~, maxid] = max(mean(squeeze(MovfNL(ids, 2, :)).^2, 2));
                %                 fNLmvctr(i, :, :) = squeeze(MovfNL(maxid, :, :));
                %                 [~, maxid] = max(NoiTab(ids, 2).^2);
                %                 fNLnsctr(i, :) = NoiTab(maxid, 1:4);
                if IsCellAlignment
                    a = Resp(i, :);
                    %                     for j = 1:nclip
                    %                         b = a(clipids'==j);
                    %                         a(clipids'==j) = b-mean(b);
                    %                     end
                    [r, lag] = xcorr(Features(findcol('RF_L'), :), a, UpSamplingFz, 'coeff');
                    r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
                    [BSlag(i, 1), mid] = max(r);
                    BSlag(i, 2) = 1-lag(mid)/UpSamplingFz;
                    BSlag(i, 1) = BSlag(i, 1)^2;
                    lockoffset = mid-100;
                    for j = 1:nclip
                        b = a(clipids'==j);
                        b = b(1:length(b)-lockoffset);
                        for k = 1:10
                            switch k
                                case 1
                                    c = Features(findcol('mresp'),clipids'==j);
                                case 2
                                    c = Features(findcol('RF_L'),clipids'==j);
                                case 3
                                    c = Features(findcol('RF_NL'),clipids'==j);
                                case 4
                                    c = Features(findcol('RF_CS'),clipids'==j);
                                case 5
                                    c = Features(findcol('MS_V'),clipids'==j);
                                case 6
                                    c = Features(findcol('TS_HP'),clipids'==j);
                                case 7
                                    c = Features(findcol('TS_LP'),clipids'==j);
                                case 8
                                    c = Features(findcol('RF_dNL'),clipids'==j);
                                case 9
                                    c = Features(findcol('MS_VA'),clipids'==j);
                                case 10
                                    c = Features(findcol('MS_VAd'),clipids'==j);
                            end
                            c = c((lockoffset+1):end);
                            r = corr(c(:), b(:));
                            switch k
                                case 1
                                    DVctr(i, j) = r^2;
                                case 2
                                    BSctr(i, j) = r^2;
                                case 3
                                    NLctr(i, j) = r^2;
                                case 4
                                    CSctr(i, j) = r^2;
                                case 5
                                    MSctr(i, j) = r^2;
                                case 6
                                    TSHctr(i, j) = r^2;
                                case 7
                                    TSLctr(i, j) = r^2;
                                case 8
                                    dNLctr(i, j) = r^2;
                                case 9
                                    MSActr(i, j) = r^2;
                                case 10
                                    MSAdctr(i, j) = r^2;
                            end
                        end
                    end
                else
                    for j = 1:nclip
                        if ~isempty(lockoffset)
                            b = Resp(i, clipids'==j);
                            b = b(1:length(b)-lockoffset);
                            for k = 1:10
                                switch k
                                    case 1
                                        a = Features(findcol('mresp'), clipids'==j);
                                    case 2
                                        a = Features(findcol('RF_L'),  clipids'==j);
                                    case 3
                                        a = Features(findcol('RF_NL'), clipids'==j);
                                    case 4
                                        a = Features(findcol('RF_CS'), clipids'==j);
                                    case 5
                                        a = Features(findcol('MS_V'),  clipids'==j);
                                    case 6
                                        a = Features(findcol('TS_HP'), clipids'==j);
                                    case 7
                                        a = Features(findcol('TS_LP'), clipids'==j);
                                    case 8
                                        a = Features(findcol('RF_dNL'), clipids'==j);
                                    case 9
                                        a = Features(findcol('MS_VA'),  clipids'==j);
                                    case 10
                                        a = Features(findcol('MS_VAd'),  clipids'==j);
                                end
                                a = a((lockoffset+1):end);
                                r = corr(a(:), b(:));
                                switch k
                                    case 1
                                        DVctr(i, j) = r^2;
                                    case 2
                                        BSctr(i, j) = r^2;
                                    case 3
                                        NLctr(i, j) = r^2;
                                    case 4
                                        CSctr(i, j) = r^2;
                                    case 5
                                        MSctr(i, j) = r^2;
                                    case 6
                                        TSHctr(i, j) = r^2;
                                    case 7
                                        TSLctr(i, j) = r^2;
                                    case 8
                                        dNLctr(i, j) = r^2;
                                    case 9
                                        MSActr(i, j) = r^2;
                                    case 10
                                        MSAdctr(i, j) = r^2;
                                end
                            end
                            
                        else
                            % mean intensity
                            [r, lag] = xcorr(Features(findcol('RF_L'), clipids'==j), Resp(i, clipids'==j), UpSamplingFz, 'coeff');
                            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
                            [maxr, mid] = max(r);
                            lagD(i, j, 1) = 1-lag(mid)/UpSamplingFz;
                            BSctr(i, j) = max(r)^2;
                            % local contrast
                            [r, lag] = xcorr(Features(findcol('RF_NL'), clipids'==j), Resp(i, clipids'==j), UpSamplingFz, 'coeff');
                            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
                            [maxr, mid] = max(r);
                            lagD(i, j, 2) = 1-lag(mid)/UpSamplingFz;
                            NLctr(i, j) = max(r)^2;
                            % center surround
                            [r, lag] = xcorr(Features(findcol('RF_CS'), clipids'==j), Resp(i, clipids'==j), UpSamplingFz, 'coeff');
                            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
                            [maxr, mid] = min(r);
                            lagD(i, j, 3) = 1-lag(mid)/UpSamplingFz;
                            CSctr(i, j) = min(r)^2;
                            % motion sensitivity
                            [r, lag] = xcorr(Features(findcol('MS_V'), clipids'==j), Resp(i, clipids'==j), UpSamplingFz, 'coeff');
                            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
                            [maxr, mid] = max(r);
                            lagD(i, j, 4) = 1-lag(mid)/UpSamplingFz;
                            MSctr(i, j) = max(r)^2;
                            % Control average response
                            [r, lag] = xcorr(Features(findcol('mresp'), clipids'==j), Resp(i, clipids'==j), UpSamplingFz, 'coeff');
                            r(lag<Txcorr(1) | lag>Txcorr(2)) = nan;
                            [maxr, mid] = max(r);
                            lagD(i, j, 5) = 1-lag(mid)/UpSamplingFz;
                            DVctr(i, j) = max(r)^2;
                        end
                    end
                end
                Plane(i, :) = median(MovTab(ids, [1 2]), 1, 'omitnan');
                typ(i) = mean(SpotTab(ids, 2));
                assert(typ(i) == SpotTab(ids(1), 2));
                nrepeats(i) = length(ids);
            end
            %                         keyboard;
            nBC = size(MovBCexp,1);
            if IsRemove5Hz
                rmids = ~ismember(Plane(:, 1), [100319 93019 92519]);
                fNLmvctr = cat(1, MovBCfNLmov(:, 1:2, :), fNLmvctr(rmids, :, :));
                fNLnsctr = [MovBCfNLnoi(:, 1:2); fNLnsctr(rmids, :)];
                Resp = [MovBCexp; Resp(rmids, :)];
                NLctr = [MovBCnl; NLctr(rmids, :)];
                BSctr = [MovBCbs; BSctr(rmids, :)];
                dNLctr = [MovBCdnl; dNLctr(rmids, :)];
                MSctr = [MovBCms; MSctr(rmids, :)];
                MSActr = [MovBCmsa; MSActr(rmids, :)];
                MSAdctr = [MovBCmsad; MSAdctr(rmids, :)];
                DVctr = [MovBCmrp; DVctr(rmids, :)];
                CSctr = [MovBCcs; CSctr(rmids, :)];
                TSHctr = [MovBCtshp; TSHctr(rmids, :)];
                TSLctr = [MovBCtslp; TSLctr(rmids, :)];
                lagD = [squeeze(MovBCspNL(:,3, :)); squeeze(lagD(rmids, :, 1))];
                BSlag = [MovBClag; BSlag(rmids, :)];
                Qual = [MovBCrep; Qual(rmids, :)];
                Rept = [MovBCrpt; rept(rmids)];
                Plane = [MovBCtab(:, 1:2); Plane(rmids, :)];
                typ = [MovBCtab(:, 3); typ(rmids, :)];
                oids = [zeros(nBC, 1); oids(rmids, :)];
                nroiids = sum(rmids);
            else
                fNLmvctr = cat(1, MovBCfNLmov(:, 1:2, :), fNLmvctr);
                fNLnsctr = [MovBCfNLnoi(:, 1:2); fNLnsctr];
                Resp = [MovBCexp; Resp];
                NLctr = [MovBCnl; NLctr];
                BSctr = [MovBCbs; BSctr];
                dNLctr = [MovBCdnl; dNLctr];
                MSctr = [MovBCms; MSctr];
                MSActr = [MovBCmsa; MSActr];
                MSAdctr = [MovBCmsad; MSAdctr];
                DVctr = [MovBCmrp; DVctr];
                CSctr = [MovBCcs; CSctr];
                TSHctr = [MovBCtshp; TSHctr];
                TSLctr = [MovBCtslp; TSLctr];
                lagD = [squeeze(MovBCspNL(:,2, :)); squeeze(lagD(:, :, 1))];
                BSlag = [MovBClag; BSlag];
                Qual = [MovBCrep; Qual];
                Rept = [MovBCrpt; rept];
                Plane = [MovBCtab(:, 1:2); Plane];
                typ = [MovBCtab(:, 3); typ];
                oids = [zeros(nBC, 1); oids];
            end
            iids = Mov2Spot;
            Pairdist = nan(nroiids+nBC, nroiids+nBC, nclip);
        case 2
            roiids = unique(SpotTab(cids, 1));
            nroiids = length(roiids);
            Resp = nan(nroiids, size(MovResp, 2));
            Qual = nan(nroiids, nclip);
            NLctr = nan(nroiids, nclip);
            BSctr = nan(nroiids, nclip);
            MSctr = nan(nroiids, nclip);
            DVctr = nan(nroiids, nclip);
            CSctr = nan(nroiids, nclip);
            typ = nan(nroiids, 1);
            Plane = nan(nroiids, 2);
            oids = roiids;
            nrepeats = nan(nroiids, 1);
            for i = 1:nroiids
                ids = find(cids & SpotTab(:, 1) == roiids(i));
                Resp(i, :) = median(MovResp(ids, :), 1, 'omitnan');
                NLctr(i, :) = median(Movnl(ids, :), 1, 'omitnan');
                BSctr(i, :) = median(squeeze(MovspNL(ids,1, :)).^2, 1, 'omitnan');
                MSctr(i, :) = median(Movms(ids, :), 1, 'omitnan');
                DVctr(i, :) = median(Movdv(ids, :), 1, 'omitnan');
                CSctr(i, :) = median(Movcs(ids, :), 1, 'omitnan');
                Qual(i, :) = median(MovQual(ids, :), 1, 'omitnan');
                Plane(i, :) = median(MovTab(ids, [1 2]), 1, 'omitnan');
                
                typ(i) = mean(SpotTab(ids, 2));
                assert(typ(i) == SpotTab(ids(1), 2));
                nrepeats(i) = length(ids);
            end
            nBC = size(MovBCexp,1);
            if IsRemove5Hz
                rmids = ~ismember(Plane(:, 1), [100319 93019 92519]);
                Resp = [MovBCexp; Resp(rmids, :)];
                NLctr = [MovBCnl; NLctr(rmids, :)];
                BSctr = [MovBCbs; BSctr(rmids, :)];
                MSctr = [MovBCms; MSctr(rmids, :)];
                DVctr = [MovBCmrp; DVctr(rmids, :)];
                CSctr = [MovBCcs; CSctr(rmids, :)];
                Qual = [MovBCrep; Qual(rmids, :)];
                Plane = [MovBCtab(:, 1:2); Plane(rmids, :)];
                typ = [MovBCtab(:, 3); typ(rmids, :)];
                oids = [zeros(nBC, 1); oids(rmids, :)];
                nroiids = sum(rmids);
                nrepeats = [zeros(nBC, 1); nrepeats(rmids)];
            else
                Resp = [MovBCexp; Resp];
                NLctr = [MovBCnl; NLctr];
                BSctr = [MovBCbs; BSctr];
                MSctr = [MovBCms; MSctr];
                DVctr = [MovBCmrp; DVctr];
                CSctr = [MovBCcs; CSctr];
                Qual = [MovBCrep; Qual];
                Plane = [MovBCtab(:, 1:2); Plane];
                typ = [MovBCtab(:, 3); typ];
                oids = [zeros(nBC, 1); oids];
                nrepeats = [zeros(nBC, 1); nrepeats];
            end
            iids = Mov2Spot;
            Pairdist = nan(nroiids+nBC, nroiids+nBC, nclip);
    end
    
    save(sprintf('./ProcessedData/FeatureReExam_%s.mat', savetimemark), 'Resp',...
        'NLctr', 'BSctr', 'MSctr', 'DVctr', 'CSctr','Qual','TSHctr', 'TSLctr',...
        'Plane', 'typ', 'oids', 'iids', 'nroiids', 'nBC', 'nrepeats', 'nclip',...
        'lagD','BSlag', 'fNLnsctr', 'fNLmvctr');
end
% run NaturalMovieVariationProjectionTest
%%
% load(sprintf('./ProcessedData/FeatureReExam_%s.mat', savetimemark));
%%
Pairdist = nan(nroiids+nBC, nroiids+nBC, nclip);
Resp = Resp-min(Resp, [], 2);
Resp = Resp./range(Resp, 2);
for i = 1:nclip
    cResp = Resp(:, clipids == i);
    cResp(Qual(:, i) < qualthr, :) = nan;
    Pairdist(:,:,i) = squareform(pdist(cResp+rand(size(cResp))*1e-3, DistanceType));
end
Pairdist = median(Pairdist, 3, 'omitnan');
%%
Pairdist(eye(size(Pairdist))==1) = nan;
rmids = ~all(isnan(Pairdist), 2);
Pairdist = Pairdist(rmids, rmids);
typ = typ(rmids);
oids = oids(rmids);
nrepeats = nrepeats(rmids);
iids = iids(rmids(1:length(iids)));
Resp = Resp(rmids, :);
Qual = Qual(rmids, :);
Rept = Rept(rmids, :);
NLctr = NLctr(rmids, :);
dNLctr = dNLctr(rmids, :);
BSctr = BSctr(rmids, :);
MSctr = MSctr(rmids, :);
MSActr = MSActr(rmids, :);
MSAdctr = MSAdctr(rmids, :);
DVctr = DVctr(rmids, :);
CSctr = CSctr(rmids, :);
TSHctr = TSHctr(rmids, :);
TSLctr = TSLctr(rmids, :);
fNLmvctr = fNLmvctr(rmids, :, :);
fNLnsctr = fNLnsctr(rmids, :, :);
lagD = lagD(rmids, :);
BSlag = BSlag(rmids, :);
Plane = Plane(rmids, :);
Pairdist(eye(size(Pairdist))==1) = 0;
%%
rmids = sum(isnan(Pairdist), 2)<1;
Pairdist = Pairdist(rmids, rmids);
typ = typ(rmids);
oids = oids(rmids);
nrepeats = nrepeats(rmids);
iids = iids(rmids(1:length(iids)));
Resp = Resp(rmids, :);
Qual = Qual(rmids, :);
Rept = Rept(rmids, :);
NLctr = NLctr(rmids, :);
dNLctr = dNLctr(rmids, :);
BSctr = BSctr(rmids, :);
MSctr = MSctr(rmids, :);
MSActr = MSActr(rmids, :);
MSAdctr = MSAdctr(rmids, :);
DVctr = DVctr(rmids, :);
CSctr = CSctr(rmids, :);
TSHctr = TSHctr(rmids, :);
TSLctr = TSLctr(rmids, :);
fNLmvctr = fNLmvctr(rmids, :, :);
fNLnsctr = fNLnsctr(rmids, :, :);
lagD = lagD(rmids, :);
BSlag = BSlag(rmids, :);
Plane = Plane(rmids, :);

%%
% figure; imagesc(isnan(Pairdist));
%%
Pairdist(eye(size(Pairdist))==1) = 0;
for i  = 1:size(Pairdist, 1)
    for j = 1:size(Pairdist, 1)
        if Pairdist(i, j) ~= Pairdist(j, i)
            Pairdist(i, j) = mean([Pairdist(i, j) Pairdist(j, i)]);
            Pairdist(j, i) = Pairdist(i, j);
        end
    end
end
Pairdist(Pairdist<0) = 1e-6;
%%
% [Y,eigvals] = cmdscale(Pairdist);
opt = statset('MaxIter', 2000);
numdim = 16;
weightmask = reshape(~isoutlier(Pairdist(:)), size(Pairdist, 1), []); % ; Pairdist<1
[Y,stress, disparity] = mdscale(Pairdist, numdim, 'criterion', 'sstress', 'weights', weightmask, 'start', 'random', 'Options', opt);

celltype = typ;
originalids = oids;
identifiedBC = iids;

violationLevel = testTriangleInequality(Pairdist);
r_recont = nan(1, size(Y, 2));
for i = 1:size(Y, 2)
    Pairdist_reconstruct = squareform(pdist(Y(:,1:i), 'euclidean'));%
    r_recont(i) = corr(Pairdist(triu(ones(size(Pairdist)), 1)==1 & weightmask),...
        Pairdist_reconstruct(triu(ones(size(Pairdist_reconstruct)), 1)==1& weightmask)).^2;
end
%%
figure; hold on
x = 1:size(Y, 2);
y = [r_recont(1), diff(r_recont)];
h = bar( x(1:numdim), y(1:numdim));
h.EdgeColor = 'w';
h.FaceColor = 0.3*ones(1, 3);
yticks(0:0.25:1);
yticklabels({'0', '', '0.5', '', '1'});
ylabel('Variance explained');

% yyaxis right
x = 1:size(Y, 2);
y = r_recont;

plot(x(1:numdim), y(1:numdim), 'b');
% yticks(0.4:0.2:1);
% yticklabels({'0.4', '', '0.8', '', '1'});
box off
xlabel('MDS components');
xticks(1:5:16);
xticklabels({'1', '6', '11', '16'});
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure5_NatMovProjection_ExplainedVariance_01', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
save(sprintf('./ProcessedData/DisplayMultidimensionalScaling_NatMov_%s.mat', savetimemark),...
    'Pairdist', 'Y', 'celltype', 'originalids', 'identifiedBC', 'Resp','clipids', 'Plane',...
    'Qual', 'NLctr', 'BSctr', 'MSctr', 'DVctr', 'CSctr', 'BSlag', 'lagD', 'TSHctr','TSLctr',...
    'MSActr', 'MSAdctr', 'fNLmvctr', 'fNLnsctr');
%%
load(sprintf('./ProcessedData/DisplayMultidimensionalScaling_NatMov_%s.mat', savetimemark));
typ = celltype;
oids = originalids;
iids = identifiedBC;

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_NatMovVarianceExplained', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
%% show distribution of each cell type in natural movies
% close all
figure;
MarkSize = 20;
ndot = size(Y, 1);
dazure = @(n) [linspace(0, 255, n); linspace(254, 0, n); linspace(255, 151, n)]'/255;
Colors = dazure(nType);
Dimids = [1 2 3];

ylabel('MDS component 2');
xlabel('MDS component 1');

for j = 1:nType
    subplot(2, 4, j); hold on
    tids = find(typ == j);
    %     scatter3(Y(:, Dimids(1)), Y(:, Dimids(2)), Y(:, Dimids(3)), MarkSize, 0.7*ones(1, 3), 'filled');
    %     scatter3(Y(tids, Dimids(1)), Y(tids, Dimids(2)), Y(tids, Dimids(3)), MarkSize, Colors(j, :), 'filled');
    scatter(Y(:, Dimids(1)), Y(:, Dimids(2)), MarkSize,  0.7*ones(1, 3), 'filled');
    scatter(Y(tids, Dimids(1)), Y(tids, Dimids(2)), MarkSize, Colors(j, :), 'filled');
    
    text(-0.2,0.35,BCTypelabels{j},'Color',Colors(j, :),'FontSize',15);
    ylabel('MDS component 2');
    xlabel('MDS component 1');
    yticks(-0.3:0.3:0.3);
    yticklabels({'-0.3', '0', '0.3'});
    xticks(-0.2:0.2:0.4);
    xticklabels({'-0.2', '0', '0.2', '0.4'});
    ylim([-0.35 0.35])
    xlim([-0.28 0.45])
end
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_NatMovCellTypeDistributioninEncodingSpace', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%% Show response trace of each BC type
normalization = @(x) (x-min(x, [], 2))./range(x, 2);
figure; hold on
t = (0:length(clipids)-1)/100;
boundaries = nan(11, 1);
for i = 1:10
    sid = find(clipids == i);
    eid = find(clipids == i+1);
    boundaries(i) = 0.5*(t(sid(end))+t(eid(1)));
end
for i = 1:nType
    for j = 1:11
        cids = find(typ == i);
        ncid = length(cids);
        %     for j = 1:ncid
        %         plot(t, 0.9*normalization(Resp(cids(j), :))-i+1, 'Color', 0.7*ones(1, 3));
        %     end
        %     plot(t, 0.9*normalization(median(Resp(cids, :), 1))-i+1, 'Color', Colors(i, :), 'LineWidth', 2);
        tids = clipids == j;
        shadePlot(t(tids), 0.9*normalization(median(Resp(cids, tids), 1))-i+1,...
            0.9*normalization(std(Resp(cids, tids), [], 1))/sqrt(ncid), Colors(i, :), 0.3, 2);
        if j == 1
            text(1, 0.8-i+1, BCTypelabels{i}, 'Color', Colors(i, :), 'FontSize', 12);
            plot(repmat(boundaries, 1, 2)', repmat([-6 1], length(boundaries), 1)', '--k');
        end
    end
end
plot([1 6], (-nType+0.5)*ones(1, 2), 'k');
plot([1 1], -nType+[0.5 1], 'k');
box off
xlim([0 max(t)]);
axis off
% yticks(1:7);
% yticklabels(BCTypelabels);
xlabel('Time (s)');
%%
% nt = t;
% for i = 1:11
%     tids = clipids == i;
%     ct = nt(tids);
%     nt(tids) = ct-min(ct);
% end
% clear ct tids
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure5_NatMovCellTypeResponses', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%%
% close all
figure;
MarkSize = 20;
ndot = size(Y, 1);
tits = {'Nonlinearlity', 'Coherent motion sensitivity', 'Control', 'Center-surround', 'Transience'};
nplot = length(tits);
mY = median(Y, 1);
Dimids = [1 2 3; 3 2 1];
nrow = size(Dimids, 1);
ChooseClip = 1:11;
cids = Qual < qualthr;
Drange = nan(nplot, 2);
ftyps = nan(ndot, nplot);
for k = 1:nrow
    for j = 1:nplot
        subplot(nrow, nplot, (k-1)*nplot+j); hold on
        switch j
            case 1
                a = squeeze(fNLmvctr(:, 1, :));
                a(squeeze(fNLmvctr(:, 2, :)).^2 < qualthr) = nan;
                [OutColor, OutIds, Drange(j, :)] = GroupScatterColor(median(log10(a(:, ChooseClip)), 2, 'omitnan'), parula(256));
               
            case 2
                a = MSAdctr;
                a(cids) = nan;
                b = BSctr;
                b(cids) = nan;
                b = b/max(b(:));
                [OutColor, OutIds, Drange(j, :)] = GroupScatterColor(log(median(1./a(:, ChooseClip), 2, 'omitnan')), parula(256));
              
            case 3
                a = BSctr;
                a(cids) = nan;
                b = Qual;
                b(cids) = nan;
                [OutColor, OutIds, Drange(j, :)] = GroupScatterColor(log(median(1./a(:, ChooseClip), 2, 'omitnan')), parula(256));
            
            case 4
                a = CSctr;
                a(cids) = nan;
                b = BSctr;
                b(cids) = nan;
                [OutColor, OutIds, Drange(j, :)] = GroupScatterColor(log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan')), parula(256));
                
            case 5
                a = TSHctr;
                a(cids) = nan;
                b = TSLctr;
                b(cids) = nan;
                [OutColor, OutIds, Drange(j, :)] = GroupScatterColor(log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan')), parula(256));
               
        end
        ftyps(OutIds, j) = typ(OutIds);
        switch k
            case 1
                for i = 1:length(OutIds)
                    plot(Y(OutIds(i), 1), Y(OutIds(i), 2), '.', 'Color', OutColor(i, :), 'MarkerSize',MarkSize);
                end
                xlabel(sprintf('MDS component %d', 1));
                xticks(-0.2:0.2:0.4);
                xticklabels({'-0.2', '0', '0.2', '0.4'});
                xlim([-0.28 0.45])
            case 2
                for i = 1:length(OutIds)
                    plot(Y(OutIds(i), 3), Y(OutIds(i), 2), '.', 'Color', OutColor(i, :), 'MarkerSize', MarkSize);
                end
                xlabel(sprintf('MDS component %d', 3));
                xticks(-0.2:0.2:0.2);
                xticklabels({'-0.2', '0', '0.2'});
                xlim([-0.22 0.29])
        end
        ylabel(sprintf('MDS component %d', 2));
        yticks(-0.3:0.3:0.3);
        yticklabels({'-0.3', '0', '0.3'});
        ylim([-0.32 0.35])
        title(tits{j});
    end
end
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% % FleNam = sprintf('%sFigure5_NatMovEncodingSpace_noise', SaveFolder);
% FleNam = sprintf('%sFigure5_NatMovEncodingSpace_natmovcorr', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);


%% By BC Types
figure;
Colors = dazure(nType);
ttbl = nan(nType, nplot);
for j = 1:nplot
    subplot(1, nplot, j); hold on
    switch j
        case 1
            a = squeeze(fNLmvctr(:, 1, :));
            ytic = [];
        case 2
            a = 1./MSAdctr;
            ytic = 1:2:5;
            yticlab = {'1', '3', '5'};
            ylm = [1 5.2];
        case 3
            a = 1./BSctr;
            ytic = [];
        case 4
            a = CSctr;
            b = BSctr;
            ytic = -5:3:1;
            yticlab = {'-5', '-2', '1'};
            ylm = [-5.5 1.5];
        case 5
            a = TSHctr;
            b = TSLctr;
            ytic = -5:4:3;
            yticlab = {'-5', '-1', '3'};
            ylm = [-5.5 3];
    end
    switch j
        case 1
             a(squeeze(fNLmvctr(:, 2, :)).^2 < qualthr) = nan;
             v = median(log10(a(:, ChooseClip)), 2, 'omitnan');
        case {2, 3}
            a(cids) = nan;
            v = log(median(a(:, ChooseClip), 2, 'omitnan'));
        case {4, 5}
            a(cids) = nan;
            b(cids) = nan;
            v = log(median(a(:, ChooseClip)./b(:, ChooseClip), 2, 'omitnan'));
    end
        
    for i = 1:nType
        tids = typ==i & ~isnan(v);
        ntid = sum(tids);
        ttbl(i, j) = ntid;
        plot(0.7*(rand(ntid, 1)-0.5)+i, v(tids), '.', 'Color', Colors(i, :));
        plot(0.7*[-0.5 0.5]+i, median(v(tids))*ones(1, 2), 'Color', Colors(i, :), 'LineWidth', 2);
    end
    if ~isempty(ytic)
        yticks(ytic);
        yticklabels(yticlab);
        ylim(ylm);
    end
    ylabel(tits{j});
    xticks(1:nType);
    xticklabels(BCTypelabels);
    keyboard;
end
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% % FleNam = sprintf('%sFigure5_NatMovEncodingSpace_noise', SaveFolder);
% FleNam = sprintf('%sFigure5_NatMovCellTypeSwamChat_natmovfeatures', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
