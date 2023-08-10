close all; clear; clc;
%%
Topic = 'Ai148_AAV-Grm6Cre';
BCTerminalPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis';
SpotPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Circle\Analysis\';
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
% savefilefolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\SpotPathDistance';
savefilefolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\SpotPathDistance_woQualtiy';

%% Load Excel Experimental Document
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
[DocNum, DocStr] = xlsread(FilNam,'ScanLocMatch', 'A:H');
ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));
%%
[Opt, Met]             = imregconfig('Multimodal');
Opt.Epsilon            = 1.5e-6;
Opt.MaximumIterations  = 300;
Opt.InitialRadius      = 6.25e-4;
%%
BCTypes = [50 51 57 58 6 7 89];
BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
%%
% ii = 8, f = 3
% ii = 10, f = 2;
SkipFiles = [8, 3; 10 2];
close all

for ii = 1:16
    clear Nodes
    switch ii
        case 1
            Day = 12019; Cel = 4; Exp = 19; % BC5o
        case 2
            Day = 53018; Cel = 2; Exp = 27; % BC5o
        case 15
            Day = 122718; Cel = 2; Exp = 19; % BC5o
        case 3
            Day = 10119; Cel = 2; Exp = 20; % BC5i
        case 4
            Day = 10119; Cel = 4; Exp = 22; % BC5i
        case 5
            Day = 11919; Cel = 1; Exp = 18; % BC5t
        case 6
            Day = 61018; Cel = 1; Exp = 28; % BC5t
        case 7
            Day = 53018; Cel = 1; Exp = 23; % XBC
        case 16
            Day = 11019; Cel = 3; Exp = 16; %  XBC
        case 8
            Day = 92418; Cel = 2; Exp = 13; % XBC
        case 9
            Day = 10519; Cel = 3; Exp = 18;% BC6
        case 10
            Day = 12019; Cel = 7; Exp = 14; % BC6
        case 11
            Day = 10919; Cel = 3; Exp = 21; % BC7
        case 12
            Day = 12119; Cel = 3; Exp = 16; % BC7
        case 13
            Day = 11819; Cel = 1; Exp = 26; % BC89
        case 14
            Day = 12219; Cel = 1; Exp = 15; % BC89
    end
    run(sprintf('%s/SkeletonGraph/Ai148_AAVGrm6Cre_SkeletonGraph_%d_%d_%d', BCTerminalPath, Day, Cel, Exp));
    if ~exist('Nodes', 'var')
        Nodes = [Nodezxy(:, [3 2]) Nodezxy(:,1)];
    end
    
    %%
    nedge = size(Edges, 1);
    npage = length(PlaneIds);
    Daystr = num2str(Day);
    if length(Daystr) == 5
        Daystr = ['0' Daystr];
    end
    %%
%     figure;
%     scatter3(Nodes(:, 1), Nodes(:, 2), Nodes(:, 3), 20, 'k', 'filled'); hold on
%     for i = 1:nedge
%         cnode = Nodes(Edges(i, :)',:);
%         plot3(cnode(:, 1), cnode(:, 2), cnode(:, 3), 'Color', 0.5*ones(1, 3));
%     end
%     view(0, -90)
    
    %% Align the scaned imaging to the 3d zstack
    PathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\CalciumZstack';
    SkeletonName = sprintf('%s/%s_%s_%d%03d_1_processed.tif', PathName, Topic, Daystr, Cel, Exp);
    
    %% Load BC checkup table
    LoadFileName = [SpotPath 'Results\CenterSurround_Reliability_individualBC_122922.mat'];
    BC = load(LoadFileName);
    Dmt = BC.Dmt;
    BCcolab = BC.Collab2num;
    Name2ColNum = BC.Name2ColNum;
    BC = BC.SummaryD;
    
    cids = BC.Table(:, 1) == Day & BC.Table(:, 2) == Cel;
    exps = unique(BC.Table(cids, 3));
    nexp = length(exps);
    skipids = SkipFiles(SkipFiles(:, 1) == ii, 2);
    for f = 1:nexp
        if ~isempty(skipids) && ismember(f, skipids)
            continue
        end
%         close all
        savefilename = sprintf('%s/%s_SpotStimLengthConstant_%d_%d_%d.mat',...
            savefilefolder, Topic, Day, Cel, exps(f));
       
        %%
        filename = sprintf('translation/Ai148_AAV-Grm6Cre_Translation_%s_%d%03d', Daystr, Cel, exps(f));
        load([SpotPath filename], 'I');
        Mov = std(I, [], 3);
        Mov = Mov/quantile(Mov(:), 0.99);
        figure; imshow(Mov);
        title(sprintf('%s-%d-%d',  Daystr, Cel, exps(f)));
        
        %         if exist(savefilename, 'file')
        %             S = load(savefilename);
        %             if isfield(S, 'dsrd') && isfield(S, 'bctype')
        %                 continue
        %             else
        %                 clear S
        %             end
        %         end
        
        % Reconstruct 3D model of the cells
        %         FilNam = sprintf('%sROIs/%s_ROIs_MorphSeg_%s_%d%03d.mat',SpotPath, Topic, Daystr, Cel, exps(f));
        %         load(FilNam, 'wROISig', 'wSlcROI', 'Cent');
        % Get the ROI and its center position
        
        %% Align the scaned imaging to the 3d zstack
        % register to the specific set of plane in zstack
        % project the ROI in two 2d to the 3D structure (nodes and edges) by the
        % minimum distance
        Pla = DocNum(DocNum(:, ColName2Ind('Day')) == Day & DocNum(:, ColName2Ind('Cel')) == Cel &...
            DocNum(:, ColName2Ind('Exp')) == exps(f), ColName2Ind('Pla'));
        rids = DocNum(:, ColName2Ind('Day')) == Day & DocNum(:, ColName2Ind('Cel')) == Cel &...
            DocNum(:, ColName2Ind('Pla')) == Pla(1);
        pids = DocNum(rids & ~isnan(DocNum(:, ColName2Ind('PlaneStart'))),...
            [ColName2Ind('PlaneStart'), ColName2Ind('PlaneEnd')]);
        if ~isempty(pids)
            PlaneIds = pids(1, 1):pids(1, 2);
        end
        Info = imfinfo(SkeletonName);
        TifLink = Tiff(SkeletonName, 'r');
        anchorImg = nan(Info(1).Height, Info(1).Width, npage);
        for i = 1:npage
            TifLink.setDirectory(PlaneIds(i));
            anchorImg(:, :, i) = TifLink.read();
        end
        TifLink.close();
        anchorImg = mean(anchorImg, 3);
        anchorImg = anchorImg./quantile(anchorImg(:), 0.995);
        Fxd = anchorImg;
        
        %%
        savename = sprintf('Ai148_AAV-Grm6Cre_CPRegistration_%s_%d%03d.mat', Daystr, Cel, exps(f));
        saveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\Analysis\ControlPointRegistration\';
        if ~exist([saveFolder savename], 'file')
            [mp,fp] = cpselect(Mov,Fxd,Wait=true);
            save([saveFolder savename], 'mp', 'fp');
        else
            load([saveFolder savename]);
        end
        t = fitgeotrans(mp,fp,"projective");
        Rfixed = imref2d(size(Fxd));
        MovRgt = imwarp(Mov,t,OutputView=Rfixed);
        
        % RefRgt = imregtform(Mov, Fxd, 'similarity', Opt, Met);
        % MovRgt = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Fxd)));
        % RefRgt = imregtform(Mov, Fxd, 'affine', Opt, Met);%'similarity'
        % MovRgt = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Fxd)));
        
        FilNam = sprintf('%sROIs/%s_ROIs_MorphSeg_%s_%d%03d.mat',SpotPath, Topic, Daystr, Cel,...
            exps(f) );
        load(FilNam, 'wROISig', 'wSlcROI');
        RmIds = any(isnan(wROISig), 2);
        wSlcROI(:, :, RmIds) = [];
        wROISig(RmIds, :) = [];
        
        nROI = size(wSlcROI, 3);
        cmap = nan(size(wSlcROI, 1), size(wSlcROI, 2));
        bbl = ones(size(wSlcROI, 1), size(wSlcROI, 2));
        for j = 1:nROI
            cIds = squeeze(wSlcROI(:, :, j));
            cmap(cIds>0) = j;
        end
        cmap = imwarp(cmap,t,'interp', 'nearest', 'OutputView',imref2d(size(Fxd)));
        cmap(cmap == 0) = nan;
        bbl = imwarp(bbl,t,'interp', 'nearest', 'OutputView',imref2d(size(Fxd)));
        bbl = bwboundaries(bbl,'noholes');
        bbl = bbl{1};
        cbbl = zeros(size(cmap, 1),size(cmap, 2), 3);
        for j = 1:size(bbl, 1)
            cbbl(bbl(j, 1), bbl(j, 2), 1) = 1;
        end
        ccmap = zeros(size(cmap, 1)*size(cmap, 2), 3);
        Colors = parula(nROI);
        % Colors = Colors(randperm(nROI), :);
        for j = 1:nROI
            for i = 1:3
                ccmap(:, i) = ccmap(:, i) + (cmap(:)==j)*Colors(j, i);
            end
        end
        ccmap = reshape(ccmap, size(cmap, 1), size(cmap, 2), 3);
        if ishandle(1)
            close(1)
        end
        figure(1);
        subplot(1, 4, 1);
        imshow(Fxd,'InitialMagnification',200);
        subplot(1, 4, 2);
        imshow(Mov,'InitialMagnification',200);
        subplot(1, 4, 3);
        C = repmat(Fxd, 1, 1, 3);
        MovRgt = repmat(MovRgt, 1, 1, 3);
        C(MovRgt>0) = C(MovRgt>0)*0.2+0.8*MovRgt(MovRgt>0);
        cC = [];
        for i = 1:3
            c = squeeze(C(:, :, i));
            c(squeeze(cbbl(:, :, 1))==1) = 0;
            cC = cat(3, cC, c);
        end
        C = cC+cbbl;
        imshow(C,'InitialMagnification',200);
        subplot(1, 4, 4);
        C = repmat(Fxd, 1, 1, 3);
        C(ccmap>0) = ccmap(ccmap>0);
        cC = [];
        for i = 1:3
            c = squeeze(C(:, :, i));
            c(squeeze(cbbl(:, :, 1))==1) = 0;
            cC = cat(3, cC, c);
        end
        C = cC+cbbl;
        imshow(C,'InitialMagnification',200);
       
        %%
        % SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
        % FleNam = sprintf('%sFigure8_Heterogeneity_registration_segmentation_%s_%d_%d', SaveFolder,...
        %      num2str(Day), Cel, Exp);
        % print('-depsc','-painters','-loose', '-r300',FleNam)
        % saveas(gcf,[FleNam '.png']);
        %%
        Cent = nan(nROI, 2);
        for i = 1:nROI
            Cent(i, :) = getROIloc(cmap == i);
        end
        tz = mean(PlaneIds);
        Cent = [Cent tz*ones(nROI, 1)];
        %% Generate temporal nodes for the edges
        ExdNodes = [Nodes (1:size(Nodes, 1))'];
        step = sqrt(0.5);
        for i = 1:nedge
            cnode = Nodes(Edges(i, :)',:);
            x = cnode(:, 1);
            y = cnode(:, 2);
            z = cnode(:, 3);
            dist = sqrt(sum(diff(cnode, [], 1).^2));
            nstep = round(dist/step);
            xq = linspace(x(1), x(2), nstep);
            yq = linspace(y(1), y(2), nstep);
            zq = linspace(z(1), z(2), nstep);
            ExdNodes = [ExdNodes; xq(2:end-1)' yq(2:end-1)' zq(2:end-1)' i*ones(length(xq(2:end-1)), 1)];
        end
        ExdNodes(:, 1:3) = ExdNodes(:, 1:3).*scales;
        Cent = Cent.*scales;
        ProjNodes = nan(nROI, 1);
        for i = 1:nROI
            [~, ProjNodes(i)] = min(sqrt(sum((ExdNodes(:, 1:3)-Cent(i, :)).^2, 2)));
        end
        %%
        IsDisplay = 0;
        if IsDisplay
            figure;
            cnode = Nodes.*scales;
            scatter3(cnode(:, 1), cnode(:, 2), cnode(:, 3), 20, 'k', 'filled'); hold on
            [X, Y] = meshgrid(min(cnode(:, 1)):max(cnode(:, 1)), min(cnode(:, 2)):max(cnode(:, 2)));
            Z = ones(size(X))*tz*scales(3);
            C = ones(size(Z, 1), size(Z, 2), 3).*reshape([235, 21, 156]/255, 1, 1, 3);
            surf(X, Y, Z, C, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            for i = 1:nedge
                cnode = Nodes(Edges(i, :)',:).*scales;
                plot3(cnode(:, 1), cnode(:, 2), cnode(:, 3), 'Color', 0.5*ones(1, 3));
            end
            for i = 1:nROI
                cCent = Cent(i, :);
                scatter3(Cent(i, 1), Cent(i, 2), Cent(i, 3), 20, 'm', 'filled');
            end
            for i = 1:nROI % for edge
                cnode = [Cent(i, :); ExdNodes(ProjNodes(i), 1:3)];
                plot3(cnode(:, 1), cnode(:, 2), cnode(:, 3), 'Color', 0.8*[1 0 0]);
            end
            view(0, -90)
        end
        % 
%         keyboard;
        %%
        % SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
        % FleNam = sprintf('%sFigure7_Heterogeneity_SkeletonPlaneIntersect_%s_%d_%d_v2',...
        %     SaveFolder, num2str(Day), Cel, Exp);
        % print('-depsc','-painters','-loose', '-r300',FleNam)
        % saveas(gcf,[FleNam '.png']);
        %% Get the distance between pair of each ROI
        % Path to the projection points
        % Projection points to the two ends
        % Use distance in 3d to replace weights
        GraphD = nan(nROI, 3);
        for i = 1:nROI
            cnode = [Cent(i, :); ExdNodes(ProjNodes(i), 1:3)];
            GraphD(i, :) = [i, nROI+i, sqrt(sum(diff(cnode, [], 1).^2))];
        end
        cNode = Nodes.*scales;
        for i = 1:nROI
            cEdge = Edges(ExdNodes(ProjNodes(i), 4), :)';
            cdist = sqrt(sum((cNode(cEdge,:)-ExdNodes(ProjNodes(i), 1:3)).^2, 2));
            cGraphD = [cEdge+nROI*2 (nROI+i)*ones(2, 1) cdist];
            GraphD = [GraphD; cGraphD];
        end
        GraphD = [GraphD; Edges+nROI*2 sqrt(sum((cNode(Edges(:, 1), :)-cNode(Edges(:, 2), :)).^2, 2))];
        
        AllNodes = [Cent; ExdNodes(ProjNodes, 1:3); cNode];
        G = graph(GraphD(:, 1), GraphD(:, 2), GraphD(:, 3));
        % calculate shortest distance for all ROI pairs
        PairD = nan(nROI, nROI);
        distD = nan(nROI, nROI);
        for i = 1:nROI
            for j = 1:nROI
                if i <= j
                    [~, PairD(i, j)] = shortestpath(G,i,j);
                    distD(i, j) = sqrt(sum((Cent(i, :)-Cent(j, :)).^2));
                else
                    PairD(i, j) = PairD(j, i);
                    distD(i, j) = distD(j, i);
                end
            end
        end
        PairD(eye(length(PairD))==1) = nan;
        distD(eye(length(distD))==1) = nan;
        %% demostrate one short path
        [P,d] = shortestpath(G,1,20);
        IsDisplay = 0;
        if IsDisplay
            for i = 1:length(P)-1
                plot3(AllNodes(P([i i+1]),1), AllNodes(P([i i+1]),2), AllNodes(P([i i+1]),3), 'Color', 0.8*[0 0 1]);
            end
        end
        %%
        addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Circle\Analysis');
        cids = find(BC.Table(:, 1) == Day & BC.Table(:, 2) == Cel & BC.Table(:, 3) == exps(f));
        Quality = BC.Table(cids, Name2ColNum('QL'));
        cids = cids(Quality.^2>0.3);
        incIds = BC.Table(cids, Name2ColNum('ROIId'));
        assert(max(incIds) <= nROI);
        dtest = BC.ClusterFeatures(cids, :);
        options.EstimationType = 'surround';
        [dsrd, ~] = estimateCenterSurround(dtest, options);
        options.EstimationType = 'linearity';
        [dlnr, ~] = estimateCenterSurround(dtest, options);
        options.EstimationType = 'transience';
        [dtrs, ~] = estimateCenterSurround(dtest, options);
        
        %%
        bctype = BC.Table(cids(1), 4);
        
        %%
        keyboard;
        continue
        %%
        X = dtest+1e-6*randn(size(dtest));
        Pairdist = squareform(pdist(X, 'correlation'));
        %%
        cPairD = PairD(incIds, incIds);
        cdistD = distD(incIds, incIds);
        BCTypes = [50 51 57 58 6 7 89];
        BCTypeLabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
        figure;
        subplot(2, 2, 1);
        Pairdist(eye(length(Pairdist))==1) = nan;
        scatter(cPairD(:), Pairdist(:), 10, 'k', 'filled'); hold on
        % gprMdl = fitrgp(cPairD(~isnan(cPairD)), Pairdist(~isnan(cPairD)),'Basis','linear','FitMethod','exact',...
        % 'PredictMethod','exact','KernelFunction','squaredexponential',...
        % 'KernelParameters', [1 0.01], 'Sigma',0.2);
        xm = linspace(quantile(cPairD(~isnan(cPairD)), 0.05), quantile(cPairD(~isnan(cPairD)), 0.95), 30);
        % % phi = [mean(std(xm));std(ym)/sqrt(2)]
        % ym = predict(gprMdl, xm');
        % plot(xm, ym, 'b');
        % shadePlot(xm, ym, diff(ed, [], 2)/2, [0 0 1], 0.5, 1);
        g = fittype('a-b*exp(-c*x)');
        x = cPairD(~isnan(cPairD));
        y = Pairdist(~isnan(cPairD));
        f0 = fit(x,y,g,'StartPoint',[mean(y) 0.1 0.1]);
        plot(xm, f0(xm), 'b');
        % [xm, ym, ed, nd] = divdotdistrib(cPairD(~isnan(cPairD)), Pairdist(~isnan(Pairdist)));
        % shadePlot(xm, ym, ed./sqrt(nd), [0 0 1], 0.5, 1);
        box off
        RealIds = ~isnan(cPairD) & ~isnan(Pairdist);
        [R, P] = corrcoef(cPairD(RealIds), Pairdist(RealIds));
        xlabel('Path distance (um)');
        ylabel('Difference (abst. unit)');
        title(sprintf('Correlation \n r=%.02G (p=%.03G)', R(2), P(2)));
        
        subplot(2, 2, 2);
        targetD = abs(dsrd-dsrd');
        targetD(eye(length(targetD))==1) = nan;
        scatter(cPairD(:), targetD(:), 10, 'k', 'filled'); hold on
        % [xm, ym, ed, nd] = divdotdistrib(cPairD(~isnan(cPairD)), targetD(~isnan(targetD)));
        % shadePlot(xm, ym, ed./sqrt(nd), [0 0 1], 0.5, 1);
        g = fittype('a-b*exp(-c*x)');
        x = cPairD(~isnan(cPairD));
        y = targetD(~isnan(cPairD));
        f0 = fit(x,y,g,'StartPoint',[mean(y) 0.1 0.1]);
        plot(xm, f0(xm), 'b');
        box off
        RealIds = ~isnan(cPairD) & ~isnan(targetD);
        [R, P] = corrcoef(cPairD(RealIds), targetD(RealIds));
        xlabel('Path distance (um)');
        ylabel('Difference (abst. unit)');
        title(sprintf('Surround strength \n r=%.02G (p=%.03G)', R(2), P(2)));
        
        subplot(2, 2, 3);
        targetD = abs(dtrs-dtrs');
        targetD(eye(length(targetD))==1) = nan;
        scatter(cPairD(:), targetD(:), 10, 'k', 'filled'); hold on
        g = fittype('a-b*exp(-c*x)');
        x = cPairD(~isnan(cPairD));
        y = targetD(~isnan(cPairD));
        f0 = fit(x,y,g,'StartPoint',[mean(y) 0.1 0.1]);
        plot(xm, f0(xm), 'b');
        % [xm, ym, ed, nd] = divdotdistrib(cPairD(~isnan(cPairD)), targetD(~isnan(targetD)));
        % shadePlot(xm, ym, ed./sqrt(nd), [0 0 1], 0.5, 1);
        box off
        RealIds = ~isnan(cPairD) & ~isnan(targetD);
        [R, P] = corrcoef(cPairD(RealIds), targetD(RealIds));
        xlabel('Path distance (um)');
        title(sprintf('Transience \n r=%.02G (p=%.03G)', R(2), P(2)));
        
        subplot(2, 2, 4);
        targetD = abs(dlnr-dlnr');
        targetD(eye(length(targetD))==1) = nan;
        scatter(cPairD(:), targetD(:), 10, 'k', 'filled'); hold on
        g = fittype('a-b*exp(-c*x)');
        x = cPairD(~isnan(cPairD));
        y = targetD(~isnan(cPairD));
        f0 = fit(x,y,g,'StartPoint',[mean(y) 0.1 0.1]);
        plot(xm, f0(xm), 'b');
        % [xm, ym, ed, nd] = divdotdistrib(cPairD(~isnan(cPairD)), targetD(~isnan(targetD)));
        % shadePlot(xm, ym, ed./sqrt(nd), [0 0 1], 0.5, 1);
        box off
        RealIds = ~isnan(cPairD) & ~isnan(targetD);
        [R, P] = corrcoef(cPairD(RealIds), targetD(RealIds));
        xlabel('Path distance (um)');
        title(sprintf('ON-OFF continuity \n r=%.02G (p=%.03G)', R(2), P(2)));
        
        sgtitle(sprintf('%s (%d-%d-%d)', BCTypeLabels{find(BC.Table(cids(1), 4)==BCTypes)},...
            Day, Cel, exps(f)));
        
        %%
        SimulationOfLengthConstant2CorrVariation
        %%
%                 keyboard;
        %%
       
        %         save(savefilename, 'LenConsts', 'np', 'Dd', 'De','Dq', 'Ds', 'Dc',...
        %             'Dcm', 'Pairdist', 'cPairD', 'bctype', 'dsrd', 'dlnr', 'dtrs' );
        save(savefilename, 'LenConsts', 'np', 'Ddp', 'Ddi', 'De','Dq', 'Da', 'Db', 'Ds','Dac',...
            'Pairdist', 'cPairD', 'bctype', 'dsrd', 'dlnr', 'dtrs', 'cdistD' );
    end
end
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure8_Heterogeneity_distancedependentvariation_%s_%d_%d', SaveFolder,...
%      num2str(Day), Cel, Exp);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
