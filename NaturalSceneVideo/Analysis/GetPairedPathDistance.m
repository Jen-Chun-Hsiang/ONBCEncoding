for q = 1:2
    filename = sprintf('/Translation/Ai148_AAV-Grm6Cre_Translation_%s_%d%03d', Daystr, Cel, expids(q));
    load([NatMovPath filename], 'I');
    Mov = std(I, [], 3);
    Mov = Mov/quantile(Mov(:), 0.99);
    savename = sprintf('Ai148_AAV-Grm6Cre_CPRegistration_%s_%d%03d.mat', Daystr, Cel, expids(q));
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
    
    FilNam = sprintf('%s/ROIs/%s_ROIs_MorphSeg_%s_%d%03d.mat',NatMovPath, Topic, Daystr, Cel,...
        expids(q) );
    load(FilNam, 'wROISig', 'wSlcROI');
    RmIds = any(isnan(wROISig), 2);
    wSlcROI(:, :, RmIds) = [];
    DAmp = wROISig;
    DAmp(RmIds, :) = [];
    DAmp = std(DAmp, [], 2);
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
    figure;
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
    %         keyboard;
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
    for i = 1:length(P)-1
        plot3(AllNodes(P([i i+1]),1), AllNodes(P([i i+1]),2), AllNodes(P([i i+1]),3), 'Color', 0.8*[0 0 1]);
    end
    switch q
        case 1
            PairD1 = PairD;
            distD1 = distD;
            DAmp1 = DAmp/max(DAmp);
        case 2
            PairD2 = PairD;
            distD2 = distD;
            DAmp2 = DAmp/max(DAmp);
    end
end