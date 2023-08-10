close all; clear all; clc;
%% Load Excel Experimental Document
FilNam = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\Recording_Document\BCTerminal_DataMatrix.xlsx';
[DocNum, DocStr] = xlsread(FilNam,'Sheet2', 'A:G');
ColName2Ind = @(String) find(strcmpi(DocStr(1,:), String));

%%
% LoadFileName = './Results/CenterSurround_Reliability_individualBC_110822.mat'; % for difference -0.2 baseline
% LoadFileName = './Results/CenterSurround_Reliability_individualBC_112222.mat'; % for deconvolution -0.5 baseline
LoadFileName = './Results/CenterSurround_Reliability_individualBC_122122.mat'; 
BC = load(LoadFileName);
Dmt = BC.Dmt;
BCcolab = BC.Collab2num;
Name2ColNum = BC.Name2ColNum;
BC = BC.SummaryD;

LoadFileName = './Results/CenterSurround_Reliability_ROI_Plexus_100722.mat';
PL = load(LoadFileName);
PLcolab = PL.Collab2num;
PL = PL.SummaryD;

%%
figure;
subplot(1, 2, 1);
scatter(BC.Table(:, 6), log10(BC.Table(:, 8)), 5, 'k', 'filled');
subplot(1, 2, 2);
scatter(PL.Table(:, 7), log10(PL.Table(:, 9)), 5, 'k', 'filled');
%%
ProbThr = 0.521;
StdThr = 6;
bc_gids = BC.CelTable(:, BCcolab('QLmed')) > ProbThr;
pl_gids = PL.Table(:, PLcolab('QL')) > ProbThr;
% bc_gids = BC.CelTable(:, BCcolab('QLmed')) > ProbThr &...
%     abs(BC.CelTable(:, BCcolab('CSmed'))) < StdThr*nanstd(BC.CelTable(:, BCcolab('CSmed'))) &...
%     abs(BC.CelTable(:, BCcolab('TSmed'))) < StdThr*nanstd(BC.CelTable(:, BCcolab('TSmed')));
% pl_gids = PL.Table(:, PLcolab('QL')) > ProbThr &...
%     abs(PL.Table(:, PLcolab('CS'))) < StdThr*nanstd(PL.Table(:, PLcolab('CS'))) &...
%     abs(PL.Table(:, PLcolab('TS'))) < StdThr*nanstd(PL.Table(:, PLcolab('TS')));
%% verify transcient
figure;
subplot(1, 2, 1);
scatter(BC.CelTable(bc_gids, BCcolab('RawAmpMeanavg')), BC.CelTable(bc_gids, BCcolab('TSavg')),5, 'k', 'filled');
subplot(1, 2, 2);
scatter(PL.Table(pl_gids, PLcolab('Rawavg')),PL.Table(pl_gids, PLcolab('TS')),  5, 'k', 'filled');

%%
% BCTypes = [5 57 58 6 7 89];
% BCTypelabels = {'BC5', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
% BCTypes = [5 50 51 57 58 6 7 89];
% BCTypelabels = {'BC5', 'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
BCTypes = [50 51 57 58 6 7 89];
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
nType = length(BCTypes);
BCb = BC.CelTable(bc_gids, :);
BCCTpr = BC.CelTprRespf(bc_gids, :);
BCCResp = BC.CelRespf(bc_gids, :);
BCCFet = BC.CelClusterFeature(bc_gids, :);
PLb = PL.Table(pl_gids, :);
PLCTpr = PL.TprRespFunc(pl_gids, :);
PLCResp = PL.RespFunc(pl_gids, :);
PLCFet = PL.ClusterFeatures(pl_gids, :);
uniBCs = unique(BCb(:, 1:3), 'rows');
uniPLs = unique(PLb(:, 1:3), 'rows');
nBC = length(uniBCs);
nPL = length(uniPLs);
BCvar = nan(nBC, 4);
PLvar = nan(nPL, 4);
BCtype = nan(nBC, 1);
Depth = nan(nPL, 1);
uniDay = unique(BCb(:, 1));
BCday = nan(nBC, 1);
BCTpr = nan(nBC, size(BC.CelTprRespf, 2));
BCResp = nan(nBC, size(BC.CelRespf, 2));
BCFet = nan(nBC, size(BC.CelClusterFeature, 2));
for i = 1:nBC
    cids = find(BCb(:, 1) == uniBCs(i, 1) & BCb(:, 2) == uniBCs(i, 2) & BCb(:, 3) == uniBCs(i, 3));
    BCvar(i, :) = [median(BCb(cids, [BCcolab('CSavg') BCcolab('TSavg')]), 1),...
        std(BCb(cids, [BCcolab('CSmed') BCcolab('TSmed')]), [], 1)];
    BCtype(i) = BCb(cids(1), BCcolab('Type'));
    BCday(i) = find(uniBCs(i, 1) == uniDay);
    BCTpr(i, :) = median(BCCTpr(cids, :), 1, 'omitnan');
    BCResp(i, :) = median(BCCResp(cids, :), 1, 'omitnan');
    BCFet(i, :) = mean(BCCFet(cids, :), 1, 'omitnan'); % maybe mean
end
uniDay = unique(PLb(:, 1));
PLday = nan(nPL, 1);
PLTpr = nan(nPL, size(PL.TprRespFunc, 2));
PLResp = nan(nPL, size(PL.RespFunc, 2));
PLFet = nan(nPL, size(PL.ClusterFeatures, 2));
for i = 1:nPL
    cids = find(PLb(:, 1) == uniPLs(i, 1) & PLb(:, 2) == uniPLs(i, 2)  & PLb(:, 3) == uniPLs(i, 3));
    PLvar(i, :) = [median(PLb(cids, [PLcolab('CS') PLcolab('TS')]), 1),...
        std(PLb(cids, [PLcolab('CS') PLcolab('TS')]), [], 1)];
    cid = DocNum(:,ColName2Ind('Day')) == uniPLs(i, 1) & DocNum(:,ColName2Ind('Cel')) == uniPLs(i, 2) &...
        DocNum(:,ColName2Ind('Exp')) == uniPLs(i, 3);
    if sum(cid) == 0
        continue
    end
    Depth(i) = DocNum(cid,ColName2Ind('Depth'));
    PLday(i) = find(uniPLs(i, 1) == uniDay);
    PLTpr(i, :) = median(PLCTpr(cids, :), 1, 'omitnan');
    PLResp(i, :) = median(PLCResp(cids, :), 1, 'omitnan');
    PLFet(i, :) = mean(PLCFet(cids, :), 1, 'omitnan');
end
bc_gids = ~isnan(BCvar(:, 1));
pl_gids = ~isnan(PLvar(:, 1));
% bc_gids = BCvar(:, 1) > 0 & BCvar(:, 2) > 0;
% pl_gids = PLvar(:, 1) > 0 & PLvar(:, 2) > 0;
BCvar = BCvar(bc_gids, :);
BCtype = BCtype(bc_gids, :);
BCday = BCday(bc_gids, :);
BCTpr = BCTpr(bc_gids, :);
BCResp = BCResp(bc_gids, :);
BCFet = BCFet(bc_gids, :);
uniBCs = uniBCs(bc_gids, :);

PLvar = PLvar(pl_gids, :);
Depth = Depth(pl_gids, :);
PLday = PLday(pl_gids, :);
PLTpr = PLTpr(pl_gids, :);
PLResp = PLResp(pl_gids, :);
PLFet = PLFet(pl_gids, :);
uniPLs = uniPLs(pl_gids, :);

nBC = size(BCvar, 1);
nPL = size(PLvar, 1);
%%
DepthBoundary = [0.55, 0.63, 0.7, 0.8];
DepthGId = nan(size(Depth));
for i = 1:length(DepthBoundary)
    if i == 1
        DepthGId(Depth < DepthBoundary(i)) = i;
        DepthGId(Depth >= DepthBoundary(i) & Depth < DepthBoundary(i+1)) = i+1;
    elseif i == length(DepthBoundary)
        DepthGId(Depth >= DepthBoundary(i)) = i+1;
    else
        DepthGId(Depth >= DepthBoundary(i) & Depth < DepthBoundary(i+1)) = i+1;
    end
end
figure; histogram(DepthGId);
%% Demostrate Response Features of each cell
close all; figure;
TargetMatrix = BCFet;
TargetMatrix = TargetMatrix./max(TargetMatrix, [], 2);
MaxV = max(TargetMatrix(:));
MinV = min(TargetMatrix(:));
GapV = (MinV-0.01*(MaxV-MinV));
nDmt = length(Dmt);
nL = 10;
ytic = 1:7;
yticlab = num2cell(Dmt);
ylab = 'Radius (um)';
nTime = size(TargetMatrix, 2)/nDmt;
rcount = 1;
SliceResp_time = [];
SliceResp_dmt = [];
SliceResp_diag = [];
SliceResp_type = [];
SliceResp_type = [];
SliceResp_full = [];
SliceResp_info = [];
DiagSlice = round(linspace(4, 21, 7)'+(1:5)) + (0:30:180)';
             
for i = 1:nType
    acids = find(BCtype == BCTypes(i));
    nloop = ceil(numel(acids)/nL);
    for k = 1:nloop
        subplot(9, 1, rcount);
        rcount = rcount + 1;
        if k == nloop
            cids = acids(1:end);
        else
            cids = acids(1:nL);
            acids(1:nL) = [];
        end
        nCell = length(cids);
        Gap = 1;
        canva = [];
        clear Mids roilabels
        for j = 1:nL
            if j > nCell
                canva = [canva; GapV*ones(nTime, nDmt); GapV*ones(Gap, nDmt)];
            else
                if j == 1
                    Mids(j, :) = [nTime/2, nTime];
                else
                    Mids(j, :) = [Mids(j-1, 1)+0.5*Mids(j-1, 2)+nTime/2+Gap, nTime];
                end
                ctile = reshape(TargetMatrix(cids(j), :), [], nDmt);
                canva = [canva; ctile; GapV*ones(Gap, nDmt)];
                roilabels{j} = sprintf('%d-%d', uniBCs(cids(j), 1), uniBCs(cids(j), 2));
                
                SliceResp_time = [SliceResp_time; mean(ctile(:, 1:4), 2)'];
                SliceResp_dmt = [SliceResp_dmt; mean(ctile(1:17, :), 1)];
                SliceResp_type = [SliceResp_type; BCTypes(i)];
                SliceResp_full = [SliceResp_full; reshape(ctile, 1, [])];
                SliceResp_info = [SliceResp_info ; uniBCs(cids(j), :)];
                SliceResp_diag = [SliceResp_diag; mean(ctile(DiagSlice), 2)'];
            end
        end
        cmap = [ones(3, 3); jet(300); ];
        imagesc(canva', [GapV MaxV]); colormap(cmap);colorbar;
        box off
        xticks(round(Mids(:, 1)));
        xticklabels(roilabels);
        yticks(ytic);
        yticklabels(yticlab);
        ylabel(ylab);
        title(BCTypelabels{i});
    end
end
keyboard;

%%
SaveDay = '40423';
SaveFileName = ['./Results/IdentifiedBCtype_Response' SaveDay '.mat'];
save(SaveFileName, 'SliceResp_info', 'SliceResp_full', 'SliceResp_type');
