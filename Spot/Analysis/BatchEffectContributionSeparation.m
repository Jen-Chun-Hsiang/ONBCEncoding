% test if distribution of two sets are different from randomly sample from
% the sets
nIter = 2000;
nSample = 10;
ndim = 2;
kldmethod = 4; % kldmethod 2: only avaiable for ndim = 2
Zt = nan(nIter, 1);
Zr = nan(nIter, nSample);
p = nan(nIter, nSample);
fullset = unique(tabd(:,1));
fullset(isnan(fullset)) = [];
sets = nchoosek(fullset, 8);
rng('shuffle');
setids = randperm(size(sets, 1));
setids = setids(1:nIter);
for i = 1:nIter
    clc
    fprintf('progress...%d/%d \n', i, nIter);
    aids = find(ismember(tabd(:,1), sets(setids(i), :)));
    oppsets = find(~ismember(fullset, sets(setids(i), :)));
    bids = find(ismember(tabd(:,1), oppsets));
    A = Y(aids+nBC, 1:ndim);
    B = Y(bids+nBC, 1:ndim);
    nA = length(aids);
%     Zt(i) = KLdivergence(A, B+(randn(size(B))-0.5)*1e-6, kldmethod);%
    D = [A; B];
    Zt(i) = KLdivergence(A, B, kldmethod);%
    
    nD = size(D, 1);
    for j = 1:nSample
        ids = randperm(nD);
        cB = D(ids(nA+1:end), :);
        %         Zr(i, j) = KLdivergence(D(ids(1:nA), :), cB+(randn(size(cB))-0.5)*1e-6, kldmethod);%
        Zr(i, j) = KLdivergence(D(ids(1:nA), :), cB, kldmethod);%
        if Zr(i, j) == Inf
            cA = D(ids(1:nA),:);
            KLdivergence(cA, cB+randn(size(cB))*1e-6, kldmethod);
            figure; scatter(cA(:, 1), cA(:, 2), 15, 'k', 'filled'); hold on
            scatter(cB(:, 1), cB(:, 2), 15, 'r', 'filled');
            keyboard;
        end
    end
    p(i, :) = Zt(i) < Zr(i, :);
end
% find the median set to represent the data
Zt_set = sets(setids, :);
[~, mids] = min(abs(Zt-median(Zt)));
mset = Zt_set(mids, :);
%%
figure; hold on
% h1 = histogram(Zr-Zt, linspace(-0.4, 0.4, 30));
h1 = histogram(Zr-Zt, linspace(-0.14, 0.1, 24));
h1.Normalization = 'probability';
h1.FaceColor = [60 120 232]/255;
h1.EdgeColor = ones(1, 3);
plot(median(Zr-Zt, 1:2)*ones(1, 2), [0 0.2], '--k');
ylabel('Probability');
xlabel('Jensen-Shannon divergence (suffle-set)');
xticks(-0.1:0.1:0.1);
xticklabels({'-0.1', '0', '0.1'});
yticks(0:0.1:0.2);
yticklabels({'0', '0.1', '0.2'});
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_BatchEffectKL', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam);
saveas(gcf,[FleNam '.png']);
%% show two example sets
ndim = 2;
kldmethod = 4; 
aids = find(ismember(tabd(:,1), [2 8 9 10 11 12 13 15]));
oppsets = find(~ismember(fullset,[2 8 9 10 11 12 13 15]));
bids = find(ismember(tabd(:,1), oppsets));
A = Y(aids+nBC, 1:ndim);
B = Y(bids+nBC, 1:ndim);
nA = length(aids);
D = [A; B];
z = KLdivergence(A, B, kldmethod);
%%
figure; hold on
pIds = find(ismember(tabd(:,1), [2 8 9 10 11 12 13 15]));
scatter3(Y(pIds, 1), Y(pIds, 2), Y(pIds, 3), 12, [91, 155, 213]/255, 'filled');
pIds = find(ismember(tabd(:,1), [1 3 4 5 6 7 14 16]));
scatter3(Y(pIds, 1), Y(pIds, 2), Y(pIds, 3), 12, [112, 179, 71]/255, 'filled');
xlim([min(Y(:, 1)) max(Y(:, 1))]);
ylim([min(Y(:, 2)) max(Y(:, 2))]);
ylabel('MDS component 2');
xlabel('MDS component 1');
title('Random split sets');
yticks(-0.4:0.4:0.8);
yticklabels({'-0.4', '0', '0.4', '0.8'});
xticks(-0.6:0.3:0.6);
xticklabels({'-0.6', '', '0', '', '0.6'});
%%
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
FleNam = sprintf('%sFigure4_SpotProjection_BatchEffectSplitDataRepresentation', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
%% Evaluate if the functional distance is larger within batch (between type)
% find the number of bc types
% (1) perform type paired distance measurement
% (2) perform within type separation estimation (ROI > 2)
% sample 
nIter = 1000;
ndim = 3;
sumdistmethod = 2;
typed = PlexusTypeD(remids, 3);
assert(length(typed) == size(tabd, 1));
p = [];
Zs = [];
Za = [];
for i = 1:16
    cids = tabd(:,1) == i;
    ctypes = unique(typed(cids));
    nctype = length(ctypes);
    if nctype < 2
        continue
    end
    ctypeperm = nchoosek(ctypes, 2);
    nperm = size(ctypeperm, 1);
    cp = nan(nperm, nIter);
    cZs = nan(nperm, 1);
    cZa = nan(nperm, nIter);
    for j = 1:nperm
        clc
        fprintf('progress...%d/%d, %d/%d \n', i, 16, j, nperm);
        aids = find(tabd(:,1) == i & typed == ctypeperm(j, 1));
        bids = find(tabd(:,1) == i & typed == ctypeperm(j, 2));
        A = Y(aids+nBC, 1:ndim);
        B = Y(bids+nBC, 1:ndim);
        nA = length(aids);
        nB = length(bids);
        cZs(j) = pairgroupdistance(A, B, sumdistmethod);
        for k = 1:nIter
            aids = find(typ == ctypeperm(j, 1));
            bids = find(typ == ctypeperm(j, 2));
            aids = randsample(aids, nA);
            bids = randsample(bids, nB);
            A = Y(aids, 1:ndim);
            B = Y(bids, 1:ndim);
            cZa(j, k) = pairgroupdistance(A, B, sumdistmethod);
        end
        cp(j, :) = cZs(j) < cZa(j, :);
    end
    Zs = [Zs; cZs i*ones(nperm, 1)];
    Za = [Za; cZa];
    p = [p; cp];
end
%%
figure; hold on
h1 = histogram(Za-Zs(:, 1), -1:0.1:1.2);
h1.Normalization = 'probability';
h1.FaceColor = [60 120 232]/255;
h1.EdgeColor = ones(1, 3);
% h2 = histogram(Za, 0:0.1:1.8);
% h2.Normalization = 'probability';
% h2.FaceColor = [232 120 60]/255;
% h2.EdgeColor = ones(1, 3);
plot(median(Za-Zs(:, 1), 1:2)*ones(1, 2), [0 0.18], '--k')
