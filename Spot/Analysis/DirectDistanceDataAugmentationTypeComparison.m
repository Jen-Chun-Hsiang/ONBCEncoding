% manually load file 'BCtype_augdata_aligned_33023.mat'
clc; %close all;
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
DistanceType = 'correlation';%'correlation''euclidean'
X = Dfull;
nP = size(X, 1);
BCTypes = [50 51 57 58 6 7 89];
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
nType = length(BCTypes);
%%
Pairdist = squareform(pdist(X, DistanceType));
Pairdist(eye(nP)== 1) = nan;
pd = fitdist(min(Pairdist, [], 2),'gamma');
distprob = 1-cdf(pd, Pairdist);
%%
getIPLProb; % a script to get IPLd and DepthTypeProb
% keyboard;
depthprob = IPL2depthprob(IPLd);
% depthprob = depthprob./sum(depthprob, 2, 'omitnan');
%% ground truth (leave one out assessment)
SliceResp_type = Dtype(Dtype ~= 81);
nBC = size(SliceResp_type, 1);
[~, dtypes] = ismember(SliceResp_type, BCTypes);
PlexusTypeGround = nan(nBC, 3);
for i = 1:nBC
    depthweight =  depthprob(dtypes(i),  dtypes);
    [PlexusTypeGround(i, 1), PlexusTypeGround(i, 2)] = max(distprob(i, 1:nBC).*depthweight,...
        [], 'omitnan');
end
PlexusTypeGround(:, 3) = dtypes(PlexusTypeGround(:, 2));
PlexusTypeGround(:, 4) = dtypes;
%% classify the cell types
PlexusTypeD = nan(nP-nBC, 3);
for i = 1:(nP-nBC)
    [PlexusTypeD(i, 1), PlexusTypeD(i, 2)] = max(distprob(i+nBC, 1:nBC).*DepthTypeProb(i, dtypes'));
end
PlexusTypeD(:, 3) = dtypes(PlexusTypeD(:, 2));


%% remove outlier
c = [PlexusTypeGround(:, 1); PlexusTypeD(:, 1)];
figure; 
subplot(1, 3, 1);
histogram(c); 
pd = fitdist(c,'gamma');
subplot(1, 3, 2);
x = 0.01:0.01:0.5;
y = pdf(pd, x);
y = y./sum(y, 'omitnan');
plot(x, y);
subplot(1, 3, 3);hold on
y = cdf(pd, x);
% y = y./sum(y, 'omitnan');
plot(x, y);
maxx = find(y>1-0.05);
maxx = x(maxx(1)-1);
plot(maxx*ones(1, 2), [0 1], '--k');
remids = find(PlexusTypeD(:, 1) <= maxx);

%%
pid_types = PlexusTypeD(:, 3);
pid_types(PlexusTypeD(:, 1) > maxx) = 0;
cids = PLTab(:, 7) > 0;
PLTab(cids, 8) = pid_types(PLTab(cids, 7));


%%
SaveDay = '030723';
SaveFileName = ['//storage1.ris.wustl.edu/kerschensteinerd/Active/Emily/NaturalSceneVideo/Analysis/Results/BCtype_augdata_aligned_assigned_' SaveDay '.mat'];
save(SaveFileName, 'PLTab', 'Dfull', 'Dtype','Ddepth', 'Dtab');

%%
%% Correlation based difference
nIter = 10000;
fea = Dfull([(1:nBC)'; remids(:)+nBC], :);
typ = [PlexusTypeGround(:, 4); PlexusTypeD(remids, 3)];
TypeDiff_val = nan(nType, nType);
TypeDiff_p = nan(nType, nType);
TypeDiff_map_p = nan(1, size(fea, 2), nType*nType);
DiffMap_ind = nan(nType*nType, 2);
BCMap = nan(nType, size(fea, 2));
BCMap_orig = nan(nType, size(fea, 2));
rng('shuffle');
Count = 1;
for i = 1:nType
    for j = 1:nType
        if j>=i
            cidis = find(typ==i);
            cidjs = find(typ==j);
            ncidi = length(cidis);
            ncidj = length(cidjs);
            if i == j
                BCMap(i, :) = mean(fea(cidis, :), 1);
                BCMap_orig(i, :) = mean(fea(Dtype==BCTypes(i), :), 1);
            end
%             diff_true = meanR(corrcoef(mean(fea(cidis, :), 1)', mean(fea(cidjs, :), 1)').^2);
            diff_true = pdist([mean(fea(cidis, :), 1); mean(fea(cidjs, :), 1)], 'correlation');
            BCMap_diff_true = mean(fea(cidis, :), 1)-mean(fea(cidjs, :), 1);
            BCMap_diff_sample = nan(nIter, size(BCMap_diff_true, 2));
            TypeDiff_val(i, j) = diff_true;
            diff_sample = nan(nIter, 1);
            for k = 1:nIter
                sid_i = randsample([cidis; cidjs], ncidi, true);
                sid_j = randsample([cidis; cidjs], ncidj, true);
%                 diff_sample(k) = meanR(corrcoef(mean(fea(sid_i, :), 1)', mean(fea(sid_j, :), 1)').^2);
                diff_sample(k) = pdist([mean(fea(sid_i, :), 1); mean(fea(sid_j, :), 1)], 'correlation');
                BCMap_diff_sample(k, :) = mean(fea(sid_i, :), 1)-mean(fea(sid_j, :), 1);
            end
            
            TypeDiff_p(i, j) = mean(diff_sample > diff_true);
            TypeDiff_map_p(:, :, Count) = 1-mean(BCMap_diff_sample.^2 < BCMap_diff_true.^2, 1);
%             if j ~= i
%                 keyboard;
%             end
            DiffMap_ind(Count, :) = [i, j];
            Count = Count + 1;
        else
            TypeDiff_val(i, j) = TypeDiff_val(j, i);
            TypeDiff_p(i, j) = TypeDiff_p(j, i);
            TypeDiff_map_p(:, :, Count) = TypeDiff_map_p(:, :, DiffMap_ind(:, 1)==j & DiffMap_ind(:, 2)==i);
            DiffMap_ind(Count, :) = [i, j];
            Count = Count + 1;
        end
        clc
        fprintf('progress... %d/%d, %d/%d (%d, %d) \n', i, nType, j, nType, ncidi, ncidj);
    end
end

%%
figure; 
subplot(1, 2, 1);
imagesc(BCMap_orig);
subplot(1, 2, 2);
imagesc(BCMap);