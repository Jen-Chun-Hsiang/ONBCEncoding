%
Resps = [];
Typs = [];
for i = 1:7
    for j = 1:5
        ids = gtyp == i & gsiz == 1 & gfrq == j;
        if j == 1
            ncel = sum(ids);
            cresp = nan(ncel, 5);
        end
        cresp(:, j) = y(ids);
    end
    Resps = [Resps; cresp];
    Typs = [Typs; i*ones(ncel, 1)];
end

keepids = isnan(Resps(:, 1));
Typs(keepids) = [];
Resps(keepids,:) = [];

%%
nIter = 2000;
nType = 7;
Frqwin = 3:5;
nFrq = length(Frqwin);
BCMap = nan(nType, nFrq);
TypeDiff_val = nan(nType, nType);
TypeDiff_p = nan(nType, nType);
DiffMap_ind = nan(nType*nType, 2);
rng('shuffle')
Count = 1;

for i = 1:nType
    for j = 1:nType
        if j >=1
            cidis = find(Typs==i);
            cidjs = find(Typs==j);
            ncidi = length(cidis);
            ncidj = length(cidjs);
            if i == j
                BCMap(i, :) = mean(Resps(cidis, Frqwin), 1);
            end
            diff_true = sqrt(sum((mean(Resps(cidis, Frqwin), 1) - mean(Resps(cidjs, Frqwin), 1)).^2));
            BCMap_diff_true = diff_true;
            TypeDiff_val(i, j) = diff_true;
            diff_sample = nan(nIter, 1);
            for k = 1:nIter
                sid_i = randsample([cidis; cidjs], ncidi, true);
                sid_j = randsample([cidis; cidjs], ncidj, true);
                diff_sample(k) = sqrt(sum((mean(Resps(sid_i, Frqwin), 1) - mean(Resps(sid_j, Frqwin), 1)).^2));
            end
            TypeDiff_p(i, j) = mean(diff_sample > diff_true);
            DiffMap_ind(Count, :) = [i, j];
            Count = Count + 1;
        else
            TypeDiff_val(i, j) = TypeDiff_val(j, i);
            TypeDiff_p(i, j) = TypeDiff_p(j, i);
            DiffMap_ind(Count, :) = [i, j];
            Count = Count + 1;
        end
        clc
        fprintf('progress... %d/%d, %d/%d (%d, %d) \n', i, nType, j, nType, ncidi, ncidj);
    end
end
%%
cTypeDiff_p = TypeDiff_p;
cTypeDiff_val = TypeDiff_val;
cTypeDiff_p(eye(size(cTypeDiff_p, 1))==1) = nan;
[~, ~, ~, pcorr] = fdr_bh(TypeDiff_p(:));
pcorr = reshape(pcorr, size(TypeDiff_p));
Targetp = pcorr;
figure;
subplot(1, 2, 1); hold on

imagesc(flipud(cTypeDiff_val)); colorbar
pv = nan(size(Targetp));
pv(Targetp<0.05) = 1;
pv(Targetp<0.01) = 2;
pv(Targetp<0.001) = 3;
for i = 1:size(Targetp, 1)
    for j = 1:size(Targetp, 2)
        if ~isnan(pv(i, j)) & j<i
            switch pv(i, j)
                case 1
                    signi = '*';
                case 2
                    signi = '**';
                case 3
                    signi = '***';
            end
            text(i, 8-j-0.1, signi, 'Color', 'r', 'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'middle', 'FontSize', 15);
        end
    end
end
xlim([0.5 7.5]);
ylim([0.5 7.5]);
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(fliplr(BCTypelabels));
colormap(gray);
title('Difference');

subplot(1, 2, 2);
ytic = 1:7;
enlargefac = 10;
c = imresize(log10(Targetp), enlargefac, 'nearest');
imagesc(c'); colorbar; hold on
visboundaries(c' < log10(0.05));
% box off
xticks((ytic-0.5)*enlargefac);
xticklabels(BCTypelabels);
yticks((ytic-0.5)*enlargefac);
yticklabels(BCTypelabels);
title('P value (log10)');

%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure2_FrequencyResponseTypeDifferenceBoostrapping', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);

%%
% save('BCtype_BlockFreq_Boostrapping.mat', 'TypeDiff_p', 'DiffMap_ind', 'TypeDiff_val', 'BCMap',...
%     'Frqwin', 'nIter');