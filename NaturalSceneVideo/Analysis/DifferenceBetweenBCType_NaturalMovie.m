%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
nIter = 5000;
%%

IsEncodingSpace = 1;
if IsEncodingSpace
    cResp = Y(:, 1:3);
else
    cResp = Resp;
end
TypeDiff_val = nan(nType, nType);
TypeDiff_p = nan(nType, nType);
rng('shuffle');
for i = 1:nType
    for j = 1:nType
        if j>=i
            aids = find(typ == i);
            bids = find(typ == j);
            na = length(aids);
            nb = length(bids);
            ra = mean(cResp(aids, :), 1, 'omitnan');
            rb = mean(cResp(bids, :), 1, 'omitnan');
            if IsEncodingSpace
                diff_true = sqrt(sum((ra-rb).^2));
            else
                diff_true = 1-median(corrVecChop(ra,rb, clipids));
            end
            TypeDiff_val(i, j) = diff_true;
            diff_sample = nan(nIter, 1);
            cids = [aids; bids];
            for k = 1:nIter
                sid_i = randsample(cids, na);
                sid_j = cids(~ismember(cids, sid_i));
                ra = mean(cResp(sid_i, :), 1, 'omitnan');
                rb = mean(cResp(sid_j, :), 1, 'omitnan');
                if IsEncodingSpace
                    diff_sample(k) = sqrt(sum((ra-rb).^2));
                else
                    diff_sample(k) = 1-median(corrVecChop(ra,rb, clipids));
                end
            end
            TypeDiff_p(i, j) = 1-mean(diff_sample < diff_true);
            clc
            fprintf('progress... %d/%d, %d/%d (%d, %d) \n', i, nType, j, nType, na, nb);
        else
            TypeDiff_val(i, j) = TypeDiff_val(j, i);
            TypeDiff_p(i, j) = TypeDiff_p(j, i);
        end
        
    end
end
%%
cTypeDiff_p = TypeDiff_p;
cTypeDiff_val = TypeDiff_val;
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource');
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
cTypeDiff_p(eye(size(cTypeDiff_p, 1))==1) = nan;
[~, ~, ~, pcorr] = fdr_bh(cTypeDiff_p(:));
pcorr = reshape(pcorr, size(cTypeDiff_p));
Targetp = pcorr;
% close all
figure;
subplot(1, 3, 1);hold on

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
subplot(1, 3, 2);
imagesc(log10(Targetp)); colorbar
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title('P value (log10)');
subplot(1, 3, 3);
imagesc(Targetp<0.05); colorbar
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title('Signficant (p <0.05)');
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_NatMovBCTypeResponsesDifference', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
% SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BipolarCellTerminal\ResultFigures\';
% FleNam = sprintf('%sFigure5_NatMovBCTypeEncodingSpaceDifference', SaveFolder);
% print('-depsc','-painters','-loose', '-r300',FleNam)
% saveas(gcf,[FleNam '.png']);
%%
map = colormap(gray);
v = TypeDiff_val;
v(eye(size(v, 1))==1) = 0;
minv = min(v(:));
maxv = max(v(:));
ncol = size(map,1);
s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
rgb_image = imresize(ind2rgb(s,map), 50, 'nearest');
rgb_image_mask = triu(rgb_image(:, :, 1), 1);

divd = sum(rgb_image, 3);
rgb_image_mask(rgb_image_mask>0 & imresize(Targetp, 50, 'nearest')<0.05) = ...
    rgb_image_mask(rgb_image_mask>0 & imresize(Targetp, 50, 'nearest')<0.05) *1.7;
rgb_image_mask = tril(rgb_image(:, :, 1), 1)+rgb_image_mask;
rgb_image(:, :, 1) = rgb_image_mask;

rgb_image = rgb_image.*divd./sum(rgb_image, 3);
close all
% figure; imagesc(rgb_image_mask);
figure; hold on;
imshow(rgb_image); 
axis on
xticks(((1:nType)-0.5)*50);
xticklabels(BCTypelabels);
yticks(((1:nType)-0.5)*50);
yticklabels(BCTypelabels);
xlim([0 350]);
ylim([0 350]);
%% direct distance
nIter = 2000;
cPairdist = Pairdist;
cPairdist(eye(size(cPairdist, 1))==1) = nan;
TypeDiff_val = nan(nType, nType);
TypeDiff_p = nan(nType, nType);
rng('shuffle');
for i = 1:nType
    for j = 1:nType
        if j>i
            aids = find(typ == i);
            bids = find(typ == j);
            na = length(aids);
            nb = length(bids);
            diff_true = median(cPairdist(aids, bids), 'all', 'omitnan');
            TypeDiff_val(i, j) = diff_true;
            diff_sample = nan(nIter, 1);
            cids = [aids; bids];
            for k = 1:nIter
                sid_i = randsample(cids, na);
                sid_j = cids(~ismember(cids, sid_i));
                diff_sample(k) = median(cPairdist(sid_i, sid_j), 'all', 'omitnan');
            end
            TypeDiff_p(i, j) = 1-mean(diff_sample < diff_true);
        else
            TypeDiff_val(i, j) = TypeDiff_val(j, i);
            TypeDiff_p(i, j) = TypeDiff_p(j, i);
        end
        clc
        fprintf('progress... %d/%d, %d/%d (%d, %d) \n', i, nType, j, nType, na, nb);
    end
end
%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource');
BCTypelabels = {'BC5o', 'BC5i', 'BC5t', 'XBC', 'BC6', 'BC7', 'BC89'};
[~, ~, ~, pcorr] = fdr_bh(TypeDiff_p(:));
pcorr = reshape(pcorr, size(TypeDiff_p));
Targetp = pcorr;
figure;
subplot(1, 3, 1);
imagesc(TypeDiff_val); colorbar
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title('Difference');
subplot(1, 3, 2);
imagesc(log10(Targetp)); colorbar
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title('P value (log10)');
subplot(1, 3, 3);
imagesc(Targetp<0.05); colorbar
xticks(1:nType);
xticklabels(BCTypelabels);
yticks(1:nType);
yticklabels(BCTypelabels);
title('Signficant (p <0.05)');