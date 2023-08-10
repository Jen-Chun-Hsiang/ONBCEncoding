function [ROImatchBinA, ROImatchAinB] = ROIMapMatching(Amap, Bmap)
    % use anchor points to map ROI from 2 different map
    % For each ROI in Bmap is caclulate with the largest overlapping ROI in
    % Amap, overlappying is excluding nan area
    % [BId, AId, CoveredA(%), CoveredB(%)]
assert(sum(size(Amap)-size(Bmap)) <eps);
aids = unique(Amap(:));
bids = unique(Bmap(:));
aids(isnan(aids)) = [];
bids(isnan(bids)) = [];
nA = length(aids);
nB = length(bids);
ROImatchBinA = nan(nB, nA);
ROImatchAinB = nan(nB, nA);
Amap = Amap(:);
Bmap = Bmap(:);
for i = 1:nB
    for j = 1:nA
        cb = Bmap == bids(i);
        ca = Amap == aids(j);
        ROImatchBinA(i,j) = sum(cb & ca)./sum(ca);
        ROImatchAinB(i,j) = sum(cb & ca)./sum(cb);
    end
end

