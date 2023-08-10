function [OutColor, OutIds, v] = GroupScatterColor(InVal, ColorMap, setrange)
switch nargin
    case 1
        nGradient = 256;
        Colors = jet(nGradient);
        setrange = [];
    case 2
        Colors = ColorMap;
        nGradient = size(ColorMap, 1);
        setrange = [];
    case 3
        Colors = ColorMap;
        nGradient = size(ColorMap, 1);
end


% MaxV = max(InVal(cids));
% MinV = min(InVal(cids));
cids = find(~isnan(InVal));
if ~isempty(setrange)
    MaxV = max(setrange);
    MinV = min(setrange);
else
    MaxV = quantile(InVal(cids), 0.9);
    MinV = quantile(InVal(cids), 0.1);
end
v = [MinV MaxV];
[~, sorted] = sort(InVal(cids));
cids = cids(sorted);
OutIds = cids;
bd = linspace(MinV, MaxV, nGradient);
vIds = nan(length(cids), 1);
for i = 1:length(cids)
    [~, sid] = min((InVal(cids(i))-bd).^2);
    vIds(i) = sid;
end

OutColor = Colors(vIds, :);
end