function rs = corrVecChop(a, b, clipids)
clips = unique(clipids);
nclip = length(clips);
rs = nan(nclip, 1);
clipids = clipids(:);
a = a(:);
b = b(:);
for i = 1:nclip
    rs(i) = corr(a(clipids == i), b(clipids == i));
end