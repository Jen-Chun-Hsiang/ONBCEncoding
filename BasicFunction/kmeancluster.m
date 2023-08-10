function [Idx, minid] = kmeancluster(X, distType)
if nargin == 1
    distType = 'sqeuclidean';
end
nIter = 20;
nX = size(X, 1);
nd = size(X, 2);
maxk = round(sqrt(nX)*10);
maxk = min([maxk 200]);
kranges = 1:1:maxk;
for k = 1:length(kranges)
    [idx, C, sumd] = kmeans(X, kranges(k), 'Distance', distType, 'Replicates', nIter); %sqeuclidean
    clear rss
    for h = 1:kranges(k)
        rss(h) = sum(sqrt(sum((X(idx==h, :)-C(h, :)).^2, 2)));
    end
    BIC(k) = nX*log(sum(sumd)/nX)+kranges(k)*log(nX);
    
end
figure; plot(BIC, 'k');hold on
plot(BIC, 'm');
[~, minid] = min(BIC);
minid = kranges(minid);
Idx = kmeans(X, minid, 'Distance', distType, 'Replicates', nIter);
