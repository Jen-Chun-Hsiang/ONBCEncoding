function [sigO, Rsquare, failflag] = AssembleClipbyQuality(sigI, clipids, method)
if nargin == 2
    method = 'top';
end

assert(size(sigI, 2) == length(clipids))
clips = unique(clipids);
nclip = length(clips);
nrepeat = size(sigI, 1);
if nrepeat == 1 
    sigO = sigI;
    Rsquare = nan;
    return
end
if nrepeat >= 14
    rng('shuffle');
end
maxSample = 1000;
topr = nan(nclip, 1);
Rsquare = nan(1, nclip);
sigO = nan(1, size(sigI, 2));
failflag = 0;
for i = 1:nclip
    % calculate repeat reliability
    cids = clipids == clips(i);
    c = sigI(:, cids);
    c = c+randn(size(c))*1e-6;
    r = corr(c').^2;
    r(eye(size(r, 1))==1) = nan;
    switch lower(method)
        case 'outmedian'
            mv = median(r, 2, 'omitnan');
            rids = ~isnan(mv);
            nrid = sum(rids);
            if nrid == 2
                sigO(cids) = median(c, 1);
                Rsquare(i) = medianR(corr(c').^2);
            else
                rids = ~isnan(mv) & ~isoutlier(mv, 'quartiles');
                nrid = sum(rids);
                c = c(rids, :);
                sigO(cids) = median(c, 1, 'omitnan');
                % get repeat reliability
                fullset = 1:nrid;
                npair = nchoosek(nrid, round(nrid*0.5));                
                if maxSample < npair
                    r = nan(maxSample, 1);
                    for j = 1:maxSample
                        aids = randsample(fullset, round(nrid*0.5));
                        bids = fullset(~ismember(fullset, aids));
                        r(j) = corr(median(c(aids, :), 1)',median(c(bids, :), 1)');
                    end
                else
                    pairs = nchoosek(fullset, round(nrid*0.5));
                    npair = size(pairs, 1);
                    r = nan(npair, 1);
                    for j = 1:npair
                        aids = pairs(j, :);
                        bids = fullset(~ismember(fullset, aids));
                        r(j) = corr(median(c(aids, :), 1)',median(c(bids, :), 1)');
                    end
                end
                Rsquare(i) = median(r.^2, 'omitnan');
            end
        case 'top'
            [topr(i), pos] = max(r.^2, [], 1:2, 'omitnan', 'linear');
            [x, y] = ind2sub(size(r), pos);
            % varify
            c = sigI([x y], cids);
            c = c+randn(size(c))*1e-6;
            sigO(cids) = mean(c, 1);
            Rsquare(i) = meanR(corr(c').^2);
            
            if abs(Rsquare(i)-topr(i))>1e-6;
                failflag = 1;
            end
    end
end