function KLd = KLdivergence(x, y, method)
% Compute the Kullback-Leibler divergence between two multivariate samples.
%   Parameters
%   ----------
%   x : 2D array (n,d)
%     Samples from distribution P, which typically represents the true
%     distribution.
%   y : 2D array (m,d)
%     Samples from distribution Q, which typically represents the approximate
%     distribution.
%   Returns
%   -------
%   out : float
%     The estimated Kullback-Leibler divergence D(P||Q).
%   References
%   ----------
%   Pérez-Cruz, F. Kullback-Leibler divergence estimation of
% continuous distributions IEEE International Symposium on Information
% Theory, 2008.
if nargin == 2
    method = 4; %JSD
end
switch method
    case 1
        x = unique(x, 'rows');
        x(any(isnan(x), 2), :) = [];
        y = unique(y, 'rows');
        y(any(isnan(y), 2), :) = [];
        xtree = KDTreeSearcher(x);
        ytree = KDTreeSearcher(y);
        n = size(x, 1);
        m = size(y, 1);
        d = size(x, 2);
        assert(d == size(y, 2));
        [~, r] = knnsearch(xtree,x,'K',2);
        r = r(:, 2);
        [~, s] = knnsearch(ytree,x,'K',1);
        KLd = d*sum(-log(r./s))/n + log(m/(n-1));
    case 3
        x = unique(x, 'rows');
        x(any(isnan(x), 2), :) = [];
        y = unique(y, 'rows');
        y(any(isnan(y), 2), :) = [];
        if size(x, 1) < size(y, 1)
            cy = x;
            x = y;
            y = cy;
        end
        n = size(x, 1);
        m = size(y, 1);
        d = size(x, 2);
        assert(d == size(y, 2));
        [~, r] = knnsearch(x,x,'K',2);
        r = r(:, 2);
        [~, s] = knnsearch(y,x,'K',1);
        KLd = d*sum(-log(r./s))/n + log(m/(n-1));
    case 2
        D = [x; y];
        ndim = ndims(D);
        spacerange = nan(ndim, 2);
        for i = 1:ndim
            spacerange(i, :) = [quantile(D(:, i), 0.05), quantile(D(:, i), 0.95)];
        end
        [X, Y] = meshgrid(linspace(spacerange(1, 1), spacerange(1, 2), 30),...
            linspace(spacerange(2, 1), spacerange(2, 2), 30));
        [fx, ~] = ksdensity(x, [X(:) Y(:)]);
        [fy, xi] = ksdensity(y, [X(:) Y(:)]);
%         rmids = fx~=0 & fy ~= 0;
%         probx = fx(rmids).*log(fx(rmids)./fy(rmids));
        probx = fx.*log((fx+eps)./(fy+eps));
        probx(isinf(probx)) = [];
%         KLd = sum(probx)/length(probx);
        KLd = mean(probx, 'omitnan');
        if isnan(KLd)
            keyboard;
        end
    case 4 % JSD
        
        D = [x; y];
        ndim = ndims(D);
        spacerange = nan(ndim, 2);
        for i = 1:ndim
            spacerange(i, :) = [quantile(D(:, i), 0.05), quantile(D(:, i), 0.95)];
        end
        [X, Y] = meshgrid(linspace(spacerange(1, 1), spacerange(1, 2), 30),...
            linspace(spacerange(2, 1), spacerange(2, 2), 30));
        fx = ksdensity(x, [X(:) Y(:)]);
        fy = ksdensity(y, [X(:) Y(:)]);
        probx = fx.*log((fx+eps)./(0.5*fy+0.5*fx+eps));
        probx(isinf(probx)) = [];
        proby = fy.*log((fy+eps)./(0.5*fy+0.5*fx+eps));
        proby(isinf(proby)) = [];
        KLd = 0.5*mean(probx, 'omitnan') +  0.5*mean(proby, 'omitnan');
        if isnan(KLd)
            keyboard;
        end
        
end