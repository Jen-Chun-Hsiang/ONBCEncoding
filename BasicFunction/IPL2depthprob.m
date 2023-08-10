function depthprob = IPL2depthprob(IPLd)
    nType = size(IPLd, 1);
    for i = 1:nType
        for j = 1:nType
            depthprob(i, j) = corr(IPLd(i, :)', IPLd(j, :)').^2;
        end
    end
end