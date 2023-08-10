function r = medianR(Rmat)
I = ones(size(Rmat));
Ids = tril(I, -1)==1;
r = median(Rmat(Ids(:)), 'omitnan');
