function r = meanR(Rmat)
I = ones(size(Rmat));
Ids = tril(I, -1)==1;
r = mean(Rmat(Ids(:)), 'omitnan');
