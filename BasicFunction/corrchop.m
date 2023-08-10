function R = corrchop(a, b, ndiv)
nL = size(a, 1);
if nargin == 2
    ndiv = ceil(nL/1200); % 60 Frame 5 Hz, 100 Sampling Hz
end
chopL = round(nL/ndiv);
for i = 1:ndiv
    sP = (i-1)*chopL+1;
    eP = min([i*chopL nL]);
    if i == 1
        dr = corr(a(sP:eP, :), b(sP:eP, :));
    else
        dr(:, :, end+1) = corr(a(sP:eP, :), b(sP:eP, :));
    end
end
R = median(dr, 3, 'omitnan');