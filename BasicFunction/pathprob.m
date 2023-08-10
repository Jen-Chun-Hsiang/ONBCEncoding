function probability = pathprob(path, pathprobmap)
nedge = length(path);
probability = nan(nedge-1, 1);
for i = 1:nedge-1
    probability(i) = pathprobmap(path(i), path(i+1));
end
probability = prod(probability);


