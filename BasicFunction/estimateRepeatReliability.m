function R2 = estimateRepeatReliability(Sigs)
% Sigs = [Repeats, Time]
if size(Sigs, 1) > size(Sigs, 2)
    Sigs = Sigs';
end
Sigs(isnan(Sigs(:, 1)), :) = [];
nRepeat = size(Sigs, 1);
comb = nchoosek(1:nRepeat, round(nRepeat/2));
ncomb = size(comb, 1);
rperm = nan(ncomb, 1);
reps = 1:nRepeat;
for i = 1:ncomb
    comb1 = comb(i, :);
    comb2 = reps(~ismember(reps, comb1));
    rperm(i) = corr(mean(Sigs(comb1, :), 1)',...
        mean(Sigs(comb2, :), 1)').^2;
end
R2 = median(rperm, 'omitnan');