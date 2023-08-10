function Edges = findedges(X, nk, distancetype)
if nargin == 1
    nk = 3;
    distancetype = 'correlation';
elseif nargin == 2
    distancetype = 'correlation';
end
nP = size(X, 1);
[Idx, pDist] = knnsearch(X, X, 'K', nk+1, 'Distance', distancetype);
Idx = Idx(:, 2:end);
pDist = pDist(:, 2:end);

% pairs = nchoosek(1:nk, 2);
% npair = size(pairs, 1);
Edges = [];
for i = 1:nk
    Edges = [Edges; (1:nP)' Idx(:, i) pDist(:, i)];
end