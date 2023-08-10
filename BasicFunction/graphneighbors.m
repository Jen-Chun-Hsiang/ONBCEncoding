function neighbormap = graphneighbors(G)
nnodes = size(G.Nodes, 1);
neighbormap = nan(nnodes, nnodes);
for i = 1:nnodes
    neighbormap(i, :) = ismember((1:nnodes), neighbors(G, i));
end
