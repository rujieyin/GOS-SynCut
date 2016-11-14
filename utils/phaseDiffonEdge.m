function EdgeDiffSum = phaseDiffonEdge(NVec, adjMat, vertPotCell, SolCell)

% phase difference on each vertex in the graph
phaseDiffVert = cellfun(@(x,y)x'*y, vertPotCell, SolCell, 'UniformOutput',0);

EdgeDiff = sparse(size(adjMat,1), size(adjMat,2));

[Ridx, Cidx] = find(triu(adjMat, 1));

for i = 1:length(Ridx)
    EdgeDiff(Ridx(i), Cidx(i)) = norm(phaseDiffVert{Ridx(i)} - phaseDiffVert{Cidx(i)}, 'fro');
end

EdgeDiffSum = zeros(length(NVec));

ClusterIdx = cumsum(NVec);
ClusterIdx = [0 ClusterIdx];

for i = 1:length(NVec)
    for j = i:length(NVec)
        EdgeDiffSum( i, j ) = sum(reshape( EdgeDiff(ClusterIdx(i)+1:ClusterIdx(i+1), ClusterIdx(j)+1:ClusterIdx(j+1)), [], 1));
    end
end
