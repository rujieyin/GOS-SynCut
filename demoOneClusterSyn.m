%% setup parameters

% clearvars
% close all
path(pathdef)
addpath(path, genpath([pwd '/utils']))

NVec = 100; % only one cluster

d = 3; % dimension of the orthogonal group
p = 1; % connection probability
% nvar = .3; % noise variance in sampling randSO3

syncRoutine = @syncSpecRelax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate random graph and edge measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = struct();

% construct binary adjacency matrix
G.adjMat = tril(rand(NVec) < p);
G.adjMat = G.adjMat + tril(G.adjMat,-1)';
G.adjMat = G.adjMat - diag(diag(G.adjMat));

% construct vertex coordinates
G.V = rand([NVec, 2]);
% shrink cluster to its centroid
G.V = G.V*.8 + repmat(mean(G.V,1),NVec,1)*.2;

% compute spectral gap
D = sum(G.adjMat);
GL = diag(1./sqrt(D))*G.adjMat*diag(1./sqrt(D));
% GL = GL - diag(diag(GL));
evals = eig(GL);
evals = sort(evals, 'descend');
G.specGap = evals(1) - evals(2);

fprintf('graph spectral gap: %f\n', G.specGap);

% assign edge measurement
vertSig = arrayfun(@(x)orth(rand(d)), 1:NVec, 'UniformOutput',0);
edgeMsr = cell(size(G.adjMat));
[rEdgeIdx, cEdgeIdx] =  find(G.adjMat);
Nedge = length(rEdgeIdx);
for j = 1 : Nedge
    edgeMsr{rEdgeIdx(j), cEdgeIdx(j)} = vertSig{rEdgeIdx(j)}*randSO3(nvar)*vertSig{cEdgeIdx(j)}';
end

fprintf('noise variance: %f \n', nvar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params = struct();
% params.vertSig = vertSig;

[wGCL, wGCL_Dvec, wGCL_W] = assembleGCL(G.adjMat, edgeMsr, d);
RelaxSol = syncRoutine(wGCL, d, wGCL_Dvec);

[v, Rmean] = var_SO3(vertSig, RelaxSol);
fprintf('Difference b/t reconstruction and ground truth signals: %f \n', v);

recMsr = cell(Nedge,1);
inputMsr = cell(Nedge,1);
cleanMsr = cell(Nedge,1);
for j = 1:Nedge
    recMsr{j} = RelaxSol{rEdgeIdx(j)}*RelaxSol{cEdgeIdx(j)}';
    inputMsr{j} = edgeMsr{rEdgeIdx(j), cEdgeIdx(j)};
    cleanMsr{j} = vertSig{rEdgeIdx(j)}*vertSig{cEdgeIdx(j)}';
end
lossrec = cellfun(@(x,y)norm(x-y,'fro'), recMsr, inputMsr);
fprintf('averaged loss on edge measurements of reconstruction: %f \n', sum(lossrec)/Nedge);
lossclean = cellfun(@(x,y)norm(x-y,'fro'), cleanMsr, inputMsr);
fprintf('averaged loss on edge measurements of ground truth: %f \n', sum(lossclean)/Nedge);

% figure;
% hist(lossclean); title('loss on each edge measurement');