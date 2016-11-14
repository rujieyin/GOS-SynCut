function [r, var, Rmeans] = varRatio(NVec, vertPotCell, SolCell)

var = zeros(length(NVec) + 1);
idx = cumsum(NVec);
idx = [0 idx];
Rmeans = cell(length(NVec)+1,1);
for i = 1:length(NVec)
    [var(i), Rmeans{i}] = var_SO3(SolCell(idx(i)+1:idx(i+1)), vertPotCell(idx(i)+1:idx(i+1)) );
end

[var(end), Rmeans{end}] = var_SO3(SolCell, vertPotCell);
r = sum(var(1:end-1))/var(end);