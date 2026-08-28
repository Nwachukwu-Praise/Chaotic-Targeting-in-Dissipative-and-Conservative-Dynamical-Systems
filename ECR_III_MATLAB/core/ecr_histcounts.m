function n = ecr_histcounts(x, edges)
%ECR_HISTCOUNTS  Bin counts that work in every MATLAB and Octave version.
%
%   N = ECR_HISTCOUNTS(X,EDGES) returns numel(EDGES)-1 counts, the k-th being
%   the number of samples with EDGES(k) <= x < EDGES(k+1); the last bin
%   includes the right edge.  Uses HISTCOUNTS where available and falls back
%   to a sort-based count otherwise (HISTC is deprecated in MATLAB).

x = x(:);
x = x(isfinite(x));
if exist('histcounts', 'file') == 2 || exist('histcounts', 'builtin') == 5
    n = histcounts(x, edges);
    n = n(:).';
    return
end
edges = edges(:).';
nb = numel(edges) - 1;
idx = zeros(size(x));
for k = 1:nb
    if k < nb
        idx = idx + (x >= edges(k) & x < edges(k+1))*k;
    else
        idx = idx + (x >= edges(k) & x <= edges(k+1))*k;
    end
end
n = zeros(1, nb);
for k = 1:nb
    n(k) = sum(idx == k);
end
end
