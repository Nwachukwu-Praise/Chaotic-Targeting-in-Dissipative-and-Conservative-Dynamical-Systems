function lab = ecr_cluster(Z, r, opt)
%ECR_CLUSTER  The "simple clustering algorithm" of ECR-III (Section 3.3).
%
%   LAB = ECR_CLUSTER(Z,R,OPT) assigns every column of Z (N x n) to a
%   cluster.  Two points are put in the same cluster whenever they lie within
%   the radius R of each other; the relation is transitive, i.e. an arbitrary
%   data point is taken as the centre of a hyper-sphere of radius R and the
%   sphere is grown over the points it captures (single-linkage chaining).
%
%   Post-processing (not specified in the paper, see docs/METHOD_NOTES.md):
%     * clusters with fewer than OPT.minClusterPts points are merged into the
%       nearest large cluster (OPT.smallCluster = 'merge') or removed
%       (LAB = 0, OPT.smallCluster = 'drop');
%     * at most OPT.maxClusters clusters are kept, the smallest ones being
%       merged into the nearest kept cluster.

if nargin < 3 || isempty(opt), opt = ecr_default_options(); end
n = size(Z,2);
lab = zeros(1,n);
if n == 0, return, end

% ---- union-find over all pairs closer than r ---------------------------
parent = 1:n;
blk = max(1, floor(2e6/max(n,1)));           % keep the distance block small
for i0 = 1:blk:n
    i1 = min(i0+blk-1, n);
    B  = Z(:, i0:i1);
    G  = B'*Z;
    sb = sum(B.^2,1)';  sz = sum(Z.^2,1);
    D2 = repmat(sb,1,n) + repmat(sz,size(B,2),1) - 2*G;
    % Only the (at most KCAP) nearest neighbours inside the radius are used
    % as union-find edges: the connected components are the same as with all
    % pairs but the edge list stays small for dense regions.
    kcap = 32;
    for t = 1:size(B,2)
        a  = i0 + t - 1;
        nb = find(D2(t,:) <= r^2);
        nb(nb == a) = [];
        if numel(nb) > kcap
            [~, o] = sort(D2(t,nb));
            nb = nb(o(1:kcap));
        end
        ra = uf_find(parent, a);
        for s = 1:numel(nb)
            rb = uf_find(parent, nb(s));
            if ra ~= rb
                parent(rb) = ra;
            end
        end
    end
end
for i = 1:n, parent(i) = uf_find(parent, i); end

% ---- compact labels -----------------------------------------------------
[~, ~, lab] = unique(parent);
lab = lab(:)';

lab = merge_small(Z, lab, opt.minClusterPts, opt.smallCluster);

% ---- cap the number of clusters ----------------------------------------
u = unique(lab(lab>0));
if numel(u) > opt.maxClusters
    cnt = zeros(1,numel(u));
    for j = 1:numel(u), cnt(j) = sum(lab == u(j)); end
    [~, ord] = sort(cnt, 'descend');
    keep = u(ord(1:opt.maxClusters));
    small = setdiff(u, keep);
    thr = 0;
    for j = 1:numel(small)
        thr = max(thr, sum(lab == small(j)));
    end
    lab = merge_small(Z, lab, thr+1, opt.smallCluster, keep);
end

% ---- renumber -----------------------------------------------------------
u = unique(lab(lab>0));
out = zeros(1,numel(lab));
for j = 1:numel(u)
    out(lab == u(j)) = j;
end
lab = out;
end

% =========================================================================
function lab = merge_small(Z, lab, minPts, mode, keepList)
u = unique(lab(lab>0));
cnt = zeros(1,numel(u));
mu  = zeros(size(Z,1), numel(u));
for j = 1:numel(u)
    m = (lab == u(j));
    cnt(j) = sum(m);
    mu(:,j) = mean(Z(:,m), 2);
end
if nargin >= 5 && ~isempty(keepList)
    big = ismember(u, keepList);
else
    big = cnt >= minPts;
end
if ~any(big)
    % No cluster reaches the minimum size (happens when the region holds
    % only a handful of points): keep the largest one and merge the rest
    % into it, rather than leaving a swarm of one-point clusters behind.
    [~, w] = max(cnt);
    big = false(size(cnt));
    big(w) = true;
end
if all(big), return, end

for j = find(~big)
    m = (lab == u(j));
    if strcmp(mode, 'drop')
        lab(m) = 0;
    else
        d = sqrt(sum((mu(:,big) - repmat(mu(:,j),1,sum(big))).^2, 1));
        ub = u(big);
        [~, w] = min(d);
        lab(m) = ub(w);
    end
end
end

% =========================================================================
function r = uf_find(parent, i)
r = i;
while parent(r) ~= r
    r = parent(r);
end
end
