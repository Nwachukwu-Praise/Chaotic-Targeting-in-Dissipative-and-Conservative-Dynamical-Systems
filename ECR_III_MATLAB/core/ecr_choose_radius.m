function [r, info] = ecr_choose_radius(Z, opt)
%ECR_CHOOSE_RADIUS  Clustering radius from the inter-data distance histogram.
%
%   [R,INFO] = ECR_CHOOSE_RADIUS(Z,OPT)
%
%   Section 3.3 of the paper: "First a radius r is chosen as the first
%   minimum of the histogram of the inter-data distances by visual analysis.
%   Data points, which lie within r-neighbourhood of each other, are assigned
%   to the same cluster."
%
%   The visual step is automated here: the histogram of all pairwise
%   distances (on a random subset of at most OPT.histSubset points) is
%   smoothed with a moving average of width OPT.histSmooth and the first
%   local minimum after the first local maximum is taken.  INFO returns the
%   histogram so that the choice can still be inspected by eye - see
%   ECR_PLOT_RADIUS.
%
%   If no interior minimum exists (a single-mode histogram, which happens for
%   compact regions) the routine falls back to twice the median
%   nearest-neighbour distance and sets INFO.fallback = true.

if nargin < 2 || isempty(opt), opt = ecr_default_options(); end

n = size(Z,2);
idx = 1:n;
if n > opt.histSubset
    idx = randperm(n, opt.histSubset);
end
X = Z(:,idx);
m = size(X,2);

% ---- pairwise distances -------------------------------------------------
G  = X'*X;
sq = diag(G);
D2 = repmat(sq,1,m) + repmat(sq',m,1) - 2*G;
D2 = max(D2, 0);
D  = sqrt(D2);
mask = triu(true(m), 1);
dv   = D(mask);
dv   = dv(isfinite(dv));

info = struct('fallback', false, 'centers', [], 'counts', [], 'nnmed', NaN);
if isempty(dv) || max(dv) <= 0
    r = 1e-6;  info.fallback = true;  return
end

% nearest neighbour distances (used for the fallback and for diagnostics)
Dn = D + diag(inf(m,1));
nn = min(Dn, [], 2);
info.nnmed = median(nn(isfinite(nn)));

% ---- histogram ----------------------------------------------------------
nb    = opt.nHistBins;
edges = linspace(0, max(dv), nb+1);
cnt   = ecr_histcounts(dv, edges);
ctr   = 0.5*(edges(1:nb) + edges(2:nb+1));

w = max(1, round(opt.histSmooth));
if w > 1
    kern = ones(1,w)/w;
    cs   = conv(cnt, kern, 'same');
else
    cs = cnt;
end
info.centers = ctr;
info.counts  = cnt;
info.smooth  = cs;

r = [];
if strcmp(opt.radiusMode, 'hist')
    % first local maximum ...
    ip = [];
    for i = 2:nb-1
        if cs(i) >= cs(i-1) && cs(i) > cs(i+1), ip = i; break, end
    end
    % ... then the first local minimum after it
    if ~isempty(ip)
        for i = ip+1:nb-1
            if cs(i) <= cs(i-1) && cs(i) < cs(i+1)
                r = ctr(i);
                break
            end
        end
    end
end

if isempty(r) || r <= 0
    r = 2*info.nnmed;
    info.fallback = true;
end
r = opt.radiusScale * r;
end
