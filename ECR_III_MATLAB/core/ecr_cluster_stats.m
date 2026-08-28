function CL = ecr_cluster_stats(Z, lab, opt)
%ECR_CLUSTER_STATS  Cluster means and normalised covariance matrices (Eq. 8).
%
%   CL = ECR_CLUSTER_STATS(Z,LAB,OPT) returns a 1 x nc struct array with
%
%     CL(j).idx    indices of the columns of Z belonging to cluster j
%     CL(j).n      number of points
%     CL(j).mu     mean of the cluster                    mu_ij  (Eq. 8)
%     CL(j).Sigma  (regularised) covariance matrix
%     CL(j).Ninv   NORMALISED covariance matrix N_ij used in Eq. (8)
%     CL(j).nmd95  95th percentile of the NMD of its own training data
%     CL(j).radius max NMD of its own training data
%
%   Eq. (8) of the paper reads
%
%       NMD_ij(z)^2 = (z - mu_ij)' * N_ij * (z - mu_ij)
%
%   with N_ij "the Normalized Covariance Matrix".  The paper does not define
%   the normalisation, so OPT.covNorm selects one of
%
%     'det'   N = inv( Sigma / det(Sigma)^(1/N) )   scale free, unit
%             determinant: clusters are compared by shape, not by size
%             (default - a big and a small cluster are equally "close" at
%             equal shape-normalised distance)
%     'plain' N = inv(Sigma)                        ordinary Mahalanobis
%     'diag'  N = diag(1./var)                      "normalised variances"
%             (ignores the cross correlations)
%     'trace' N = inv( Sigma / (trace(Sigma)/N) )
%
%   See docs/METHOD_NOTES.md for the discussion.

if nargin < 3 || isempty(opt), opt = ecr_default_options(); end
N  = size(Z,1);
u  = unique(lab(lab>0));
CL = repmat(struct('idx',[],'n',0,'mu',zeros(N,1),'Sigma',eye(N), ...
                   'Ninv',eye(N),'nmd95',0,'radius',0), 1, numel(u));

for j = 1:numel(u)
    idx = find(lab == u(j));
    X   = Z(:,idx);
    n   = numel(idx);
    mu  = mean(X,2);
    Xc  = X - repmat(mu,1,n);

    if n > 1
        Sigma = (Xc*Xc')/(n-1);
    else
        Sigma = zeros(N);
    end

    % isotropic fallback for degenerate clusters
    s2 = mean(sum(Xc.^2,1))/max(N,1);
    if ~isfinite(s2) || s2 <= 0, s2 = 1; end
    if n <= N || ~all(isfinite(Sigma(:))) || rcond(Sigma + eps*eye(N)) < 1e-10
        Sigma = s2*eye(N);
    end
    Sigma = Sigma + max(opt.covReg*trace(Sigma)/N, 1e-12)*eye(N);
    if rcond(Sigma) < 1e-10                    % last resort: isotropic
        Sigma = max(s2, 1e-12)*eye(N);
    end

    switch lower(opt.covNorm)
        case 'det'
            dt = det(Sigma);
            if ~isfinite(dt) || dt <= 0
                Sn = Sigma;
            else
                Sn = Sigma / dt^(1/N);
            end
            Ninv = inv(Sn);
        case 'plain'
            Ninv = inv(Sigma);
        case 'diag'
            v = diag(Sigma);
            v(v <= 0 | ~isfinite(v)) = s2;
            Ninv = diag(1./v);
        case 'trace'
            Sn   = Sigma / (trace(Sigma)/N);
            Ninv = inv(Sn);
        otherwise
            error('ecr_cluster_stats:covNorm', ...
                  'unknown covNorm "%s"', opt.covNorm);
    end
    Ninv = (Ninv + Ninv')/2;

    CL(j).idx   = idx;
    CL(j).n     = n;
    CL(j).mu    = mu;
    CL(j).Sigma = Sigma;
    CL(j).Ninv  = Ninv;

    d = sqrt(max(sum(Xc .* (Ninv*Xc), 1), 0));
    CL(j).nmd95  = prctile_simple(d, 95);
    CL(j).radius = max(d);
end
end

% =========================================================================
function q = prctile_simple(x, p)
x = sort(x(isfinite(x)));
if isempty(x), q = 0; return, end
if numel(x) == 1, q = x; return, end
pos = 1 + (p/100)*(numel(x)-1);
lo  = floor(pos); hi = ceil(pos);
q   = x(lo) + (pos-lo)*(x(hi)-x(lo));
end
