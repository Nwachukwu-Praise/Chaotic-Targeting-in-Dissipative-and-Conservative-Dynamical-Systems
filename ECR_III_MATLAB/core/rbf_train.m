function net = rbf_train(X, Y, opt)
%RBF_TRAIN  Train a radial basis function network  Y ~ W*[phi(X); 1].
%
%   NET = RBF_TRAIN(X,Y,OPT) fits a Gaussian RBF network mapping the columns
%   of X (d x n, states z_n) to the columns of Y (q x n, normalised parameter
%   perturbations).  Section 3 of the paper: "RBF-based neural networks have
%   been preferred due to their better local approximation capability".
%
%   Training is the standard two-stage scheme, which is what makes ECR cheap
%   to train compared with back-propagation:
%     1. centres      - k-means on the (standardised) inputs
%     2. widths       - mean distance to the nearest centres times OPT.rbfSpread
%     3. weights      - regularised linear least squares (closed form)
%
%   NET fields: mu, sd (input standardisation), C (centres), sig (widths),
%   W (weights), n (training size), const (fallback constant output).

if nargin < 3 || isempty(opt), opt = ecr_default_options(); end
[d, n] = size(X);
q = size(Y,1);

net = struct('mu',zeros(d,1),'sd',ones(d,1),'C',[],'sig',[],'W',[], ...
             'n',n,'const',zeros(q,1),'q',q,'d',d);

if n == 0
    return
end
net.const = mean(Y, 2);
if n < opt.rbfMinPts
    return                                  % too little data: constant model
end

net.mu = mean(X, 2);
sd     = std(X, 0, 2);
sd(sd < 1e-12 | ~isfinite(sd)) = 1;
net.sd = sd;
Xn = (X - repmat(net.mu,1,n)) ./ repmat(net.sd,1,n);

k = max(1, min(opt.nCenters, floor(n/2)));
C = kmeans_simple(Xn, k);
k = size(C,2);

% ---- widths -------------------------------------------------------------
sig = zeros(1,k);
if k == 1
    sig(1) = max(mean(sqrt(sum((Xn - repmat(C,1,n)).^2,1))), 1e-6);
else
    for j = 1:k
        dj = sqrt(sum((C - repmat(C(:,j),1,k)).^2, 1));
        dj(j) = Inf;
        dj = sort(dj);
        sig(j) = mean(dj(1:min(3,k-1)));
    end
    sig(~isfinite(sig) | sig < 1e-6) = 1e-6;
end
sig = opt.rbfSpread * sig;

% ---- weights ------------------------------------------------------------
Phi = rbf_design(Xn, C, sig);
A   = Phi'*Phi + opt.rbfLambda*eye(size(Phi,2));
W   = A \ (Phi' * Y');                       % (k+1) x q

net.C   = C;
net.sig = sig;
net.W   = W;
end
