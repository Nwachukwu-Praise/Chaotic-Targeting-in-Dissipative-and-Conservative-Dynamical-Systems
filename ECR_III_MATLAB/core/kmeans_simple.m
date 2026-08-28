function [C, lab] = kmeans_simple(X, k, iters)
%KMEANS_SIMPLE  Minimal k-means (no toolbox required), k-means++ seeding.
%
%   [C,LAB] = KMEANS_SIMPLE(X,K,ITERS) clusters the columns of X (d x n) into
%   K groups.  C is d x K, LAB is 1 x n.  Used only to place RBF centres.

if nargin < 3 || isempty(iters), iters = 25; end
[d, n] = size(X);
k = min(k, n);

% ---- k-means++ initialisation ------------------------------------------
C = zeros(d, k);
C(:,1) = X(:, randi(n));
if k > 1
    d2 = sum((X - repmat(C(:,1),1,n)).^2, 1);
    for j = 2:k
        if sum(d2) <= 0
            C(:,j) = X(:, randi(n));
        else
            c = cumsum(d2/sum(d2));
            C(:,j) = X(:, find(rand <= c, 1, 'first'));
        end
        d2 = min(d2, sum((X - repmat(C(:,j),1,n)).^2, 1));
    end
end

% ---- Lloyd iterations ---------------------------------------------------
lab = ones(1,n);
for it = 1:iters
    Dm = zeros(k, n);
    for j = 1:k
        Dm(j,:) = sum((X - repmat(C(:,j),1,n)).^2, 1);
    end
    [~, newlab] = min(Dm, [], 1);
    if it > 1 && all(newlab == lab), break, end
    lab = newlab;
    for j = 1:k
        m = (lab == j);
        if any(m)
            C(:,j) = mean(X(:,m), 2);
        else
            C(:,j) = X(:, randi(n));       % restart an empty cluster
        end
    end
end
end
