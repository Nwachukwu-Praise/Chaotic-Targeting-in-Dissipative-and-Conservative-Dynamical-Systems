function nfail = test_rbf()
%TEST_RBF  The RBF network must interpolate a smooth 2-D function well and
%   must fall back to a constant when starved of data.
nfail = 0;
ecr_seed(7);

f = @(X) [sin(X(1,:)) + 0.3*X(2,:).^2; cos(0.5*X(1,:)).*X(2,:)];
X = [4*rand(1,600)-2; 4*rand(1,600)-2];
Y = f(X);

opt = ecr_default_options();
opt.nCenters = 40;
net = rbf_train(X, Y, opt);

Xt = [4*rand(1,300)-2; 4*rand(1,300)-2];
Yt = f(Xt);
E  = rbf_eval(net, Xt) - Yt;
rmse = sqrt(mean(E(:).^2));
nfail = nfail + check(rmse < 0.05, sprintf('RBF test RMSE = %.4f (want < 0.05)', rmse));

% a network trained on too few points returns the mean of the targets
net2 = rbf_train(X(:,1:3), Y(:,1:3), opt);
y2   = rbf_eval(net2, Xt(:,1:5));
nfail = nfail + check(max(max(abs(y2 - repmat(mean(Y(:,1:3),2),1,5)))) < 1e-12, ...
                      'RBF constant fallback for tiny data sets');

% k-means returns the requested number of centres and is deterministic given
% a seed
ecr_seed(3);  C1 = kmeans_simple(X, 5);
ecr_seed(3);  C2 = kmeans_simple(X, 5);
nfail = nfail + check(size(C1,2) == 5 && isequal(C1,C2), 'kmeans_simple size/repeatability');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
