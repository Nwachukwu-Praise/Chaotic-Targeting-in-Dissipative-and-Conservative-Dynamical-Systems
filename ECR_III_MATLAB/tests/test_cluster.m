function nfail = test_cluster()
%TEST_CLUSTER  Radius selection and the simple r-neighbourhood clustering.
nfail = 0;
ecr_seed(11);

% three well separated blobs
Z = [0.1*randn(2,200), 0.1*randn(2,200) + [5;0], 0.1*randn(2,200) + [0;5]];

opt = ecr_default_options();
opt.minClusterPts = 10;

[r, info] = ecr_choose_radius(Z, opt);
nfail = nfail + check(r > 0.2 && r < 4.6, ...
    sprintf('radius %.3f should separate the blobs (0.2 < r < 4.6)', r));
nfail = nfail + check(~info.fallback, 'histogram minimum found for a multi-modal set');

lab = ecr_cluster(Z, r, opt);
nfail = nfail + check(numel(unique(lab(lab>0))) == 3, ...
    sprintf('found %d clusters, expected 3', numel(unique(lab(lab>0)))));

% chaining: a dense bridge between two blobs must merge them
Zb = [0.1*randn(2,200), 0.1*randn(2,200)+[3;0], [linspace(0,3,60); zeros(1,60)]];
lab2 = ecr_cluster(Zb, 0.4, opt);
nfail = nfail + check(numel(unique(lab2(lab2>0))) == 1, ...
    'r-neighbourhood chaining merges bridged blobs');

% small clusters are merged, not left dangling
Zs = [0.1*randn(2,200), 0.1*randn(2,3)+[9;9]];
lab3 = ecr_cluster(Zs, 0.5, opt);
nfail = nfail + check(all(lab3 > 0) && numel(unique(lab3)) == 1, ...
    'clusters below minClusterPts are merged');

% cluster statistics
CL = ecr_cluster_stats(Z, lab, opt);
nfail = nfail + check(numel(CL) == 3 && abs(sum([CL.n]) - size(Z,2)) == 0, ...
    'cluster statistics cover every point');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
