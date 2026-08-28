function nfail = test_nmd()
%TEST_NMD  Properties of the Normalised Mahalanobis Distance of Eq. (8).
nfail = 0;
ecr_seed(5);

A = [2 0.6; 0.6 0.5];
Z = A*randn(2,4000) + [3; -1];
lab = ones(1, size(Z,2));

opt = ecr_default_options();

% --- distance to its own mean is zero, distances are non-negative --------
CL = ecr_cluster_stats(Z, lab, opt);
d  = ecr_nmd(CL, [CL(1).mu, Z(:,1:10)]);
nfail = nfail + check(abs(d(1)) < 1e-10, 'NMD of the cluster mean is zero');
nfail = nfail + check(all(d >= 0), 'NMD is non-negative');

% --- 'plain' reproduces the textbook Mahalanobis distance ----------------
o2 = opt; o2.covNorm = 'plain'; o2.covReg = 0;
CL2 = ecr_cluster_stats(Z, lab, o2);
x   = Z(:,1:50);
dx  = x - repmat(CL2(1).mu,1,50);
ref = sqrt(sum(dx .* (CL2(1).Sigma\dx), 1));
nfail = nfail + check(max(abs(ecr_nmd(CL2,x) - ref)) < 1e-8, ...
                      'covNorm=plain equals the Mahalanobis distance');

% --- 'plain' is invariant when the whole cluster is scaled --------------
CLb = ecr_cluster_stats(10*Z, lab, o2);
da  = ecr_nmd(CL2, Z(:,1:50));
db  = ecr_nmd(CLb, 10*Z(:,1:50));
nfail = nfail + check(max(abs(da - db)./max(da,1e-9)) < 1e-6, ...
                      'covNorm=plain is scale invariant');

% --- 'det' removes only the size, so for an isotropic cluster it reduces
%     to the plain Euclidean distance ------------------------------------
o3  = opt; o3.covNorm = 'det'; o3.covReg = 0;
Zi  = 3*randn(2,4000) + [1;2];
CLi = ecr_cluster_stats(Zi, ones(1,4000), o3);
x   = Zi(:,1:50);
de  = sqrt(sum((x - repmat(CLi(1).mu,1,50)).^2, 1));
dn  = ecr_nmd(CLi, x);
nfail = nfail + check(max(abs(dn - de)./de) < 0.05, ...
                      'covNorm=det equals the Euclidean distance for isotropic clusters');

% --- the 95th percentile bookkeeping used by the gate is sane -----------
nfail = nfail + check(CL(1).nmd95 > 0 && CL(1).nmd95 <= CL(1).radius, ...
                      'nmd95 <= max training NMD');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
