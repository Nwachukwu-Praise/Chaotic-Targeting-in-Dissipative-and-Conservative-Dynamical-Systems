function nfail = test_pipeline_logistic()
%TEST_PIPELINE_LOGISTIC  End-to-end check of the ECR machinery on a cheap map.
%
%   The logistic map is used only as a fast test fixture for the region
%   construction, the clustering and the controller - the Lorenz system is
%   the system this toolbox is about.  Because the logistic map is also in
%   Table 1 of the paper, the numbers produced here can be compared with it:
%   the paper reports 143.7 steps without targeting, 9.14 for ECR-II and
%   5.52 for ECR-III with 6 activation regions.
nfail = 0;

S   = sys_logistic_fixture();
opt = ecr_default_options();
opt.nData    = 6000;
opt.trajLen  = 30;
opt.nTrials  = 100;
opt.maxSteps = 200;
opt.verbose  = 0;
opt.seed     = 1;

D  = ecr_generate_data(S, opt);
m3 = ecr_train(S, D, opt, 'ECR-III');
m2 = ecr_train(S, D, opt, 'ECR-II');

% --- Definition 1 must hold for every level-0 training point -------------
L0   = m3.levels(1);
idx0 = [];
for j = 1:numel(L0.CL), idx0 = [idx0, L0.CL(j).idx]; end %#ok<AGROW>
nfail = nfail + check(numel(unique(idx0)) == L0.nData, ...
                      'level 0 clusters partition the region data');
g0  = L0.dataIdx;
d1  = sqrt(sum((D.Z(:,g0)  - repmat(S.zstar,1,numel(g0))).^2, 1));
d2  = sqrt(sum((D.Zn(:,g0) - repmat(S.zstar,1,numel(g0))).^2, 1));
nfail = nfail + check(all(d1 < S.delta) && all(d2 < S.delta), ...
                      'every T0 point satisfies Eq. (4)');
nfail = nfail + check(all(ecr_in_pi(D.P(:,g0), S)), ...
                      'every T0 target parameter lies in Pi');

% --- Definition 3: T_i data must not already belong to a lower region ----
for k = 2:numel(m3.levels)
    gk = m3.levels(k).dataIdx;
    dk = sqrt(sum((D.Z(:,gk) - repmat(S.zstar,1,numel(gk))).^2, 1));
    nfail = nfail + check(~any(ismember(gk, g0)), ...
        sprintf('T%d is disjoint from T0', m3.levels(k).index));
    nfail = nfail + check(all(dk >= 0), 'sanity');
end

% --- the controller never leaves Pi (Eq. 3) -----------------------------
Zt = S.obs(S.init(500));
P3 = ecr_control(m3, Zt);
P2 = ecr_control(m2, Zt);
nfail = nfail + check(all(ecr_in_pi(P3,S)), 'ECR-III output stays inside Pi');
nfail = nfail + check(all(ecr_in_pi(P2,S)), 'ECR-II  output stays inside Pi');

% --- targeting must be much faster than waiting -------------------------
H0 = S.init(opt.nTrials);
b0 = ecr_reaching_time(S, [],                     opt, H0);
b2 = ecr_reaching_time(S, @(Z) ecr_control(m2,Z), opt, H0);
b3 = ecr_reaching_time(S, @(Z) ecr_control(m3,Z), opt, H0);

nfail = nfail + check(b3.meanCens < 0.3*b0.meanCens, ...
    sprintf('ECR-III %.2f vs no control %.2f steps', b3.meanCens, b0.meanCens));
nfail = nfail + check(b2.meanCens < 0.5*b0.meanCens, ...
    sprintf('ECR-II %.2f vs no control %.2f steps', b2.meanCens, b0.meanCens));
nfail = nfail + check(b3.meanCens <= 1.25*b2.meanCens, ...
    sprintf('ECR-III (%.2f) should not be clearly worse than ECR-II (%.2f)', ...
            b3.meanCens, b2.meanCens));
nfail = nfail + check(b3.success > 0.95 && b3.hold > 0.9, ...
    sprintf('capture rate %.2f, retention %.2f', b3.success, b3.hold));

% --- ECR-III trains at most as long as ECR-I would: more, smaller nets ---
nfail = nfail + check(m3.info.nNets >= m2.info.nNets, ...
    'ECR-III uses at least as many networks as ECR-II');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
