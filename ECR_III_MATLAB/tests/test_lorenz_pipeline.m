function nfail = test_lorenz_pipeline()
%TEST_LORENZ_PIPELINE  Small end-to-end ECR-III run on the Lorenz system.
%
%   Deliberately small (a few thousand map steps) so that the suite stays
%   runnable; the full experiment is DEMO_LORENZ_ECR3.
nfail = 0;

S   = sys_lorenz_real();
opt = ecr_default_options();
opt.nData    = 2500;
opt.trajLen  = 25;
opt.Kmax     = 4;
opt.nTrials  = 40;
opt.maxSteps = 120;
opt.holdSteps= 10;
opt.verbose  = 0;
opt.seed     = 2;

D  = ecr_generate_data(S, opt);
nfail = nfail + check(size(D.Z,2) > 0.9*opt.nData, 'most generated triples are valid');

m3 = ecr_train(S, D, opt, 'ECR-III');
nfail = nfail + check(numel(m3.levels) >= 2, ...
    sprintf('at least T0 and T1 were built (got %d regions)', numel(m3.levels)));

% the controller must respect Eq. (3) everywhere, including far from the data
Zfar = repmat(S.zstar,1,200) + 50*randn(2,200);
Pf   = ecr_control(m3, Zfar);
nfail = nfail + check(all(ecr_in_pi(Pf,S)), 'controller stays in Pi far from the data');

Zt = S.obs(S.init(100));
Pt = ecr_control(m3, Zt);
nfail = nfail + check(all(ecr_in_pi(Pt,S)), 'controller stays in Pi on the attractor');

H0 = S.init(opt.nTrials);
b0 = ecr_reaching_time(S, [],                     opt, H0);
b3 = ecr_reaching_time(S, @(Z) ecr_control(m3,Z), opt, H0);

nfail = nfail + check(b0.meanCens > 8 && b0.meanCens < 30, ...
    sprintf('uncontrolled reaching time %.2f (paper: 15.09)', b0.meanCens));
nfail = nfail + check(b3.meanCens < 0.6*b0.meanCens, ...
    sprintf('ECR-III %.2f vs uncontrolled %.2f', b3.meanCens, b0.meanCens));
nfail = nfail + check(b3.success > 0.95, sprintf('capture rate %.2f', b3.success));
nfail = nfail + check(b3.hold > 0.8, sprintf('retention after capture %.2f', b3.hold));

% noise must not break the controller
optn = opt; optn.noise = 0.0008;
bn = ecr_reaching_time(S, @(Z) ecr_control(m3,Z), optn, H0);
nfail = nfail + check(bn.meanCens < 0.8*b0.meanCens, ...
    sprintf('ECR-III with 0.08%% noise %.2f', bn.meanCens));
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
