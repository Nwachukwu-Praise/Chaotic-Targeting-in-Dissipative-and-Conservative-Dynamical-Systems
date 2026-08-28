function nfail = test_control()
%TEST_CONTROL  Controller plumbing: histogram helper, selection rules, escape.
nfail = 0;
ecr_seed(4);

% --- the histogram helper agrees with a manual count ---------------------
x = [0 0.5 1 1.5 2 2.5 3];
e = 0:1:3;
n = ecr_histcounts(x, e);
nfail = nfail + check(isequal(n, [2 2 3]), ...
    sprintf('ecr_histcounts gave [%s], expected [2 2 3]', num2str(n)));

% --- a trained model on the cheap fixture --------------------------------
S   = sys_logistic_fixture();
opt = ecr_default_options();
opt.nData = 4000;  opt.trajLen = 25;  opt.verbose = 0;  opt.Kmax = 4;
D = ecr_generate_data(S, opt);
m = ecr_train(S, D, opt, 'ECR-III');

Z = S.obs(S.init(300));

% --- both selection rules stay inside Pi and agree about "no targeting" --
o1 = opt;  o1.selection = 'nmd';    m1 = m;  m1.opt = o1;
o2 = opt;  o2.selection = 'level';  m2 = m;  m2.opt = o2;
[P1, i1] = ecr_control(m1, Z);
[P2, i2] = ecr_control(m2, Z);
nfail = nfail + check(all(ecr_in_pi(P1,S)) && all(ecr_in_pi(P2,S)), ...
                      'both selection rules respect Pi');
nfail = nfail + check(all(i2.region(i2.region>=0) >= 0), 'region indices are sane');
nfail = nfail + check(mean(i2.region <= i1.region) > 0.5, ...
                      'selection=level does not prefer higher regions');

% --- states far away from every cluster get the nominal parameters -------
Zfar = 1e3*ones(S.N, 20);
[Pf, iff] = ecr_control(m, Zfar);
nfail = nfail + check(all(iff.region == -1) && ...
                      max(max(abs(Pf - repmat(S.pnom,1,20)))) == 0, ...
                      'no targeting far outside the attractor');

% --- the NMD gate really does reject far clusters -----------------------
o3 = opt; o3.nmdGate = 1e-6;  m3 = m;  m3.opt = o3;
[~, i3] = ecr_control(m3, Z);
nfail = nfail + check(mean(i3.region >= 0) < mean(i1.region >= 0), ...
                      'a tighter NMD gate claims fewer states');

% --- escape rule: a controller stuck outside S0 must fall back to p_nom --
opt2 = opt; opt2.escapeAfter = 3; opt2.maxSteps = 12; opt2.holdSteps = 0;
stuck = @(z) deal(S.pnom + S.dpmax, struct('region', 1, 'cluster', 1, 'nmd', 0));
R = ecr_simulate(S, stuck, S.init(1), opt2);
nfail = nfail + check(any(all(abs(R.P - repmat(S.pnom,1,size(R.P,2))) < 1e-12, 1)), ...
                      'escapeAfter applies the nominal parameters');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
