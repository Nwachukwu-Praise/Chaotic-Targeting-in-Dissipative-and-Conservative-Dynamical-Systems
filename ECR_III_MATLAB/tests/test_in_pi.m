function nfail = test_in_pi()
%TEST_IN_PI  Eq. (3) membership test.
nfail = 0;
S = sys_lorenz_real(struct('refine', false));

P = [S.pnom, S.pnom + S.dpmax*0.99, S.pnom - S.dpmax*0.99, ...
     S.pnom + S.dpmax*1.01, [NaN;0;0]];
tf = ecr_in_pi(P, S);
nfail = nfail + check(isequal(tf, [true true true false false]), ...
                      'ecr_in_pi boundary behaviour');

% each parameter is tested separately
P2 = repmat(S.pnom,1,3);
P2(1,1) = S.pnom(1) + 2*S.dpmax(1);
P2(2,2) = S.pnom(2) - 2*S.dpmax(2);
P2(3,3) = S.pnom(3) + 2*S.dpmax(3);
nfail = nfail + check(~any(ecr_in_pi(P2,S)), 'ecr_in_pi per-parameter test');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
