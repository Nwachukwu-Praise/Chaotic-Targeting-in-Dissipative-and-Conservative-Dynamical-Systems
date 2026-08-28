function nfail = test_lorenz_system()
%TEST_LORENZ_SYSTEM  The Poincare map, its fixed point and the delay embedding.
nfail = 0;
ecr_seed(1);

S = sys_lorenz_real();

% --- the target really is a fixed point of the Poincare map --------------
[~, Pz] = S.step(S.hfromz(S.zstar), S.pnom);
nfail = nfail + check(norm(Pz - S.zstar) < 1e-6, ...
    sprintf('||P(z*)-z*|| = %.2e', norm(Pz - S.zstar)));

% --- and it agrees with the value quoted in Table 2 ----------------------
nfail = nfail + check(norm(S.zstar - [14.2387; 39.7934]) < 0.05, ...
    sprintf('z* = [%.4f %.4f] vs Table 2 [14.2387 39.7934]', S.zstar(1), S.zstar(2)));

% --- section value equals sqrt(beta*(rho-1)) = 8.4853 --------------------
nfail = nfail + check(abs(S.cfg.sec_val - 8.4853) < 1e-3, 'y_PSS = 8.4853');

% --- the map is deterministic and obs/hfromz are consistent -------------
H = S.init(20);
[H1, Z1] = S.step(H, repmat(S.pnom,1,20));
[H2, Z2] = S.step(H, repmat(S.pnom,1,20));
nfail = nfail + check(isequal(H1,H2) && isequal(Z1,Z2), 'Poincare map is deterministic');
nfail = nfail + check(max(max(abs(S.obs(S.hfromz(Z1)) - Z1))) < 1e-12, ...
                      'obs(hfromz(z)) == z');

% --- the fixed point is a saddle: one unstable, one strongly stable dir.
h = 1e-4;  J = zeros(2);
for j = 1:2
    dz = zeros(2,1); dz(j) = h;
    [~, Pp] = S.step(S.hfromz(S.zstar+dz), S.pnom);
    [~, Pm] = S.step(S.hfromz(S.zstar-dz), S.pnom);
    J(:,j) = (Pp - Pm)/(2*h);
end
e = sort(abs(eig(J)), 'descend');
nfail = nfail + check(e(1) > 1 && e(2) < 1, ...
    sprintf('eigenvalues of the Poincare map |l| = [%.3g %.3g]', e(1), e(2)));

% --- the parameters really do move the map ------------------------------
sens = zeros(1,3);
for j = 1:3
    dp = zeros(3,1); dp(j) = S.dpmax(j);
    [~, Pp] = S.step(S.hfromz(S.zstar), S.pnom+dp);
    sens(j) = norm(Pp - S.zstar);
end
% (sigma is the weakest actuator: on this section it moves the image by far
%  less than delta, rho and beta are the effective ones)
nfail = nfail + check(sum(sens > S.delta) >= 2 && all(sens > 0), ...
    sprintf('dp_max moves the image by [%.2f %.2f %.2f], delta = %.2f', ...
            sens(1), sens(2), sens(3), S.delta));

% --- attractor statistics match Table 2 (RMS of the measured state) -----
Z = S.obs(S.init(300));
rx = sqrt(mean(Z(1,:).^2));
nfail = nfail + check(abs(rx - 13.43) < 1.0, ...
    sprintf('RMS(x) on the section = %.2f (Table 2: 13.43)', rx));

% --- delay coordinates: z* is a fixed point of the delay map -------------
SD = sys_lorenz_delay();
R  = sys_lorenz_real();
hd = [R.hfromz(R.zstar); zeros(SD.cfg.L,1)];
for k = 1:2, hd = SD.step(hd, SD.pnom); end
[~, zd] = SD.step(hd, SD.pnom);
nfail = nfail + check(norm(zd - SD.zstar) < 1e-4, ...
    sprintf('delay-coordinate fixed point residual %.2e', norm(zd - SD.zstar)));
nfail = nfail + check(abs(SD.zstar(1) - R.zstar(1)) < 1e-6, ...
    'first delay coordinate equals the measured x on the section');
end

function n = check(c, msg)
n = 0;
if ~c, n = 1; fprintf('   ! %s\n', msg); end
end
