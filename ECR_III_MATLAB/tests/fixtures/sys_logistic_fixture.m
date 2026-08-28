function S = sys_logistic_fixture()
%SYS_LOGISTIC_FIXTURE  Logistic map  x_{n+1} = p*x_n*(1-x_n)   (Table 2 of the paper)
%
%   z* = 0.7435, p_nom = 3.9, dp_max = 0.10, x_rms = 0.6608, delta = 0.01.
%
%   Iplikci & Denizhan, "An improved neural network based targeting method for
%   chaotic dynamics", Chaos Solitons & Fractals 17 (2003) 523-529, Table 2.

S = ecr_system_template();

S.name   = 'logistic';
S.N      = 1;
S.nh     = 1;
S.r      = 1;
S.pnom   = 3.9;
S.dpmax  = 0.10;
S.delta  = 0.01;
S.rms    = 0.6608;
S.pnames = {'p'};
S.type   = 'map';

% Unstable fixed point of the nominal map: x* = 1 - 1/p_nom = 0.74359
S.zstar  = 1 - 1/S.pnom;

S.step   = @stepfun;
S.obs    = @(H) H;
S.hfromz = @(Z) Z;
S.init   = @initfun;

S.notes  = sprintf(['Logistic map, 1 state / 1 control parameter.  ' ...
                    'Paper quotes x* = 0.7435; the exact fixed point ' ...
                    '1-1/p_nom = %.5f is used here.'], S.zstar);
end

% -------------------------------------------------------------------------
function [Hn, Zn] = stepfun(H, P)
% One iteration of the logistic map, vectorised over columns.
Hn = P(1,:) .* H(1,:) .* (1 - H(1,:));
Hn = min(max(Hn, -10), 10);          % guard against escape to infinity
Zn = Hn;
end

% -------------------------------------------------------------------------
function H = initfun(m)
% m random points on the attractor (transient discarded).
H = 0.05 + 0.9*rand(1, m);
for k = 1:200
    H = 3.9 * H .* (1 - H);
end
end
