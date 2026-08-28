function S = sys_lorenz_delay(opts)
%SYS_LORENZ_DELAY  Lorenz system observed in delay coordinates on a Poincare
%                  surface of section (Table 2 of the paper).
%
%   S = SYS_LORENZ_DELAY()
%   S = SYS_LORENZ_DELAY(OPTS)
%
%   Only a single scalar measurement s(t) of the Lorenz flow is available.
%   At the n-th crossing of the surface of section y = y_PSS = 8.4853 the
%   controller observes the delay vector
%
%       z_n = [ s(t_n) ; s(t_n - T) ; s(t_n - 2T) ] ,   T = 100 ms
%
%   which is the "delay coord." case of Table 2.  Because the delay vector
%   does not determine the underlying flow state, the ECR machinery works on
%   a *hidden* state
%
%       h = [ x ; y ; z ; s(t_n-ds) ; s(t_n-2ds) ; ... ; s(t_n-L*ds) ]
%
%   with ds = T/nsub and L = nsub*(m-1), so that the delayed samples needed
%   for the observation are stored exactly on the sample grid.  Only S.obs(h)
%   is ever shown to the controller; h is used solely to propagate the system
%   and to allow the region tests of Definitions 1-3 to be evaluated.
%
%   Control parameters:  p = [sigma; rho; beta],  p_nom = [10; 28; 8/3],
%   dp_max = [0.30; 0.84; 0.08]  (Table 2).
%
%   OPTS fields (all optional)
%     meas      handle @(u)->s, default @(u) u(1)   (the x coordinate)
%     T         delay time, default 0.1  (= 100 ms)
%     m         embedding dimension, default 3
%     nsub      samples per delay, default 10
%     dt        RK4 step, default 1e-3
%     sec_dir   crossing direction, default -1 (see SYS_LORENZ_REAL)
%     delta     S0 radius, default 0.30
%     interp    interpolation for off-grid samples, default 'linear'
%
%   NOTE ON THE TARGET.  Table 2 quotes z* = [-3.0988 -3.2911 -3.8340] and
%   an RMS of 2.675 for the delay-coordinate case.  Those numbers cannot be
%   reproduced with s = x (or y, or z) on the y = 8.4853 section - the paper
%   does not state which signal was measured, nor how it was scaled.  This
%   file therefore computes the target itself: it takes the period-1 unstable
%   periodic orbit found in real coordinates (which *does* reproduce the
%   quoted [14.2387 39.7934] to three decimals) and evaluates its delay
%   vector.  See docs/METHOD_NOTES.md.

if nargin < 1, opts = struct(); end
opts = setdef(opts, 'meas',    @(U) U(1,:));   % vectorised: 3 x M -> 1 x M
opts = setdef(opts, 'T',       0.1);
opts = setdef(opts, 'm',       3);
opts = setdef(opts, 'nsub',    10);
opts = setdef(opts, 'dt',      1e-3);
opts = setdef(opts, 'tmax',    30);
opts = setdef(opts, 'tskip',   0.05);
opts = setdef(opts, 'sec_dir', -1);
opts = setdef(opts, 'delta',   0.30);
opts = setdef(opts, 'interp',  'linear');

S = ecr_system_template();

S.name   = 'lorenz_delay';
S.N      = opts.m;
S.r      = 3;
S.pnom   = [10; 28; 8/3];
S.dpmax  = [0.30; 0.84; 0.08];
S.delta  = opts.delta;
S.pnames = {'sigma','rho','beta'};
S.type   = 'delay';

cfg.dt      = opts.dt;
cfg.tmax    = opts.tmax;
cfg.tskip   = opts.tskip;
cfg.sec_idx = 2;
cfg.sec_val = sqrt(S.pnom(3)*(S.pnom(2)-1));    % 8.4853
cfg.sec_dir = opts.sec_dir;
cfg.meas    = opts.meas;
cfg.nref    = 8;
cfg.T       = opts.T;
cfg.m       = opts.m;
cfg.nsub    = opts.nsub;
cfg.ds      = opts.T/opts.nsub;
cfg.L       = opts.nsub*(opts.m-1);
cfg.interp  = opts.interp;
S.cfg       = cfg;
S.nh        = 3 + cfg.L;

S.step   = @(H,P) step_delay(H, P, cfg);
S.obs    = @(H) obs_delay(H, cfg);
S.hfromz = [];                       % delay vector cannot be lifted exactly
S.init   = @(m) init_delay(m, S.pnom, cfg);

% ---- target: delay vector of the period-1 UPO ---------------------------
ropts.sec_dir = cfg.sec_dir;  ropts.dt = cfg.dt;  ropts.tmax = cfg.tmax;
ropts.tskip   = cfg.tskip;
R  = sys_lorenz_real(ropts);
h  = [R.hfromz(R.zstar); zeros(cfg.L,1)];
for k = 1:2                                     % two steps fill the buffer
    h = step_delay(h, S.pnom, cfg);
end
S.zstar = obs_delay(h, cfg);

% ---- RMS of the measured signal on the attractor ------------------------
H     = S.init(200);
Zs    = S.obs(H);
S.rms = sqrt(mean(Zs(1,:).^2));

S.notes = sprintf(['Lorenz, delay coordinates (m = %d, T = %g, section ' ...
    'y = %.4f dir %+d).  z* = [%s] computed from the period-1 UPO; the ' ...
    'paper quotes [-3.0988 -3.2911 -3.8340] for an unspecified measured ' ...
    'signal.  Measured-signal RMS here = %.3f (paper: 2.675).'], ...
    cfg.m, cfg.T, cfg.sec_val, cfg.sec_dir, ...
    strtrim(sprintf('%.4f ', S.zstar)), S.rms);
end

% =========================================================================
function Z = obs_delay(H, cfg)
%OBS_DELAY  z_n = [s(t_n) s(t_n-T) ... s(t_n-(m-1)T)] read off the buffer.
M = size(H,2);
Z = zeros(cfg.m, M);
Z(1,:) = cfg.meas(H(1:3,:));
for j = 1:cfg.m-1
    Z(j+1,:) = H(3 + j*cfg.nsub, :);
end
end

% =========================================================================
function [Hn, Zn] = step_delay(H, P, cfg)
%STEP_DELAY  Advance to the next section crossing and rebuild the delay
%            buffer by interpolating the recorded measurement history.
U0  = H(1:3,:);
BUF = H(4:end,:);
[U1, ~, ok, BUFn] = flow_poincare_step_vec(@lorenz_ode_vec, U0, P, cfg, BUF);
Hn = [U1; BUFn];
Hn(:, ~ok) = NaN;
Zn = obs_delay(Hn, cfg);
end

% =========================================================================
function H = init_delay(m, pnom, cfg)
%INIT_DELAY  m independent hidden states (flow state + delay buffer) taken
%   from the attractor.  Two section crossings are run first so that every
%   sample in the delay buffer comes from real recorded history.
P = repmat(pnom, 1, m);
U = repmat([1;1;20], 1, m) + 5*randn(3,m);
for k = 1:round(10/cfg.dt)                     % burn-in
    U = flow_rk4_step_vec(@lorenz_ode_vec, U, P, cfg.dt);
end
H = [U; zeros(cfg.L, m)];
for k = 1:2
    H = step_delay(H, P, cfg);
    bad = any(~isfinite(H),1);
    if any(bad)
        Ub = repmat([1;1;20],1,sum(bad)) + 5*randn(3,sum(bad));
        for j = 1:round(10/cfg.dt)
            Ub = flow_rk4_step_vec(@lorenz_ode_vec, Ub, P(:,bad), cfg.dt);
        end
        H(:,bad) = [Ub; zeros(cfg.L, sum(bad))];
    end
end
bad = any(~isfinite(H),1);
if any(bad)
    good = find(~bad, 1);
    H(:,bad) = repmat(H(:,good), 1, sum(bad));
end
end

% =========================================================================
function o = setdef(o, f, v)
if ~isfield(o, f) || isempty(o.(f)), o.(f) = v; end
end
