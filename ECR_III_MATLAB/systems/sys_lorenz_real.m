function S = sys_lorenz_real(opts)
%SYS_LORENZ_REAL  Lorenz system observed in real coordinates on a Poincare
%                 surface of section (Table 2 of the paper).
%
%   S = SYS_LORENZ_REAL()
%   S = SYS_LORENZ_REAL(OPTS)
%
%   The continuous Lorenz system
%
%       dx/dt = sigma*(y-x),  dy/dt = rho*x - y - x*z,  dz/dt = x*y - beta*z
%
%   is reduced to the discrete map of Eq. (1) by taking the Poincare surface
%   of section  y = y_PSS = 8.4853  (the y-coordinate of the C+ equilibrium,
%   sqrt(beta*(rho-1)) ).  The observed state is
%
%       z_n = [x_n ; z_n]   at the n-th crossing of the section,
%
%   the hidden state is the full [x;y;z] on the section, so the observed
%   state can be lifted back exactly (S.hfromz is available).
%
%   Control parameters:  p = [sigma; rho; beta]
%       p_nom  = [10 ; 28 ; 8/3]
%       dp_max = [0.30; 0.84; 0.08]      (Table 2)
%   Target (fixed point of the Poincare map, i.e. a period-1 UPO of the flow):
%       z* = [14.2387 ; 39.7934]         (Table 2, refined numerically here)
%   OGY / S0 region radius: delta = 0.30.
%
%   OPTS fields (all optional)
%     dt        RK4 step                                   (default 2e-3)
%     tmax      max flight time between crossings          (default 30)
%     tskip     dead time after a crossing                 (default 0.05)
%     sec_dir   -1 (default) = crossings with dy/dt < 0.  This is the branch
%               that reproduces the target quoted in Table 2 and an x-RMS of
%               13.5 (paper: 13.43).  With +1 the surface of section contains
%               the C+ equilibrium itself and the fixed-point search
%               degenerates onto it, so use -1 unless you know better.
%     zstar     override the target                        (default: refined)
%     refine    refine z* by Newton on the Poincare map    (default true)
%     delta     S0 radius                                  (default 0.30)

if nargin < 1, opts = struct(); end
opts = setdef(opts, 'dt',      2e-3);
opts = setdef(opts, 'tmax',    30);
opts = setdef(opts, 'tskip',   0.05);
opts = setdef(opts, 'sec_dir', -1);
opts = setdef(opts, 'refine',  true);
opts = setdef(opts, 'delta',   0.30);
opts = setdef(opts, 'zstar',   [14.2387; 39.7934]);

S = ecr_system_template();

S.name   = 'lorenz_real';
S.N      = 2;
S.nh     = 3;
S.r      = 3;
S.pnom   = [10; 28; 8/3];
S.dpmax  = [0.30; 0.84; 0.08];
S.delta  = opts.delta;
S.rms    = 13.43;                      % Table 2 (used to scale noise)
S.pnames = {'sigma','rho','beta'};
S.type   = 'poincare';

cfg.dt      = opts.dt;
cfg.tmax    = opts.tmax;
cfg.tskip   = opts.tskip;
cfg.sec_idx = 2;                       % section on the y-coordinate
cfg.sec_val = sqrt(S.pnom(3)*(S.pnom(2)-1));   % = 8.4853
cfg.sec_dir = opts.sec_dir;
cfg.meas    = [];
cfg.nref    = 10;
S.cfg       = cfg;

S.step   = @(H,P) step_lorenz(H, P, cfg);
S.obs    = @(H) H([1 3], :);
S.hfromz = @(Z) [Z(1,:); cfg.sec_val*ones(1,size(Z,2)); Z(2,:)];
S.init   = @(m) init_lorenz(m, S.pnom, cfg);

S.zstar  = opts.zstar(:);
if opts.refine
    zs = refine_fixed_point(S, opts.zstar(:));
    if ~isempty(zs), S.zstar = zs; end
end

S.notes = sprintf(['Lorenz, real coordinates.  Section y = %.4f (dir %+d), ' ...
                   'RK4 dt = %g.  z* = [%.4f %.4f] (paper: [14.2387 39.7934]).'], ...
                   cfg.sec_val, cfg.sec_dir, cfg.dt, S.zstar(1), S.zstar(2));
end

% =========================================================================
function [Hn, Zn] = step_lorenz(H, P, cfg)
[U1, ~, ok] = flow_poincare_step_vec(@lorenz_ode_vec, H, P, cfg);
Hn = U1;
Hn(:, ~ok) = NaN;                       % lost the section -> invalid sample
Zn = Hn([1 3], :);
end

% =========================================================================
function H = init_lorenz(m, pnom, cfg)
%INIT_LORENZ  m independent states on the attractor, taken as the first
%   section crossing after a burn-in from random initial conditions.
P = repmat(pnom, 1, m);
U = repmat([1;1;20], 1, m) + 5*randn(3,m);
for k = 1:round(10/cfg.dt)                      % burn-in, 10 time units
    U = flow_rk4_step_vec(@lorenz_ode_vec, U, P, cfg.dt);
end
[H, ~, ok] = flow_poincare_step_vec(@lorenz_ode_vec, U, P, cfg);
tries = 0;
while ~all(ok) && tries < 5                     % retry the odd failure
    bad = find(~ok);
    Ub  = repmat([1;1;20], 1, numel(bad)) + 5*randn(3,numel(bad));
    Pb  = repmat(pnom, 1, numel(bad));
    for k = 1:round(10/cfg.dt)
        Ub = flow_rk4_step_vec(@lorenz_ode_vec, Ub, Pb, cfg.dt);
    end
    [Hb, ~, okb] = flow_poincare_step_vec(@lorenz_ode_vec, Ub, Pb, cfg);
    H(:,bad) = Hb;  ok(bad) = okb;
    tries = tries + 1;
end
H(:, ~ok) = repmat(H(:, find(ok, 1)), 1, sum(~ok));
end

% =========================================================================
function zs = refine_fixed_point(S, z0)
%REFINE_FIXED_POINT  Newton iteration on  F(z) = P(z) - z  using finite
%   differences of the Poincare map at nominal parameters.
zs = [];
z  = z0;
for it = 1:12
    [~, Pz] = S.step(S.hfromz(z), S.pnom);
    if any(~isfinite(Pz)), return, end
    F = Pz - z;
    if norm(F) < 1e-10, break, end
    h = 1e-5;
    J = zeros(2,2);
    for j = 1:2
        dz = zeros(2,1); dz(j) = h;
        [~, Pp] = S.step(S.hfromz(z+dz), S.pnom);
        [~, Pm] = S.step(S.hfromz(z-dz), S.pnom);
        if any(~isfinite(Pp)) || any(~isfinite(Pm)), return, end
        J(:,j) = (Pp - Pm)/(2*h);
    end
    A = (J - eye(2));
    if rcond(A) < 1e-12, return, end
    step = -A\F;
    if norm(step) > 5, step = 5*step/norm(step); end
    z = z + step;
end
if norm(F) < 1e-6
    zs = z;
end
end

% =========================================================================
function o = setdef(o, f, v)
if ~isfield(o, f) || isempty(o.(f)), o.(f) = v; end
end
