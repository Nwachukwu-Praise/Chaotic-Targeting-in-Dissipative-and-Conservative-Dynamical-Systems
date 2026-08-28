function [U1, T1, ok, BUFn] = flow_poincare_step_vec(odefun, U0, P, cfg, BUF)
%FLOW_POINCARE_STEP_VEC  Poincare map of many trajectories in lock-step.
%
%   [U1,T1,OK]      = FLOW_POINCARE_STEP_VEC(ODEFUN,U0,P,CFG)
%   [U1,T1,OK,BUFN] = FLOW_POINCARE_STEP_VEC(ODEFUN,U0,P,CFG,BUF)
%
%   Integrates M copies of the flow (columns of U0, each with its own
%   parameter column of P) until each one crosses the surface of section
%
%       u(CFG.sec_idx) = CFG.sec_val   in direction CFG.sec_dir,
%
%   using one vectorised RK4 step for all copies at a time.  A copy is
%   dropped from the active set as soon as it has crossed, so the cost is set
%   by the *longest* flight rather than by their sum.  This is what makes the
%   Lorenz experiments of the paper tractable in MATLAB/Octave.
%
%   The crossing time is located by Newton's method on the section function
%   (a bisection is not accurate enough: it leaves a crossing-time error that
%   swamps finite-difference Jacobians of the Poincare map).
%
%   DELAY COORDINATES.  If CFG.meas is a (vectorised) measurement handle and
%   BUF (L x M) holds the samples s(t_prev - j*CFG.ds), j = 1..L, the routine
%   also returns BUFN, the same buffer re-sampled at the new crossing time,
%   obtained by interpolating the measurement history recorded along the
%   flight.  CFG.L is L and CFG.interp the interpolation method.
%
%   Outputs: U1 (nu x M) crossing states, T1 (1 x M) flight times, OK (1 x M).

M   = size(U0,2);
nu  = size(U0,1);
dt  = cfg.dt;
useMeas = isfield(cfg,'meas') && ~isempty(cfg.meas) && nargin >= 5;

U1 = NaN(nu, M);
T1 = NaN(1, M);
ok = false(1, M);
if useMeas
    L    = cfg.L;
    BUFn = NaN(L, M);
    W    = ceil(L*cfg.ds/dt) + 4;
    Rt   = NaN(1, W);        % ring buffer of sample times (shared)
    Rs   = NaN(W, M);        % ring buffer of measurements (per column)
    ptr  = 0;
    [Rt, Rs, ptr] = ring_push(Rt, Rs, ptr, 0, cfg.meas(U0));
else
    BUFn = [];
end

act  = 1:M;                  % original indices of the active columns
U    = U0;
Pa   = P;
t    = 0;
gprev = cfg.sec_dir*(U(cfg.sec_idx,:) - cfg.sec_val);
nmax  = ceil(cfg.tmax/dt);

for n = 1:nmax
    Uprev = U;
    U = flow_rk4_step_vec(odefun, U, Pa, dt);
    t = t + dt;

    bad = any(~isfinite(U),1) | (max(abs(U),[],1) > 1e6);
    g   = cfg.sec_dir*(U(cfg.sec_idx,:) - cfg.sec_val);

    if useMeas
        [Rt, Rs, ptr] = ring_push(Rt, Rs, ptr, t, cfg.meas(U));
    end

    cross = (t > cfg.tskip) & (gprev < 0) & (g >= 0) & ~bad;

    if any(cross)
        c  = find(cross);
        Up = Uprev(:,c);  Pp = Pa(:,c);
        gp = gprev(c);    gn = g(c);
        h  = dt*(-gp)./(gn - gp);
        h  = min(max(h, 0), dt);
        for it = 1:max(6, cfg.nref)
            Uc  = flow_rk4_step_vec(odefun, Up, Pp, h);
            gc  = cfg.sec_dir*(Uc(cfg.sec_idx,:) - cfg.sec_val);
            dUc = odefun(Uc, Pp);
            dc  = cfg.sec_dir*dUc(cfg.sec_idx,:);
            hn  = h - gc./dc;
            good = isfinite(hn) & hn >= 0 & hn <= dt;
            h(good) = hn(good);
            if max(abs(gc)) < 1e-13, break, end
        end
        Uc = flow_rk4_step_vec(odefun, Up, Pp, h);
        Uc(cfg.sec_idx,:) = cfg.sec_val;
        tc = (t - dt) + h;

        U1(:,act(c)) = Uc;
        T1(act(c))   = tc;
        ok(act(c))   = true;

        if useMeas
            [Rtc, ord] = ring_times(Rt, ptr);
            sc = cfg.meas(Uc);
            for q = 1:numel(c)
                col  = c(q);
                told = -cfg.ds*(L:-1:1);
                tall = [told, Rtc, tc(q)];
                sall = [BUF(L:-1:1, act(col)).', Rs(ord, col).', sc(q)];
                keep = isfinite(tall) & isfinite(sall);
                tall = tall(keep);  sall = sall(keep);
                [tall, iu] = unique(tall);
                sall = sall(iu);
                tq   = tc(q) - cfg.ds*(1:L);
                BUFn(:,act(col)) = interp1(tall, sall, tq, cfg.interp, 'extrap').';
            end
        end
    end

    drop = cross | bad;
    if any(drop)
        keep = ~drop;
        U = U(:,keep);  Pa = Pa(:,keep);  gprev = g(keep);  act = act(keep);
        if useMeas, Rs = Rs(:,keep); end
        if isempty(act), return, end
    else
        gprev = g;
    end
end
end

% =========================================================================
function [Rt, Rs, ptr] = ring_push(Rt, Rs, ptr, t, s)
W   = numel(Rt);
ptr = mod(ptr, W) + 1;
Rt(ptr)   = t;
Rs(ptr,:) = s;
end

% =========================================================================
function [Rtc, ord] = ring_times(Rt, ptr)
W   = numel(Rt);
ord = [ptr+1:W, 1:ptr];
Rtc = Rt(ord);
m   = isfinite(Rtc);
ord = ord(m);
Rtc = Rtc(m);
end
