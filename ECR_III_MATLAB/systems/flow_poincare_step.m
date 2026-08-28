function [u1, t1, trec, srec, ok] = flow_poincare_step(odefun, u0, p, cfg)
%FLOW_POINCARE_STEP  Integrate a flow up to its next Poincare section crossing.
%
%   [U1,T1,TREC,SREC,OK] = FLOW_POINCARE_STEP(ODEFUN,U0,P,CFG)
%
%   Integrates du/dt = ODEFUN(u,P) from U0 (which normally lies *on* the
%   surface of section) until the next crossing of the surface
%
%        u(CFG.sec_idx) = CFG.sec_val ,  crossed in direction CFG.sec_dir
%
%   Inputs
%     CFG.dt        integration step (fixed step RK4)
%     CFG.tmax      give up after this much time (OK = false)
%     CFG.tskip     minimum time before a crossing is accepted (prevents the
%                   detection of the crossing we start on)
%     CFG.meas      handle @(u) -> scalar measurement s(t) (may be [])
%     CFG.nref      number of bisection refinements of the crossing time (8)
%
%   Outputs
%     U1     state at the crossing            (nu x 1)
%     T1     elapsed time to the crossing
%     TREC   1 x K vector of recorded times (0 ... T1)
%     SREC   1 x K vector of recorded measurements (empty if CFG.meas is [])
%     OK     true if a crossing was found within CFG.tmax and the trajectory
%            stayed bounded
%
%   The record TREC/SREC is what makes delay-coordinate observation possible:
%   the delay vector at the crossing is obtained by interpolating SREC.

if ~isfield(cfg,'nref') || isempty(cfg.nref), cfg.nref = 8; end
haveMeas = isfield(cfg,'meas') && ~isempty(cfg.meas);

dt   = cfg.dt;
nmax = ceil(cfg.tmax/dt) + 1;

trec = zeros(1, nmax);
if haveMeas, srec = zeros(1, nmax); else, srec = []; end

u  = u0(:);
t  = 0;
k  = 1;
trec(1) = 0;
if haveMeas, srec(1) = cfg.meas(u); end

gprev = cfg.sec_dir * (u(cfg.sec_idx) - cfg.sec_val);
ok    = false;
u1    = u;  t1 = 0;

for n = 1:nmax-1
    uprev = u;
    u = flow_rk4_step(odefun, u, p, dt);
    t = t + dt;

    if any(~isfinite(u)) || norm(u(1:min(3,numel(u)))) > 1e6
        ok = false;
        trec = trec(1:k); if haveMeas, srec = srec(1:k); end
        return
    end

    k = k + 1;
    trec(k) = t;
    if haveMeas, srec(k) = cfg.meas(u); end

    g = cfg.sec_dir * (u(cfg.sec_idx) - cfg.sec_val);
    if t > cfg.tskip && gprev < 0 && g >= 0
        % ---- refine the crossing time -----------------------------------
        % Newton iteration on  g(h) = dir*(u_idx(uprev flowed by h) - val),
        % started from the linear interpolation of the bracketing step.
        % (Newton is used rather than bisection because the finite-difference
        %  Jacobians of the Poincare map are otherwise swamped by the
        %  crossing-time resolution error.)
        h  = dt * (-gprev)/(g - gprev);
        h  = min(max(h, 0), dt);
        uc = uprev;
        for j = 1:max(cfg.nref, 6)
            uc = flow_rk4_step(odefun, uprev, p, h);
            gc = cfg.sec_dir * (uc(cfg.sec_idx) - cfg.sec_val);
            duc = odefun(uc, p);
            dc  = cfg.sec_dir * duc(cfg.sec_idx);
            if abs(gc) < 1e-13 || abs(dc) < eps, break, end
            hn = h - gc/dc;
            if ~isfinite(hn) || hn < 0 || hn > dt + 1e-12, break, end
            if abs(hn - h) < 1e-15, h = hn; break, end
            h = hn;
        end
        uc = flow_rk4_step(odefun, uprev, p, h);
        tc = (t - dt) + h;
        uc(cfg.sec_idx) = cfg.sec_val;      % project exactly onto the section

        k = k + 1;
        trec(k) = tc;
        if haveMeas, srec(k) = cfg.meas(uc); end

        u1 = uc;  t1 = tc;  ok = true;
        trec = trec(1:k); if haveMeas, srec = srec(1:k); end
        return
    end
    gprev = g;
end

trec = trec(1:k); if haveMeas, srec = srec(1:k); end
end
