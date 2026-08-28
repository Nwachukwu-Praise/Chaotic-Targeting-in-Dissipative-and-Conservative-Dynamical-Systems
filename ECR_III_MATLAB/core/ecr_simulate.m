function R = ecr_simulate(S, ctrl, h0, opt)
%ECR_SIMULATE  Closed-loop run of one trajectory (Fig. 4 control loop).
%
%   R = ECR_SIMULATE(S,CTRL,H0,OPT)
%
%   CTRL is a function handle @(z) -> p (use [] for the uncontrolled system,
%   which is the reference used to judge how much targeting gains).  The
%   controller only ever sees the *measured* state, i.e. the observed state
%   plus OPT.noise * S.rms * randn - the noisy case of Table 1.
%
%   R.Z      N x n   observed states
%   R.Zm     N x n   measured states (what the controller saw)
%   R.P      r x n   applied parameters
%   R.reg    1 x n   control region used at each step (-1 = none)
%   R.reach  number of steps needed to enter the delta-ball of z* (NaN if the
%            trajectory never reached it within OPT.maxSteps)
%   R.hold   fraction of the OPT.holdSteps steps after capture that stayed
%            inside the delta-ball

if nargin < 4 || isempty(opt), opt = ecr_default_options(); end

h  = h0(:);
n  = 0;
nmax = opt.maxSteps + opt.holdSteps;
Z  = zeros(S.N, nmax);  Zm = Z;  P = zeros(S.r, nmax);  reg = -ones(1,nmax);
inS0 = false(1,nmax);
reach = NaN;
ctrlRun = 0;

for k = 1:nmax
    z  = S.obs(h);
    zm = z + opt.noise*S.rms*randn(S.N,1);
    n  = k;
    Z(:,k) = z;  Zm(:,k) = zm;
    inS0(k) = norm(zm - S.zstar) < S.delta;

    if isnan(reach) && inS0(k)
        reach = k - 1;                       % steps taken to get there
    end
    if ~isnan(reach) && k >= reach + 1 + opt.holdSteps
        break
    end
    if isnan(reach) && k >= opt.maxSteps
        break
    end

    if isempty(ctrl)
        p = S.pnom;
    else
        [p, inf1] = ctrl(zm);
        if isstruct(inf1) && isfield(inf1,'region'), reg(k) = inf1.region(1); end
        % escape from a controlled periodic orbit (see OPT.escapeAfter)
        if opt.escapeAfter > 0 && isnan(reach)
            if reg(k) >= 0, ctrlRun = ctrlRun + 1; else, ctrlRun = 0; end
            if ctrlRun >= opt.escapeAfter
                p = S.pnom;  reg(k) = -1;  ctrlRun = 0;
            end
        end
    end
    P(:,k) = p;

    h = S.step(h, p);
    if any(~isfinite(h))
        break
    end
end

R.Z = Z(:,1:n);  R.Zm = Zm(:,1:n);  R.P = P(:,1:n);  R.reg = reg(1:n);
R.inS0 = inS0(1:n);
R.reach = reach;
if ~isnan(reach)
    seg = inS0(min(reach+2,n):n);
    if isempty(seg), R.hold = NaN; else, R.hold = mean(seg); end
else
    R.hold = NaN;
end
end
