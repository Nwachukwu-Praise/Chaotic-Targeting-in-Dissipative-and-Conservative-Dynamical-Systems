function B = ecr_reaching_time(S, ctrl, opt, H0)
%ECR_REACHING_TIME  Average reaching time benchmark (Table 1 of the paper).
%
%   B = ECR_REACHING_TIME(S,CTRL,OPT)
%   B = ECR_REACHING_TIME(S,CTRL,OPT,H0)
%
%   "The average reaching time is calculated as the average of the time it
%   takes the system to enter the OGY region starting from many different
%   initial conditions" (Section 4).  OPT.nTrials initial conditions are used
%   (pass H0 to reuse exactly the same set for every method, which is what
%   makes the comparison fair).  All trials are propagated simultaneously.
%
%   B.reach      1 x nTrials  steps to first entry (NaN = never within maxSteps)
%   B.mean       mean over the successful trials
%   B.meanCens   mean with unsuccessful trials counted as OPT.maxSteps
%   B.median, B.std, B.success, B.hold  (retention after capture)
%   B.pUsed      fraction of steps at which targeting was actually applied
%   B.regionHist histogram of the control regions used

if nargin < 3 || isempty(opt), opt = ecr_default_options(); end
if nargin < 4 || isempty(H0),  H0 = S.init(opt.nTrials); end

nT   = size(H0,2);
H    = H0;
reach = NaN(1,nT);
holdN = zeros(1,nT);  holdIn = zeros(1,nT);
alive = true(1,nT);
ctrlRun = zeros(1,nT);          % consecutive controlled steps without capture
nAct  = 0;  nStep = 0;
regionHist = [];

for k = 1:opt.maxSteps + opt.holdSteps
    Z  = S.obs(H);
    Zm = Z + opt.noise*S.rms*randn(S.N, nT);
    d  = sqrt(sum((Zm - repmat(S.zstar,1,nT)).^2, 1));
    inS0 = d < S.delta;

    new = alive & inS0 & isnan(reach);
    reach(new) = k - 1;

    % retention bookkeeping for trials already captured
    cap = alive & ~isnan(reach) & ~new;
    holdN(cap)  = holdN(cap) + 1;
    holdIn(cap) = holdIn(cap) + inS0(cap);

    done = ~isnan(reach) & (holdN >= opt.holdSteps);
    alive = alive & ~done;
    if k >= opt.maxSteps
        alive = alive & ~isnan(reach);         % give up on the stragglers
    end
    if ~any(alive), break, end

    cols = find(alive);
    if isempty(ctrl)
        P = repmat(S.pnom, 1, numel(cols));
    else
        [P, inf1] = ctrl(Zm(:,cols));
        if isstruct(inf1) && isfield(inf1,'region')
            regionHist = [regionHist, inf1.region]; %#ok<AGROW>
            nAct  = nAct  + sum(inf1.region >= 0);
            % escape from a controlled periodic orbit (see OPT.escapeAfter)
            if opt.escapeAfter > 0
                act = inf1.region >= 0;
                nc  = ~isnan(reach(cols));
                ctrlRun(cols(act & ~nc)) = ctrlRun(cols(act & ~nc)) + 1;
                ctrlRun(cols(~act | nc)) = 0;
                esc = cols(ctrlRun(cols) >= opt.escapeAfter);
                if ~isempty(esc)
                    [~, loc] = ismember(esc, cols);
                    P(:,loc) = repmat(S.pnom, 1, numel(esc));
                    ctrlRun(esc) = 0;
                end
            end
        end
    end
    nStep = nStep + numel(cols);

    Hn = S.step(H(:,cols), P);
    bad = any(~isfinite(Hn),1);
    if any(bad)
        alive(cols(bad)) = false;
    end
    H(:,cols) = Hn;
end

ok = ~isnan(reach);
B.reach    = reach;
B.mean     = mean(reach(ok));
r2 = reach; r2(~ok) = opt.maxSteps;          % censored at maxSteps
B.meanCens = mean(r2);
B.median  = median(reach(ok));
B.std     = std(reach(ok));
B.success = mean(ok);
B.hold    = mean(holdIn(holdN>0) ./ holdN(holdN>0));
B.pUsed   = nAct/max(nStep,1);
B.nTrials = nT;
if ~isempty(regionHist)
    u = unique(regionHist);
    c = zeros(1, numel(u));
    for k = 1:numel(u), c(k) = sum(regionHist == u(k)); end
    B.regionHist = [u; c];
else
    B.regionHist = [];
end
end
