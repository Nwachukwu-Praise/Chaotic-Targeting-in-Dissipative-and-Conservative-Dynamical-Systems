function [Pout, ok, info] = ecr_level_eval(L, S, Z, opt)
%ECR_LEVEL_EVAL  Response of one control region to a batch of states.
%
%   [P,OK,INFO] = ECR_LEVEL_EVAL(L,S,Z,OPT)
%
%   For every column of Z the cluster of region L with the smallest NMD
%   (Eq. 8) is selected, its network NN_ij is evaluated and the resulting
%   parameter vector is accepted if it lies inside Pi (Eq. 3).  OK(k) is
%   therefore the practical membership test "z_k belongs to this control
%   region" used by Definition 3 / Fig. 3 of the paper.
%
%   OPT.nFallback > 1 lets the search continue with the next-nearest cluster
%   when the nearest one returns a parameter outside Pi (an addition to the
%   paper - set OPT.nFallback = 1 for the literal algorithm).
%   OPT.nmdGate = g rejects a cluster whose NMD exceeds g times the 95th
%   percentile NMD of its own training data (Inf disables the gate).
%
%   INFO.cluster  index of the cluster used (0 = none)
%   INFO.nmd      NMD to that cluster (Inf if none)

if nargin < 4 || isempty(opt), opt = ecr_default_options(); end

M    = size(Z,2);
nc   = numel(L.CL);
Pout = repmat(S.pnom, 1, M);
ok   = false(1, M);
info = struct('cluster', zeros(1,M), 'nmd', inf(1,M));
if nc == 0 || M == 0, return, end

Dm = ecr_nmd(L.CL, Z);
if nc > 1
    [Ds, Ord] = sort(Dm, 1);
else
    Ds = Dm;  Ord = ones(1, M);
end

nTry    = min(max(1, opt.nFallback), nc);
if ~L.useNMD, nTry = 1; end
pending = true(1, M);

for t = 1:nTry
    if ~any(pending), break, end
    jrow = Ord(t, :);
    for j = unique(jrow(pending))
        cols = find(pending & (jrow == j));
        if isempty(cols), continue, end
        if L.useNMD && isfinite(opt.nmdGate)
            lim  = opt.nmdGate * max(L.CL(j).nmd95, 1e-12);
            cols = cols(Ds(t,cols) <= lim);
            if isempty(cols), continue, end
        end
        dp   = rbf_eval(L.nets{j}, Z(:,cols));
        good = all(abs(dp) <= 1 + 1e-9, 1) & all(isfinite(dp), 1);
        if any(good)
            gc = cols(good);
            Pout(:,gc) = repmat(S.pnom,1,numel(gc)) + ...
                         repmat(S.dpmax,1,numel(gc)) .* dp(:,good);
            ok(gc) = true;
            info.cluster(gc) = j;
            info.nmd(gc)     = Ds(t,gc);
            pending(gc) = false;
        end
    end
end
end
