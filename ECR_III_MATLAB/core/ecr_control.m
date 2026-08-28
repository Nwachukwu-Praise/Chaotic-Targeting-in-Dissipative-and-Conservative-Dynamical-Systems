function [Pout, info] = ecr_control(model, Z)
%ECR_CONTROL  On-line ECR controller (Fig. 4 of the paper).
%
%   [P,INFO] = ECR_CONTROL(MODEL,Z) returns the parameter vector applied to
%   every column of Z.  States that are not recognised by any control region
%   get the nominal parameters ("If NN_ij produces a p_n not in Pi, no
%   targeting is applied").
%
%   ECR-II   the state is fed sequentially to NN_0, NN_1, ... and the region
%            with the smallest index whose network returns an acceptable
%            parameter change is regarded as the region of the state
%            (Section 3.2).
%   ECR-III  the region is determined analytically: the cluster C_ij with the
%            minimum Normalised Mahalanobis Distance (Eq. 8) over all regions
%            is taken as the region of the state and its network NN_ij is
%            activated (Section 3.3).
%
%   INFO.region   region index used (-1 = none, nominal parameters applied)
%   INFO.cluster  cluster index inside that region
%   INFO.nmd      NMD to the selected cluster

S   = model.S;
opt = model.opt;
M   = size(Z,2);

Pout = repmat(S.pnom, 1, M);
info = struct('region', -ones(1,M), 'cluster', zeros(1,M), 'nmd', inf(1,M));

if strcmpi(model.variant, 'ECR-II')
    pending = true(1,M);
    for k = 1:numel(model.levels)
        if ~any(pending), break, end
        cols = find(pending);
        [Pk, okk, ik] = ecr_level_eval(model.levels(k), S, Z(:,cols), opt);
        if any(okk)
            gc = cols(okk);
            Pout(:,gc)       = Pk(:,okk);
            info.region(gc)  = model.levels(k).index;
            info.cluster(gc) = ik.cluster(okk);
            info.nmd(gc)     = ik.nmd(okk);
            pending(gc) = false;
        end
    end
    return
end

% ------------------------------- ECR-III ---------------------------------
F  = model.flat;
nc = numel(F.CL);
if nc == 0 || M == 0, return, end

Dm = ecr_nmd(F.CL, Z);                 % Eq. (8), always the true NMD
key = Dm;                              % what the search is ordered by
if isfield(opt,'selection') && strcmpi(opt.selection,'level')
    % order by region index first, by NMD inside a region: the state is
    % steered by the *lowest* region that recognises it (ECR-II ordering)
    % while the cluster inside that region is still chosen by Eq. (8).
    key = repmat(F.level(:), 1, M)*1e6 + min(Dm, 1e5);
elseif opt.levelPriority
    % break near-ties in favour of the lower region index
    key = Dm + 1e-9*repmat(F.level(:), 1, M);
end
if nc > 1
    [~, Ord] = sort(key, 1);
else
    Ord = ones(1,M);
end

nTry    = min(max(1, opt.nFallback), nc);
pending = true(1, M);
for t = 1:nTry
    if ~any(pending), break, end
    jrow = Ord(t,:);
    for j = unique(jrow(pending))
        cols = find(pending & (jrow == j));
        if isempty(cols), continue, end
        if isfinite(opt.nmdGate)
            lim  = opt.nmdGate * max(F.CL(j).nmd95, 1e-12);
            cols = cols(Dm(j,cols) <= lim);
            if isempty(cols), continue, end
        end
        dp   = rbf_eval(F.nets{j}, Z(:,cols));
        good = all(abs(dp) <= 1 + 1e-9, 1) & all(isfinite(dp), 1);
        if any(good)
            gc = cols(good);
            Pout(:,gc) = repmat(S.pnom,1,numel(gc)) + ...
                         repmat(S.dpmax,1,numel(gc)) .* dp(:,good);
            info.region(gc)  = F.level(j);
            info.cluster(gc) = F.jj(j);
            info.nmd(gc)     = Dm(j,gc);
            pending(gc) = false;
        end
    end
end
end
