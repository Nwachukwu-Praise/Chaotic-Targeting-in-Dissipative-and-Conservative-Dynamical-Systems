function model = ecr_train(S, D, opt, variant)
%ECR_TRAIN  Build the extended control regions and their neural networks.
%
%   MODEL = ECR_TRAIN(S,D,OPT,VARIANT) with VARIANT = 'ECR-III' (default) or
%   'ECR-II'.
%
%   The control regions are built exactly in the order of Definitions 1 and 3
%   of the paper, using nothing but the recorded triples (z_n, p_n, z_{n+1}):
%
%   Level 0 (Definition 1, Eq. 4)
%       z_n in S0 = T0   <=>   ||z_n - z*|| < delta  and  ||G(z_n,p_n) - z*|| < delta
%       The recorded p_n is the evidence that such a p in Pi exists, and it is
%       used as the training target: NN_0 is the local (fine) controller.
%
%   Level i (Definition 3, Eq. 7)
%       z_n in T_i   <=>   z_n not in T_0..T_{i-1}   and   G(z_n,p_n) in T_{i-1}
%       Membership of the successor in T_{i-1} is decided exactly as in Fig. 3
%       of the paper: the state is fed to the already trained level-(i-1)
%       model and accepted if that model returns a parameter vector inside Pi.
%       No knowledge of the system equations is used anywhere.
%
%   ECR-II  : one RBF network per control region, regions searched in
%             ascending index order on line (Section 3.2).
%   ECR-III : each control region is first split into clusters C_ij by the
%             simple radius clustering of Section 3.3; every cluster gets its
%             own network NN_ij and is described by its mean mu_ij and its
%             normalised covariance matrix N_ij, so that the region of a
%             current state can be found analytically from the NMD of Eq. (8).
%
%   MODEL fields
%     variant, S, opt, zstar, delta
%     levels(k)   .index   region index i (0-based)
%                 .radius  clustering radius r (NaN for ECR-II)
%                 .CL      cluster statistics (ECR_CLUSTER_STATS)
%                 .nets    1 x nc cell array of RBF networks
%                 .useNMD  true for ECR-III
%                 .nData   number of training points of the region
%                 .rmse    per cluster training error (normalised units)
%     info        counts, training time, network count

if nargin < 3 || isempty(opt),     opt = ecr_default_options(); end
if nargin < 4 || isempty(variant), variant = 'ECR-III'; end
useCluster = strcmpi(variant, 'ECR-III');

t0 = tic;
Z  = D.Z;  Zn = D.Zn;  P = D.P;
M  = size(Z,2);

% normalised parameter perturbation used as network target: dp in [-1,1]^r
Y  = (P - repmat(S.pnom,1,M)) ./ repmat(S.dpmax,1,M);

dz  = sqrt(sum((Z  - repmat(S.zstar,1,M)).^2, 1));
dzn = sqrt(sum((Zn - repmat(S.zstar,1,M)).^2, 1));

model = struct('variant', variant, 'S', S, 'opt', opt, ...
               'zstar', S.zstar, 'delta', S.delta, 'levels', [], 'info', []);

% ---------------------------------------------------------------- level 0
sel = (dz < S.delta) & (dzn < S.delta);              % Definition 1 (Eq. 4)
if sum(sel) < opt.rbfMinPts
    error('ecr_train:noS0', ...
        ['Only %d training points fall inside S0.  Increase opt.nData, ' ...
         'opt.enrichFactor or delta.'], sum(sel));
end
% Training pairs for the local controller NN_0.  Definition 1 only asks for
% the successor to stay inside the delta-ball; among those pairs the ones
% that land *well* inside it teach a much sharper local controller, so they
% are preferred whenever there are enough of them (OPT.targetTighten, an
% addition to the paper - set it to 1 for the literal reading).
selT = sel & (dzn < opt.targetTighten*S.delta);
if sum(selT) >= max(opt.rbfMinPts, 0.25*sum(sel))
    sel0 = selT;
else
    sel0 = sel;
end
L = build_level(Z(:,sel0), Y(:,sel0), opt, useCluster, 0, find(sel0));
levels = L;
claimed = false(1,M);
[~, okm] = ecr_level_eval(levels(1), S, Z, opt);
claimed = claimed | (okm & (dz < S.delta));

if opt.verbose
    fprintf('[train %s] T0: %d points, %d cluster(s), claims %d states\n', ...
            variant, sum(sel), numel(levels(1).CL), sum(claimed));
end

% ------------------------------------------------------------- levels 1..K
for i = 1:opt.Kmax
    [~, prevOK] = ecr_level_eval(levels(end), S, Zn, opt);   % Fig. 3
    if i == 1
        % the successor must be in T0: also honour the delta test of Def. 1
        prevOK = prevOK & (dzn < S.delta);
    end
    cand = prevOK & ~claimed;
    if sum(cand) < max(opt.minLevelPts, opt.rbfMinPts)
        if opt.verbose
            fprintf('[train %s] stopping at K = %d (only %d candidates for T%d)\n', ...
                    variant, i-1, sum(cand), i);
        end
        break
    end
    L = build_level(Z(:,cand), Y(:,cand), opt, useCluster, i, find(cand));
    levels(end+1) = L; %#ok<AGROW>

    [~, okm] = ecr_level_eval(levels(end), S, Z, opt);
    claimed  = claimed | okm;

    if opt.verbose
        fprintf('[train %s] T%d: %d points, %d cluster(s), r = %.4g, claimed total %d/%d\n', ...
                variant, i, sum(cand), numel(L.CL), L.radius, sum(claimed), M);
    end
end

model.levels = levels;
model.flat   = ecr_flatten(levels);

nNets = 0;
for k = 1:numel(levels), nNets = nNets + numel(levels(k).nets); end
model.info = struct('trainTime', toc(t0), 'nLevels', numel(levels), ...
                    'nNets', nNets, 'nData', M, ...
                    'claimed', sum(claimed)/M);

if opt.verbose
    fprintf('[train %s] %d regions, %d networks, %.1f s\n', ...
            variant, numel(levels), nNets, model.info.trainTime);
end
end

% =========================================================================
function L = build_level(Zi, Yi, opt, useCluster, index, dataIdx)
%BUILD_LEVEL  Cluster one control region and train one network per cluster.
n = size(Zi,2);
if useCluster
    [r, rinfo] = ecr_choose_radius(Zi, opt);
    lab = ecr_cluster(Zi, r, opt);
    if all(lab == 0), lab = ones(1,n); end
else
    r = NaN;  rinfo = struct(); lab = ones(1,n);
end
CL = ecr_cluster_stats(Zi, lab, opt);

nets = cell(1, numel(CL));
rmse = zeros(1, numel(CL));
for j = 1:numel(CL)
    idx = CL(j).idx;
    nets{j} = rbf_train(Zi(:,idx), Yi(:,idx), opt);
    Yh = rbf_eval(nets{j}, Zi(:,idx));
    rmse(j) = sqrt(mean(mean((Yh - Yi(:,idx)).^2)));
end

L = struct('index', index, 'radius', r, 'radiusInfo', rinfo, 'CL', CL, ...
           'nets', {nets}, 'useNMD', useCluster, 'nData', n, 'rmse', rmse, ...
           'dataIdx', dataIdx);
end
