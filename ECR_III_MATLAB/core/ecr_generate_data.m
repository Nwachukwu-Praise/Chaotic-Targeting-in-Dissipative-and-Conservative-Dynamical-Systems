function D = ecr_generate_data(S, opt)
%ECR_GENERATE_DATA  Collect the (z_n, p_n, z_{n+1}) triples used to train ECR.
%
%   D = ECR_GENERATE_DATA(S,OPT)
%
%   The system is run from many initial conditions on the attractor while a
%   *random* parameter vector p_n, uniformly distributed inside the allowable
%   set Pi of Eq. (3), is applied at every step.  This is the "training data
%   gathered with random parameter variations within Pi" of Section 3.1 - the
%   only information the ECR method ever uses about the plant.
%
%   Output fields (M = number of valid triples)
%     D.H    nh x M   hidden state before the step (bookkeeping only)
%     D.Z     N x M   observed state z_n
%     D.P     r x M   applied parameter vector p_n
%     D.Hn   nh x M   hidden state after the step
%     D.Zn    N x M   observed state z_{n+1} = G(z_n,p_n)
%     D.opt           options used
%
%   Note on ENRICHMENT (OPT.enrich): a chaotic trajectory visits the tiny
%   target region S0 rarely, so the level-0 network would be starved of data.
%   States that come within OPT.enrichRadius*delta of z* are therefore
%   re-simulated OPT.enrichFactor times with fresh random parameters, and -
%   when the system allows an observed state to be lifted back to a full
%   state (S.hfromz, i.e. real-coordinate observation) - a few states are
%   seeded directly around z*.  Both are data-collection conveniences; they
%   use no model knowledge beyond the ability to re-run the experiment.

if nargin < 2 || isempty(opt), opt = ecr_default_options(); end
ecr_seed(opt.seed);

nTraj = max(4, ceil(opt.nData/opt.trajLen));
if opt.verbose
    fprintf('[data] %s: %d trajectories x %d steps\n', S.name, nTraj, opt.trajLen);
end

H = S.init(nTraj);

Hc = cell(1, opt.trajLen);  Zc = Hc;  Pc = Hc;  Hnc = Hc;  Znc = Hc;
for k = 1:opt.trajLen
    P = repmat(S.pnom,1,nTraj) + repmat(S.dpmax,1,nTraj).*(2*rand(S.r,nTraj)-1);
    Z = S.obs(H);
    [Hn, Zn] = S.step(H, P);

    Hc{k} = H;  Zc{k} = Z;  Pc{k} = P;  Hnc{k} = Hn;  Znc{k} = Zn;

    bad = any(~isfinite(Hn), 1);
    if any(bad)
        Hn(:,bad) = S.init(sum(bad));       % re-seed lost trajectories
    end
    H = Hn;
    if opt.verbose && mod(k, max(1,round(opt.trajLen/5))) == 0
        fprintf('       step %d/%d\n', k, opt.trajLen);
    end
end

D.H  = [Hc{:}];   D.Z  = [Zc{:}];   D.P = [Pc{:}];
D.Hn = [Hnc{:}];  D.Zn = [Znc{:}];

% ---- enrichment near the target ----------------------------------------
if opt.enrich
    Hs = [];
    d  = sqrt(sum((D.Z - repmat(S.zstar,1,size(D.Z,2))).^2, 1));

    % (a) states in the neighbourhood of the target: try several parameters
    near = find(d < opt.enrichRadius*S.delta);
    if ~isempty(near)
        near = near(randperm(numel(near)));
        near = near(1:min(numel(near), ceil(opt.nData/10)));
        Hs = repmat(D.H(:,near), 1, opt.enrichFactor);
    end
    % (b) when the observation can be lifted back to a full state, seed
    %     directly inside the delta-ball and in its neighbourhood
    if ~isempty(S.hfromz)
        nSeed = ceil(opt.nData/20);
        Zb = randn(S.N, nSeed);
        Zb = Zb ./ repmat(sqrt(sum(Zb.^2,1)), S.N, 1);
        Zb = repmat(S.zstar,1,nSeed) + ...
             Zb .* repmat(S.delta*rand(1,nSeed).^(1/S.N), S.N, 1);   % in the ball
        Zo = repmat(S.zstar,1,nSeed) + ...
             S.delta*opt.enrichRadius*(2*rand(S.N,nSeed)-1);          % around it
        Hs = [Hs, S.hfromz([Zb, Zo])];
    end
    if ~isempty(Hs)
        Ps = repmat(S.pnom,1,size(Hs,2)) + ...
             repmat(S.dpmax,1,size(Hs,2)).*(2*rand(S.r,size(Hs,2))-1);
        [Hns, Zns] = S.step(Hs, Ps);
        D.H  = [D.H,  Hs];   D.Z  = [D.Z,  S.obs(Hs)];  D.P = [D.P, Ps];
        D.Hn = [D.Hn, Hns];  D.Zn = [D.Zn, Zns];
        if opt.verbose
            fprintf('       enrichment: +%d triples near the target\n', size(Hs,2));
        end
    end
end

% ---- drop invalid triples ----------------------------------------------
good = all(isfinite(D.Z),1) & all(isfinite(D.Zn),1) & all(isfinite(D.Hn),1);
D.H = D.H(:,good);  D.Z = D.Z(:,good);  D.P = D.P(:,good);
D.Hn = D.Hn(:,good); D.Zn = D.Zn(:,good);
D.opt = opt;

if opt.verbose
    dz = sqrt(sum((D.Z - repmat(S.zstar,1,size(D.Z,2))).^2,1));
    fprintf('[data] %d triples, %d of them inside the delta-ball of z*\n', ...
            size(D.Z,2), sum(dz < S.delta));
end
end
