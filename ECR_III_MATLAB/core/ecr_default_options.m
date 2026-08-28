function opt = ecr_default_options()
%ECR_DEFAULT_OPTIONS  Default options of the ECR-II / ECR-III implementation.
%
%   opt = ECR_DEFAULT_OPTIONS() returns a struct with every tunable knob of
%   the toolbox.  Override the fields you care about and pass the struct to
%   ECR_GENERATE_DATA / ECR_TRAIN / ECR_REACHING_TIME.
%
%   Where a choice is *not* fixed by the paper the default is flagged
%   "[interpretation]" and is discussed in docs/METHOD_NOTES.md.

opt = struct();

% ---- data generation ----------------------------------------------------
opt.nData        = 20000;   % number of (z_n, p_n, z_{n+1}) training triples
opt.trajLen      = 40;      % steps per trajectory before re-seeding
opt.enrich       = true;    % re-simulate states close to z* to populate S0
opt.enrichFactor = 6;       % branches per near-target state       [interpretation]
opt.enrichRadius = 4;       % "close to z*" = within this many delta
opt.seed         = 1;       % RNG seed

% ---- region construction (Definitions 1-3) ------------------------------
opt.Kmax         = 5;       % maximum number of control regions T_1..T_K
opt.minLevelPts  = 40;      % stop adding regions below this much data
opt.targetTighten= 0.5;     % NN_0 is trained on the S0 pairs whose successor
                            % lands within this fraction of delta (1 = the
                            % literal Definition 1)               [interpretation]

% ---- ECR-III clustering -------------------------------------------------
opt.radiusMode   = 'hist';  % 'hist' = first minimum of the inter-distance
                            % histogram (paper), 'nn' = k-th nearest neighbour
opt.nHistBins    = 60;      % bins of that histogram
opt.histSubset   = 1200;    % points used to build the histogram
opt.histSmooth   = 3;       % moving-average width applied before the search
opt.radiusScale  = 1.0;     % multiplies the radius found by the histogram
opt.minClusterPts= 12;      % clusters smaller than this are merged / dropped
opt.smallCluster = 'merge'; % 'merge' into nearest cluster, or 'drop'
opt.maxClusters  = 40;      % safety cap per control region

% ---- normalised Mahalanobis distance (Eq. 8) ----------------------------
opt.covNorm      = 'det';   % 'det' | 'plain' | 'corr' | 'trace'  [interpretation]
opt.covReg       = 1e-3;    % ridge added to the covariance (times its trace)
opt.nmdGate      = 2.0;     % accept a cluster only if NMD <= gate * NMD95 of
                            % its own training data (Inf = literal paper) [interpretation]

% ---- RBF networks -------------------------------------------------------
opt.nCenters     = 20;      % RBF centres per network
opt.rbfSpread    = 1.2;     % width factor
opt.rbfLambda    = 1e-6;    % ridge regularisation of the output weights
opt.rbfMinPts    = 6;       % below this a network degenerates to the mean p

% ---- on-line control ----------------------------------------------------
opt.nFallback    = 3;       % clusters tried in NMD order before giving up
                            % (1 = literal paper reading)          [interpretation]
opt.levelPriority= false;   % true = prefer low-index regions on exact ties
opt.selection    = 'nmd';   % ECR-III region selection:
                            %   'nmd'   nearest cluster over ALL regions
                            %           (the literal Section 3.3 rule)
                            %   'level' lowest region index first, nearest
                            %           cluster inside that region  [interpretation]

% ---- simulation / benchmarking ------------------------------------------
opt.noise        = 0;       % measurement noise, fraction of S.rms
opt.escapeAfter  = 12;      % if targeting has been applied for this many
                            % consecutive steps without capturing the target,
                            % apply p_nom for one step.  Deterministic
                            % feedback can otherwise lock a trajectory into a
                            % controlled periodic orbit that never enters S0
                            % (0 disables the escape)              [interpretation]
opt.maxSteps     = 300;     % give up after this many map steps
opt.nTrials      = 100;     % initial conditions per benchmark
opt.holdSteps    = 25;      % steps of local control checked after capture
opt.verbose      = 1;
end
