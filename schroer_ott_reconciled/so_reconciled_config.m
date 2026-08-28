function cfg = so_reconciled_config()
%SO_RECONCILED_CONFIG Base configuration for the reconciled build.
%
% Numerically identical to the fundamental-folder default configuration.
% Everything added here is presentation, validation or bookkeeping; nothing
% in this file changes how a control is selected.
%
% The default geometry is the DIAGONAL case: the source sits low-right in
% the chaotic sea and the target high-left, so the route brackets two
% direct-hyperbolic resonances (omega = 3/5 and 2/3).

cfg.k = 1.25;
cfg.controlAmplitude = 0.003;

% ---- geometry (diagonal case; overridden by so_case_config) -------------
cfg.caseName = "diagonal";
cfg.caseLabel = "Diagonal transfer (source low-right, target high-left)";
cfg.sourceRectangle = struct('xMin', 0.64, 'xMax', 0.69, ...
    'yMin', 0.48, 'yMax', 0.52, 'id', 'source');
cfg.targetRectangle = struct('xMin', 0.21, 'xMax', 0.25, ...
    'yMin', 0.71, 'yMax', 0.75, 'id', 'final-target');

cfg.transport.allowEquivalentTargetCells = false;
cfg.transport.targetLiftShift = 0;
cfg.transport.liftSelectionPolicy = 'declared-fixed-target-cell';

% An empty rotation bracket is a legitimate route (see so_construct_route).
cfg.route.allowEmptyBracket = true;

% ---- periodic-orbit search ---------------------------------------------
cfg.orbit.minPeriod = 2;
cfg.orbit.maxPeriod = 5;
cfg.orbit.residualTolerance = 1e-10;
cfg.orbit.deduplicationTolerance = 2e-8;
cfg.orbit.traceTolerance = 1e-7;
cfg.orbit.newtonMaxIterations = 80;
cfg.orbit.newtonTolerance = 1e-12;
cfg.orbit.seedXCount = 49;
cfg.orbit.seedYHalfWidth = 0.45;
cfg.orbit.seedYOffsets = linspace(-0.45, 0.45, 13);

% ---- targeting budgets --------------------------------------------------
cfg.proxyTargetRadius = 0.030;
cfg.maxForwardIterations = 11;
cfg.maxBackwardIterations = 8;
cfg.maxTotalTransferTime = 19;
cfg.maxStages = 16;

% ---- adaptive curve resolution -----------------------------------------
cfg.curve.initialControlSamples = 17;
cfg.curve.initialBoundarySamples = 64;
cfg.curve.maxGapFraction = 0.12;
cfg.curve.midpointToleranceFraction = 0.06;
cfg.curve.maxPoints = 1e5;
cfg.curve.maxSubdivisionDepth = 15;

cfg.intersectionTolerance = 1e-9;
cfg.containmentTolerance = 1e-9;
cfg.propagationConsistencyTolerance = 1e-8;

cfg.equationOne.L = 1;
cfg.equationOne.discrepancyFactor = 10;

% ---- phase-portrait background (figures only) ---------------------------
% The background is the uncontrolled standard map itself: a grid of seeds
% iterated forward and plotted as small dots, so islands appear as closed
% curves and the connected chaotic region appears as stipple, exactly as in
% Figure 3 of Schroer and Ott.
cfg.background.enable = true;
cfg.background.xSeedCount = 60;
cfg.background.ySeedCount = 34;
cfg.background.iterations = 260;
cfg.background.yPadding = 0.14;
cfg.background.markerSize = 1;
cfg.background.color = [0.62 0.62 0.62];
cfg.background.cacheFile = 'phase_portrait_cache.mat';

% ---- 50-source ensemble validation -------------------------------------
cfg.ensemble.trialCount = 50;
cfg.ensemble.randomSeed = 1729;
cfg.ensemble.runtimeLimitSecondsPerTrial = 240;
cfg.ensemble.checkpointEvery = 1;

% ---- noise validation ---------------------------------------------------
% Schroer and Ott report additive noise up to 1e-2 * delta, where delta is
% the full length of the admissible control segment.
cfg.noise.sigmaFactors = [0, 1e-3, 1e-2, 1e-1];
cfg.noise.headlineFactor = 1e-2;
cfg.noise.realisationsPerLevel = 40;
cfg.noise.randomSeed = 314159;

cfg.randomSeed = 1;
cfg.firstLightMode = true;

cfg.outputDirectory = fullfile(pwd, 'outputs');
cfg.figureDirectory = fullfile(cfg.outputDirectory, 'figures');
cfg.saveFigures = true;
cfg.figureFormats = {'png'};
cfg.verbose = true;
end
