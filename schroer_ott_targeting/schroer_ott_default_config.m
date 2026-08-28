function cfg = schroer_ott_default_config()
%SCHROER_OTT_DEFAULT_CONFIG Central first-light configuration.
%
% The defaults are deliberately modest.  They target the k = 1.25
% normalized standard-map demonstration on the cylinder S^1 x R.  The
% source and target rectangles are reproducible demonstration rectangles,
% not digitized paper coordinates.  The source is deliberately outside all
% retained proxy balls and outside every admissible zero-time proxy kick.

cfg.k = 1.25;
cfg.controlAmplitude = 0.003;

cfg.sourceRectangle = struct('xMin', 0.859409, 'xMax', 0.909409, ...
    'yMin', 0.486716, 'yMax', 0.526716, 'id', 'source');
cfg.targetRectangle = struct('xMin', 0.21, 'xMax', 0.25, ...
    'yMin', 0.72, 'yMax', 0.76, 'id', 'final-target');

cfg.transport.allowEquivalentTargetCells = false;
cfg.transport.targetLiftShift = 0;
cfg.transport.liftSelectionPolicy = 'declared-fixed-target-cell';

cfg.routeBracket.margin = 1e-10;
cfg.routeBracket.lowerBoundary = 'open';
cfg.routeBracket.upperBoundary = 'open';

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

cfg.proxyTargetRadius = 0.030;
cfg.maxForwardIterations = 8;
cfg.maxBackwardIterations = 6;
cfg.maxTotalTransferTime = 14;
cfg.maxStages = 16;

cfg.curve.initialControlSamples = 17;
cfg.curve.initialBoundarySamples = 64;
cfg.curve.maxGapFraction = 0.12;
cfg.curve.midpointToleranceFraction = 0.06;
cfg.curve.maxPoints = 1e5;
cfg.curve.maxSubdivisionDepth = 15;

cfg.intersectionBackend = 'indexed';
cfg.intersectionTolerance = 1e-9;
cfg.containmentTolerance = 1e-9;
cfg.propagationConsistencyTolerance = 1e-8;

cfg.equationOne.L = 1;
cfg.equationOne.LConvention = 'declared order-one phase-space scale';
cfg.equationOne.discrepancyFactor = 10;
cfg.equationOne.targetSizeConvention = 'diameter';

cfg.randomSeed = 1;
cfg.firstLightMode = true;

cfg.outputDirectory = fullfile(pwd, 'outputs');
cfg.figureDirectory = fullfile(cfg.outputDirectory, 'figures');
cfg.saveFigures = true;
cfg.verbose = true;
end
