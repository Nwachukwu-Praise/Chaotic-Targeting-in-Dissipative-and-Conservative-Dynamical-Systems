function study = so_run_resolution_study(saveOutputs)
%SO_RUN_RESOLUTION_STUDY Compact three-level deterministic resolution check.
if nargin < 1
    saveOutputs = true;
end
baseCfg = schroer_ott_default_config();
baseCfg.saveFigures = false;
if saveOutputs && ~exist(baseCfg.outputDirectory, 'dir')
    mkdir(baseCfg.outputDirectory);
end
levels = resolution_levels(baseCfg);
rows = {};
results = cell(numel(levels), 1);
timer = tic;
for i = 1:numel(levels)
    cfg = levels(i).cfg;
    catalogue = so_enumerate_periodic_orbits(cfg);
    route = so_construct_first_light_route(catalogue, cfg);
    trialTimer = tic;
    try
        result = so_multistage_targeting(cfg, catalogue, route);
        results{i} = result;
        rows(end + 1, :) = result_row(levels(i).name, cfg, result, toc(trialTimer)); %#ok<AGROW>
    catch err
        results{i} = err;
        rows(end + 1, :) = failure_row(levels(i).name, cfg, err, toc(trialTimer)); %#ok<AGROW>
    end
end
names = {'level','initialControlSamples','initialBoundarySamples','maxGapFraction', ...
    'midpointToleranceFraction','maxSubdivisionDepth','maxPoints','success', ...
    'failureCategory','selectedRoute','controls','splits','switchingIndices', ...
    'totalIterations','maxResidual','maxReplayError','targetContainment', ...
    'unresolvedSplits','minimumCertified','peakCurvePointCount','runtimeSeconds'};
comparison = cell2table(rows, 'VariableNames', names);
study.levels = levels;
study.comparison = comparison;
study.results = results;
study.runtimeSeconds = toc(timer);
if saveOutputs
    save(fullfile(baseCfg.outputDirectory, 'resolution_study.mat'), 'study', '-v7.3');
    writetable(comparison, fullfile(baseCfg.outputDirectory, 'resolution_study.csv'));
end
end

function levels = resolution_levels(baseCfg)
levels = struct('name', {}, 'cfg', {});
levels(1).name = "reduced";
levels(1).cfg = baseCfg;
levels(1).cfg.curve.initialControlSamples = 13;
levels(1).cfg.curve.initialBoundarySamples = 48;
levels(1).cfg.curve.maxGapFraction = 0.14;
levels(1).cfg.curve.midpointToleranceFraction = 0.07;
levels(1).cfg.curve.maxPoints = 50000;

levels(2).name = "baseline";
levels(2).cfg = baseCfg;

levels(3).name = "increased";
levels(3).cfg = baseCfg;
levels(3).cfg.curve.initialControlSamples = 21;
levels(3).cfg.curve.initialBoundarySamples = 80;
levels(3).cfg.curve.maxGapFraction = 0.10;
levels(3).cfg.curve.midpointToleranceFraction = 0.05;
levels(3).cfg.curve.maxPoints = 120000;
end

function row = result_row(levelName, cfg, result, runtimeSeconds)
splits = arrayfun(@(s) mat2str(s.provisionalConnection.forwardBackwardSplit), ...
    result.stagePlans, 'UniformOutput', false);
switches = arrayfun(@(s) s.selectedSwitchIndex, result.stagePlans);
maxResidual = max(arrayfun(@(s) s.provisionalConnection.intersectionResidual, result.stagePlans));
maxReplay = max(result.executedControls.maxPathwiseReplayError);
certified = all(arrayfun(@(s) s.provisionalConnection.timeMinimumCertified, result.stagePlans));
row = {levelName, cfg.curve.initialControlSamples, cfg.curve.initialBoundarySamples, ...
    cfg.curve.maxGapFraction, cfg.curve.midpointToleranceFraction, ...
    cfg.curve.maxSubdivisionDepth, cfg.curve.maxPoints, result.targetReached, "", ...
    mat2str(result.route.rotationNumbers, 6), result.numberOfControls, ...
    strjoin(string(splits), '; '), mat2str(switches), result.totalExecutedIterations, ...
    maxResidual, maxReplay, result.targetContained, result.performanceProfile.unresolvedSplits, ...
    certified, result.performanceProfile.peakCurvePointCount, runtimeSeconds};
end

function row = failure_row(levelName, cfg, err, runtimeSeconds)
row = {levelName, cfg.curve.initialControlSamples, cfg.curve.initialBoundarySamples, ...
    cfg.curve.maxGapFraction, cfg.curve.midpointToleranceFraction, ...
    cfg.curve.maxSubdivisionDepth, cfg.curve.maxPoints, false, string(err.identifier), ...
    "", 0, "", "", Inf, NaN, NaN, false, NaN, false, NaN, runtimeSeconds};
end
