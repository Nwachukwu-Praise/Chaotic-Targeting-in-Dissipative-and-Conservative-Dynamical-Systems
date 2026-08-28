function ensemble = so_run_paper_ensemble(saveOutputs)
%SO_RUN_PAPER_ENSEMBLE Run 50 estimated Figure-3 source points.
if nargin < 1
    saveOutputs = true;
end
[baseCfg, paper] = so_paper_comparison_config();
baseCfg.saveFigures = false;
baseCfg.runtimeLimitSecondsPerTrial = 120;
if saveOutputs
    ensure_output_dirs(baseCfg);
end
rng(1729);
sourceStates = sample_sources(paper.source, 50);
catalogue = so_enumerate_periodic_orbits(baseCfg);

checkpointFile = fullfile(baseCfg.outputDirectory, 'paper_50_source_ensemble_checkpoint.mat');
rows = {};
results = cell(50, 1);
startTrial = 1;
if saveOutputs && isfile(checkpointFile)
    checkpoint = load(checkpointFile, 'rows', 'results', 'sourceStates', 'lastCompletedTrial');
    rows = checkpoint.rows;
    results = checkpoint.results;
    sourceStates = checkpoint.sourceStates;
    startTrial = checkpoint.lastCompletedTrial + 1;
end
timer = tic;
for trial = startTrial:50
    cfg = baseCfg;
    z0 = sourceStates(:, trial);
    cfg.sourceRectangle = point_rectangle(z0, sprintf('paper-source-%02d', trial));
    trialTimer = tic;
    cfg.runtime.startTime = trialTimer;
    cfg.runtime.limitSeconds = baseCfg.runtimeLimitSecondsPerTrial;
    route = so_construct_omega_route(catalogue, cfg, paper.expectedRotationNumbers, ...
        'estimated Figure-3 seven-resonance route');
    try
        result = so_multistage_targeting(cfg, catalogue, route);
        results{trial} = result;
        [success, failureCategory] = trial_status(result);
        rows(end + 1, :) = trial_row(trial, z0, success, failureCategory, result, toc(trialTimer)); %#ok<AGROW>
    catch err
        results{trial} = err;
        rows(end + 1, :) = {trial, z0(1), z0(2), false, string(err.identifier), ...
            "", 0, Inf, NaN, NaN, NaN, NaN, NaN, NaN, false, 0, false, toc(trialTimer)}; %#ok<AGROW>
    end
    if saveOutputs
        lastCompletedTrial = trial; %#ok<NASGU>
        save(checkpointFile, 'rows', 'results', 'sourceStates', 'lastCompletedTrial', ...
            'paper', 'baseCfg', '-v7.3');
    end
end

names = {'trial','sourceX','sourceY','success','failureCategory','selectedRoute', ...
    'numberOfControls','totalControlledIterations','maxAbsControl','totalAbsControl', ...
    'maxIntersectionResidual','maxPathwiseReplayError','finalX','finalY', ...
    'finalContainment','unresolvedSplits','minimumCertified','runtimeSeconds'};
trialTable = cell2table(rows, 'VariableNames', names);
summary = ensemble_summary(trialTable, paper.publishedControlledRange);

ensemble.configuration = baseCfg;
ensemble.paperMetadata = paper;
ensemble.randomSeed = 1729;
ensemble.sourceStates = sourceStates;
ensemble.trials = trialTable;
ensemble.summary = summary;
ensemble.results = results;
ensemble.runtimeSeconds = toc(timer);
if saveOutputs
    save(fullfile(baseCfg.outputDirectory, 'paper_50_source_ensemble.mat'), 'ensemble', '-v7.3');
    writetable(trialTable, fullfile(baseCfg.outputDirectory, 'paper_50_source_trials.csv'));
    writetable(summary, fullfile(baseCfg.outputDirectory, 'paper_50_source_summary.csv'));
end
end

function sourceStates = sample_sources(rect, n)
sourceStates = [rect.xMin + (rect.xMax - rect.xMin) .* rand(1, n); ...
    rect.yMin + (rect.yMax - rect.yMin) .* rand(1, n)];
end

function rect = point_rectangle(z, id)
rect = struct('xMin', z(1), 'xMax', z(1), 'yMin', z(2), 'yMax', z(2), 'id', id);
end

function [success, failureCategory] = trial_status(result)
success = isfield(result, 'targetReached') && result.targetReached;
if success
    failureCategory = "";
elseif isfield(result, 'failureDiagnostics') && ~isempty(result.failureDiagnostics)
    failureCategory = string(result.failureDiagnostics(end).failureCategory);
else
    failureCategory = "not_reached";
end
end

function row = trial_row(trial, z0, success, failureCategory, result, runtimeSeconds)
if success
    controls = result.executedControls.controlY;
    maxAbsControl = max(abs(controls));
    totalAbsControl = sum(abs(controls));
    maxResidual = max(arrayfun(@(s) s.provisionalConnection.intersectionResidual, result.stagePlans));
    maxReplay = max(result.executedControls.maxPathwiseReplayError);
    finalState = result.finalState;
    unresolved = result.performanceProfile.unresolvedSplits;
    certified = all(arrayfun(@(s) s.provisionalConnection.timeMinimumCertified, result.stagePlans));
else
    maxAbsControl = NaN;
    totalAbsControl = NaN;
    maxResidual = NaN;
    maxReplay = NaN;
    finalState = [NaN; NaN];
    unresolved = result.performanceProfile.unresolvedSplits;
    certified = false;
end
row = {trial, z0(1), z0(2), success, string(failureCategory), ...
    mat2str(result.route.rotationNumbers, 6), result.numberOfControls, ...
    result.totalExecutedIterations, maxAbsControl, totalAbsControl, maxResidual, ...
    maxReplay, finalState(1), finalState(2), result.targetContained, unresolved, ...
    certified, runtimeSeconds};
end

function summary = ensemble_summary(trialTable, publishedRange)
finiteMask = trialTable.success & isfinite(trialTable.totalControlledIterations);
counts = trialTable.totalControlledIterations(finiteMask);
if isempty(counts)
    stats = {height(trialTable), 0, 0, NaN, NaN, NaN, NaN, NaN, NaN, "", 0, 0, 0};
else
    stats = {height(trialTable), sum(finiteMask), sum(finiteMask) / height(trialTable), ...
        min(counts), max(counts), mean(counts), median(counts), std(counts), ...
        local_iqr(counts), mat2str(counts.'), sum(counts < publishedRange(1)), ...
        sum(counts >= publishedRange(1) & counts <= publishedRange(2)), ...
        sum(counts > publishedRange(2))};
end
summary = cell2table(stats, 'VariableNames', {'attempts','successes','successFraction', ...
    'minimum','maximum','mean','median','standardDeviation','interquartileRange', ...
    'finiteIterationCounts','numberBelow125','numberWithin125To132','numberAbove132'});
end

function v = local_iqr(x)
%LOCAL_IQR Interquartile range without the Statistics Toolbox.
%
%   The previous code called iqr, which lives in the Statistics and Machine
%   Learning Toolbox.  It was never reached while the ensemble returned no
%   successes, so the dependency was latent.  Linear-interpolated quartiles
%   are used here, matching the default prctile convention.
x = sort(x(:));
n = numel(x);
if n == 0
    v = NaN;
    return;
elseif n == 1
    v = 0;
    return;
end
v = local_percentile(x, 75) - local_percentile(x, 25);
end

function q = local_percentile(sortedX, pct)
n = numel(sortedX);
% Sample points sit at (i - 0.5)/n; outside that range the extremes are used.
pos = pct / 100 * n + 0.5;
if pos <= 1
    q = sortedX(1);
elseif pos >= n
    q = sortedX(n);
else
    lo = floor(pos);
    frac = pos - lo;
    q = (1 - frac) * sortedX(lo) + frac * sortedX(lo + 1);
end
end

function ensure_output_dirs(cfg)
if ~exist(cfg.outputDirectory, 'dir')
    mkdir(cfg.outputDirectory);
end
if ~exist(cfg.figureDirectory, 'dir')
    mkdir(cfg.figureDirectory);
end
end
