function ensemble = so_run_source_ensemble(varargin)
%SO_RUN_SOURCE_ENSEMBLE 50 random source points in the source rectangle.
%
%   ensemble = SO_RUN_SOURCE_ENSEMBLE()
%   ensemble = SO_RUN_SOURCE_ENSEMBLE('Case','vertical','TrialCount',10)
%
% This is the validation Schroer and Ott report for their Figure 3: rather
% than one favourable initial condition, draw many at random from the source
% region and report the spread of controlled transfer times.  Their sentence
% is "For 50 source points chosen at random in the gray source rectangle,
% the controlled trajectories reached the target rectangle after 125 to 132
% steps, while the uncontrolled transport times varied greatly between 1119
% and 3.77e6 steps."
%
% Two things are measured per trial and both matter:
%   * the CONTROLLED iteration count, whose spread should be narrow --
%     that is the claim that targeting time is set by the phase-space
%     structure, not by the particular initial condition;
%   * the UNCONTROLLED first-passage time from the same point, which should
%     be enormous and wildly variable -- that is what the control buys.
%
% Name/value options:
%   'Case'          case name                            (default 'diagonal')
%   'TrialCount'    number of source points              (default 50)
%   'Seed'          RNG seed for the source draw         (default 1729)
%   'RuntimeLimit'  seconds per trial before abort       (default from cfg)
%   'Uncontrolled'  also measure natural transport times (default true)
%   'UncontrolledCap' iteration cap for that measurement (default 4e6)
%   'Resume'        continue from a checkpoint if present(default true)
%   'Save'          write MAT/CSV                        (default true)
%   'Plot'          draw the ensemble figure             (default true)
%   'Verbose'       per-trial progress line              (default true)

opts = so_parse_options(varargin, struct( ...
    'Case', 'diagonal', 'TrialCount', [], 'Seed', [], 'RuntimeLimit', [], ...
    'Uncontrolled', true, 'UncontrolledCap', 4e6, 'Resume', true, ...
    'Budget', [20, 18, 38], ...
    'Save', true, 'Plot', true, 'Verbose', true));

[baseCfg, meta] = so_case_config(opts.Case);
baseCfg.saveFigures = false;

% ---- the ensemble needs a bigger budget than the case it is drawn from ---
%
% A case configuration is tuned for the transfer from the CENTRE of its
% source rectangle.  The ensemble draws from the whole rectangle, and for
% the diagonal geometry that straddles y = 1/2: sources below it pick up the
% omega = 1/2 resonance and take a three-resonance route, which the
% two-resonance budget cannot resolve.
%
% Measured on this geometry, 50 sources, same seed:
%
%     budget nF/nB/tau      sources reaching the target
%     11/8/19  (diagonal)          24 / 50
%     16/14/30 (vertical)          43 / 50
%     20/18/38                     49 / 50
%
% Every failure at 11/8/19 was a three-resonance source, and re-running
% those failures individually at 20/18/38 rescued all but one -- which is
% what identifies the cause as the budget rather than the geometry.  The
% cost is real: the failing trials currently bail out early, so raising the
% budget takes the sweep from roughly ten minutes to the better part of an
% hour.  That is the price of the trials actually finishing.
%
% Pass 'Budget', [] to inherit the case budgets unchanged.
if ~isempty(opts.Budget)
    if numel(opts.Budget) ~= 3
        error('SchroerOtt:BadEnsembleBudget', ...
            'Budget must be [maxForward, maxBackward, maxTotalTransferTime].');
    end
    baseCfg.maxForwardIterations = opts.Budget(1);
    baseCfg.maxBackwardIterations = opts.Budget(2);
    baseCfg.maxTotalTransferTime = opts.Budget(3);
end
if isempty(opts.TrialCount),   opts.TrialCount = baseCfg.ensemble.trialCount; end
if isempty(opts.Seed),         opts.Seed = baseCfg.ensemble.randomSeed; end
if isempty(opts.RuntimeLimit), opts.RuntimeLimit = baseCfg.ensemble.runtimeLimitSecondsPerTrial; end

outDir = fullfile(baseCfg.outputDirectory, 'ensemble');
figDir = fullfile(outDir, 'figures');
if opts.Save
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(figDir, 'dir'), mkdir(figDir); end
end
checkpointFile = fullfile(outDir, 'source_ensemble_checkpoint.mat');

rng(opts.Seed);
rect = baseCfg.sourceRectangle;
sourceStates = [rect.xMin + (rect.xMax - rect.xMin) * rand(1, opts.TrialCount); ...
                rect.yMin + (rect.yMax - rect.yMin) * rand(1, opts.TrialCount)];

if opts.Verbose
    fprintf('\nEnumerating periodic orbits for the ensemble...\n');
end
catalogue = so_enumerate_periodic_orbits(baseCfg);

rows = {};
results = cell(opts.TrialCount, 1);
startTrial = 1;
if opts.Resume && opts.Save && isfile(checkpointFile)
    chk = load(checkpointFile);
    if isfield(chk, 'sourceStates') && isequal(size(chk.sourceStates), size(sourceStates)) ...
            && max(abs(chk.sourceStates(:) - sourceStates(:))) < 1e-12
        rows = chk.rows;
        results = chk.results;
        startTrial = chk.lastCompletedTrial + 1;
        if opts.Verbose && startTrial > 1
            fprintf('Resuming from checkpoint at trial %d.\n', startTrial);
        end
    end
end

timer = tic;
for trial = startTrial:opts.TrialCount
    z0 = sourceStates(:, trial);
    cfg = baseCfg;
    cfg.sourceRectangle = struct('xMin', z0(1), 'xMax', z0(1), ...
        'yMin', z0(2), 'yMax', z0(2), 'id', sprintf('ensemble-source-%02d', trial));
    trialTimer = tic;
    cfg.runtime.startTime = trialTimer;
    cfg.runtime.limitSeconds = opts.RuntimeLimit;

    try
        route = so_construct_route(catalogue, cfg);
        result = so_multistage_targeting(cfg, catalogue, route);
        results{trial} = strip_result(result);
        row = trial_row(trial, z0, result, toc(trialTimer));
    catch err
        results{trial} = err;
        row = failed_row(trial, z0, err, toc(trialTimer));
    end

    rows(end + 1, :) = row; %#ok<AGROW>

    if opts.Verbose
        fprintf('  trial %2d/%2d  source (%.4f, %.4f)  %-28s  controlled iterations = %s  [%.1f s]\n', ...
            trial, opts.TrialCount, z0(1), z0(2), status_text(row{4}, row{5}), ...
            num_text(row{8}), row{19});
    end
    if opts.Save && mod(trial, baseCfg.ensemble.checkpointEvery) == 0
        lastCompletedTrial = trial; %#ok<NASGU>
        save(checkpointFile, 'rows', 'results', 'sourceStates', ...
            'lastCompletedTrial', 'baseCfg', 'meta', '-v7.3');
    end
end

names = {'trial','sourceX','sourceY','success','failureCategory','routeOmegas', ...
    'proxyCount','totalControlledIterations','numberOfControls','maxAbsControl', ...
    'totalAbsControl','maxIntersectionResidual','maxConsistencyError','finalX', ...
    'finalY','finalContainment','allMinimaCertified','unresolvedSplits', ...
    'runtimeSeconds','uncontrolledTransportTime'};
trials = cell2table(rows, 'VariableNames', names);

% The uncontrolled first-passage times are measured for every source at
% once: 50 independent orbits stepped together is a vector operation, while
% one orbit at a time would be millions of scalar iterations per trial.
if opts.Uncontrolled
    if opts.Verbose
        fprintf('\nMeasuring uncontrolled transport times (cap %g iterations)...\n', ...
            opts.UncontrolledCap);
    end
    trials.uncontrolledTransportTime = so_uncontrolled_transport_time( ...
        sourceStates, baseCfg, opts.UncontrolledCap).';
end
summary = ensemble_summary(trials);

ensemble.caseName = baseCfg.caseName;
ensemble.caseMetadata = meta;
ensemble.configuration = baseCfg;
ensemble.randomSeed = opts.Seed;
ensemble.sourceStates = sourceStates;
ensemble.trials = trials;
ensemble.summary = summary;
ensemble.results = results;
ensemble.runtimeSeconds = toc(timer);

if opts.Verbose
    fprintf('\n================ ENSEMBLE SUMMARY (%s) ================\n', char(baseCfg.caseName));
    disp(summary);
end
if opts.Save
    save(fullfile(outDir, 'source_ensemble.mat'), 'ensemble', '-v7.3');
    writetable(trials, fullfile(outDir, 'source_ensemble_trials.csv'));
    writetable(summary, fullfile(outDir, 'source_ensemble_summary.csv'));
end
if opts.Plot
    so_plot_ensemble(ensemble, opts.Save);
end
end

% ======================================================================= %

function r = strip_result(result)
%STRIP_RESULT Keep the ensemble MAT small: drop the per-trial curve data.
%
% A full result carries the orbit catalogue and every switch probe, each of
% which holds adaptively refined curves running to thousands of points.  Kept
% for 50 trials that is gigabytes for data identical across trials, so the
% shared and the recomputable parts are dropped here.
drop = {'orbitCatalogue', 'switchEvaluations', 'diagnosticManifolds'};
r = result;
for i = 1:numel(drop)
    if isfield(r, drop{i})
        r = rmfield(r, drop{i});
    end
end
for i = 1:numel(r.stagePlans)
    r.stagePlans(i).switchProbeConnections = [];
    r.stagePlans(i).selectedNextConnection = so_empty_connection();
end
end

function row = trial_row(trial, z0, result, runtimeSeconds)
success = result.targetReached;
if success
    failure = "";
else
    if ~isempty(result.failureDiagnostics)
        failure = string(result.failureDiagnostics(end).failureCategory);
    else
        failure = "not_reached";
    end
end
if result.numberOfControls > 0
    maxControl = max(abs(result.executedControls.controlY));
    totalControl = sum(abs(result.executedControls.controlY));
    maxConsistency = max(result.executedControls.propagationConsistencyError);
else
    maxControl = 0; totalControl = 0; maxConsistency = 0;
end
resid = -Inf;
certified = true;
for s = 1:numel(result.stagePlans)
    c = result.stagePlans(s).provisionalConnection;
    if c.success
        resid = max(resid, c.intersectionResidual);
        certified = certified && c.timeMinimumCertified;
    else
        certified = false;
    end
end
if ~isfinite(resid), resid = NaN; end
if success
    iterations = result.totalExecutedIterations;
else
    iterations = NaN;
end
if isempty(result.route.rotationNumbers)
    omegaText = "(empty bracket)";
else
    omegaText = string(mat2str(result.route.rotationNumbers, 5));
end
row = {trial, z0(1), z0(2), success, failure, omegaText, ...
    numel(result.route.rotationNumbers), iterations, result.numberOfControls, ...
    maxControl, totalControl, resid, maxConsistency, ...
    result.finalState(1), result.finalState(2), result.targetContained, ...
    certified, height(result.resolutionFailures), runtimeSeconds, NaN};
end

function row = failed_row(trial, z0, err, runtimeSeconds)
row = {trial, z0(1), z0(2), false, string(err.identifier), "", 0, NaN, 0, ...
    NaN, NaN, NaN, NaN, NaN, NaN, false, false, 0, runtimeSeconds, NaN};
end

function summary = ensemble_summary(trials)
mask = trials.success & isfinite(trials.totalControlledIterations);
counts = trials.totalControlledIterations(mask);
unc = trials.uncontrolledTransportTime(isfinite(trials.uncontrolledTransportTime));
if isempty(counts)
    stats = {height(trials), 0, 0, NaN, NaN, NaN, NaN, NaN, NaN, "", ...
        NaN, NaN, NaN, NaN};
else
    stats = {height(trials), sum(mask), sum(mask) / height(trials), ...
        min(counts), max(counts), mean(counts), median(counts), std(counts), ...
        so_interquartile_range(counts), string(mat2str(sort(counts).')), ...
        emptyfun(@min, unc), emptyfun(@max, unc), emptyfun(@median, unc), ...
        median_ratio(unc, counts)};
end
summary = cell2table(stats, 'VariableNames', {'attempts','successes', ...
    'successFraction','minIterations','maxIterations','meanIterations', ...
    'medianIterations','stdIterations','iqrIterations','allIterationCounts', ...
    'minUncontrolled','maxUncontrolled','medianUncontrolled','medianSpeedupFactor'});
end

function v = emptyfun(fn, x)
if isempty(x), v = NaN; else, v = fn(x); end
end

function v = median_ratio(unc, counts)
if isempty(unc) || isempty(counts)
    v = NaN;
else
    v = median(unc) / median(counts);
end
end

function t = status_text(success, failure)
if success
    t = 'target reached';
else
    t = ['FAILED: ', char(failure)];
end
end

function t = num_text(v)
if isnan(v)
    t = '   --';
else
    t = sprintf('%5g', v);
end
end
