function ensemble = so_run_paper_ensemble_double_runtime(paperEnsembleAction, varargin)
%SO_RUN_PAPER_ENSEMBLE_DOUBLE_RUNTIME Versioned 50-source doubled-budget run.
%
%   The default action is 'estimate'.  It audits the existing 120 s
%   paper-comparison ensemble, verifies the saved source states against the
%   original seed, and reports the versioned 240 s budget without launching
%   the full rerun.
%
%   To start or resume the complete rerun, call
%
%       so_run_paper_ensemble_double_runtime('resumeDoubleRuntime')
%
%   Name-value options used by tests and pilots include SourceIndices,
%   RuntimeBudgetSeconds, UseParallel, NumWorkers, OutputDirectory,
%   OutputStem, SaveOutputs, Quiet, and Force.

if nargin < 1 || isempty(paperEnsembleAction)
    paperEnsembleAction = "estimate";
end
paperEnsembleAction = string(paperEnsembleAction);
opts = parse_options(varargin{:});
paths = output_paths(opts);
audit = runtime_audit(opts, paths);

switch lower(paperEnsembleAction)
    case "estimate"
        ensemble = audit;
        if opts.saveOutputs
            ensure_output_dir(paths.outputDirectory);
            comparison = runtime_comparison_table(audit.previousTerminations, table(), ...
                audit.originalRuntimeBudgetSeconds, audit.newRuntimeBudgetSeconds);
            writetable(comparison, paths.comparisonCsv);
            log_message(paths, opts, sprintf(['estimate only: original %.6g s, ', ...
                'new %.6g s, previous timeouts %d/50'], ...
                audit.originalRuntimeBudgetSeconds, audit.newRuntimeBudgetSeconds, ...
                audit.previousTerminationCounts.time_budget_exhausted));
        end
    case "load"
        if isfile(paths.finalMat)
            loaded = load(paths.finalMat, 'ensemble');
            ensemble = loaded.ensemble;
        else
            ensemble = audit;
        end
    case "pilot"
        if isempty(opts.sourceIndices)
            opts.sourceIndices = first_previous_timeout(audit.previousTerminations);
        end
        if strcmp(opts.outputStem, default_output_stem())
            opts.outputStem = default_output_stem() + "_pilot";
        end
        opts.useParallel = false;
        paths = output_paths(opts);
        ensemble = run_requested_sources(audit, opts, paths, "pilot");
    case "resumedoubleruntime"
        ensemble = run_requested_sources(audit, opts, paths, "resumeDoubleRuntime");
    otherwise
        error('SchroerOtt:UnknownPaperEnsembleAction', ...
            'Unknown paperEnsembleAction "%s".', paperEnsembleAction);
end
end

function opts = parse_options(varargin)
opts.outputDirectory = fullfile(pwd, 'outputs', 'paper_comparison');
opts.outputStem = default_output_stem();
opts.sourceIndices = [];
opts.useParallel = "auto";
opts.numWorkers = [];
opts.runtimeBudgetSeconds = NaN;
opts.saveOutputs = true;
opts.quiet = false;
opts.force = false;
opts.originalFile = "";
if mod(numel(varargin), 2) ~= 0
    error('SchroerOtt:InvalidNameValue', 'Options must be supplied as name-value pairs.');
end
for i = 1:2:numel(varargin)
    name = lower(string(varargin{i}));
    value = varargin{i + 1};
    switch name
        case "outputdirectory"
            opts.outputDirectory = char(value);
        case "outputstem"
            opts.outputStem = string(value);
        case "sourceindices"
            opts.sourceIndices = value(:).';
        case "useparallel"
            if islogical(value)
                opts.useParallel = value;
            else
                opts.useParallel = string(value);
            end
        case "numworkers"
            opts.numWorkers = value;
        case "runtimebudgetseconds"
            opts.runtimeBudgetSeconds = value;
        case "saveoutputs"
            opts.saveOutputs = logical(value);
        case "quiet"
            opts.quiet = logical(value);
        case "force"
            opts.force = logical(value);
        case "originalfile"
            opts.originalFile = char(value);
        otherwise
            error('SchroerOtt:UnknownOption', 'Unknown option "%s".', name);
    end
end
end

function stem = default_output_stem()
stem = "schroer_ott_paper_ensemble_double_runtime";
end

function paths = output_paths(opts)
paths.outputDirectory = char(opts.outputDirectory);
stem = string(opts.outputStem);
paths.outputStem = char(stem);
paths.finalMat = fullfile(paths.outputDirectory, char(stem + ".mat"));
paths.checkpointMat = fullfile(paths.outputDirectory, char(stem + "_checkpoint.mat"));
paths.trialsCsv = fullfile(paths.outputDirectory, char(stem + ".csv"));
paths.comparisonCsv = fullfile(paths.outputDirectory, char(stem + "_runtime_comparison.csv"));
if strcmp(stem, default_output_stem())
    paths.comparisonCsv = fullfile(paths.outputDirectory, ...
        'schroer_ott_paper_ensemble_runtime_comparison.csv');
end
paths.log = fullfile(paths.outputDirectory, paths.outputStem + ".log");
end

function audit = runtime_audit(opts, paths)
[baseCfg, paper] = so_paper_comparison_config();
baseCfg.saveFigures = false;
if strlength(string(opts.originalFile)) > 0
    originalFile = opts.originalFile;
else
    originalFile = fullfile(baseCfg.outputDirectory, 'paper_50_source_ensemble.mat');
end
if ~isfile(originalFile)
    error('SchroerOtt:MissingOriginalPaperEnsemble', ...
        'Cannot find the original ensemble file: %s', originalFile);
end
loaded = load(originalFile, 'ensemble');
original = loaded.ensemble;
originalBudget = original.configuration.runtimeLimitSecondsPerTrial;
newBudget = 2 * originalBudget;

sourceStates = original.sourceStates;
seed = original.randomSeed;
regeneratedSourceStates = sample_sources(paper.source, size(sourceStates, 2), seed);
sourceStatesMatchOriginal = isequaln(sourceStates, regeneratedSourceStates);
if ~sourceStatesMatchOriginal
    maxMismatch = max(abs(sourceStates(:) - regeneratedSourceStates(:)));
    error('SchroerOtt:PaperSourceStatesMismatch', ...
        ['The saved 50 source states do not match regeneration from seed %g. ', ...
        'Maximum absolute mismatch is %.17g.'], seed, maxMismatch);
end

newCfg = baseCfg;
newCfg.runtimeLimitSecondsPerTrial = newBudget;
originalScientific = remove_runtime_metadata(original.configuration);
newScientific = remove_runtime_metadata(newCfg);
scientificSettingsPreserved = isequaln(originalScientific, newScientific);
fingerprintPayload = struct();
fingerprintPayload.scientificConfiguration = newScientific;
fingerprintPayload.paperMetadata = paper;
fingerprintPayload.originalRuntimeBudgetSeconds = originalBudget;
fingerprintPayload.newRuntimeBudgetSeconds = newBudget;
configurationFingerprint = so_configuration_fingerprint(fingerprintPayload);

previousTerminations = previous_termination_table(original, originalBudget);
previousTerminationCounts = termination_count_struct(previousTerminations.terminationReason);

audit.action = "estimate";
audit.originalFile = originalFile;
audit.originalOutputFiles = original_output_files(baseCfg.outputDirectory);
audit.outputPaths = paths;
audit.originalRuntimeBudgetSeconds = originalBudget;
audit.newRuntimeBudgetSeconds = newBudget;
audit.runtimeMultiplier = 2;
audit.budgetAppliesTo = "one complete source trial: route construction plus multistage targeting";
audit.timeoutEnforcement = ["so_runtime_exceeded(cfg)"; ...
    "so_multistage_targeting stage, switch-probe and final checks"; ...
    "so_resolve_connection total-time and split checks"; ...
    "timeout_connection(connection)"];
audit.independentLimits = independent_limits(newCfg, paper);
audit.originalControlAmplitude = original.configuration.controlAmplitude;
audit.newControlAmplitude = newCfg.controlAmplitude;
audit.baseConfiguration = newCfg;
audit.paperMetadata = paper;
audit.randomSeed = seed;
audit.sourceStates = sourceStates;
audit.regeneratedSourceStates = regeneratedSourceStates;
audit.sourceStatesMatchOriginal = sourceStatesMatchOriginal;
audit.scientificSettingsPreserved = scientificSettingsPreserved;
audit.configurationFingerprint = configurationFingerprint;
audit.previousTerminations = previousTerminations;
audit.previousTerminationCounts = previousTerminationCounts;
audit.availableWorkers = available_worker_count();
audit.selectedWorkersForFullRun = max(1, min(audit.availableWorkers, 50));
audit.fullRunCommand = "so_run_paper_ensemble_double_runtime('resumeDoubleRuntime')";
audit.latexStatusStatement = "None of the $50$ sampled sources reached the target within the runtime budget imposed in the initial paper-comparison ensemble. The outcome is therefore classified as computationally unresolved at the declared budget, rather than as evidence that the Schroer--Ott targeting construction fails for these sources.";
end

function sourceStates = sample_sources(rect, n, seed)
rng(seed);
sourceStates = [rect.xMin + (rect.xMax - rect.xMin) .* rand(1, n); ...
    rect.yMin + (rect.yMax - rect.yMin) .* rand(1, n)];
end

function cfg = remove_runtime_metadata(cfg)
if isfield(cfg, 'runtimeLimitSecondsPerTrial')
    cfg = rmfield(cfg, 'runtimeLimitSecondsPerTrial');
end
if isfield(cfg, 'runtime')
    cfg = rmfield(cfg, 'runtime');
end
end

function files = original_output_files(outputDirectory)
names = ["paper_50_source_ensemble.mat", "paper_50_source_ensemble_checkpoint.mat", ...
    "paper_50_source_trials.csv", "paper_50_source_summary.csv"];
files = strings(numel(names), 1);
for i = 1:numel(names)
    files(i) = string(fullfile(outputDirectory, names(i)));
end
end

function limits = independent_limits(cfg, paper)
limits.runtimeLimitSecondsPerSource = cfg.runtimeLimitSecondsPerTrial;
limits.maxForwardIterations = cfg.maxForwardIterations;
limits.maxBackwardIterations = cfg.maxBackwardIterations;
limits.maxTotalTransferTime = cfg.maxTotalTransferTime;
limits.maxStages = cfg.maxStages;
limits.maxControl = cfg.controlAmplitude;
limits.orbitMinPeriod = cfg.orbit.minPeriod;
limits.orbitMaxPeriod = cfg.orbit.maxPeriod;
limits.orbitNewtonMaxIterations = cfg.orbit.newtonMaxIterations;
limits.orbitNewtonTolerance = cfg.orbit.newtonTolerance;
limits.orbitSeedXCount = cfg.orbit.seedXCount;
limits.orbitSeedYOffsets = cfg.orbit.seedYOffsets;
limits.curveInitialControlSamples = cfg.curve.initialControlSamples;
limits.curveInitialBoundarySamples = cfg.curve.initialBoundarySamples;
limits.curveMaxPoints = cfg.curve.maxPoints;
limits.curveMaxSubdivisionDepth = cfg.curve.maxSubdivisionDepth;
limits.curveMaxGapFraction = cfg.curve.maxGapFraction;
limits.curveMidpointToleranceFraction = cfg.curve.midpointToleranceFraction;
limits.intersectionTolerance = cfg.intersectionTolerance;
limits.containmentTolerance = cfg.containmentTolerance;
limits.propagationConsistencyTolerance = cfg.propagationConsistencyTolerance;
limits.routeRotationNumbers = paper.expectedRotationNumbers;
limits.switchProbeRule = "j = 0:min(provisional.totalIterations, cfg.maxTotalTransferTime)";
limits.switchProbePruning = "so_should_prune_probe(j, Jbest)";
limits.publishedIterationRangeIsDiagnosticOnly = paper.publishedControlledRange;
end

function tbl = previous_termination_table(original, originalBudget)
old = original.trials;
n = height(old);
sourceState = strings(n, 1);
terminationReason = strings(n, 1);
legacyFailureCategory = strings(n, 1);
for i = 1:n
    sourceState(i) = mat2str([old.sourceX(i); old.sourceY(i)], 17);
    if old.success(i)
        legacyFailureCategory(i) = "";
        terminationReason(i) = "success";
    else
        legacyFailureCategory(i) = string(old.failureCategory(i));
        terminationReason(i) = map_legacy_failure(legacyFailureCategory(i));
    end
end
tbl = table(old.trial, sourceState, old.sourceX, old.sourceY, old.success, ...
    legacyFailureCategory, terminationReason, repmat(originalBudget, n, 1), ...
    old.runtimeSeconds, old.totalControlledIterations, old.numberOfControls, ...
    'VariableNames', {'sourceIndex','sourceState','sourceX','sourceY','found', ...
    'legacyFailureCategory','terminationReason','runtimeBudgetSeconds', ...
    'runtimeSeconds','iterations','controls'});
end

function reason = map_legacy_failure(legacy)
switch string(legacy)
    case "runtime_limit_exceeded"
        reason = "time_budget_exhausted";
    case {"all_switch_probes_failed", "switch_probe_budget_exhausted"}
        reason = "candidate_budget_exhausted";
    case {"provisional_connection_failed", "final_connection_failed", "not_reached"}
        reason = "no_eligible_route";
    case {"propagation_consistency_failed", "final_independent_containment_failed", ...
            "control_bound_violation"}
        reason = "invalid_replay";
    otherwise
        reason = "numerical_failure";
end
end

function counts = termination_count_struct(reasons)
names = ["success", "time_budget_exhausted", "iteration_limit_reached", ...
    "candidate_budget_exhausted", "no_eligible_route", "numerical_failure", ...
    "invalid_replay"];
counts = struct();
for i = 1:numel(names)
    counts.(char(names(i))) = sum(string(reasons) == names(i));
end
end

function idx = first_previous_timeout(previousTerminations)
mask = previousTerminations.terminationReason == "time_budget_exhausted";
idx = previousTerminations.sourceIndex(find(mask, 1, 'first'));
end

function ensemble = run_requested_sources(audit, opts, paths, action)
ensure_output_dir(paths.outputDirectory);
sourceIndices = opts.sourceIndices;
if isempty(sourceIndices)
    sourceIndices = 1:size(audit.sourceStates, 2);
end
sourceIndices = unique(sourceIndices, 'stable');
if any(sourceIndices < 1) || any(sourceIndices > size(audit.sourceStates, 2))
    error('SchroerOtt:InvalidSourceIndex', 'Source indices must be in the range 1:50.');
end
runtimeBudgetSeconds = audit.newRuntimeBudgetSeconds;
if isfinite(opts.runtimeBudgetSeconds)
    runtimeBudgetSeconds = opts.runtimeBudgetSeconds;
end
if strcmp(action, "resumeDoubleRuntime") && numel(sourceIndices) == 50 && ...
        abs(runtimeBudgetSeconds - audit.newRuntimeBudgetSeconds) > 10 * eps(audit.newRuntimeBudgetSeconds)
    error('SchroerOtt:FullRunBudgetOverride', ...
        'The complete doubled-runtime ensemble must use %.17g seconds per source.', ...
        audit.newRuntimeBudgetSeconds);
end

metadata = run_metadata(audit, paths, runtimeBudgetSeconds, sourceIndices, action);
[records, results, completedMask] = load_or_initialize_checkpoint(paths, audit, metadata, opts);
pending = sourceIndices(~completedMask(sourceIndices));
skippedCompleted = sourceIndices(completedMask(sourceIndices));

baseCfg = audit.baseConfiguration;
baseCfg.runtimeLimitSecondsPerTrial = runtimeBudgetSeconds;
baseCfg.outputDirectory = paths.outputDirectory;
baseCfg.figureDirectory = fullfile(paths.outputDirectory, 'figures');
baseCfg.saveFigures = false;

log_message(paths, opts, sprintf('%s: %d requested, %d already complete, %d pending', ...
    action, numel(sourceIndices), numel(skippedCompleted), numel(pending)));

runTimer = tic;
catalogue = so_enumerate_periodic_orbits(baseCfg);
workerInfo = choose_workers(opts, numel(pending));

if workerInfo.useParallel
    try
        [records, results, completedMask] = run_parallel(pending, records, results, ...
            completedMask, audit, baseCfg, catalogue, runtimeBudgetSeconds, metadata, ...
            paths, opts, runTimer, workerInfo);
    catch err
        log_message(paths, opts, sprintf('parallel fallback: %s', err.message));
        workerInfo.useParallel = false;
        workerInfo.selectedWorkers = 1;
        [records, results, completedMask] = run_serial(pending, records, results, ...
            completedMask, audit, baseCfg, catalogue, runtimeBudgetSeconds, metadata, ...
            paths, opts, runTimer, workerInfo);
    end
else
    [records, results, completedMask] = run_serial(pending, records, results, ...
        completedMask, audit, baseCfg, catalogue, runtimeBudgetSeconds, metadata, ...
        paths, opts, runTimer, workerInfo);
end

newlyCompleted = pending(completedMask(pending));
ensemble = finalize_ensemble(records, results, completedMask, audit, metadata, paths, ...
    sourceIndices, skippedCompleted, newlyCompleted, workerInfo, toc(runTimer));
if opts.saveOutputs
    save_final_outputs(paths, ensemble);
end
end

function metadata = run_metadata(audit, paths, runtimeBudgetSeconds, sourceIndices, action)
metadata.action = string(action);
metadata.createdAt = string(datetime('now'));
metadata.outputPaths = paths;
metadata.originalRuntimeBudgetSeconds = audit.originalRuntimeBudgetSeconds;
metadata.newRuntimeBudgetSeconds = audit.newRuntimeBudgetSeconds;
metadata.activeRuntimeBudgetSeconds = runtimeBudgetSeconds;
metadata.runtimeMultiplier = audit.runtimeMultiplier;
metadata.randomSeed = audit.randomSeed;
metadata.sourceIndices = sourceIndices;
metadata.configurationFingerprint = audit.configurationFingerprint;
metadata.scientificSettingsPreserved = audit.scientificSettingsPreserved;
metadata.sourceStatesMatchOriginal = audit.sourceStatesMatchOriginal;
metadata.fullRunCommand = audit.fullRunCommand;
end

function [records, results, completedMask] = load_or_initialize_checkpoint(paths, audit, metadata, opts)
recordTemplate = so_empty_paper_ensemble_record();
records = repmat(recordTemplate, size(audit.sourceStates, 2), 1);
results = cell(size(audit.sourceStates, 2), 1);
completedMask = false(size(audit.sourceStates, 2), 1);
if opts.saveOutputs && isfile(paths.checkpointMat) && ~opts.force
    checkpoint = load(paths.checkpointMat, 'records', 'results', 'completedMask', ...
        'sourceStates', 'metadata');
    if ~isequaln(checkpoint.sourceStates, audit.sourceStates)
        error('SchroerOtt:CheckpointSourceMismatch', ...
            'Checkpoint source states do not match the preserved original source states.');
    end
    if isfield(checkpoint, 'metadata') && ...
            abs(checkpoint.metadata.activeRuntimeBudgetSeconds - metadata.activeRuntimeBudgetSeconds) > ...
            10 * eps(metadata.activeRuntimeBudgetSeconds)
        error('SchroerOtt:CheckpointBudgetMismatch', ...
            'Checkpoint runtime budget does not match this run. Use Force=true to start fresh.');
    end
    records = checkpoint.records;
    results = checkpoint.results;
    completedMask = checkpoint.completedMask;
elseif opts.saveOutputs && isfile(paths.checkpointMat) && opts.force
    log_message(paths, opts, 'Force=true: existing checkpoint ignored for this output stem');
end
end

function info = choose_workers(opts, pendingCount)
available = available_worker_count();
info.availableWorkers = available;
info.selectedWorkers = 1;
info.useParallel = false;
if pendingCount < 2
    return;
end
requested = opts.useParallel;
if islogical(requested)
    wantParallel = requested;
else
    wantParallel = any(lower(string(requested)) == ["auto", "true", "yes"]);
end
if wantParallel && available > 1
    info.useParallel = true;
    info.selectedWorkers = min(available, pendingCount);
    if ~isempty(opts.numWorkers)
        info.selectedWorkers = max(1, min([opts.numWorkers, available, pendingCount]));
    end
end
end

function n = available_worker_count()
n = 0;
try
    if license('test', 'Distrib_Computing_Toolbox')
        cluster = parcluster('local');
        n = cluster.NumWorkers;
    end
catch
    n = 0;
end
end

function [records, results, completedMask] = run_serial(pending, records, results, ...
    completedMask, audit, baseCfg, catalogue, runtimeBudgetSeconds, metadata, paths, opts, runTimer, workerInfo)
for i = 1:numel(pending)
    idx = pending(i);
    [record, result] = so_paper_ensemble_double_runtime_trial(idx, ...
        audit.sourceStates(:, idx), baseCfg, audit.paperMetadata, catalogue, ...
        runtimeBudgetSeconds, audit.randomSeed, audit.configurationFingerprint);
    records(idx) = record;
    results{idx} = result;
    completedMask(idx) = true;
    save_checkpoint_atomic(paths.checkpointMat, records, results, completedMask, ...
        audit.sourceStates, metadata, opts);
    report_progress(records, completedMask, idx, paths, opts, runTimer, workerInfo);
end
end

function [records, results, completedMask] = run_parallel(pending, records, results, ...
    completedMask, audit, baseCfg, catalogue, runtimeBudgetSeconds, metadata, paths, opts, runTimer, workerInfo)
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', workerInfo.selectedWorkers);
elseif pool.NumWorkers < workerInfo.selectedWorkers
    workerInfo.selectedWorkers = pool.NumWorkers;
end
futures(numel(pending), 1) = parallel.FevalFuture;
for i = 1:numel(pending)
    idx = pending(i);
    futures(i) = parfeval(pool, @so_paper_ensemble_double_runtime_trial, 2, idx, ...
        audit.sourceStates(:, idx), baseCfg, audit.paperMetadata, catalogue, ...
        runtimeBudgetSeconds, audit.randomSeed, audit.configurationFingerprint);
end
remaining = numel(futures);
while remaining > 0
    [completedFutureIndex, record, result] = fetchNext(futures);
    idx = pending(completedFutureIndex);
    records(idx) = record;
    results{idx} = result;
    completedMask(idx) = true;
    save_checkpoint_atomic(paths.checkpointMat, records, results, completedMask, ...
        audit.sourceStates, metadata, opts);
    report_progress(records, completedMask, idx, paths, opts, runTimer, workerInfo);
    remaining = remaining - 1;
end
end

function save_checkpoint_atomic(checkpointFile, records, results, completedMask, sourceStates, metadata, opts)
if ~opts.saveOutputs
    return;
end
tmp = [tempname(fileparts(checkpointFile)), '.mat'];
lastUpdated = string(datetime('now')); %#ok<NASGU>
save(tmp, 'records', 'results', 'completedMask', 'sourceStates', 'metadata', ...
    'lastUpdated', '-v7.3');
movefile(tmp, checkpointFile, 'f');
end

function report_progress(records, completedMask, idx, paths, opts, runTimer, workerInfo)
completed = sum(completedMask);
successCount = sum([records(completedMask).found]);
reasons = string({records(completedMask).terminationReason});
timeoutCount = sum(reasons == "time_budget_exhausted");
elapsed = toc(runTimer);
remaining = numel(completedMask) - completed;
if workerInfo.useParallel
    estimatedRemaining = elapsed / max(completed, 1) * remaining / max(workerInfo.selectedWorkers, 1);
else
    estimatedRemaining = elapsed / max(completed, 1) * remaining;
end
msg = sprintf(['completed %d/50 | source %d | %s | elapsed %.1f s | ', ...
    'successes %d | timeouts %d | estimated remaining %.1f s'], ...
    completed, idx, records(idx).terminationReason, elapsed, successCount, ...
    timeoutCount, estimatedRemaining);
log_message(paths, opts, msg);
end

function ensemble = finalize_ensemble(records, results, completedMask, audit, metadata, paths, ...
    sourceIndices, skippedCompleted, newlyCompleted, workerInfo, runtimeSeconds)
trialTable = records_to_table(records(completedMask));
summary = ensemble_summary(trialTable, sourceIndices);
runtimeComparison = runtime_comparison_table(audit.previousTerminations, trialTable, ...
    audit.originalRuntimeBudgetSeconds, metadata.activeRuntimeBudgetSeconds);
ensemble.metadata = metadata;
ensemble.outputPaths = paths;
ensemble.configuration = audit.baseConfiguration;
ensemble.paperMetadata = audit.paperMetadata;
ensemble.randomSeed = audit.randomSeed;
ensemble.sourceStates = audit.sourceStates;
ensemble.sourceStatesMatchOriginal = audit.sourceStatesMatchOriginal;
ensemble.scientificSettingsPreserved = audit.scientificSettingsPreserved;
ensemble.configurationFingerprint = audit.configurationFingerprint;
ensemble.records = records;
ensemble.results = results;
ensemble.completedMask = completedMask;
ensemble.requestedSourceIndices = sourceIndices;
ensemble.skippedCompletedIndices = skippedCompleted;
ensemble.newlyCompletedIndices = newlyCompleted;
ensemble.trials = trialTable;
ensemble.summary = summary;
ensemble.terminationCounts = termination_count_table(trialTable);
ensemble.previousTerminations = audit.previousTerminations;
ensemble.previousTerminationCounts = audit.previousTerminationCounts;
ensemble.runtimeComparison = runtimeComparison;
ensemble.independentLimits = audit.independentLimits;
ensemble.workerInfo = workerInfo;
ensemble.runtimeSeconds = runtimeSeconds;
ensemble.fullRunCommand = audit.fullRunCommand;
end

function tbl = records_to_table(records)
names = {'sourceIndex','sourceState','sourceX','sourceY','found','terminationReason', ...
    'legacyFailureCategory','runtimeBudgetSeconds','runtimeSeconds','iterations', ...
    'controls','controlMagnitudes','selectedRoute','switchPoints','finalX','finalY', ...
    'finalTargetDistance','seed','configurationFingerprint','maxAbsControl', ...
    'totalAbsControl','maxIntersectionResidual','maxPathwiseReplayError', ...
    'targetContained','minimumCertified','unresolvedSplits','errorIdentifier','errorMessage'};
if isempty(records)
    tbl = cell2table(cell(0, numel(names)), 'VariableNames', names);
    return;
end
rows = cell(numel(records), numel(names));
for i = 1:numel(records)
    r = records(i);
    rows(i, :) = {r.sourceIndex, mat2str(r.sourceState, 17), r.sourceX, r.sourceY, ...
        r.found, string(r.terminationReason), string(r.legacyFailureCategory), ...
        r.runtimeBudgetSeconds, r.runtimeSeconds, r.iterations, r.controls, ...
        mat2str(r.controlMagnitudes, 17), mat2str(r.selectedRoute, 17), ...
        mat2str(r.switchPoints, 17), r.finalX, r.finalY, r.finalTargetDistance, ...
        r.seed, string(r.configurationFingerprint), r.maxAbsControl, ...
        r.totalAbsControl, r.maxIntersectionResidual, r.maxPathwiseReplayError, ...
        r.targetContained, r.minimumCertified, r.unresolvedSplits, ...
        string(r.errorIdentifier), string(r.errorMessage)};
end
tbl = cell2table(rows, 'VariableNames', names);
end

function summary = ensemble_summary(trialTable, sourceIndices)
if isempty(trialTable) || height(trialTable) == 0
    stats = {numel(sourceIndices), 0, 0, 0, 0, 0, NaN, NaN, NaN, NaN};
else
    successes = sum(trialTable.found);
    timeouts = sum(trialTable.terminationReason == "time_budget_exhausted");
    finiteCounts = trialTable.iterations(trialTable.found & isfinite(trialTable.iterations));
    if isempty(finiteCounts)
        minIter = NaN;
        maxIter = NaN;
        meanIter = NaN;
        medianIter = NaN;
    else
        minIter = min(finiteCounts);
        maxIter = max(finiteCounts);
        meanIter = mean(finiteCounts);
        medianIter = median(finiteCounts);
    end
    stats = {numel(sourceIndices), height(trialTable), successes, ...
        successes / max(height(trialTable), 1), timeouts, ...
        sum(trialTable.terminationReason ~= "success" & ...
        trialTable.terminationReason ~= "time_budget_exhausted"), ...
        minIter, maxIter, meanIter, medianIter};
end
summary = cell2table(stats, 'VariableNames', {'requestedSources','completedSources', ...
    'successes','successFractionAmongCompleted','timeBudgetExhausted', ...
    'nonTimeoutTerminations','minimumIterations','maximumIterations', ...
    'meanIterations','medianIterations'});
end

function counts = termination_count_table(trialTable)
reasons = ["success"; "time_budget_exhausted"; "iteration_limit_reached"; ...
    "candidate_budget_exhausted"; "no_eligible_route"; "numerical_failure"; ...
    "invalid_replay"];
count = zeros(numel(reasons), 1);
for i = 1:numel(reasons)
    if isempty(trialTable) || height(trialTable) == 0
        count(i) = 0;
    else
        count(i) = sum(trialTable.terminationReason == reasons(i));
    end
end
counts = table(reasons, count, 'VariableNames', {'terminationReason','count'});
end

function comparison = runtime_comparison_table(previous, current, originalBudget, newBudget)
n = height(previous);
rows = cell(n, 11);
for i = 1:n
    idx = previous.sourceIndex(i);
    newMask = false;
    if ~isempty(current) && height(current) > 0
        newMask = current.sourceIndex == idx;
    end
    if any(newMask)
        row = current(find(newMask, 1, 'first'), :);
        newFound = row.found;
        newReason = row.terminationReason;
        newRuntime = row.runtimeSeconds;
    else
        newFound = false;
        newReason = "";
        newRuntime = NaN;
    end
    rows(i, :) = {idx, previous.sourceState(i), previous.found(i), ...
        previous.terminationReason(i), originalBudget, previous.runtimeSeconds(i), ...
        newFound, newReason, newBudget, newRuntime, newRuntime - previous.runtimeSeconds(i)};
end
comparison = cell2table(rows, 'VariableNames', {'sourceIndex','sourceState', ...
    'originalFound','originalTerminationReason','originalRuntimeBudgetSeconds', ...
    'originalRuntimeSeconds','newFound','newTerminationReason', ...
    'newRuntimeBudgetSeconds','newRuntimeSeconds','runtimeSecondsDelta'});
end

function save_final_outputs(paths, ensemble)
save(paths.finalMat, 'ensemble', '-v7.3');
writetable(ensemble.trials, paths.trialsCsv);
writetable(ensemble.runtimeComparison, paths.comparisonCsv);
log_message(paths, struct('quiet', true, 'saveOutputs', true), ...
    sprintf('saved %s with %d completed sources', paths.finalMat, height(ensemble.trials)));
end

function log_message(paths, opts, msg)
if ~isfield(opts, 'quiet') || ~opts.quiet
    fprintf('%s\n', msg);
end
if isfield(opts, 'saveOutputs') && opts.saveOutputs
    ensure_output_dir(fileparts(paths.log));
    fid = fopen(paths.log, 'a');
    if fid >= 0
        cleaner = onCleanup(@() fclose(fid));
        fprintf(fid, '[%s] %s\n', char(datetime('now')), msg);
        clear cleaner;
    end
end
end

function ensure_output_dir(folder)
if ~exist(folder, 'dir')
    mkdir(folder);
end
end
