function results = run_noise_sweep( ...
    sourceState, target, params, noiseControl, noise, sweep, initialSearch)
%RUN_NOISE_SWEEP Completion-driven bisection-only noise sweep.
%
% The final experiment contains 110 predetermined records.  A pilot run is
% not a disposable smoke sweep: it runs trial 1 at each configured noise
% level and writes those eleven trials into the genuine final checkpoint.

verifiedIdentifier = verified_bisection_identifier();
validate_sweep_inputs(noiseControl, noise, sweep, verifiedIdentifier);
if nargin < 7 || isempty(initialSearch)
    initialSearch = search_parameter_to_target( ...
        sourceState, target, params, noiseControl, verifiedIdentifier);
end
assert_verified_bisection_search( ...
    initialSearch, verifiedIdentifier, 'noise sweep initial search', target);

config = build_bisection_noise_config( ...
    sourceState, target, params, noiseControl, noise, sweep);
sigmaNoiseValues = sweep.sigmaNoiseValues(:).';
modeList = sweep.modes(:).';
numTrials = sweep.numTrials;
totalTrials = config.totalTrials;
trialPlan = make_trial_plan(sweep.baseSeed, modeList, ...
    sigmaNoiseValues, numTrials);
runMode = sweep_run_mode(sweep);
selectedRecordIndices = selected_records_for_run(runMode, trialPlan);

parallelInfo = sweep_parallel_info(sweep);
[parallelInfo, pool] = ensure_process_pool(parallelInfo); %#ok<ASGLU>
execution = parallelInfo;
execution.scheduler = 'parfeval/fetchNext completion-driven';
execution.runMode = runMode;
execution.selectedRecordIndices = selectedRecordIndices;
execution.clientOnlyWrites = true;

fprintf(['Noise sweep execution mode: %s | selected workers: %d | ', ...
    'configured trials: %d | requested this run: %d\n'], ...
    parallelInfo.executionMode, parallelInfo.selectedWorkers, ...
    totalTrials, numel(selectedRecordIndices));

[trials, completed, runStartedAt] = initialize_or_resume( ...
    sourceState, target, sweep, config, totalTrials);
pendingRecords = selectedRecordIndices(~completed(selectedRecordIndices));
runTimer = tic;

if isempty(pendingRecords)
    fprintf('All requested records for %s are already complete.\n', runMode);
else
    if strcmp(parallelInfo.executionMode, 'parallel')
        [trials, completed] = run_parallel_records(pendingRecords, ...
            trialPlan, trials, completed, sourceState, target, params, ...
            noiseControl, noise, initialSearch, config, sweep, ...
            runStartedAt, execution, runTimer);
    else
        [trials, completed] = run_serial_records(pendingRecords, ...
            trialPlan, trials, completed, sourceState, target, params, ...
            noiseControl, noise, initialSearch, config, sweep, ...
            runStartedAt, execution, runTimer);
    end
end

results = assemble_results(trials, completed, config, sweep, ...
    initialSearch, sourceState, target, noise, noiseControl, execution);

checkpoint = make_checkpoint(trials, completed, config, sweep, ...
    runStartedAt, all(completed), execution);
save_checkpoint(checkpoint, sweep.checkpointMatFile);

if all(completed)
    save_results(results, sweep);
else
    fprintf(['Partial final checkpoint saved with %d/%d completed records. ', ...
        'No final MAT/CSV result was written.\n'], ...
        sum(completed), numel(completed));
end
end

function validate_sweep_inputs(noiseControl, noise, sweep, verifiedIdentifier)
if ~isfield(noiseControl, 'searchMethod') || ...
        ~strcmp(noiseControl.searchMethod, verifiedIdentifier)
    error('noiseControl.searchMethod must be %s.', verifiedIdentifier);
end
if ~isfield(sweep, 'searchMethod') || ...
        ~strcmp(sweep.searchMethod, verifiedIdentifier)
    error('sweep.searchMethod must be %s.', verifiedIdentifier);
end
expectedSigmaNoiseValues = [0, 0.01, 0.03, 0.05, 0.075, 0.1, ...
    0.15, 0.2, 0.3, 0.5, 1.0];
if ~isfield(sweep, 'sigmaNoiseValues') || ...
        ~isequal(sweep.sigmaNoiseValues(:).', expectedSigmaNoiseValues)
    error(['run_noise_sweep requires the final eleven sigma_noise ', ...
        'values: ', mat2str(expectedSigmaNoiseValues), '.']);
end
if ~isfield(sweep, 'numTrials') || sweep.numTrials ~= 10
    error('run_noise_sweep requires sweep.numTrials == 10.');
end
if ~isfield(sweep, 'modes') || ...
        ~isequal(sweep.modes(:).', {'paperCycleCadenceRetargetMode'})
    error('The repaired final sweep uses only paperCycleCadenceRetargetMode.');
end
if sweep.numTrials * numel(sweep.sigmaNoiseValues) * ...
        numel(sweep.modes) ~= 110
    error('run_noise_sweep requires exactly 110 configured trials.');
end
if ~isfield(noise, 'maxControlledCrossings') || ...
        noise.maxControlledCrossings ~= 60
    error('run_noise_sweep requires noise.maxControlledCrossings == 60.');
end
if ~isfield(noise, 'retargetCadencePolicy') || ...
        ~strcmp(noise.retargetCadencePolicy, 'meanReturnCycleMatched')
    error('Cycle-cadence final sweep requires meanReturnCycleMatched policy.');
end
if ~isfield(noise, 'retargetStepInterval') || ...
        ~isfinite(noise.retargetStepInterval) || ...
        noise.retargetStepInterval < 1
    error('Cycle-cadence final sweep requires a finite retargetStepInterval.');
end
if ~isfield(noise, 'cadenceRatio') || abs(noise.cadenceRatio - 1) > 0.05
    error('Cycle-cadence final sweep requires abs(cadenceRatio - 1) <= 0.05.');
end
if ~isfield(noise, 'maxRk4Steps') || ~isfinite(noise.maxRk4Steps)
    error('noise.maxRk4Steps must be one finite full-trial RK4-step budget.');
end
if isfield(noise, 'requiredRk4Steps') && ...
        noise.maxRk4Steps < noise.requiredRk4Steps
    error('noise.maxRk4Steps is below the measured required RK4-step budget.');
end
if ~isfield(sweep, 'checkpointMatFile') || ...
        ~isfield(sweep, 'outputMatFile') || ...
        ~isfield(sweep, 'summaryCsvFile') || ...
        ~isfield(sweep, 'trialsCsvFile')
    error('sweep output filenames are required.');
end
if ~isfield(sweep, 'isSmokeTest') || sweep.isSmokeTest
    error('run_noise_sweep only accepts the full configuration: sweep.isSmokeTest must be false.');
end
end

function runMode = sweep_run_mode(sweep)
runMode = 'resumeFull';
if isfield(sweep, 'runMode') && ~isempty(sweep.runMode)
    runMode = char(string(sweep.runMode));
end
valid = {'pilot', 'resumeFull', 'regenerateFull'};
if ~any(strcmp(runMode, valid))
    error('Unknown sweep.runMode: %s', runMode);
end
end

function selected = selected_records_for_run(runMode, trialPlan)
switch runMode
    case 'pilot'
        selected = find([trialPlan.trialIndex] == 1);
    otherwise
        selected = 1:numel(trialPlan);
end
end

function [trials, completed] = run_serial_records(pendingRecords, ...
    trialPlan, trials, completed, sourceState, target, params, ...
    noiseControl, noise, initialSearch, config, sweep, runStartedAt, ...
    execution, runTimer)
for k = 1:numel(pendingRecords)
    recordIndex = pendingRecords(k);
    plan = trialPlan(recordIndex);
    trial = run_trial_from_plan(plan, sourceState, target, params, ...
        noiseControl, noise, initialSearch);
    [trials, completed] = store_completed_trial(trials, completed, ...
        recordIndex, trial, plan, config, sweep, runStartedAt, ...
        execution, runTimer);
end
end

function [trials, completed] = run_parallel_records(pendingRecords, ...
    trialPlan, trials, completed, sourceState, target, params, ...
    noiseControl, noise, initialSearch, config, sweep, runStartedAt, ...
    execution, runTimer)

pool = gcp('nocreate');
if isempty(pool)
    error('Parallel execution was requested but no process pool is active.');
end

maxActive = max(1, execution.selectedWorkers);
nextToSubmit = 1;
futures = parallel.FevalFuture.empty(0, 1);
futureRecordIndices = zeros(0, 1);

while nextToSubmit <= numel(pendingRecords) && ...
        numel(futures) < min(maxActive, numel(pendingRecords))
    [futures, futureRecordIndices, nextToSubmit] = submit_next_future( ...
        futures, futureRecordIndices, nextToSubmit, pendingRecords, ...
        trialPlan, sourceState, target, params, noiseControl, noise, ...
        initialSearch);
end

fetched = 0;
while fetched < numel(pendingRecords)
    [futureIndex, trial] = fetchNext(futures);
    fetched = fetched + 1;
    recordIndex = futureRecordIndices(futureIndex);
    plan = trialPlan(recordIndex);
    [trials, completed] = store_completed_trial(trials, completed, ...
        recordIndex, trial, plan, config, sweep, runStartedAt, ...
        execution, runTimer);

    if nextToSubmit <= numel(pendingRecords)
        [futures, futureRecordIndices, nextToSubmit] = submit_next_future( ...
            futures, futureRecordIndices, nextToSubmit, pendingRecords, ...
            trialPlan, sourceState, target, params, noiseControl, noise, ...
            initialSearch);
    end
end
end

function [futures, futureRecordIndices, nextToSubmit] = submit_next_future( ...
    futures, futureRecordIndices, nextToSubmit, pendingRecords, ...
    trialPlan, sourceState, target, params, noiseControl, noise, ...
    initialSearch)
recordIndex = pendingRecords(nextToSubmit);
plan = trialPlan(recordIndex);
futures(end + 1, 1) = parfeval(@run_noisy_targeting_trial, 1, ...
    sourceState, target, params, noiseControl, noise, plan.sigmaNoise, ...
    plan.trialIndex, plan.seed, plan.mode, initialSearch, false);
futureRecordIndices(end + 1, 1) = recordIndex;
nextToSubmit = nextToSubmit + 1;
end

function [trials, completed] = store_completed_trial(trials, completed, ...
    recordIndex, trial, plan, config, sweep, runStartedAt, execution, ...
    runTimer)
validate_trial_against_plan(trial, plan);
validate_bisection_trial(trial, config);
if ~should_retain_full_trajectory(plan, config.sigmaNoiseValues)
    trial = strip_full_trajectory(trial);
end
trials(recordIndex) = trial;
completed(recordIndex) = true;

checkpoint = make_checkpoint(trials, completed, config, sweep, ...
    runStartedAt, false, execution);
save_checkpoint(checkpoint, sweep.checkpointMatFile);
line = progress_line(completed, numel(completed), plan, trial, runTimer);
fprintf('%s\n', line);
append_progress_log(sweep, line);
end

function trial = run_trial_from_plan(plan, sourceState, target, params, ...
    noiseControl, noise, initialSearch)
trial = run_noisy_targeting_trial( ...
    sourceState, target, params, noiseControl, noise, ...
    plan.sigmaNoise, plan.trialIndex, plan.seed, plan.mode, ...
    initialSearch, false);
end

function validate_trial_against_plan(trial, plan)
if trial.trialNumber ~= plan.trialIndex
    error('Completed trialNumber does not match the predetermined plan.');
end
if trial.seed ~= plan.seed
    error('Completed trial seed does not match the predetermined plan.');
end
if abs(trial.sigmaNoise - plan.sigmaNoise) > 10*eps
    error('Completed trial sigmaNoise does not match the predetermined plan.');
end
if ~strcmp(trial.mode, plan.mode)
    error('Completed trial mode does not match the predetermined plan.');
end
end

function tf = should_retain_full_trajectory(plan, sigmaNoiseValues)
middleIndex = ceil((numel(sigmaNoiseValues) + 1) / 2);
retainedSigmaIndices = unique([1, middleIndex, numel(sigmaNoiseValues)]);
tf = plan.trialIndex == 1 && ...
    any(plan.sigmaIndex == retainedSigmaIndices);
end

function trial = strip_full_trajectory(trial)
trial.t = [];
trial.X = [];
end

function [trials, completed, runStartedAt] = initialize_or_resume( ...
    sourceState, target, sweep, config, totalTrials)
runStartedAt = datetime('now');
template = empty_trial_template(sourceState, target);
trials(1, totalTrials) = template;
completed = false(1, totalTrials);
if ~isfield(sweep, 'resumeFromCheckpoint') || ~sweep.resumeFromCheckpoint || ...
        ~isfile(sweep.checkpointMatFile)
    return;
end

loaded = load(sweep.checkpointMatFile);
if ~isfield(loaded, 'checkpoint')
    return;
end
checkpoint = loaded.checkpoint;
if ~isfield(checkpoint, 'configuration') || ...
        ~isfield(checkpoint, 'trials') || ~isfield(checkpoint, 'completed')
    return;
end
[valid, differences] = validate_bisection_noise_config( ...
    checkpoint.configuration, config);
if ~valid
    error('Checkpoint configuration mismatch: %s', ...
        strjoin(cellstr(differences), '; '));
end
if numel(checkpoint.trials) ~= totalTrials || ...
        numel(checkpoint.completed) ~= totalTrials
    error('Checkpoint trial dimensions do not match the current sweep.');
end
trials = checkpoint.trials;
completed = checkpoint.completed;
for i = find(completed)
    validate_bisection_trial(trials(i), config);
end
if isfield(checkpoint, 'runStartedAt')
    runStartedAt = checkpoint.runStartedAt;
end
end

function recordIndex = linear_record_index( ...
    modeIndex, sigmaIndex, trialIndex, numSigmas, numTrials)
recordIndex = (modeIndex - 1) * numSigmas * numTrials + ...
    (sigmaIndex - 1) * numTrials + trialIndex;
end

function trialPlan = make_trial_plan(baseSeed, modeList, ...
    sigmaNoiseValues, numTrials)
totalTrials = numel(modeList) * numel(sigmaNoiseValues) * numTrials;
template = struct('recordIndex', NaN, 'modeIndex', NaN, ...
    'sigmaIndex', NaN, 'trialIndex', NaN, 'mode', '', ...
    'sigmaNoise', NaN, 'seed', NaN);
trialPlan(1, totalTrials) = template;
for modeIndex = 1:numel(modeList)
    for sigmaIndex = 1:numel(sigmaNoiseValues)
        for trialIndex = 1:numTrials
            recordIndex = linear_record_index(modeIndex, sigmaIndex, ...
                trialIndex, numel(sigmaNoiseValues), numTrials);
            trialPlan(recordIndex).recordIndex = recordIndex;
            trialPlan(recordIndex).modeIndex = modeIndex;
            trialPlan(recordIndex).sigmaIndex = sigmaIndex;
            trialPlan(recordIndex).trialIndex = trialIndex;
            trialPlan(recordIndex).mode = modeList{modeIndex};
            trialPlan(recordIndex).sigmaNoise = sigmaNoiseValues(sigmaIndex);
            trialPlan(recordIndex).seed = trial_seed( ...
                baseSeed, modeIndex, sigmaIndex, trialIndex);
        end
    end
end
end

function seed = trial_seed(baseSeed, modeIndex, sigmaIndex, trialIndex)
seed = baseSeed + 100000 * modeIndex + 1000 * sigmaIndex + trialIndex;
end

function parallelInfo = sweep_parallel_info(sweep)
if isfield(sweep, 'parallel') && isstruct(sweep.parallel) && ...
        isfield(sweep.parallel, 'selectedWorkers')
    parallelInfo = sweep.parallel;
else
    requestedWorkers = 'auto';
    if isfield(sweep, 'parallel') && isstruct(sweep.parallel) && ...
            isfield(sweep.parallel, 'requestedWorkers')
        requestedWorkers = sweep.parallel.requestedWorkers;
    end
    parallelInfo = select_noise_parallel_workers(requestedWorkers);
end
end

function [parallelInfo, pool] = ensure_process_pool(parallelInfo)
pool = [];
if ~isfield(parallelInfo, 'hasParallelToolbox') || ...
        ~parallelInfo.hasParallelToolbox || ...
        parallelInfo.selectedWorkers <= 1
    parallelInfo.executionMode = 'serial';
    parallelInfo.selectedWorkers = 1;
    return;
end

try
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers ~= parallelInfo.selectedWorkers
        delete(pool);
        pool = [];
    end
    if isempty(pool)
        pool = parpool('Processes', parallelInfo.selectedWorkers);
    end
    parallelInfo.executionMode = 'parallel';
catch problem
    warning(['Process-pool creation failed; falling back to serial ', ...
        'noise trials. Reason: %s'], problem.message);
    parallelInfo.executionMode = 'serial';
    parallelInfo.selectedWorkers = 1;
    parallelInfo.poolCreationWarning = problem.message;
end
end

function checkpoint = make_checkpoint( ...
    trials, completed, config, sweep, runStartedAt, complete, execution)
checkpoint.trials = trials;
checkpoint.completed = completed;
checkpoint.configuration = config;
checkpoint.configurationSignature = config.signatureText;
checkpoint.sweep = sweep;
checkpoint.verifiedSearchMethod = config.verifiedSearchMethod;
checkpoint.retargetingMode = config.retargetingMode;
checkpoint.completedTrialCount = sum(completed);
checkpoint.totalTrials = numel(completed);
checkpoint.nextIncompleteTrial = find(~completed, 1, 'first');
checkpoint.complete = complete && all(completed);
checkpoint.isSmokeTest = false;
checkpoint.isPilotOnly = ~checkpoint.complete;
checkpoint.trajectoryRetentionPolicy = sweep.trajectoryRetentionPolicy;
checkpoint.execution = execution;
checkpoint.runStartedAt = runStartedAt;
checkpoint.lastSavedAt = datetime('now');
end

function save_checkpoint(checkpoint, checkpointFile)
temporaryFile = [checkpointFile, '.tmp'];
save(temporaryFile, 'checkpoint', '-v7.3');
movefile(temporaryFile, checkpointFile, 'f');
end

function append_progress_log(sweep, line)
if ~isfield(sweep, 'progressLogFile') || isempty(sweep.progressLogFile)
    return;
end
fid = fopen(sweep.progressLogFile, 'a');
if fid < 0
    warning('Could not append to progress log: %s', sweep.progressLogFile);
    return;
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s | %s\n', char(datetime('now')), line);
end

function line = progress_line(completed, totalTrials, plan, trial, runTimer)
completedCount = sum(completed);
elapsedMinutes = toc(runTimer) / 60;
meanMinutes = elapsedMinutes / max(completedCount, 1);
remainingMinutes = meanMinutes * (totalTrials - completedCount);
line = sprintf(['Trial %d/%d | record %d | sigma_noise = %.3g | ', ...
    'realization %d/10 | success = %d | crossings = %g | ', ...
    'steps = %d | retargets = %d | seconds = %.1f | ', ...
    'elapsed = %.1f min | ETA = %.1f min'], ...
    completedCount, totalTrials, plan.recordIndex, plan.sigmaNoise, ...
    plan.trialIndex, trial.hit, trial.numCrossingsToTarget, ...
    trial.numRk4Steps, trial.numRetargetingEvents, ...
    trial.computationalTime, elapsedMinutes, remainingMinutes);
end

function results = assemble_results(trials, completed, config, sweep, ...
    initialSearch, sourceState, target, noise, noiseControl, execution)
completedTrials = trials(completed);
if isempty(completedTrials)
    trialsTable = table();
    summaryTable = table();
else
    trialsTable = make_trials_table(completedTrials);
    summaryTable = make_summary_table(completedTrials, ...
        config.sigmaNoiseValues, config.modes);
end

results.trials = trials;
results.completed = all(completed);
results.completedTrialCount = sum(completed);
results.completedRecordIndices = find(completed);
results.isPilotOnly = ~results.completed;
results.trialsTable = trialsTable;
results.summaryTable = summaryTable;
results.sigmaNoiseValues = config.sigmaNoiseValues;
results.numTrials = config.numTrials;
results.modes = config.modes;
results.mode = strjoin(config.modes, ', ');
results.initialSearch = initialSearch;
results.sourceState = sourceState;
results.target = target;
results.noise = noise;
results.control = noiseControl;
results.searchMethod = config.verifiedSearchMethod;
results.verifiedBisectionIdentifier = config.verifiedSearchMethod;
results.configuration = config;
results.configurationSignature = config.signatureText;
results.isSmokeTest = false;
results.totalTrials = config.totalTrials;
results.outputType = 'full bisection noise results';
results.execution = execution;
results.trajectoryRetentionPolicy = sweep.trajectoryRetentionPolicy;
results.retainedTrajectoryRecordIndices = ...
    find(arrayfun(@(trial) ~isempty(trial.X), trials));
results.noiseAxisDefinition = ...
    'per-step per-coordinate Gaussian standard deviation';
results.directlyComparableToPaperFigure6 = false;
end

function save_results(results, sweep)
noiseResults = results;
temporaryFile = [sweep.outputMatFile, '.tmp'];
save(temporaryFile, 'noiseResults', '-v7.3');
movefile(temporaryFile, sweep.outputMatFile, 'f');
writetable(results.summaryTable, sweep.summaryCsvFile);
writetable(results.trialsTable, sweep.trialsCsvFile);
end

function validate_bisection_trial(trial, config)
verifiedIdentifier = config.verifiedSearchMethod;
if ~strcmp(trial.searchMethod, verifiedIdentifier)
    error('Trial searchMethod is not the verified bisection identifier.');
end
if ~isfield(trial, 'initialSearchPerformed') || ~trial.initialSearchPerformed
    error('Trial lacks initial search provenance.');
end
if ~strcmp(trial.initialSearchMethod, verifiedIdentifier)
    error('Trial initial search did not use verified bisection.');
end
if trial.initialSearchGridEvaluations ~= 0
    error('Trial initial search contains grid evaluations.');
end
if isfield(trial, 'numRetargetingEvents') && ...
        trial.numRetargetingEvents ~= numel(trial.retargetingStepIndices)
    error('Trial retargeting event count does not match stored step indices.');
end
if any(strcmp(trial.retargetingSearchMethods, 'hybridGridBisection'))
    error('Trial contains hybridGridBisection provenance.');
end
if ~isempty(trial.retargetingSearchMethods) && ...
        ~all(strcmp(trial.retargetingSearchMethods, verifiedIdentifier))
    error('A trial retargeting event did not use verified bisection.');
end
if any(trial.retargetingGridEvaluations ~= 0)
    error('A trial retargeting search contains hybrid grid evaluations.');
end
if ~isempty(trial.retargetingCrossingOrderTestUsed) && ...
        ~all(trial.retargetingCrossingOrderTestUsed)
    error('A trial retargeting search lacks crossing-order metadata.');
end
if any(contains(trial.retargetingSelectionMethods, 'grid-first'))
    error('A trial used hybrid grid-first candidate selection.');
end
if ~strcmp(trial.mode, 'paperCycleCadenceRetargetMode')
    error('Final noise sweep trial did not use paperCycleCadenceRetargetMode.');
end
if ~strcmp(trial.retargetCadencePolicy, config.retargetCadencePolicy)
    error('Cycle-cadence trial policy mismatch.');
end
if trial.paperReferenceStepCount ~= config.paperReferenceStepCount
    error('Cycle-cadence trial paperReferenceStepCount mismatch.');
end
if trial.retargetStepInterval ~= config.retargetStepInterval
    error('Cycle-cadence trial retargetStepInterval mismatch.');
end
if abs(trial.cadenceRatio - config.cadenceRatio) > 100*eps
    error('Cycle-cadence trial cadenceRatio mismatch.');
end
if trial.numRetargetingEvents > floor(trial.numRk4Steps / ...
        trial.retargetStepInterval)
    error('Cycle-cadence retarget count exceeds scheduled cadence.');
end
stepIndices = trial.retargetingStepIndices;
if ~isempty(stepIndices) && ...
        (~all(mod(stepIndices, trial.retargetStepInterval) == 0) || ...
        (numel(stepIndices) > 1 && ...
        any(diff(stepIndices) ~= trial.retargetStepInterval)))
    error('Cycle-cadence trial did not retarget at scheduled RK4 steps.');
end
end

function trial = empty_trial_template(sourceState, target)
trial = initialize_noisy_trial_record(sourceState, target, NaN, NaN, NaN, ...
    'paperCycleCadenceRetargetMode', verified_bisection_identifier(), ...
    false);
end

function trialsTable = make_trials_table(trials)
numRows = numel(trials);
SigmaNoise = NaN(numRows, 1);
TrialNumber = NaN(numRows, 1);
Seed = NaN(numRows, 1);
Mode = cell(numRows, 1);
SearchMethod = cell(numRows, 1);
Hit = false(numRows, 1);
NumCrossingsToTarget = NaN(numRows, 1);
NumAcceptedCrossings = NaN(numRows, 1);
NumRejectedCrossings = NaN(numRows, 1);
PhysicalTime = NaN(numRows, 1);
NumRk4Steps = NaN(numRows, 1);
SelectedPInitial = NaN(numRows, 1);
SelectedPFinal = NaN(numRows, 1);
SelectedPValues = cell(numRows, 1);
InitialSearchPerformed = false(numRows, 1);
InitialSearchFound = false(numRows, 1);
InitialSearchMethod = cell(numRows, 1);
InitialSearchSelectedP = NaN(numRows, 1);
InitialSearchHorizon = NaN(numRows, 1);
NumSearchEvaluations = NaN(numRows, 1);
NumCompleteBisectionSearches = NaN(numRows, 1);
NumRetargetingEvents = NaN(numRows, 1);
NumRetargetingFailures = NaN(numRows, 1);
RetargetingCrossingIndices = cell(numRows, 1);
RetargetingStepIndices = cell(numRows, 1);
RetargetingTimes = cell(numRows, 1);
RetargetingSearchMethods = cell(numRows, 1);
RetargetingSelectedP = cell(numRows, 1);
RetargetStepInterval = NaN(numRows, 1);
RetargetCadencePolicy = cell(numRows, 1);
CadenceRatio = NaN(numRows, 1);
FinalTargetError = NaN(numRows, 1);
FinalStateX = NaN(numRows, 1);
FinalStateY = NaN(numRows, 1);
FinalStateZ = NaN(numRows, 1);
IntegrationFailed = false(numRows, 1);
StepBudgetExhausted = false(numRows, 1);
AcceptedCrossingLimitExhausted = false(numRows, 1);
RejectedCrossingLimitExhausted = false(numRows, 1);
FailureReason = cell(numRows, 1);
ComputationalTime = NaN(numRows, 1);

for i = 1:numRows
    SigmaNoise(i) = trials(i).sigmaNoise;
    TrialNumber(i) = trials(i).trialNumber;
    Seed(i) = trials(i).seed;
    Mode{i} = trials(i).mode;
    SearchMethod{i} = trials(i).searchMethod;
    Hit(i) = trials(i).hit;
    NumCrossingsToTarget(i) = get_numeric_field(trials(i), ...
        'numCrossingsToTarget', NaN);
    NumAcceptedCrossings(i) = get_numeric_field(trials(i), ...
        'numAcceptedCrossings', NaN);
    NumRejectedCrossings(i) = get_numeric_field(trials(i), ...
        'numRejectedCrossings', NaN);
    PhysicalTime(i) = trials(i).physicalTime;
    NumRk4Steps(i) = trials(i).numRk4Steps;
    SelectedPInitial(i) = trials(i).selectedPInitial;
    SelectedPFinal(i) = trials(i).selectedPFinal;
    SelectedPValues{i} = numeric_values_to_text(trials(i).selectedPValues);
    InitialSearchPerformed(i) = logical_field(trials(i), ...
        'initialSearchPerformed');
    InitialSearchFound(i) = logical_field(trials(i), 'initialSearchFound');
    InitialSearchMethod{i} = trials(i).initialSearchMethod;
    InitialSearchSelectedP(i) = get_numeric_field(trials(i), ...
        'initialSearchSelectedP', NaN);
    InitialSearchHorizon(i) = get_numeric_field(trials(i), ...
        'initialSearchHorizon', NaN);
    NumSearchEvaluations(i) = get_numeric_field(trials(i), ...
        'numSearchEvaluations', NaN);
    NumCompleteBisectionSearches(i) = get_numeric_field(trials(i), ...
        'numCompleteBisectionSearches', NaN);
    NumRetargetingEvents(i) = trials(i).numRetargetingEvents;
    NumRetargetingFailures(i) = trials(i).numRetargetingFailures;
    RetargetingCrossingIndices{i} = ...
        numeric_values_to_text(trials(i).retargetingCrossingIndices);
    RetargetingStepIndices{i} = ...
        numeric_values_to_text(trials(i).retargetingStepIndices);
    RetargetingTimes{i} = numeric_values_to_text(trials(i).retargetingTimes);
    RetargetingSearchMethods{i} = ...
        strjoin(trials(i).retargetingSearchMethods, ';');
    RetargetingSelectedP{i} = ...
        numeric_values_to_text(trials(i).retargetingSelectedP);
    RetargetStepInterval(i) = trials(i).retargetStepInterval;
    RetargetCadencePolicy{i} = trials(i).retargetCadencePolicy;
    CadenceRatio(i) = trials(i).cadenceRatio;
    FinalTargetError(i) = trials(i).finalDistance;
    FinalStateX(i) = trials(i).finalState(1);
    FinalStateY(i) = trials(i).finalState(2);
    FinalStateZ(i) = trials(i).finalState(3);
    IntegrationFailed(i) = trials(i).integrationFailed;
    StepBudgetExhausted(i) = logical_field(trials(i), ...
        'stepBudgetExhausted');
    AcceptedCrossingLimitExhausted(i) = logical_field(trials(i), ...
        'acceptedCrossingLimitExhausted');
    RejectedCrossingLimitExhausted(i) = logical_field(trials(i), ...
        'rejectedCrossingLimitExhausted');
    FailureReason{i} = trials(i).failureReason;
    ComputationalTime(i) = trials(i).computationalTime;
end

trialsTable = table(SigmaNoise, TrialNumber, Seed, Mode, SearchMethod, ...
    Hit, NumCrossingsToTarget, NumAcceptedCrossings, ...
    NumRejectedCrossings, PhysicalTime, NumRk4Steps, ...
    SelectedPInitial, SelectedPFinal, SelectedPValues, ...
    InitialSearchPerformed, InitialSearchFound, InitialSearchMethod, ...
    InitialSearchSelectedP, InitialSearchHorizon, ...
    NumSearchEvaluations, NumCompleteBisectionSearches, ...
    NumRetargetingEvents, NumRetargetingFailures, ...
    RetargetingCrossingIndices, RetargetingStepIndices, ...
    RetargetingTimes, RetargetingSearchMethods, RetargetingSelectedP, ...
    RetargetStepInterval, RetargetCadencePolicy, CadenceRatio, ...
    FinalTargetError, FinalStateX, FinalStateY, FinalStateZ, ...
    IntegrationFailed, StepBudgetExhausted, ...
    AcceptedCrossingLimitExhausted, RejectedCrossingLimitExhausted, ...
    FailureReason, ComputationalTime);
end

function summaryTable = make_summary_table(trials, sigmaNoiseValues, modeList)
numRows = numel(sigmaNoiseValues) * numel(modeList);
Mode = cell(numRows, 1);
SigmaNoise = NaN(numRows, 1);
Successes = NaN(numRows, 1);
TotalTrials = NaN(numRows, 1);
SuccessFraction = NaN(numRows, 1);
WilsonLower95 = NaN(numRows, 1);
WilsonUpper95 = NaN(numRows, 1);
MeanCrossingsToTargetAmongSuccesses = NaN(numRows, 1);
MedianCrossingsToTargetAmongSuccesses = NaN(numRows, 1);
MeanAcceptedCrossingsAllTrials = NaN(numRows, 1);
MeanFinalTargetError = NaN(numRows, 1);
MedianFinalTargetError = NaN(numRows, 1);
StdFinalTargetError = NaN(numRows, 1);
MeanRk4Steps = NaN(numRows, 1);
MeanPhysicalTime = NaN(numRows, 1);
MeanRetargetingEvents = NaN(numRows, 1);
MedianRetargetingEvents = NaN(numRows, 1);
MeanCompleteBisectionSearches = NaN(numRows, 1);
MeanRetargetingFailures = NaN(numRows, 1);
NumIntegrationFailures = NaN(numRows, 1);
NumStepBudgetExhaustions = NaN(numRows, 1);
NumAcceptedCrossingLimitExhaustions = NaN(numRows, 1);
NumRejectedCrossingLimitExhaustions = NaN(numRows, 1);
MeanComputationalTime = NaN(numRows, 1);

row = 0;
trialModes = {trials.mode};
trialSigmas = [trials.sigmaNoise];
for modeIndex = 1:numel(modeList)
    thisMode = modeList{modeIndex};
    for sigmaIndex = 1:numel(sigmaNoiseValues)
        row = row + 1;
        Mode{row} = thisMode;
        SigmaNoise(row) = sigmaNoiseValues(sigmaIndex);
        indices = find(strcmp(trialModes, thisMode) & ...
            abs(trialSigmas - sigmaNoiseValues(sigmaIndex)) < 10*eps);
        if isempty(indices)
            continue;
        end
        hits = [trials(indices).hit];
        Successes(row) = sum(hits);
        TotalTrials(row) = numel(indices);
        SuccessFraction(row) = mean(hits);
        [WilsonLower95(row), WilsonUpper95(row)] = ...
            wilson_interval(Successes(row), TotalTrials(row));
        successIndices = indices(hits);
        crossingsToTarget = [trials(successIndices).numCrossingsToTarget];
        MeanCrossingsToTargetAmongSuccesses(row) = ...
            mean_finite(crossingsToTarget);
        MedianCrossingsToTargetAmongSuccesses(row) = ...
            median_finite(crossingsToTarget);
        MeanAcceptedCrossingsAllTrials(row) = ...
            mean_finite([trials(indices).numAcceptedCrossings]);
        finalErrors = [trials(indices).finalDistance];
        MeanFinalTargetError(row) = mean_finite(finalErrors);
        MedianFinalTargetError(row) = median_finite(finalErrors);
        StdFinalTargetError(row) = std_finite(finalErrors);
        MeanRk4Steps(row) = mean_finite([trials(indices).numRk4Steps]);
        MeanPhysicalTime(row) = mean_finite([trials(indices).physicalTime]);
        MeanRetargetingEvents(row) = ...
            mean_finite([trials(indices).numRetargetingEvents]);
        MedianRetargetingEvents(row) = ...
            median_finite([trials(indices).numRetargetingEvents]);
        MeanCompleteBisectionSearches(row) = ...
            mean_finite([trials(indices).numCompleteBisectionSearches]);
        MeanRetargetingFailures(row) = ...
            mean_finite([trials(indices).numRetargetingFailures]);
        NumIntegrationFailures(row) = ...
            sum([trials(indices).integrationFailed]);
        NumStepBudgetExhaustions(row) = ...
            sum([trials(indices).stepBudgetExhausted]);
        NumAcceptedCrossingLimitExhaustions(row) = ...
            sum([trials(indices).acceptedCrossingLimitExhausted]);
        NumRejectedCrossingLimitExhaustions(row) = ...
            sum([trials(indices).rejectedCrossingLimitExhausted]);
        MeanComputationalTime(row) = ...
            mean_finite([trials(indices).computationalTime]);
    end
end

summaryTable = table(Mode, SigmaNoise, Successes, TotalTrials, ...
    SuccessFraction, WilsonLower95, WilsonUpper95, ...
    MeanCrossingsToTargetAmongSuccesses, ...
    MedianCrossingsToTargetAmongSuccesses, ...
    MeanAcceptedCrossingsAllTrials, MeanFinalTargetError, ...
    MedianFinalTargetError, StdFinalTargetError, MeanRk4Steps, ...
    MeanPhysicalTime, MeanRetargetingEvents, ...
    MedianRetargetingEvents, MeanCompleteBisectionSearches, ...
    MeanRetargetingFailures, NumIntegrationFailures, ...
    NumStepBudgetExhaustions, NumAcceptedCrossingLimitExhaustions, ...
    NumRejectedCrossingLimitExhaustions, MeanComputationalTime);
end

function [lower, upper] = wilson_interval(successes, total)
if total <= 0
    lower = NaN;
    upper = NaN;
    return;
end
z = 1.95996398454005;
p = successes / total;
denominator = 1 + z^2 / total;
centre = p + z^2 / (2 * total);
spread = z * sqrt((p * (1 - p) + z^2 / (4 * total)) / total);
lower = max(0, (centre - spread) / denominator);
upper = min(1, (centre + spread) / denominator);
end

function textValue = numeric_values_to_text(values)
if isempty(values)
    textValue = '';
    return;
end
parts = cell(1, numel(values));
for i = 1:numel(values)
    parts{i} = sprintf('%.12g', values(i));
end
textValue = strjoin(parts, ';');
end

function value = get_numeric_field(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function tf = logical_field(s, fieldName)
tf = isfield(s, fieldName) && logical(s.(fieldName));
end

function value = mean_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = mean(x);
end
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function value = std_finite(x)
x = x(isfinite(x));
if numel(x) <= 1
    value = NaN;
else
    value = std(x);
end
end
