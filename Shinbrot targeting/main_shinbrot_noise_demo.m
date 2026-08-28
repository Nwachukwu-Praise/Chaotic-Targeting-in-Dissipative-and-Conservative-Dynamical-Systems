function workflow = main_shinbrot_noise_demo(noiseAction)
%MAIN_SHINBROT_NOISE_DEMO Paper-cycle noise workflow controller.
%
% noiseAction:
%   'skip'           return without running the noise workflow
%   'validate'       run/load deterministic and zero-noise gates only
%   'loadFull'       load a complete, configuration-matched full result
%   'pilot'          run trial 1 at each final sigma into final checkpoint
%   'resumeFull'     resume the final 110-trial sweep from checkpoint
%   'regenerateFull' explicit full-run alias; also resumes checkpoint
%
% The final noise experiment uses paperCycleCadenceRetargetMode only.  The
% literal paper40StepRetargetMode and sectionRetargetMode remain callable
% for audit/extension purposes but are not final scientific output.

if nargin < 1 || isempty(noiseAction)
    noiseAction = 'loadFull';
end
if isstring(noiseAction)
    noiseAction = char(noiseAction);
end
validActions = {'skip', 'validate', 'loadFull', ...
    'pilot', 'resumeFull', 'regenerateFull'};
if ~any(strcmp(noiseAction, validActions))
    error('Unknown noiseAction: %s', noiseAction);
end

workflow = struct();
workflow.noiseAction = noiseAction;
workflow.noiseResults = [];
workflow.validation = struct();
workflow.readyForFullSweep = false;
workflow.readyText = 'NOT READY FOR FULL SWEEP';

if strcmp(noiseAction, 'skip')
    fprintf('Noise workflow skipped by noiseAction.%s', newline);
    return;
end

fprintf('========== Paper-Cycle Bisection Noisy Shinbrot Lorenz Targeting ==========%s', newline);
verifiedBisectionIdentifier = verified_bisection_identifier();

[params, target, deterministicControl, noise] = base_problem_settings();
noiseControl = deterministicControl;
noiseControl.searchMethod = verifiedBisectionIdentifier;
assert(strcmp(noiseControl.searchMethod, verifiedBisectionIdentifier));

fprintf('Loading or generating deterministic source state X_150...%s', newline);
[mapData, sourceState, sourceIndex] = nominal_source_state(params, target);
fprintf('Source X_s = %.12g, target X_t = %.12g%s', ...
    sourceState(1), target.x, newline);

noise = set_measured_noise_cadence_and_budget(noise, mapData, 100);
cadence = shinbrot_noise_cadence_diagnostics(mapData, noise, 100);
fprintf(['Mean return time Tbar = %.6g, steps per mean return = %.1f, ', ...
    'cycle-matched interval = %d steps, cadence ratio = %.6g%s'], ...
    cadence.Tbar, cadence.stepsPerMeanReturn, ...
    noise.retargetStepInterval, cadence.cadenceRatio, newline);
fprintf(['Literal 40-step time at dt = %.4g is %.6g, cadence ratio %.6g; ', ...
    'that is about %.1fx too frequent.%s'], ...
    noise.dt, cadence.literalFortyStepPhysicalTime, ...
    cadence.literalFortyStepCadenceRatio, ...
    1 / cadence.literalFortyStepCadenceRatio, newline);
fprintf(['Computed full-trial RK4-step budget = %d ', ...
    '(safety factor %.3g).%s'], ...
    noise.maxRk4Steps, noise.budgetSafetyFactor, newline);

sweep = bisection_sweep_settings(verifiedBisectionIdentifier, sourceIndex);
assert(strcmp(sweep.searchMethod, verifiedBisectionIdentifier));
parallelInfo = select_noise_parallel_workers(sweep.parallel.requestedWorkers);
sweep.parallel = parallelInfo;

fprintf('Finding initial deterministic bisection perturbation...%s', newline);
initialSearch = search_parameter_to_target( ...
    sourceState, target, params, noiseControl, verifiedBisectionIdentifier);
assert_verified_bisection_search(initialSearch, verifiedBisectionIdentifier, ...
    'main noise initial search', target);
if ~initialSearch.found
    error('The initial deterministic bisection search did not find a hit.');
end
fprintf('Initial pControl = %.12g, horizon = %d%s', ...
    initialSearch.selectedP, initialSearch.horizon, newline);

currentConfig = build_bisection_noise_config( ...
    sourceState, target, params, noiseControl, noise, sweep);

[validation, validationLoaded, validationDifferences] = ...
    load_or_run_validation_cache(sweep.validationCacheFile, ...
    currentConfig, sourceState, target, params, noise, ...
    sweep, initialSearch, sourceIndex);
if validationLoaded
    fprintf('Loaded valid cached zero-noise validation from %s.%s', ...
        sweep.validationCacheFile, newline);
elseif ~isempty(validationDifferences)
    fprintf('Rebuilt validation cache because: %s%s', ...
        strjoin(cellstr(validationDifferences), '; '), newline);
end

zeroNoiseGate = validation.zeroNoiseGate;
runtimeEstimate = estimate_noise_runtime( ...
    zeroNoiseGate, currentConfig.totalTrials, parallelInfo);

assert(noise.maxControlledCrossings == 60);
assert(noise.maxRk4Steps >= noise.requiredRk4Steps);
assert(isequal(sweep.sigmaNoiseValues, final_noise_sigma_values()));
assert(sweep.numTrials == 10);
assert(sweep.totalTrials == 110);
assert(~sweep.isSmokeTest);
assert(zeroNoiseGate.allPassed);

workflow.validation.params = params;
workflow.validation.target = target;
workflow.validation.noise = noise;
workflow.validation.sweep = sweep;
workflow.validation.sourceState = sourceState;
workflow.validation.sourceIndex = sourceIndex;
workflow.validation.cadence = cadence;
workflow.validation.initialSearch = initialSearch;
workflow.validation.openLoopComparison = validation.openLoopComparison;
workflow.validation.openLoopConvergence = validation.openLoopConvergence;
workflow.validation.zeroNoiseGate = zeroNoiseGate;
workflow.validation.currentConfig = currentConfig;
workflow.validation.parallel = parallelInfo;
workflow.validation.runtimeEstimate = runtimeEstimate;
workflow.validation.validationCacheLoaded = validationLoaded;
workflow.readyForFullSweep = true;
workflow.readyText = 'READY FOR FULL 110-TRIAL, 60-CROSSING NOISE SWEEP';

switch noiseAction
    case 'validate'
        print_ready_message(workflow, sweep, noise);
        return;

    case 'loadFull'
        [noiseResults, savedNoiseLoaded, differences] = ...
            load_valid_bisection_noise_results( ...
            sweep.outputMatFile, currentConfig);
        if savedNoiseLoaded
            workflow.noiseResults = noiseResults;
            fprintf('Loaded valid saved full noise results from %s.%s', ...
                sweep.outputMatFile, newline);
            report_full_results(noiseResults, target, params, noise);
        else
            warning(['No valid saved full paper-cycle noise result is ', ...
                'available. The deterministic investigations will continue. ', ...
                'Validation differences: %s'], ...
                strjoin(cellstr(differences), '; '));
            print_ready_message(workflow, sweep, noise);
        end
        return;

    case 'pilot'
        sweep.runMode = 'pilot';
        fprintf('%sExperiment type: FINAL-CHECKPOINT PILOT%s', newline, newline);
        fprintf(['Pilot records: trial 1 at each of the %d configured ', ...
            'sigma_noise levels.%s'], numel(sweep.sigmaNoiseValues), newline);
        fprintf('Checkpoint file: %s%s', sweep.checkpointMatFile, newline);
        noiseResults = run_noise_sweep( ...
            sourceState, target, params, noiseControl, noise, ...
            sweep, initialSearch);
        noiseResults.cadenceDiagnostics = cadence;
        noiseResults.validation = validation;
        workflow.noiseResults = noiseResults;
        workflow.pilotReport = report_pilot_results( ...
            noiseResults, parallelInfo, noise);
        return;

    case {'resumeFull', 'regenerateFull'}
        sweep.runMode = noiseAction;
        print_full_sweep_header(sweep, noise, parallelInfo, currentConfig);
        assert(zeroNoiseGate.allPassed);
        noiseResults = run_noise_sweep( ...
            sourceState, target, params, noiseControl, noise, ...
            sweep, initialSearch);
        noiseResults.cadenceDiagnostics = cadence;
        noiseResults.validation = validation;
        workflow.noiseResults = noiseResults;
        if noiseResults.completed
            report_full_results(noiseResults, target, params, noise);
        else
            fprintf('Full sweep remains incomplete: %d/%d records complete.%s', ...
                noiseResults.completedTrialCount, noiseResults.totalTrials, newline);
        end
end
end

function [params, target, control, noise] = base_problem_settings()
params.sigma = 10;
params.rho = 28;
params.beta = 8/3;
params.zSection = 26.921;
params.xSectionMin = 8.0;
params.crossingDirection = +1;
params.odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
params.eventDisableTime = 1e-6;
params.maxStepOffSection = 1e-5;
params.sectionTolerance = 1e-9;
params.maxTimeToNextValidSection = 100;
params.maxRejectedCrossings = 100;

target.x = 13.729;
target.state = [13.729; 19.585; params.zSection];
target.tolerance = 0.008;

control.deltaP = 0.1;
control.maxSearchCrossings = 12;
control.numParameterSamples = 81;
control.bisectionIterations = 24;
control.maxDiscontinuityIsolationDepth = 24;
control.maxDiscontinuityIntervals = 256;

noise.dt = 0.001;
noise.maxTimeToNextValidSection = 80;
noise.maxControlledCrossings = 60;
noise.maxRejectedCrossings = 250;
noise.maxStateNorm = 1e5;
noise.axisDefinition = ...
    'per-step per-coordinate Gaussian standard deviation';
noise.directlyComparableToPaperFigure6 = false;
noise.budgetSafetyFactor = 2.5;
noise.paperReferenceStepCount = 40;
noise.retargetCadencePolicy = 'meanReturnCycleMatched';
end

function [mapData, sourceState, sourceIndex] = nominal_source_state(params, target)
cacheFile = fullfile(fileparts(mfilename('fullpath')), 'shinbrot_nominal_map.mat');
if isfile(cacheFile)
    cached = load(cacheFile, 'mapData');
    mapData = cached.mapData;
else
    X0 = [1; 1; 1];
    numMapCrossings = 1200;
    mapData = generate_return_map(X0, numMapCrossings, params, 0);
    save(cacheFile, 'mapData', '-v7.3');
end
sourceIndex = 150;
if size(mapData.states, 2) < sourceIndex
    error('Only %d valid section crossings were generated.', ...
        size(mapData.states, 2));
end
sourceState = mapData.states(:, sourceIndex);
if abs(sourceState(1) - target.x) <= target.tolerance
    error('Configured source X_150 lies inside the target interval.');
end
end

function noise = set_measured_noise_cadence_and_budget(noise, mapData, burnIn)
eventTimes = mapData.eventTimes(:);
if numel(eventTimes) < burnIn + 2
    error('Not enough section returns to estimate the RK4 step budget.');
end
intervals = diff(eventTimes((burnIn + 1):end));
intervals = intervals(isfinite(intervals) & intervals > 0);
if isempty(intervals)
    error('No positive return intervals were available for budget estimation.');
end
noise.budgetBurnIn = burnIn;
noise.meanReturnTimeForBudget = mean(intervals);
noise.estimatedStepsPerCrossing = noise.meanReturnTimeForBudget / noise.dt;
noise.retargetStepInterval = max(1, round( ...
    noise.meanReturnTimeForBudget / noise.dt));
noise.retargetPhysicalTime = noise.retargetStepInterval * noise.dt;
noise.cadenceRatio = ...
    noise.retargetPhysicalTime / noise.meanReturnTimeForBudget;
noise.requiredRk4Steps = ceil(noise.budgetSafetyFactor * ...
    noise.maxControlledCrossings * noise.estimatedStepsPerCrossing);
noise.maxRk4Steps = noise.requiredRk4Steps;
end

function sweep = bisection_sweep_settings(verifiedBisectionIdentifier, sourceIndex)
sweep.sigmaNoiseValues = final_noise_sigma_values();
sweep.numTrials = 10;
sweep.baseSeed = 91000;
sweep.modes = {'paperCycleCadenceRetargetMode'};
sweep.totalTrials = numel(sweep.sigmaNoiseValues) * ...
    sweep.numTrials * numel(sweep.modes);
sweep.searchMethod = verifiedBisectionIdentifier;
sweep.sourceIndex = sourceIndex;
sweep.outputMatFile = 'noise_results_paper_cycle_retarget_11sigma.mat';
sweep.checkpointMatFile = ...
    'noise_results_paper_cycle_retarget_11sigma_checkpoint.mat';
sweep.summaryCsvFile = 'noise_summary_paper_cycle_retarget_11sigma.csv';
sweep.trialsCsvFile = 'noise_trials_paper_cycle_retarget_11sigma.csv';
sweep.progressLogFile = 'noise_progress_paper_cycle_retarget_11sigma.log';
sweep.validationCacheFile = 'noise_validation_paper_cycle_retarget_11sigma.mat';
sweep.resumeFromCheckpoint = true;
sweep.isSmokeTest = false;
sweep.trajectoryRetentionPolicy = ...
    ['retain full RK4 trajectories only for trial 1 at the lowest, ', ...
    'middle, and highest configured sigma_noise levels'];
sweep.parallel.requestedWorkers = 'auto';
end

function values = final_noise_sigma_values()
values = [0, 0.01, 0.03, 0.05, 0.075, 0.1, ...
    0.15, 0.2, 0.3, 0.5, 1.0];
end

function [validation, loaded, differences] = load_or_run_validation_cache( ...
    cacheFile, currentConfig, sourceState, target, params, ...
    noise, sweep, initialSearch, sourceIndex)
loaded = false;
differences = strings(0, 1);
if isfile(cacheFile)
    saved = load(cacheFile);
    if isfield(saved, 'validation') && isfield(saved.validation, 'configuration')
        [valid, differences] = validate_bisection_noise_config( ...
            saved.validation.configuration, currentConfig);
        if valid && isfield(saved.validation, 'zeroNoiseGate') && ...
                saved.validation.zeroNoiseGate.allPassed
            validation = saved.validation;
            loaded = true;
            return;
        end
    else
        differences(end + 1, 1) = "validation cache missing validation struct";
    end
else
    differences(end + 1, 1) = "validation cache file not found";
end

fprintf('%s========== Open-Loop Zero-Noise Solver Comparison ==========%s', ...
    newline, newline);
openLoopComparison = fixed_step_open_loop_replay( ...
    sourceState, target, params, noise, ...
    initialSearch.selectedP, initialSearch.horizon);
disp(openLoopComparison.table);
if ~openLoopComparison.adaptiveHit
    error(['Open-loop adaptive replay did not reach the target with ', ...
        'the selected constant pControl.']);
end

openLoopConvergence = [];
if ~openLoopComparison.fixedStepHit
    warning(['Configured-dt fixed-step open-loop replay did not hit ', ...
        'the target with the same constant pControl. Running a reduced ', ...
        'timestep convergence diagnostic before any full sweep.']);
    openLoopConvergence = fixed_step_open_loop_convergence( ...
        sourceState, target, params, noise, ...
        initialSearch.selectedP, initialSearch.horizon, ...
        noise.dt ./ [1, 2, 4]);
    disp(openLoopConvergence(:, {'horizon', 'selectedP', ...
        'fixedStepFinalX', 'fixedStepTargetError', ...
        'absoluteSolverDiscrepancy', ...
        'discrepancyFractionOfTolerance', 'rk4StepCount', ...
        'acceptedCrossingCount', 'rejectedCrossingCount', ...
        'fixedStepHit'}));
    if ~any(openLoopConvergence.fixedStepHit)
        error(['Open-loop fixed-step convergence diagnostic did not ', ...
            'recover the deterministic target hit.']);
    end
end

zeroGateOptions.sourceIndex = sourceIndex;
zeroGateOptions.noise = noise;
zeroGateOptions.verbose = true;
zeroNoiseGate = test_zero_noise_gate(sweep.modes, zeroGateOptions);

validation.configuration = currentConfig;
validation.sourceIndex = sourceIndex;
validation.initialSearch = initialSearch;
validation.openLoopComparison = openLoopComparison;
validation.openLoopConvergence = openLoopConvergence;
validation.zeroNoiseGate = zeroNoiseGate;
validation.savedAt = datetime('now');
temporaryFile = [cacheFile, '.tmp'];
save(temporaryFile, 'validation', '-v7.3');
movefile(temporaryFile, cacheFile, 'f');
end

function estimate = estimate_noise_runtime( ...
    zeroNoiseGate, totalTrials, parallelInfo)
seconds = zeroNoiseGate.table.seconds;
seconds = seconds(isfinite(seconds) & seconds > 0);
if isempty(seconds)
    measuredTrialSeconds = NaN;
else
    measuredTrialSeconds = median(seconds);
end
estimate.measuredValidationTrialSeconds = measuredTrialSeconds;
estimate.totalTrials = totalTrials;
estimate.serialSeconds = measuredTrialSeconds * totalTrials;
if strcmp(parallelInfo.executionMode, 'parallel')
    estimate.parallelEfficiencyAssumption = ...
        parallelInfo.parallelEfficiencyAssumption;
    estimate.parallelSeconds = estimate.serialSeconds / ...
        (parallelInfo.selectedWorkers * ...
        estimate.parallelEfficiencyAssumption);
else
    estimate.parallelEfficiencyAssumption = 1;
    estimate.parallelSeconds = estimate.serialSeconds;
end
estimate.note = ['Estimate is based on the sigma_noise = 0 validation ', ...
    'trial time; noisy failures may take longer.'];
end

function report = report_pilot_results(noiseResults, parallelInfo, noise)
completed = noiseResults.completedRecordIndices(:);
trials = noiseResults.trials(completed);
RecordIndex = completed;
SigmaNoise = [trials.sigmaNoise].';
TrialNumber = [trials.trialNumber].';
Seed = [trials.seed].';
Hit = [trials.hit].';
CrossingsToTarget = [trials.numCrossingsToTarget].';
AcceptedCrossings = [trials.numAcceptedCrossings].';
Rk4Steps = [trials.numRk4Steps].';
RetargetingEvents = [trials.numRetargetingEvents].';
CompleteBisectionSearches = [trials.numCompleteBisectionSearches].';
Seconds = [trials.computationalTime].';
report.table = table(RecordIndex, SigmaNoise, TrialNumber, Seed, Hit, ...
    CrossingsToTarget, AcceptedCrossings, Rk4Steps, RetargetingEvents, ...
    CompleteBisectionSearches, Seconds);

completedCount = numel(completed);
remainingCount = noiseResults.totalTrials - completedCount;
meanSeconds = mean(Seconds(isfinite(Seconds) & Seconds > 0));
medianSeconds = median(Seconds(isfinite(Seconds) & Seconds > 0));
maxSeconds = max(Seconds(isfinite(Seconds) & Seconds > 0));
if isempty(meanSeconds)
    meanSeconds = NaN;
    medianSeconds = NaN;
    maxSeconds = NaN;
end
efficiency = parallelInfo.parallelEfficiencyAssumption;
workers = max(1, parallelInfo.selectedWorkers);

report.completedCount = completedCount;
report.remainingCount = remainingCount;
report.meanPilotSeconds = meanSeconds;
report.medianPilotSeconds = medianSeconds;
report.slowestPilotSeconds = maxSeconds;
report.projectedRemainingSerialSeconds = meanSeconds * remainingCount;
report.projectedRemainingParallelSeconds = ...
    report.projectedRemainingSerialSeconds / (workers * efficiency);
report.conservativeRemainingSerialSeconds = maxSeconds * remainingCount;
report.conservativeRemainingParallelSeconds = ...
    report.conservativeRemainingSerialSeconds / (workers * efficiency);
report.cadenceRatio = noise.cadenceRatio;

fprintf('%s========== Final Pilot Trials (%d Sigma Levels) ==========%s', ...
    newline, numel(unique(SigmaNoise)), newline);
disp(report.table);
fprintf('Completed final checkpoint records: %d/%d%s', ...
    completedCount, noiseResults.totalTrials, newline);
fprintf('Remaining records for resumeFull: %d%s', remainingCount, newline);
fprintf('Measured cadence ratio: %.6g%s', report.cadenceRatio, newline);
fprintf('Mean pilot duration: %s%s', ...
    duration_text(meanSeconds), newline);
fprintf('Median pilot duration: %s%s', ...
    duration_text(medianSeconds), newline);
fprintf('Slowest pilot duration: %s%s', ...
    duration_text(maxSeconds), newline);
fprintf('Projected remaining serial duration: %s%s', ...
    duration_text(report.projectedRemainingSerialSeconds), newline);
fprintf('Projected remaining parallel duration: %s%s', ...
    duration_text(report.projectedRemainingParallelSeconds), newline);
fprintf('Conservative remaining serial duration: %s%s', ...
    duration_text(report.conservativeRemainingSerialSeconds), newline);
fprintf('Conservative remaining parallel duration: %s%s', ...
    duration_text(report.conservativeRemainingParallelSeconds), newline);
if report.conservativeRemainingParallelSeconds > 3600
    fprintf(['Projection gate: conservative parallel estimate exceeds ', ...
        '60 minutes, so the remaining %d trials were not started.%s'], ...
        report.remainingCount, newline);
else
    fprintf(['Projection gate: conservative parallel estimate is within ', ...
        '60 minutes, but the remaining %d trials were still not started ', ...
        'without explicit instruction.%s'], report.remainingCount, newline);
end
fprintf(['Resume command after approval: ', ...
    'main_shinbrot_noise_demo(''resumeFull'')%s'], newline);
end

function print_ready_message(workflow, sweep, noise)
parallelInfo = workflow.validation.parallel;
estimate = workflow.validation.runtimeEstimate;
fprintf('%s%s%s%s', newline, workflow.readyText, newline, newline);
fprintf('Noise levels:%s', newline);
disp(sweep.sigmaNoiseValues);
fprintf('Trials per level: %d%s', sweep.numTrials, newline);
fprintf('Total trials: %d%s', sweep.totalTrials, newline);
fprintf('Retargeting mode: %s%s', sweep.modes{1}, newline);
fprintf('Cycle-matched retarget interval: %d RK4 steps%s', ...
    noise.retargetStepInterval, newline);
fprintf('Cadence ratio: %.6g%s', noise.cadenceRatio, newline);
fprintf('Maximum accepted crossings per trial: %d%s', ...
    noise.maxControlledCrossings, newline);
fprintf('Global RK4-step budget per trial: %d%s', ...
    noise.maxRk4Steps, newline);
fprintf('Available process workers: %d%s', ...
    parallelInfo.availableWorkers, newline);
fprintf('Requested workers: %s%s', ...
    worker_request_text(parallelInfo.requestedWorkers), newline);
fprintf('Selected CPU workers: %d%s', ...
    parallelInfo.selectedWorkers, newline);
fprintf('Execution mode: %s%s', parallelInfo.executionMode, newline);
fprintf('Measured validation-trial time: %.3g s%s', ...
    estimate.measuredValidationTrialSeconds, newline);
fprintf('Estimated serial duration: %s%s', ...
    duration_text(estimate.serialSeconds), newline);
fprintf('Estimated parallel duration: %s%s', ...
    duration_text(estimate.parallelSeconds), newline);
fprintf('Parallel efficiency assumption: %.2f%s', ...
    estimate.parallelEfficiencyAssumption, newline);
fprintf('Checkpoint filename: %s%s%s', ...
    sweep.checkpointMatFile, newline, newline);
fprintf(['Set noiseAction = ''pilot'' to run trial 1 at each configured ', ...
    'sigma_noise level.%s'], newline);
fprintf('Set noiseAction = ''resumeFull'' to continue the remaining records.%s', newline);
end

function print_full_sweep_header(sweep, noise, parallelInfo, currentConfig)
fprintf('%sExperiment type: FULL NOISE SWEEP%s', newline, newline);
fprintf('Maximum accepted crossings per trial: %d%s', ...
    noise.maxControlledCrossings, newline);
fprintf('Global RK4-step budget per trial: %d%s', ...
    noise.maxRk4Steps, newline);
fprintf('Noise levels: %d%s', numel(sweep.sigmaNoiseValues), newline);
fprintf('Trials per level: %d%s', sweep.numTrials, newline);
fprintf('Total trials: %d%s', currentConfig.totalTrials, newline);
fprintf('Retargeting mode: %s%s', sweep.modes{1}, newline);
fprintf('Retarget step interval: %d%s', noise.retargetStepInterval, newline);
fprintf('Cadence ratio: %.6g%s', noise.cadenceRatio, newline);
fprintf('Available process workers: %d%s', ...
    parallelInfo.availableWorkers, newline);
fprintf('Requested workers: %s%s', ...
    worker_request_text(parallelInfo.requestedWorkers), newline);
fprintf('Selected CPU workers: %d%s', ...
    parallelInfo.selectedWorkers, newline);
fprintf('Execution mode: %s%s', parallelInfo.executionMode, newline);
fprintf('Checkpoint file: %s%s', sweep.checkpointMatFile, newline);
end

function textValue = duration_text(seconds)
if ~isfinite(seconds)
    textValue = 'unavailable';
elseif seconds < 90
    textValue = sprintf('%.1f s', seconds);
elseif seconds < 7200
    textValue = sprintf('%.1f min', seconds / 60);
else
    textValue = sprintf('%.2f h', seconds / 3600);
end
end

function textValue = worker_request_text(requestedWorkers)
if isnumeric(requestedWorkers)
    textValue = sprintf('%d', requestedWorkers);
else
    textValue = char(string(requestedWorkers));
end
end

function report_full_results(noiseResults, target, params, noise)
fprintf('%s========== Full Paper-Cycle Gaussian-Noise Summary ==========%s', ...
    newline, newline);
disp(noiseResults.summaryTable);

fprintf('%s========== Full Paper-Cycle Gaussian-Noise Trials ==========%s', ...
    newline, newline);
disp(noiseResults.trialsTable);

plot_noise_results(noiseResults, target, params, noise);

fprintf('Verified bisection identifier: %s%s', ...
    noiseResults.verifiedBisectionIdentifier, newline);
fprintf('Noise retargeting mode: %s%s', noiseResults.modes{1}, newline);
fprintf('Total configured trials: %d%s', noiseResults.totalTrials, newline);
fprintf('Result completed: %s%s', logical_text_noise(noiseResults.completed), newline);
fprintf('Directly comparable to Shinbrot Fig. 6 axis: false%s', newline);
end

function textValue = logical_text_noise(value)
if value
    textValue = 'true';
else
    textValue = 'false';
end
end
