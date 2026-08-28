function results = test_zero_noise_gate(modes, options)
%TEST_ZERO_NOISE_GATE Acceptance gate for the noisy trial at sigma = 0.
%
%   results = TEST_ZERO_NOISE_GATE(modes, options)
%
%   At sigma_noise = 0 the noisy trial reduces to a deterministic fixed-step
%   integration of the same problem the deterministic search already solved.
%   It must therefore reach the target.  If it does not, nothing measured at
%   sigma > 0 can be attributed to noise, because the procedure has already
%   failed without any.
%
%   This gate exists to make a sigma = 0 failure impossible to mistake for
%   a noise-robustness result.  A failure here is a numerical-workflow
%   problem, not evidence about stochastic robustness.
%
%   The gate is deliberately cheap: one trial per mode.  Run it before any
%   sweep.
%
%   WHAT IS ASSERTED
%     * the deterministic search from the source succeeds;
%     * the sigma = 0 trial reaches the target neighbourhood;
%     * it does so within the declared crossing limit;
%     * at least one accepted crossing is completed, so the trial is not
%       terminated by a step budget before the dynamics are exercised;
%     * no integration failure occurs;
%     * no retargeting search fails to return a candidate.
%
%   WHAT IS REPORTED BUT NOT ASSERTED
%     * the number of accepted crossings, which need not equal the
%       deterministic horizon: a mode that retargets before the horizon is
%       reached replaces the control mid-plan and may arrive by a different
%       route;
%     * the discrepancy between the deterministic adaptive-solver final
%       error and the fixed-step sigma = 0 final error from the same source
%       and the same perturbation.  Subsection 3.1.9 promises this
%       measurement so that the discretisation change can be distinguished
%       from the effect of the noise itself.  It is a measurement, not a
%       pass criterion.
%
%   modes   : cellstr or string array, default both implemented modes.
%   options : optional struct
%             .sourceIndex   default 150
%             .seed          default 20260815
%             .noise         a full noise struct, overriding the defaults
%             .verbose       default true

if nargin < 1 || isempty(modes)
    modes = {'paperCycleCadenceRetargetMode'};
end
modes = cellstr(modes);
if nargin < 2
    options = struct();
end
verbose = local_opt(options, 'verbose', true);
seed = local_opt(options, 'seed', 20260815);
sourceIndex = local_opt(options, 'sourceIndex', 150);

params = local_params();
target = local_target(params);
control = local_control();
noise = local_opt(options, 'noise', local_noise());
if any(strcmp(modes, 'paper40StepRetargetMode')) && ...
        ~isfield(options, 'noise')
    noise.retargetCadencePolicy = 'literalFortyStepHistorical';
    noise.retargetStepInterval = 40;
    noise.retargetPhysicalTime = 40 * noise.dt;
    noise.cadenceRatio = NaN;
end

verifiedIdentifier = verified_bisection_identifier();
noiseControl = control;
noiseControl.searchMethod = verifiedIdentifier;

mapData = local_nominal_map(params);
sourceState = mapData.states(:, sourceIndex);
while abs(sourceState(1) - target.x) <= target.tolerance
    sourceIndex = sourceIndex + 1;
    sourceState = mapData.states(:, sourceIndex);
end
noise = ensure_global_step_budget(noise, mapData);

if verbose
    fprintf('\n=== Zero-noise acceptance gate ===\n');
    fprintf('implementation      : %s\n', verifiedIdentifier);
    fprintf('source index %d, X_s = %.8g, X_t = %.8g, eps_t = %.4g\n', ...
        sourceIndex, sourceState(1), target.x, target.tolerance);
    fprintf('dt = %g, maxControlledCrossings = %d, maxTimeToNextValidSection = %g\n', ...
        noise.dt, noise.maxControlledCrossings, noise.maxTimeToNextValidSection);
    fprintf('mean return time of the nominal map is about %.4g, so one\n', ...
        local_mean_return_time(mapData));
    fprintf('accepted crossing costs roughly %.0f fixed steps.\n\n', ...
        local_mean_return_time(mapData) / noise.dt);
end

% ---- deterministic reference -------------------------------------------
deterministicSearch = search_parameter_to_target( ...
    sourceState, target, params, noiseControl, verifiedIdentifier);
assert(deterministicSearch.found, ...
    'Zero-noise gate: the deterministic search did not find a candidate.');

deterministicReplay = local_adaptive_replay(sourceState, target, params, ...
    deterministicSearch.selectedP, deterministicSearch.horizon);
assert(deterministicReplay.hit, ...
    'Zero-noise gate: the deterministic candidate does not replay into the target.');

if verbose
    fprintf('deterministic reference (adaptive solver)\n');
    fprintf('  p                 : %.12g\n', deterministicSearch.selectedP);
    fprintf('  horizon           : %d\n', deterministicSearch.horizon);
    fprintf('  replay error      : %.6g\n\n', deterministicReplay.targetError);
end

% ---- one sigma = 0 trial per mode --------------------------------------
rows = {};
failures = strings(0, 1);

for m = 1:numel(modes)
    mode = modes{m};
    trial = run_noisy_targeting_trial(sourceState, target, params, ...
        noiseControl, noise, 0, 1, seed, mode, deterministicSearch, false);

    retargetsFired = trial.numRetargetingEvents;
    controlUnchanged = isequaln(trial.selectedPFinal, trial.selectedPInitial);

    problems = strings(0, 1);
    if ~trial.hit
        problems(end + 1, 1) = "did not reach the target"; %#ok<AGROW>
    end
    if trial.integrationFailed
        problems(end + 1, 1) = "integration failed: " + string(trial.failureReason); %#ok<AGROW>
    end
    if ~(trial.numAcceptedCrossings >= 1)
        problems(end + 1, 1) = ...
            "no accepted crossing completed, so the trial ended before the dynamics were exercised"; %#ok<AGROW>
    end
    if trial.numAcceptedCrossings > noise.maxControlledCrossings
        problems(end + 1, 1) = "exceeded the declared crossing limit"; %#ok<AGROW>
    end
    if trial.numRetargetingFailures > 0
        problems(end + 1, 1) = sprintf("%d retargeting searches returned no candidate", ...
            trial.numRetargetingFailures); %#ok<AGROW>
    end

    passed = isempty(problems);
    if ~passed
        failures(end + 1, 1) = string(mode) + ": " + strjoin(problems, "; "); %#ok<AGROW>
    end

    discretisationDifference = abs(trial.finalDistance - deterministicReplay.targetError);

    rows(end + 1, :) = {string(mode), passed, trial.hit, ...
        trial.numCrossingsToTarget, trial.numAcceptedCrossings, ...
        deterministicSearch.horizon, trial.finalDistance, ...
        deterministicReplay.targetError, discretisationDifference, ...
        discretisationDifference / target.tolerance, retargetsFired, ...
        trial.numSearchEvaluations, controlUnchanged, trial.numRk4Steps, ...
        trial.stepBudgetExhausted, trial.integrationFailed, ...
        string(trial.failureReason), trial.computationalTime}; %#ok<AGROW>

    if verbose
        fprintf('%s\n', mode);
        fprintf('  reached target    : %d\n', trial.hit);
        fprintf('  crossings to hit  : %g   (deterministic horizon %d)\n', ...
            trial.numCrossingsToTarget, deterministicSearch.horizon);
        fprintf('  accepted crossings: %d\n', trial.numAcceptedCrossings);
        fprintf('  final error       : %.6g\n', trial.finalDistance);
        fprintf('  retargets fired   : %d\n', retargetsFired);
        fprintf('  control unchanged : %d\n', controlUnchanged);
        fprintf('  RK4 steps         : %d\n', trial.numRk4Steps);
        if ~passed
            fprintf('  FAILED            : %s\n', strjoin(problems, '; '));
        else
            fprintf('  passed\n');
        end
        fprintf('\n');
    end
end

results.table = cell2table(rows, 'VariableNames', {'mode','passed','hit', ...
    'crossingsToTarget','acceptedCrossings','deterministicHorizon', ...
    'fixedStepFinalError', ...
    'adaptiveFinalError','discretisationDifference', ...
    'discretisationFractionOfTolerance','retargetsFired', ...
    'searchEvaluations','controlUnchanged', 'rk4Steps', ...
    'stepBudgetExhausted','integrationFailed','failureReason','seconds'});
results.deterministicSearch = deterministicSearch;
results.deterministicReplay = deterministicReplay;
results.sourceState = sourceState;
results.sourceIndex = sourceIndex;
results.noise = noise;
results.failures = failures;
results.allPassed = isempty(failures);

if verbose
    disp(results.table);
    fprintf('modes passed : %d of %d\n', sum(results.table.passed), height(results.table));
end

if ~results.allPassed
    error('SchroerOtt:ZeroNoiseGate', ...
        ['The zero-noise gate failed and no sweep should be run:\n  %s\n', ...
        'A failure here is not a noise result. It means the trial does not ', ...
        'reproduce the deterministic solution when no noise is applied.'], ...
        strjoin(failures, sprintf('\n  ')));
end

if verbose
    fprintf(['\nGate passed. The discretisation difference above is the adaptive\n', ...
        'versus fixed-step discrepancy that Subsection 3.1.9 promises to report,\n', ...
        'measured from the same source and the same perturbation.\n']);
    fprintf('==================================\n');
end
end

% -------------------------------------------------------------------------

function replay = local_adaptive_replay(sourceState, target, params, pControl, horizon)
%LOCAL_ADAPTIVE_REPLAY Deterministic reference using the adaptive solver.
replay.hit = false;
replay.targetError = NaN;
replay.numCrossings = 0;
if ~isfinite(pControl) || ~isfinite(horizon) || horizon < 1
    return;
end
currentState = sourceState(:);
for crossing = 1:horizon
    segment = next_valid_section_crossing(currentState, params, pControl);
    if ~segment.success
        return;
    end
    currentState = segment.eventState;
    replay.numCrossings = crossing;
end
replay.targetError = abs(currentState(1) - target.x);
replay.hit = replay.targetError <= target.tolerance;
end

function value = local_mean_return_time(mapData)
if isfield(mapData, 'eventTimes') && numel(mapData.eventTimes) > 101
    intervals = diff(mapData.eventTimes(101:end));
    value = mean(intervals(isfinite(intervals) & intervals > 0));
else
    value = NaN;
end
end

function value = local_opt(options, name, defaultValue)
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

function params = local_params()
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
end

function target = local_target(params)
target.x = 13.729;
target.state = [13.729; 19.585; params.zSection];
target.tolerance = 0.008;
end

function control = local_control()
control.deltaP = 0.1;
control.maxSearchCrossings = 12;
control.numParameterSamples = 81;
control.bisectionIterations = 24;
control.maxDiscontinuityIsolationDepth = 24;
control.maxDiscontinuityIntervals = 256;
end

function noise = local_noise()
noise.dt = 0.001;
noise.maxTimeToNextValidSection = 80;
noise.maxControlledCrossings = 60;
noise.maxRejectedCrossings = 250;
noise.maxStateNorm = 1e5;
noise.maxRk4Steps = 12000;
noise.paperReferenceStepCount = 40;
noise.retargetCadencePolicy = 'meanReturnCycleMatched';
noise.axisDefinition = 'per-step per-coordinate Gaussian standard deviation';
noise.directlyComparableToPaperFigure6 = false;
end

function mapData = local_nominal_map(params)
cacheFile = fullfile(fileparts(mfilename('fullpath')), 'shinbrot_nominal_map.mat');
if isfile(cacheFile)
    cached = load(cacheFile, 'mapData');
    mapData = cached.mapData;
    return;
end
fprintf('Generating the nominal return map (1200 crossings)...\n');
mapData = generate_return_map([1; 1; 1], 1200, params, 0);
save(cacheFile, 'mapData', '-v7.3');
end

function noise = ensure_global_step_budget(noise, mapData)
meanReturnTime = local_mean_return_time(mapData);
if ~isfinite(meanReturnTime) || meanReturnTime <= 0
    return;
end
safetyFactor = 2.5;
requiredSteps = noise.maxControlledCrossings * meanReturnTime / noise.dt;
noise.maxRk4Steps = max(noise.maxRk4Steps, ceil(safetyFactor * requiredSteps));
noise.meanReturnTimeForBudget = meanReturnTime;
noise.estimatedStepsPerCrossing = meanReturnTime / noise.dt;
noise.budgetSafetyFactor = safetyFactor;
noise.requiredRk4Steps = ceil(safetyFactor * requiredSteps);
if isfield(noise, 'retargetCadencePolicy') && ...
        strcmp(noise.retargetCadencePolicy, 'meanReturnCycleMatched')
    noise.retargetStepInterval = max(1, round(meanReturnTime / noise.dt));
    noise.retargetPhysicalTime = noise.retargetStepInterval * noise.dt;
    noise.cadenceRatio = noise.retargetPhysicalTime / meanReturnTime;
end
end
