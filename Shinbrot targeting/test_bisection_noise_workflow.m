function testResults = test_bisection_noise_workflow()
%TEST_BISECTION_NOISE_WORKFLOW Non-generating checks for full noise config.
%
% This test deliberately does not run a reduced Monte Carlo sweep and does
% not write noise MAT/CSV files.  The numerical zero-noise propagation gate
% is covered by test_zero_noise_gate.

close all;
verifiedIdentifier = verified_bisection_identifier();

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
control.searchMethod = verifiedIdentifier;

mapData = generate_return_map([1; 1; 1], 200, params, 0);
sourceIndex = 150;
sourceState = mapData.states(:, sourceIndex);

bisectionSearch = search_parameter_to_target( ...
    sourceState, target, params, control, verifiedIdentifier);
assert(bisectionSearch.found, 'Verified bisection did not find a hit.');
assert_verified_bisection_search( ...
    bisectionSearch, verifiedIdentifier, 'targeted deterministic bisection', target);

noise.dt = 0.001;
noise.maxTimeToNextValidSection = 80;
noise.maxControlledCrossings = 60;
noise.maxRejectedCrossings = 250;
noise.maxStateNorm = 1e5;
noise.budgetSafetyFactor = 2.5;
noise.paperReferenceStepCount = 40;
noise.retargetCadencePolicy = 'meanReturnCycleMatched';
noise.axisDefinition = 'per-step per-coordinate Gaussian standard deviation';
noise.directlyComparableToPaperFigure6 = false;
noise = set_measured_noise_step_budget(noise, mapData, 100);

sweep.sigmaNoiseValues = ...
    [0, 0.01, 0.03, 0.05, 0.075, 0.1, ...
     0.15, 0.2, 0.3, 0.5, 1.0];
sweep.numTrials = 10;
sweep.baseSeed = 91000;
sweep.modes = {'paperCycleCadenceRetargetMode'};
sweep.totalTrials = numel(sweep.sigmaNoiseValues) * ...
    sweep.numTrials * numel(sweep.modes);
sweep.searchMethod = verifiedIdentifier;
sweep.sourceIndex = sourceIndex;
sweep.outputMatFile = 'noise_results_paper_cycle_retarget_11sigma.mat';
sweep.checkpointMatFile = 'noise_results_paper_cycle_retarget_11sigma_checkpoint.mat';
sweep.summaryCsvFile = 'noise_summary_paper_cycle_retarget_11sigma.csv';
sweep.trialsCsvFile = 'noise_trials_paper_cycle_retarget_11sigma.csv';
sweep.progressLogFile = 'noise_progress_paper_cycle_retarget_11sigma.log';
sweep.resumeFromCheckpoint = true;
sweep.isSmokeTest = false;
sweep.trajectoryRetentionPolicy = ...
    ['retain full RK4 trajectories only for trial 1 at the lowest, ', ...
    'middle, and highest configured sigma_noise levels'];

config = build_bisection_noise_config( ...
    sourceState, target, params, control, noise, sweep);
assert(config.totalTrials == 110);
assert(noise.maxControlledCrossings == 60);
assert(noise.maxRk4Steps >= noise.requiredRk4Steps);
assert(isequal(sweep.sigmaNoiseValues, ...
    [0, 0.01, 0.03, 0.05, 0.075, 0.1, ...
     0.15, 0.2, 0.3, 0.5, 1.0]));
assert(sweep.numTrials == 10);
assert(sweep.totalTrials == 110);
assert(~sweep.isSmokeTest);
assert(strcmp(sweep.modes{1}, 'paperCycleCadenceRetargetMode'));
assert(abs(noise.cadenceRatio - 1) <= 0.05);
assert(abs(noise.retargetStepInterval - round(noise.meanReturnTimeForBudget / noise.dt)) <= 0);
assert(noise.retargetStepInterval > 1000);

segment = next_valid_section_crossing_noisy_rk4( ...
    sourceState, params, set_sigma_noise(noise, 0), bisectionSearch.selectedP);
assert(segment.crossingUsesIntegratedSegmentEndpoints);
assert(segment.instantaneousNoiseJumpCrossingsExcluded);

badModeRejected = false;
try
    run_noisy_targeting_trial(sourceState, target, params, control, ...
        noise, 0, 1, 1234, 'paperRetargetMode', bisectionSearch, false);
catch
    badModeRejected = true;
end
assert(badModeRejected, 'Obsolete paperRetargetMode was not rejected.');

mismatchConfig = config;
mismatchConfig.deltaP = 999;
[validConfig, ~] = validate_bisection_noise_config(config, mismatchConfig);
assert(~validConfig, 'Configuration mismatch was not rejected.');

oldCheckpointRejected = false;
if isfile('noise_results_paper40_step_retarget_checkpoint.mat')
    old = load('noise_results_paper40_step_retarget_checkpoint.mat');
    if isfield(old, 'checkpoint') && isfield(old.checkpoint, 'configuration')
        [oldValid, oldDifferences] = validate_bisection_noise_config( ...
            old.checkpoint.configuration, config);
        oldCheckpointRejected = ~oldValid && ...
            any(contains(oldDifferences, 'retargetingMode') | ...
            contains(oldDifferences, 'modes') | ...
            contains(oldDifferences, 'retargetStepInterval'));
    end
else
    oldCheckpointRejected = true;
end
assert(oldCheckpointRejected, ...
    'The old literal-40-step checkpoint was not rejected.');

[~, missingAccepted, missingDifferences] = load_valid_bisection_noise_results( ...
    'definitely_missing_full_noise_result.mat', config);
assert(~missingAccepted && any(contains(missingDifferences, "file not found")), ...
    'Missing full result was not reported accurately.');

testResults.verifiedIdentifier = verifiedIdentifier;
testResults.bisectionSelectedP = bisectionSearch.selectedP;
testResults.bisectionHorizon = bisectionSearch.horizon;
testResults.totalConfiguredTrials = config.totalTrials;
testResults.badModeRejected = badModeRejected;
testResults.missingFullResultRejected = ~missingAccepted;
testResults.oldPaper40CheckpointRejected = oldCheckpointRejected;
testResults.retargetStepInterval = noise.retargetStepInterval;
testResults.cadenceRatio = noise.cadenceRatio;
testResults.noiseMaxControlledCrossings = noise.maxControlledCrossings;
testResults.noiseMaxRk4Steps = noise.maxRk4Steps;

disp(struct2table(testResults));
fprintf('Bisection-only full-noise workflow validation checks passed.\n');
end

function noiseOut = set_sigma_noise(noiseIn, sigmaNoise)
noiseOut = noiseIn;
noiseOut.sigmaNoise = sigmaNoise;
end

function noise = set_measured_noise_step_budget(noise, mapData, burnIn)
eventTimes = mapData.eventTimes(:);
intervals = diff(eventTimes((burnIn + 1):end));
intervals = intervals(isfinite(intervals) & intervals > 0);
if isempty(intervals)
    error('No valid return intervals were available for the RK4 budget.');
end
noise.budgetBurnIn = burnIn;
noise.meanReturnTimeForBudget = mean(intervals);
noise.estimatedStepsPerCrossing = noise.meanReturnTimeForBudget / noise.dt;
noise.requiredRk4Steps = ceil(noise.budgetSafetyFactor * ...
    noise.maxControlledCrossings * noise.estimatedStepsPerCrossing);
noise.maxRk4Steps = noise.requiredRk4Steps;
noise.retargetStepInterval = max(1, round(noise.meanReturnTimeForBudget / noise.dt));
noise.retargetPhysicalTime = noise.retargetStepInterval * noise.dt;
noise.cadenceRatio = noise.retargetPhysicalTime / noise.meanReturnTimeForBudget;
end
