function results = test_shinbrot_paper_repair()
%TEST_SHINBROT_PAPER_REPAIR Focused guards for the repaired paper path.
%
% The test deliberately avoids the full Monte Carlo sweep.  It checks the
% paper bisection implementation and configuration fingerprints.  No
% reduced noise sweep or reduced-result artifact is generated here.

projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);
rng(7);

[params, target, control, noise] = local_settings();
sourceIndex = 150;
mapData = generate_return_map([1; 1; 1], 220, params, 0);
sourceState = mapData.states(:, sourceIndex);
noise = local_step_budget(noise, mapData, 100);

results = struct();
results.paperFunction = which('run_shinbrot_paper_bisection');
assert(~isempty(results.paperFunction), ...
    'run_shinbrot_paper_bisection was not found on the MATLAB path.');
assert(strcmp(verified_bisection_identifier(), 'shinbrotPaperBisection'), ...
    'verified_bisection_identifier does not name the paper path.');

paperSource = fileread(results.paperFunction);
forbiddenTokens = {'Lipschitz', 'minCertifiedDepth', ...
    'linspace', 'grid-first', 'hybridGridBisection', ...
    'selectSmallestParameter'};
for k = 1:numel(forbiddenTokens)
    assert(~contains(paperSource, forbiddenTokens{k}), ...
        'Forbidden token present in paper path: %s', forbiddenTokens{k});
end

paperSearch = search_parameter_to_target( ...
    sourceState, target, params, control, 'shinbrotPaperBisection');
assert(strcmp(paperSearch.method, 'shinbrotPaperBisection'));
assert(paperSearch.gridEvaluations == 0, ...
    'Paper bisection performed grid evaluations.');
assert(paperSearch.crossingOrderTestUsed, ...
    'Paper bisection did not use crossing-order signatures.');
assert(isequal(paperSearch.horizonsExamined(:).', ...
    1:numel(paperSearch.horizonsExamined)), ...
    'Horizons were not examined in ascending order.');
assert(paperSearch.parameterEvaluations > 0, ...
    'No parameter evaluations were recorded.');
assert(paperSearch.crossingSignatureEvaluations == ...
    paperSearch.parameterEvaluations, ...
    'Every paper evaluation should record a crossing signature.');

if paperSearch.found
    assert(abs(paperSearch.selectedP) <= control.deltaP + 10*eps, ...
        'Selected p lies outside the paper control bound.');
    assert(paperSearch.paperReplayVerified, ...
        'Selected paper p was not replay verified.');
    assert(paperSearch.finalTargetError <= target.tolerance, ...
        'Replay target error exceeds tolerance.');
    if paperSearch.refinementsCompleted > 0
        assert(paperSearch.refinementsCompleted == control.bisectionIterations, ...
            'Final bracket did not complete the declared 24 refinements.');
    end
end

dispatcherSource = fileread(which('search_parameter_to_target'));
assert(contains(dispatcherSource, 'case ''shinbrotPaperBisection''') && ...
    contains(dispatcherSource, ...
    'case ''shinbrotDiscontinuityAwareBisection'''), ...
    'Paper and historical extension dispatcher cases are not distinct.');

currentConfig = build_bisection_noise_config( ...
    sourceState, target, params, control, noise, local_sweep());
oldConfig = currentConfig;
oldConfig.modes = {'paper40StepRetargetMode'};
oldConfig.retargetingMode = 'paper40StepRetargetMode';
oldConfig.retargetStepInterval = 40;
[configValid, differences] = validate_bisection_noise_config( ...
    oldConfig, currentConfig);
assert(~configValid && any(contains(differences, 'modes')), ...
    'A stale paper40 configuration was not rejected.');

rng(1234);
drawsA = randn(3, 2000);
rng(1234);
drawsB = randn(3, 2000);
assert(isequal(drawsA, drawsB), ...
    'Gaussian draws are not reproducible from the seed.');
drawCorrelation = corrcoef(drawsA.');
offDiagonal = drawCorrelation(~eye(3));
assert(max(abs(offDiagonal)) < 0.08, ...
    'Gaussian coordinate draws are unexpectedly correlated.');

results.paperSearch = paperSearch;
results.configurationDifferencesForOldMode = differences;
fprintf('test_shinbrot_paper_repair passed.\n');
end

function [params, target, control, noise] = local_settings()
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
control.searchMethod = 'shinbrotPaperBisection';

noise.dt = 0.001;
noise.maxTimeToNextValidSection = 80;
noise.maxControlledCrossings = 60;
noise.maxRejectedCrossings = 250;
noise.maxStateNorm = 1e5;
noise.budgetSafetyFactor = 2.5;
noise.paperReferenceStepCount = 40;
noise.retargetCadencePolicy = 'meanReturnCycleMatched';
end

function sweep = local_sweep()
sweep.sigmaNoiseValues = ...
    [0, 0.01, 0.03, 0.05, 0.075, 0.1, ...
     0.15, 0.2, 0.3, 0.5, 1.0];
sweep.numTrials = 10;
sweep.baseSeed = 91000;
sweep.modes = {'paperCycleCadenceRetargetMode'};
sweep.totalTrials = numel(sweep.sigmaNoiseValues) * ...
    sweep.numTrials * numel(sweep.modes);
sweep.searchMethod = 'shinbrotPaperBisection';
sweep.sourceIndex = 150;
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
end

function noise = local_step_budget(noise, mapData, burnIn)
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
