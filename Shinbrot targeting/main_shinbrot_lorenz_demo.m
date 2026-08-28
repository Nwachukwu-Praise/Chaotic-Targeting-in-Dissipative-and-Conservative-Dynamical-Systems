%MAIN_SHINBROT_LORENZ_DEMO Parameter-perturbation targeting for Lorenz.
%
% Two numerical searches are compared: the independent Shinbrot-style
% full-interval method and the modified grid-first hybrid method. In both,
% one constant parameter perturbation is applied from the original section
% source through the selected targeting horizon, then removed after the hit.

clear; clc; close all;
rng(7);

%% Lorenz system and accepted Poincare half-plane
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

%% Fixed-point target and bounded search
target.x = 13.729;
target.state = [13.729; 19.585; params.zSection];
target.tolerance = 0.008;

control.deltaP = 0.1;
control.maxSearchCrossings = 12;
control.numParameterSamples = 81;
control.bisectionIterations = 24;
control.maxDiscontinuityIsolationDepth = 24;
control.maxDiscontinuityIntervals = 256;
control.searchMethod = 'hybridGridBisection';

benchmarkSettings.repetitions = 5;
runSourceEnsembleStudy = false;
sourceEnsembleSettings.candidateIndices = 100:300;
sourceEnsembleSettings.numSources = 60;
sourceEnsembleSettings.sourceSeparation = 1.5;
sourceEnsembleSettings.outputMatFile = 'shinbrot_source_ensemble.mat';
sourceEnsembleSettings.methodResultsCsv = ...
    'shinbrot_source_ensemble_method_results.csv';
sourceEnsembleSettings.pairedResultsCsv = ...
    'shinbrot_source_ensemble_paired_results.csv';

%% Generate the nominal return map and choose one source
fprintf('Generating the uncontrolled Poincare return map...\n');
X0 = [1; 1; 1];
numMapCrossings = 1200;
mapData = generate_return_map(X0, numMapCrossings, params, 0);
if size(mapData.states, 2) < numMapCrossings
    error('Only %d valid section crossings were generated.', ...
        size(mapData.states, 2));
end

sourceIndex = 150;
sourceState = mapData.states(:, sourceIndex);
while abs(sourceState(1) - target.x) <= target.tolerance
    sourceIndex = sourceIndex + 1;
    sourceState = mapData.states(:, sourceIndex);
end
sourceX = sourceState(1);

%% Deterministic targeting, recurrence, RMS, and timing benchmark
fprintf('Running the uncontrolled recurrence baseline...\n');
uncontrolled = uncontrolled_hitting_time_1d( ...
    sourceState, target, params, 2000);

fprintf('Running the grid-first hybrid targeting search...\n');
controlled = shinbrot_parameter_targeting( ...
    sourceState, target, params, control);
search = controlled.search;

fprintf('Computing the post-transient time-weighted RHS RMS...\n');
rhsRms = compute_rhs_rms(mapData, params, control, 100);

fprintf(['Benchmarking the Shinbrot discontinuity-aware bisection and ', ...
    'grid-first hybrid methods...\n']);
searchBenchmark = benchmark_shinbrot_search_methods( ...
    sourceState, target, params, control, benchmarkSettings);

%% Saved or explicitly requested sixty-source benchmark
sourceEnsemble = [];
if runSourceEnsembleStudy
    fprintf('Running the sixty-source deterministic benchmark once...\n');
    sourceEnsemble = benchmark_shinbrot_source_ensemble( ...
        mapData, target, params, control, sourceEnsembleSettings);
elseif isfile(sourceEnsembleSettings.outputMatFile)
    savedEnsemble = load(sourceEnsembleSettings.outputMatFile);
    if is_valid_source_ensemble(savedEnsemble, target, control)
        fprintf('Loading the saved sixty-source deterministic benchmark...\n');
        sourceEnsemble = savedEnsemble;
    else
        fprintf(['Saved ensemble results do not match the current ', ...
            'two-method settings; no ensemble search was run.\n']);
    end
else
    fprintf(['No saved sixty-source benchmark found. Set ', ...
        'runSourceEnsembleStudy = true to generate it.\n']);
end

%% Three benchmark/report views
fprintf('\n========== Benchmark 1: Single-Source Median Timing ==========\n');
disp(searchBenchmark.comparisonTable);
if ~isempty(sourceEnsemble)
    fprintf('\n========== Benchmark 2: Sixty-Source Method Summary ==========\n');
    disp(sourceEnsemble.methodSummary);
    fprintf('\n========== Full 120-Row Method Results ==========\n');
    disp(sourceEnsemble.methodResultsTable);
    fprintf('\n========== Benchmark 3: Sixty-Source Paired Summary ==========\n');
    disp(sourceEnsemble.summary);
    fprintf('\n========== Full 60-Row Paired Results ==========\n');
    disp(sourceEnsemble.pairedResultsTable);
    plot_shinbrot_source_ensemble(sourceEnsemble.pairedResultsTable);
end

%% Deterministic result summary
fprintf('\n========== Shinbrot Lorenz Targeting ==========\n');
fprintf('Target X_t:                               %.6f\n', target.x);
fprintf('Source index and X_s:                     %d, %.6f\n', ...
    sourceIndex, sourceX);
fprintf('Selected p and horizon:                   %.12g, %g\n', ...
    controlled.selectedParameter, search.horizon);
fprintf('Selection detail:                        %s\n', ...
    search.selectionMethod);
fprintf('Complete grid horizons evaluated:        %g\n', ...
    search.gridHorizonsEvaluated);
fprintf('Grid evaluations:                        %g\n', ...
    search.gridEvaluations);
if isfinite(search.horizon)
    fprintf('Direct hits at selected horizon:         %g\n', ...
        search.directGridHitsByHorizon(search.horizon));
    fprintf('Adjacent brackets at selected horizon:   %g\n', ...
        search.adjacentBracketsByHorizon(search.horizon));
    fprintf('Verified candidates at selected horizon: %g\n', ...
        search.verifiedCandidatesByHorizon(search.horizon));
end
fprintf('Final/total bisection iterations:         %g / %g\n', ...
    search.finalBracketBisectionIterations, ...
    search.totalBisectionIterations);
fprintf('Final target error:                      %.12g\n', ...
    search.finalTargetError);
fprintf('Controlled/uncontrolled crossings:       %g / %g\n', ...
    controlled.numCrossings, uncontrolled.numCrossings);
fprintf('Time-weighted RMS (f_x, f_y, f_z):        %.6g, %.6g, %.6g\n', ...
    rhsRms.rmsFx, rhsRms.rmsFy, rhsRms.rmsFz);
fprintf('DeltaP / RMS f_y:                        %.6g\n', ...
    rhsRms.parameterToRmsFyRatio);
fprintf('================================================\n');

%% Visual results
plot_shinbrot_results( ...
    mapData, uncontrolled, controlled, target, params, control);

function rhsRms = compute_rhs_rms( ...
    mapData, params, control, transientCrossingIndex)
%COMPUTE_RHS_RMS Time-weighted RMS over a nominal post-transient trajectory.

if transientCrossingIndex < 1 || ...
        transientCrossingIndex > numel(mapData.eventTimes)
    error('Transient crossing index %d is unavailable.', ...
        transientCrossingIndex);
end
transientTime = mapData.eventTimes(transientCrossingIndex);
keep = mapData.t >= transientTime;
t = mapData.t(keep);
X = mapData.X(keep, :);
if numel(t) < 2
    error('At least two post-transient samples are required for RMS.');
end
duration = t(end) - t(1);
if ~isfinite(duration) || duration <= 0
    error('RMS averaging duration must be finite and positive.');
end

F = NaN(size(X));
for i = 1:size(X, 1)
    F(i, :) = lorenz_rhs(t(i), X(i, :).', params, 0).';
end
rhsRms.rmsFx = sqrt(trapz(t, F(:, 1).^2) / duration);
rhsRms.rmsFy = sqrt(trapz(t, F(:, 2).^2) / duration);
rhsRms.rmsFz = sqrt(trapz(t, F(:, 3).^2) / duration);
rhsRms.transientCrossingIndex = transientCrossingIndex;
rhsRms.transientTime = transientTime;
rhsRms.averagingDuration = duration;
rhsRms.maxParameterEffect = control.deltaP;
rhsRms.parameterToRmsFyRatio = control.deltaP / rhsRms.rmsFy;
rhsRms.paperApproximateRmsFy = 100;
rhsRms.rmsFyDifferenceFromPaperApproximation = ...
    rhsRms.rmsFy - rhsRms.paperApproximateRmsFy;
end

function valid = is_valid_source_ensemble(saved, target, control)
%IS_VALID_SOURCE_ENSEMBLE Validate saved results without recomputation.

required = {'sourceIndices', 'methodResultsTable', ...
    'pairedResultsTable', 'summary', 'methodSummary', ...
    'studySettings', 'target', 'control'};
valid = all(isfield(saved, required));
if ~valid
    return;
end

methods = saved.methodResultsTable.Method;
methodVariables = saved.methodResultsTable.Properties.VariableNames;
pairedVariables = saved.pairedResultsTable.Properties.VariableNames;
valid = istable(saved.methodResultsTable) && ...
    istable(saved.pairedResultsTable) && istable(saved.summary) && ...
    istable(saved.methodSummary) && ...
    numel(saved.sourceIndices) == 60 && ...
    numel(unique(saved.sourceIndices)) == 60 && ...
    all(saved.sourceIndices >= 100 & saved.sourceIndices <= 300) && ...
    height(saved.methodResultsTable) == 120 && ...
    height(saved.pairedResultsTable) == 60 && ...
    sum(strcmp(methods, verified_bisection_identifier())) == 60 && ...
    sum(strcmp(methods, 'hybridGridBisection')) == 60 && ...
    all(saved.methodResultsTable.GridEvaluations(strcmp(methods, ...
        verified_bisection_identifier())) == 0) && ...
    all(ismember({'FailureCode', 'FailureReason'}, methodVariables)) && ...
    all(ismember({'ShinbrotFinalTargetError', ...
        'HybridFinalTargetError', 'ShinbrotReplayTargetError', ...
        'HybridReplayTargetError'}, pairedVariables)) && ...
    saved.studySettings.numSources == 60 && ...
    saved.studySettings.sourceSeparation == 1.5 && ...
    isequal(saved.studySettings.candidateIndices(:).', 100:300) && ...
    saved.control.maxSearchCrossings == control.maxSearchCrossings && ...
    saved.control.deltaP == control.deltaP && ...
    saved.target.x == target.x && ...
    saved.target.tolerance == target.tolerance;
end
