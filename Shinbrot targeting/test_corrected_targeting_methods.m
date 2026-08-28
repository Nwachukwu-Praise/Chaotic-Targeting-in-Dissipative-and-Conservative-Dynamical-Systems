function report = test_corrected_targeting_methods()
%TEST_CORRECTED_TARGETING_METHODS Required grid logic and Lorenz replay tests.

[params, target, control] = settings();
mapData = generate_return_map([1; 1; 1], 300, params, 0);
sourceState = mapData.states(:, 150);

shinbrot = search_parameter_to_target(sourceState, target, params, ...
    control, verified_bisection_identifier());
hybrid = search_parameter_to_target(sourceState, target, params, ...
    control, 'hybridGridBisection');

names = {
    'grid runs despite equal-sign outer endpoints'
    'direct interior grid hit is detected'
    'multiple direct hits are retained before selection'
    'all adjacent brackets are identified and refined'
    'hybrid bisection stays inside one grid interval'
    'sign change alone is not a verified candidate'
    'smallest successful horizon is selected'
    'smallest absolute p is selected'
    'target error is the final tie-breaker'
    'failure is structured and numerical'
    'selected p replays from original 3D source'
    'both methods share dynamical and numerical settings'
    'coarse-grid-only method is not callable'};
passed = false(13, 1);

features = find_hybrid_grid_features([-0.1, 0, 0.1], [1, -1, 1], 0.01);
assert(isempty(features.directHitIndices));
assert(size(features.bracketIndexPairs, 1) == 2);
passed(1) = true;

features = find_hybrid_grid_features([-0.1, 0, 0.1], [2, 0.001, 3], 0.008);
assert(isequal(features.directHitIndices, 2));
passed(2) = true;

features = find_hybrid_grid_features( ...
    [-0.1, -0.02, 0.01, 0.1], [1, 0.004, -0.006, 1], 0.008);
assert(isequal(features.directHitIndices, [2, 3]));
candidates = candidate_array([-0.02, 0.01], [0.004, 0.006]);
[selected, retained] = select_verified_hybrid_candidate(candidates);
assert(numel(retained) == 2 && selected.p == 0.01);
passed(3) = true;

features = find_hybrid_grid_features( ...
    [-0.1, -0.05, 0, 0.05, 0.1], [1, -1, 1, -1, 1], 0.008);
assert(size(features.bracketIndexPairs, 1) == 4);
assert(hybrid.totalBisectionIterations >= ...
    sum(hybrid.adjacentBracketsByHorizon));
passed(4) = true;

gridSpacing = 2 * control.deltaP / (control.numParameterSamples - 1);
assert(hybrid.finalBracketBisectionIterations == 0 || ...
    abs((hybrid.bracketHigh - hybrid.bracketLow) - gridSpacing) < 1e-12);
assert(all(abs(diff(features.bracketIntervals, 1, 2) - 0.05) < 1e-12));
passed(5) = true;

discontinuousFeatures = find_hybrid_grid_features( ...
    [-0.01, 0.01], [-1, 1], 0.008);
assert(size(discontinuousFeatures.bracketIndexPairs, 1) == 1);
selectionFailed = false;
try
    select_verified_hybrid_candidate(struct( ...
        'p', {}, 'x', {}, 'error', {}, 'origin', {}, ...
        'bracketLow', {}, 'bracketHigh', {}, 'iterations', {}));
catch
    selectionFailed = true;
end
assert(selectionFailed);
passed(6) = true;

assert(hybrid.found && hybrid.gridEvaluations == ...
    control.numParameterSamples * hybrid.horizon);
assert(all(hybrid.verifiedCandidatesByHorizon(1:hybrid.horizon-1) == 0));
passed(7) = true;

candidates = candidate_array([-0.03, 0.01, 0.02], [1e-4, 7e-3, 2e-3]);
selected = select_verified_hybrid_candidate(candidates);
assert(selected.p == 0.01);
passed(8) = true;

candidates = candidate_array([-0.02, 0.02], [6e-3, 2e-3]);
selected = select_verified_hybrid_candidate(candidates);
assert(selected.p == 0.02);
passed(9) = true;

unreachableTarget = target;
unreachableTarget.x = 1e6;
unreachableTarget.tolerance = 1e-12;
failureControl = control;
failureControl.maxSearchCrossings = 1;
failed = search_parameter_to_target(sourceState, unreachableTarget, ...
    params, failureControl, 'hybridGridBisection');
assert(~failed.found && strcmp(failed.failureCode, 'noVerifiedCandidate'));
assert(~isempty(failed.failureReason) && isnan(failed.selectedP));
passed(10) = true;

[replayHit, replayError] = replay_search( ...
    sourceState, target, params, hybrid);
assert(replayHit && replayError <= target.tolerance);
passed(11) = true;

assert(shinbrot.gridEvaluations == 0);
assert(isequal(shinbrot.pGrid, hybrid.pGrid));
assert(size(shinbrot.xByHorizon, 1) == control.maxSearchCrossings);
assert(size(hybrid.xByHorizon, 1) == control.maxSearchCrossings);
assert(shinbrot.found && hybrid.found);
passed(12) = true;

coarseRejected = false;
try
    search_parameter_to_target(sourceState, target, params, ...
        control, 'coarseGridOnlyBaseline');
catch exception
    coarseRejected = contains(exception.message, 'Unknown search method');
end
assert(coarseRejected);
passed(13) = true;

report = table((1:13).', names, passed, ...
    'VariableNames', {'TestNumber', 'Requirement', 'Passed'});
disp(report);
assert(all(passed));
end

function candidates = candidate_array(pValues, errors)
%Candidate records for deterministic selection tests.

template = struct('p', NaN, 'x', NaN, 'error', NaN, ...
    'origin', 'synthetic verified candidate', ...
    'bracketLow', NaN, 'bracketHigh', NaN, 'iterations', 0);
candidates = repmat(template, 1, numel(pValues));
for i = 1:numel(pValues)
    candidates(i).p = pValues(i);
    candidates(i).x = errors(i);
    candidates(i).error = errors(i);
    candidates(i).bracketLow = pValues(i);
    candidates(i).bracketHigh = pValues(i);
end
end

function [hit, targetError] = replay_search(sourceState, target, params, search)
%Directly replay a selected parameter from the original full source state.

currentState = sourceState(:);
for crossing = 1:search.horizon
    segment = next_valid_section_crossing( ...
        currentState, params, search.selectedP);
    assert(segment.success);
    currentState = segment.eventState;
end
targetError = abs(currentState(1) - target.x);
hit = targetError <= target.tolerance;
end

function [params, target, control] = settings()
%Shared dynamical and numerical settings.

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
end
