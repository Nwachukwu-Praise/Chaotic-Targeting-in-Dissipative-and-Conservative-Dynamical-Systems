function benchmark = benchmark_shinbrot_search_methods( ...
    sourceState, target, params, control, benchmarkSettings)
%BENCHMARK_SHINBROT_SEARCH_METHODS Compare the two retained search methods.
%
% The benchmark starts after source selection and times only parameter search
% and controlled replay. Return-map generation, plotting, and table assembly
% are outside the timing boundary.

if nargin < 5 || isempty(benchmarkSettings)
    benchmarkSettings.repetitions = 5;
end
if ~isfield(benchmarkSettings, 'repetitions')
    benchmarkSettings.repetitions = 5;
end

methods = {verified_bisection_identifier(), 'hybridGridBisection'};
numMethods = numel(methods);

Method = methods(:);
HitFound = false(numMethods, 1);
SelectedP = NaN(numMethods, 1);
Horizon = NaN(numMethods, 1);
FinalError = NaN(numMethods, 1);
ParameterEvaluations = NaN(numMethods, 1);
CrossingPropagations = NaN(numMethods, 1);
FinalBracketBisectionIterations = NaN(numMethods, 1);
TotalBisectionIterations = NaN(numMethods, 1);
GridEvaluations = NaN(numMethods, 1);
MedianSearchTime = NaN(numMethods, 1);
MedianReplayTime = NaN(numMethods, 1);
MedianSearchAndReplayTime = NaN(numMethods, 1);
ControlledReplayCrossings = NaN(numMethods, 1);
ReplayHit = false(numMethods, 1);
SelectionMethod = cell(numMethods, 1);
FailureCode = cell(numMethods, 1);
FailureReason = cell(numMethods, 1);
searchResults = cell(numMethods, 1);
replayResults = cell(numMethods, 1);

for methodIndex = 1:numMethods
    methodControl = control;
    methodControl.searchMethod = methods{methodIndex};

    warmSearch = search_parameter_to_target( ...
        sourceState, target, params, methodControl);
    replay_controlled_search(sourceState, target, params, warmSearch);

    searchTimes = NaN(benchmarkSettings.repetitions, 1);
    replayTimes = NaN(benchmarkSettings.repetitions, 1);
    totalTimes = NaN(benchmarkSettings.repetitions, 1);
    timedSearch = warmSearch;
    timedReplay = replay_controlled_search( ...
        sourceState, target, params, warmSearch);

    for rep = 1:benchmarkSettings.repetitions
        totalTimer = tic;
        searchTimer = tic;
        timedSearch = search_parameter_to_target( ...
            sourceState, target, params, methodControl);
        searchTimes(rep) = toc(searchTimer);

        replayTimer = tic;
        timedReplay = replay_controlled_search( ...
            sourceState, target, params, timedSearch);
        replayTimes(rep) = toc(replayTimer);
        totalTimes(rep) = toc(totalTimer);
    end

    searchResults{methodIndex} = timedSearch;
    replayResults{methodIndex} = timedReplay;
    HitFound(methodIndex) = timedSearch.found;
    SelectedP(methodIndex) = timedSearch.selectedP;
    if isfinite(timedSearch.horizon)
        Horizon(methodIndex) = timedSearch.horizon;
    end
    FinalError(methodIndex) = timedSearch.finalTargetError;
    ParameterEvaluations(methodIndex) = timedSearch.parameterEvaluations;
    CrossingPropagations(methodIndex) = timedSearch.crossingPropagations;
    FinalBracketBisectionIterations(methodIndex) = ...
        timedSearch.finalBracketBisectionIterations;
    TotalBisectionIterations(methodIndex) = ...
        timedSearch.totalBisectionIterations;
    GridEvaluations(methodIndex) = timedSearch.gridEvaluations;
    MedianSearchTime(methodIndex) = median(searchTimes);
    MedianReplayTime(methodIndex) = median(replayTimes);
    MedianSearchAndReplayTime(methodIndex) = median(totalTimes);
    if isfinite(timedReplay.numCrossings)
        ControlledReplayCrossings(methodIndex) = timedReplay.numCrossings;
    end
    ReplayHit(methodIndex) = timedReplay.hit;
    SelectionMethod{methodIndex} = timedSearch.selectionMethod;
    FailureCode{methodIndex} = timedSearch.failureCode;
    FailureReason{methodIndex} = timedSearch.failureReason;
end

comparisonTable = table(Method, HitFound, SelectedP, Horizon, ...
    FinalError, ParameterEvaluations, CrossingPropagations, ...
    FinalBracketBisectionIterations, TotalBisectionIterations, ...
    GridEvaluations, MedianSearchTime, ControlledReplayCrossings, ...
    ReplayHit, MedianReplayTime, MedianSearchAndReplayTime, ...
    SelectionMethod, FailureCode, FailureReason);

benchmark.methods = methods;
benchmark.repetitions = benchmarkSettings.repetitions;
benchmark.comparisonTable = comparisonTable;
benchmark.searchResults = searchResults;
benchmark.replayResults = replayResults;
end

function replay = replay_controlled_search(sourceState, target, params, search)
%REPLAY_CONTROLLED_SEARCH Verify the selected p from the original 3D source.

replay.hit = false;
replay.numCrossings = Inf;
replay.finalError = Inf;
replay.failed = false;
if ~search.found || ~isfinite(search.selectedP) || ~isfinite(search.horizon)
    replay.failed = true;
    return;
end

currentState = sourceState(:);
for crossing = 1:search.horizon
    segment = next_valid_section_crossing( ...
        currentState, params, search.selectedP);
    if ~segment.success
        replay.failed = true;
        return;
    end
    currentState = segment.eventState;
end

replay.finalError = abs(currentState(1) - target.x);
replay.hit = replay.finalError <= target.tolerance;
replay.numCrossings = search.horizon;
end
