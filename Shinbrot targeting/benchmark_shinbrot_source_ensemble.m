function ensemble = benchmark_shinbrot_source_ensemble( ...
    mapData, target, params, control, studySettings)
%BENCHMARK_SHINBROT_SOURCE_ENSEMBLE Compare two methods over sixty sources.
%
% This diagnostic runs each existing targeting method exactly once per source.
% SearchTime is therefore a single-run descriptive time, not a replacement
% for the five-repetition computing-time benchmark.

if nargin < 5 || isempty(studySettings)
    studySettings = struct();
end
studySettings = apply_study_defaults(studySettings);

if control.maxSearchCrossings ~= 12
    error('The source ensemble requires control.maxSearchCrossings = 12.');
end

candidateIndices = studySettings.candidateIndices(:).';
if any(candidateIndices < 1) || ...
        any(candidateIndices > size(mapData.states, 2))
    error(['Candidate source indices must exist in mapData.states. ', ...
        'Available states: %d; requested range: %d:%d.'], ...
        size(mapData.states, 2), candidateIndices(1), candidateIndices(end));
end

sourceXByCandidate = mapData.states(1, candidateIndices);
validMask = abs(sourceXByCandidate - target.x) >= studySettings.sourceSeparation;
validIndices = candidateIndices(validMask);
if numel(validIndices) < studySettings.numSources
    error(['Only %d source points satisfying the source-target margin exist in ', ...
        'the requested range; %d are required.'], ...
        numel(validIndices), studySettings.numSources);
end

positions = round(linspace(1, numel(validIndices), ...
    studySettings.numSources));
sourceIndices = validIndices(positions);
if numel(sourceIndices) ~= studySettings.numSources || ...
        numel(unique(sourceIndices)) ~= studySettings.numSources
    error(['Deterministic spacing did not produce %d unique source ', ...
        'indices.'], studySettings.numSources);
end
if any(abs(mapData.states(1, sourceIndices) - target.x) < ...
        studySettings.sourceSeparation)
    error('A selected source violates the scalar source-target margin.');
end

methods = { ...
    verified_bisection_identifier(), ...
    'hybridGridBisection'};
numSources = numel(sourceIndices);
numMethods = numel(methods);
numRows = numSources * numMethods;

SourceIndex = NaN(numRows, 1);
SourceX = NaN(numRows, 1);
Method = cell(numRows, 1);
HitFound = false(numRows, 1);
ReplayHit = false(numRows, 1);
SelectedP = NaN(numRows, 1);
AbsoluteP = NaN(numRows, 1);
Horizon = NaN(numRows, 1);
FinalTargetError = NaN(numRows, 1);
ReplayTargetError = NaN(numRows, 1);
ParameterEvaluations = NaN(numRows, 1);
CrossingPropagations = NaN(numRows, 1);
FinalBracketBisectionIterations = NaN(numRows, 1);
TotalBisectionIterations = NaN(numRows, 1);
GridEvaluations = NaN(numRows, 1);
SearchTime = NaN(numRows, 1);
SelectionMethod = cell(numRows, 1);
FailureCode = cell(numRows, 1);
FailureReason = cell(numRows, 1);
% Discontinuity-awareness diagnostics.  Subsection 3.1.5 states that the
% counts of detected discontinuities and of reductions performed are
% reported with the results; these columns are what make that possible at
% the ensemble level.  The grid method performs no discontinuity isolation,
% so zeros are the expected and informative value for it.
DiscontinuitiesDetected = NaN(numRows, 1);
DiscontinuityReductions = NaN(numRows, 1);
ContinuousIntervalsExamined = NaN(numRows, 1);
CrossingSignatureEvaluations = NaN(numRows, 1);
CertifiedParameterResolution = NaN(numRows, 1);
ParameterResolutionFloorReached = false(numRows, 1);
IntervalBudgetExhausted = false(numRows, 1);
% Independence of the replay check.  A divergence of exactly zero means the
% search error and the replay error were produced by the same propagation
% and the replay verified nothing; see test_shinbrot_replay_negative_control.
ReplayDivergence = NaN(numRows, 1);
searchResults = cell(numSources, numMethods);
replayResults = cell(numSources, numMethods);

row = 0;
for sourceNumber = 1:numSources
    sourceIndex = sourceIndices(sourceNumber);
    sourceState = mapData.states(:, sourceIndex);

    for methodNumber = 1:numMethods
        row = row + 1;
        methodName = methods{methodNumber};
        methodControl = control;
        methodControl.searchMethod = methodName;

        search = search_parameter_to_target( ...
            sourceState, target, params, methodControl, methodName);
        replay = replay_selected_parameter( ...
            sourceState, target, params, search);

        searchResults{sourceNumber, methodNumber} = search;
        replayResults{sourceNumber, methodNumber} = replay;

        SourceIndex(row) = sourceIndex;
        SourceX(row) = sourceState(1);
        Method{row} = methodName;
        HitFound(row) = search.found;
        ReplayHit(row) = replay.hit;
        SelectedP(row) = search.selectedP;
        if isfinite(search.selectedP)
            AbsoluteP(row) = abs(search.selectedP);
        end
        if isfinite(search.horizon)
            Horizon(row) = search.horizon;
        end
        FinalTargetError(row) = search.finalTargetError;
        ReplayTargetError(row) = replay.targetError;
        ParameterEvaluations(row) = search.parameterEvaluations;
        CrossingPropagations(row) = search.crossingPropagations;
        FinalBracketBisectionIterations(row) = ...
            search.finalBracketBisectionIterations;
        TotalBisectionIterations(row) = search.totalBisectionIterations;
        GridEvaluations(row) = search.gridEvaluations;
        SearchTime(row) = search.searchTime;
        SelectionMethod{row} = search.selectionMethod;
        FailureCode{row} = search.failureCode;
        FailureReason{row} = search.failureReason;

        DiscontinuitiesDetected(row) = ...
            field_or_default(search, 'discontinuitiesDetected', NaN);
        DiscontinuityReductions(row) = ...
            field_or_default(search, 'discontinuityReductions', NaN);
        ContinuousIntervalsExamined(row) = ...
            field_or_default(search, 'continuousIntervalsExamined', NaN);
        CrossingSignatureEvaluations(row) = ...
            field_or_default(search, 'crossingSignatureEvaluations', NaN);
        CertifiedParameterResolution(row) = ...
            field_or_default(search, 'certifiedParameterResolution', NaN);
        ParameterResolutionFloorReached(row) = ...
            logical(field_or_default(search, 'parameterResolutionFloorReached', false));
        IntervalBudgetExhausted(row) = ...
            logical(field_or_default(search, 'intervalBudgetExhausted', false));
        ReplayDivergence(row) = ...
            abs(search.finalTargetError - replay.targetError);
    end
end

methodResultsTable = table(SourceIndex, SourceX, Method, HitFound, ...
    ReplayHit, SelectedP, AbsoluteP, Horizon, FinalTargetError, ...
    ReplayTargetError, ParameterEvaluations, CrossingPropagations, ...
    FinalBracketBisectionIterations, TotalBisectionIterations, ...
    GridEvaluations, SearchTime, SelectionMethod, FailureCode, ...
    FailureReason, DiscontinuitiesDetected, DiscontinuityReductions, ...
    ContinuousIntervalsExamined, CrossingSignatureEvaluations, ...
    CertifiedParameterResolution, ParameterResolutionFloorReached, ...
    IntervalBudgetExhausted, ReplayDivergence);

shinbrotMask = strcmp(methodResultsTable.Method, ...
    verified_bisection_identifier());
if any(methodResultsTable.GridEvaluations(shinbrotMask) ~= 0)
    error('The Shinbrot ensemble results unexpectedly used grid evaluations.');
end

pairedResultsTable = build_paired_table( ...
    methodResultsTable, sourceIndices);
summary = summarize_paired_results(pairedResultsTable);
methodSummary = summarize_methods(methodResultsTable);

ensemble.sourceIndices = sourceIndices;
ensemble.methodResultsTable = methodResultsTable;
ensemble.pairedResultsTable = pairedResultsTable;
ensemble.summary = summary;
ensemble.methodSummary = methodSummary;
ensemble.studySettings = studySettings;
ensemble.target = target;
ensemble.control = control;
ensemble.searchResults = searchResults;
ensemble.replayResults = replayResults;

save_ensemble_results(ensemble);
end

function settings = apply_study_defaults(settings)
%APPLY_STUDY_DEFAULTS Fill deterministic study and output settings.

if ~isfield(settings, 'candidateIndices')
    settings.candidateIndices = 100:300;
end
if ~isfield(settings, 'numSources')
    settings.numSources = 60;
end
if ~isfield(settings, 'sourceSeparation')
    settings.sourceSeparation = 1.5;
end
if ~isfield(settings, 'outputMatFile')
    settings.outputMatFile = 'shinbrot_source_ensemble.mat';
end
if ~isfield(settings, 'methodResultsCsv')
    settings.methodResultsCsv = ...
        'shinbrot_source_ensemble_method_results.csv';
end
if ~isfield(settings, 'pairedResultsCsv')
    settings.pairedResultsCsv = ...
        'shinbrot_source_ensemble_paired_results.csv';
end
settings.singleRunDescriptiveTiming = true;
end

function replay = replay_selected_parameter( ...
    sourceState, target, params, search)
%REPLAY_SELECTED_PARAMETER Independently verify the reported p and horizon.

replay.hit = false;
replay.targetError = Inf;
replay.failed = false;
replay.numCrossings = 0;

if ~search.found || ~isfinite(search.selectedP) || ...
        ~isfinite(search.horizon)
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
    replay.numCrossings = crossing;
end

replay.targetError = abs(currentState(1) - target.x);
replay.hit = replay.targetError <= target.tolerance;
end

function value = field_or_default(s, fieldName, defaultValue)
%FIELD_OR_DEFAULT Read a diagnostic field that not every search path sets.
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function paired = build_paired_table(methodResults, sourceIndices)
%BUILD_PAIRED_TABLE Put both method results on one row per source.

numSources = numel(sourceIndices);
SourceIndex = sourceIndices(:);
SourceX = NaN(numSources, 1);
ShinbrotHit = false(numSources, 1);
HybridHit = false(numSources, 1);
ShinbrotReplayHit = false(numSources, 1);
HybridReplayHit = false(numSources, 1);
ShinbrotP = NaN(numSources, 1);
HybridP = NaN(numSources, 1);
ShinbrotAbsoluteP = NaN(numSources, 1);
HybridAbsoluteP = NaN(numSources, 1);
AbsolutePDifference = NaN(numSources, 1);
ShinbrotHorizon = NaN(numSources, 1);
HybridHorizon = NaN(numSources, 1);
ShinbrotFinalTargetError = NaN(numSources, 1);
HybridFinalTargetError = NaN(numSources, 1);
ShinbrotReplayTargetError = NaN(numSources, 1);
HybridReplayTargetError = NaN(numSources, 1);
SameHorizon = false(numSources, 1);
BothReplayHit = false(numSources, 1);
SmallerMagnitudeMethod = repmat({'not comparable'}, numSources, 1);

for sourceNumber = 1:numSources
    sourceIndex = sourceIndices(sourceNumber);
    sourceRows = methodResults.SourceIndex == sourceIndex;
    shinbrotRow = find(sourceRows & strcmp(methodResults.Method, ...
        verified_bisection_identifier()));
    hybridRow = find(sourceRows & strcmp(methodResults.Method, ...
        'hybridGridBisection'));

    if numel(shinbrotRow) ~= 1 || numel(hybridRow) ~= 1
        error('Each source must have exactly one result from each method.');
    end

    SourceX(sourceNumber) = methodResults.SourceX(shinbrotRow);
    ShinbrotHit(sourceNumber) = methodResults.HitFound(shinbrotRow);
    HybridHit(sourceNumber) = methodResults.HitFound(hybridRow);
    ShinbrotReplayHit(sourceNumber) = methodResults.ReplayHit(shinbrotRow);
    HybridReplayHit(sourceNumber) = methodResults.ReplayHit(hybridRow);
    ShinbrotP(sourceNumber) = methodResults.SelectedP(shinbrotRow);
    HybridP(sourceNumber) = methodResults.SelectedP(hybridRow);
    ShinbrotAbsoluteP(sourceNumber) = methodResults.AbsoluteP(shinbrotRow);
    HybridAbsoluteP(sourceNumber) = methodResults.AbsoluteP(hybridRow);
    ShinbrotHorizon(sourceNumber) = methodResults.Horizon(shinbrotRow);
    HybridHorizon(sourceNumber) = methodResults.Horizon(hybridRow);
    ShinbrotFinalTargetError(sourceNumber) = ...
        methodResults.FinalTargetError(shinbrotRow);
    HybridFinalTargetError(sourceNumber) = ...
        methodResults.FinalTargetError(hybridRow);
    ShinbrotReplayTargetError(sourceNumber) = ...
        methodResults.ReplayTargetError(shinbrotRow);
    HybridReplayTargetError(sourceNumber) = ...
        methodResults.ReplayTargetError(hybridRow);
    SameHorizon(sourceNumber) = ...
        isfinite(ShinbrotHorizon(sourceNumber)) && ...
        ShinbrotHorizon(sourceNumber) == HybridHorizon(sourceNumber);
    BothReplayHit(sourceNumber) = ...
        ShinbrotReplayHit(sourceNumber) && HybridReplayHit(sourceNumber);

    if BothReplayHit(sourceNumber)
        AbsolutePDifference(sourceNumber) = ...
            HybridAbsoluteP(sourceNumber) - ...
            ShinbrotAbsoluteP(sourceNumber);
        if AbsolutePDifference(sourceNumber) < 0
            SmallerMagnitudeMethod{sourceNumber} = 'hybrid';
        elseif AbsolutePDifference(sourceNumber) > 0
            SmallerMagnitudeMethod{sourceNumber} = 'Shinbrot';
        else
            SmallerMagnitudeMethod{sourceNumber} = 'equal';
        end
    end
end

paired = table(SourceIndex, SourceX, ShinbrotHit, HybridHit, ...
    ShinbrotReplayHit, HybridReplayHit, ShinbrotP, HybridP, ...
    ShinbrotAbsoluteP, HybridAbsoluteP, AbsolutePDifference, ...
    ShinbrotHorizon, HybridHorizon, ShinbrotFinalTargetError, ...
    HybridFinalTargetError, ShinbrotReplayTargetError, ...
    HybridReplayTargetError, SameHorizon, BothReplayHit, ...
    SmallerMagnitudeMethod);
end

function summary = summarize_paired_results(paired)
%SUMMARIZE_PAIRED_RESULTS Summarize verified matched-horizon comparisons.

matched = paired.BothReplayHit & paired.SameHorizon;
different = paired.BothReplayHit & ~paired.SameHorizon;
difference = paired.AbsolutePDifference(matched);
shinbrotMagnitude = paired.ShinbrotAbsoluteP(matched);
hybridMagnitude = paired.HybridAbsoluteP(matched);

NumberSourcesTested = height(paired);
NumberShinbrotReplayHits = sum(paired.ShinbrotReplayHit);
NumberHybridReplayHits = sum(paired.HybridReplayHit);
NumberBothReplayHit = sum(paired.BothReplayHit);
NumberSameSuccessfulHorizon = sum(matched);
NumberDifferentSuccessfulHorizons = sum(different);
MatchedHybridSmaller = sum(difference < 0);
MatchedShinbrotSmaller = sum(difference > 0);
MatchedTies = sum(difference == 0);

shinbrotReplayErrors = paired.ShinbrotReplayTargetError( ...
    paired.ShinbrotReplayHit & isfinite(paired.ShinbrotReplayTargetError));
hybridReplayErrors = paired.HybridReplayTargetError( ...
    paired.HybridReplayHit & isfinite(paired.HybridReplayTargetError));
[MeanShinbrotReplayTargetError, MedianShinbrotReplayTargetError, ...
    MaximumShinbrotReplayTargetError] = summarize_replay_errors( ...
    shinbrotReplayErrors);
[MeanHybridReplayTargetError, MedianHybridReplayTargetError, ...
    MaximumHybridReplayTargetError] = summarize_replay_errors( ...
    hybridReplayErrors);

HybridSmallerSourceIndices = {paired.SourceIndex( ...
    matched & paired.AbsolutePDifference < 0).'};
ShinbrotSmallerSourceIndices = {paired.SourceIndex( ...
    matched & paired.AbsolutePDifference > 0).'};
TieSourceIndices = {paired.SourceIndex( ...
    matched & paired.AbsolutePDifference == 0).'};

if isempty(difference)
    MeanShinbrotAbsoluteP = NaN;
    MedianShinbrotAbsoluteP = NaN;
    MeanHybridAbsoluteP = NaN;
    MedianHybridAbsoluteP = NaN;
    MeanAbsolutePDifference = NaN;
    MedianAbsolutePDifference = NaN;
    TotalShinbrotAbsoluteP = NaN;
    TotalHybridAbsoluteP = NaN;
    NetAbsolutePDifference = NaN;
    AbsoluteNetDifference = NaN;
    SumAbsolutePairedDifference = NaN;
    MinimumShinbrotAbsoluteP = NaN;
    MinimumShinbrotSourceIndex = NaN;
    MaximumShinbrotAbsoluteP = NaN;
    MaximumShinbrotSourceIndex = NaN;
    MinimumHybridAbsoluteP = NaN;
    MinimumHybridSourceIndex = NaN;
    MaximumHybridAbsoluteP = NaN;
    MaximumHybridSourceIndex = NaN;
else
    matchedSourceIndices = paired.SourceIndex(matched);
    MeanShinbrotAbsoluteP = mean(shinbrotMagnitude);
    MedianShinbrotAbsoluteP = median(shinbrotMagnitude);
    MeanHybridAbsoluteP = mean(hybridMagnitude);
    MedianHybridAbsoluteP = median(hybridMagnitude);
    MeanAbsolutePDifference = mean(difference);
    MedianAbsolutePDifference = median(difference);
    TotalShinbrotAbsoluteP = sum(shinbrotMagnitude);
    TotalHybridAbsoluteP = sum(hybridMagnitude);
    NetAbsolutePDifference = sum(difference);
    AbsoluteNetDifference = abs(NetAbsolutePDifference);
    SumAbsolutePairedDifference = sum(abs(difference));
    [MinimumShinbrotAbsoluteP, location] = min(shinbrotMagnitude);
    MinimumShinbrotSourceIndex = matchedSourceIndices(location);
    [MaximumShinbrotAbsoluteP, location] = max(shinbrotMagnitude);
    MaximumShinbrotSourceIndex = matchedSourceIndices(location);
    [MinimumHybridAbsoluteP, location] = min(hybridMagnitude);
    MinimumHybridSourceIndex = matchedSourceIndices(location);
    [MaximumHybridAbsoluteP, location] = max(hybridMagnitude);
    MaximumHybridSourceIndex = matchedSourceIndices(location);
end

ShinbrotSuccessRate = NumberShinbrotReplayHits / NumberSourcesTested;
HybridSuccessRate = NumberHybridReplayHits / NumberSourcesTested;

summary = table(NumberSourcesTested, NumberShinbrotReplayHits, ...
    NumberHybridReplayHits, ShinbrotSuccessRate, HybridSuccessRate, ...
    NumberBothReplayHit, ...
    NumberSameSuccessfulHorizon, NumberDifferentSuccessfulHorizons, ...
    MatchedHybridSmaller, MatchedShinbrotSmaller, MatchedTies, ...
    HybridSmallerSourceIndices, ShinbrotSmallerSourceIndices, ...
    TieSourceIndices, MeanShinbrotAbsoluteP, ...
    MedianShinbrotAbsoluteP, MeanHybridAbsoluteP, ...
    MedianHybridAbsoluteP, MeanAbsolutePDifference, ...
    MedianAbsolutePDifference, TotalShinbrotAbsoluteP, ...
    TotalHybridAbsoluteP, NetAbsolutePDifference, ...
    AbsoluteNetDifference, SumAbsolutePairedDifference, ...
    MinimumShinbrotAbsoluteP, MinimumShinbrotSourceIndex, ...
    MaximumShinbrotAbsoluteP, MaximumShinbrotSourceIndex, ...
    MinimumHybridAbsoluteP, MinimumHybridSourceIndex, ...
    MaximumHybridAbsoluteP, MaximumHybridSourceIndex, ...
    MeanShinbrotReplayTargetError, MedianShinbrotReplayTargetError, ...
    MaximumShinbrotReplayTargetError, MeanHybridReplayTargetError, ...
    MedianHybridReplayTargetError, MaximumHybridReplayTargetError);
end

function methodSummary = summarize_methods(methodResults)
%SUMMARIZE_METHODS Aggregate the sixty-source two-method benchmark.

methods = {verified_bisection_identifier(); 'hybridGridBisection'};
numMethods = numel(methods);
Method = methods;
NumberSources = zeros(numMethods, 1);
NumberReplayHits = zeros(numMethods, 1);
SuccessRate = zeros(numMethods, 1);
MeanHorizon = NaN(numMethods, 1);
MedianHorizon = NaN(numMethods, 1);
MeanAbsoluteP = NaN(numMethods, 1);
MedianAbsoluteP = NaN(numMethods, 1);
MeanReplayTargetError = NaN(numMethods, 1);
MedianReplayTargetError = NaN(numMethods, 1);
MaximumReplayTargetError = NaN(numMethods, 1);
MeanSearchTime = NaN(numMethods, 1);
MedianSearchTime = NaN(numMethods, 1);

for methodNumber = 1:numMethods
    rows = strcmp(methodResults.Method, methods{methodNumber});
    verified = rows & methodResults.ReplayHit;
    NumberSources(methodNumber) = sum(rows);
    NumberReplayHits(methodNumber) = sum(verified);
    SuccessRate(methodNumber) = ...
        NumberReplayHits(methodNumber) / NumberSources(methodNumber);
    if any(verified)
        MeanHorizon(methodNumber) = mean(methodResults.Horizon(verified));
        MedianHorizon(methodNumber) = median(methodResults.Horizon(verified));
        MeanAbsoluteP(methodNumber) = mean(methodResults.AbsoluteP(verified));
        MedianAbsoluteP(methodNumber) = median(methodResults.AbsoluteP(verified));
        MeanReplayTargetError(methodNumber) = ...
            mean(methodResults.ReplayTargetError(verified));
        MedianReplayTargetError(methodNumber) = ...
            median(methodResults.ReplayTargetError(verified));
        MaximumReplayTargetError(methodNumber) = ...
            max(methodResults.ReplayTargetError(verified));
    end
    MeanSearchTime(methodNumber) = mean(methodResults.SearchTime(rows));
    MedianSearchTime(methodNumber) = median(methodResults.SearchTime(rows));
end

methodSummary = table(Method, NumberSources, NumberReplayHits, SuccessRate, ...
    MeanHorizon, MedianHorizon, MeanAbsoluteP, MedianAbsoluteP, ...
    MeanReplayTargetError, MedianReplayTargetError, ...
    MaximumReplayTargetError, MeanSearchTime, MedianSearchTime);
end

function [meanError, medianError, maximumError] = ...
    summarize_replay_errors(replayErrors)
%SUMMARIZE_REPLAY_ERRORS Summarize finite errors from verified replays.

if isempty(replayErrors)
    meanError = NaN;
    medianError = NaN;
    maximumError = NaN;
else
    meanError = mean(replayErrors);
    medianError = median(replayErrors);
    maximumError = max(replayErrors);
end
end

function save_ensemble_results(ensemble)
%SAVE_ENSEMBLE_RESULTS Save only ensemble-specific result files.

sourceIndices = ensemble.sourceIndices;
methodResultsTable = ensemble.methodResultsTable;
pairedResultsTable = ensemble.pairedResultsTable;
summary = ensemble.summary;
methodSummary = ensemble.methodSummary;
studySettings = ensemble.studySettings;
target = ensemble.target;
control = ensemble.control;

save(studySettings.outputMatFile, 'sourceIndices', ...
    'methodResultsTable', 'pairedResultsTable', 'summary', ...
    'methodSummary', ...
    'studySettings', 'target', 'control');
writetable(methodResultsTable, studySettings.methodResultsCsv);
writetable(pairedResultsTable, studySettings.pairedResultsCsv);
end