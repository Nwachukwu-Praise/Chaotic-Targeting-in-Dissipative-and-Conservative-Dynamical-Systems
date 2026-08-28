function trial = record_noisy_search_event( ...
    trial, eventKind, crossingIndex, stepIndex, eventTime, search, accepted)
%RECORD_NOISY_SEARCH_EVENT Store initial and retargeting search provenance.
%
% eventKind is 'initial' or 'retarget'.  The initial search contributes to
% numSearchEvaluations but is not a retargeting event.

trial.numSearchEvaluations = trial.numSearchEvaluations + 1;
trial.numCompleteBisectionSearches = trial.numCompleteBisectionSearches + 1;

switch eventKind
    case 'initial'
        trial.initialSearch = compact_search_record( ...
            eventKind, crossingIndex, stepIndex, eventTime, search, accepted);
        trial.initialSearchPerformed = true;
        trial.initialSearchMethod = search.method;
        trial.initialSearchSelectionMethod = search.selectionMethod;
        trial.initialSearchGridEvaluations = search.gridEvaluations;
        trial.initialSearchCrossingOrderTestUsed = search.crossingOrderTestUsed;
        trial.initialSearchFound = search.found;
        trial.initialSearchSelectedP = search.selectedP;
        trial.initialSearchHorizon = search.horizon;
        trial.initialSearchFinalTargetError = search.finalTargetError;

    case 'retarget'
        trial.numRetargetingEvents = trial.numRetargetingEvents + 1;
        trial.retargetingCrossingIndices(end + 1) = crossingIndex;
        trial.retargetingStepIndices(end + 1) = stepIndex;
        trial.retargetingTimes(end + 1) = eventTime;
        trial.retargetingSearchMethods{end + 1} = search.method;
        trial.retargetingSelectedP(end + 1) = search.selectedP;
        trial.retargetingSearchFound(end + 1) = search.found;
        trial.retargetingSelectionMethods{end + 1} = search.selectionMethod;
        trial.retargetingGridEvaluations(end + 1) = search.gridEvaluations;
        trial.retargetingCrossingOrderTestUsed(end + 1) = ...
            search.crossingOrderTestUsed;
        if ~accepted
            trial.numRetargetingFailures = trial.numRetargetingFailures + 1;
        end

    otherwise
        error('Unknown noisy search event kind: %s', eventKind);
end
end

function record = compact_search_record( ...
    eventKind, crossingIndex, stepIndex, eventTime, search, accepted)
record.eventKind = eventKind;
record.crossingIndex = crossingIndex;
record.stepIndex = stepIndex;
record.eventTime = eventTime;
record.accepted = accepted;
record.method = search.method;
record.selectionMethod = search.selectionMethod;
record.found = search.found;
record.selectedP = search.selectedP;
record.horizon = search.horizon;
record.finalTargetError = search.finalTargetError;
record.parameterEvaluations = search.parameterEvaluations;
record.gridEvaluations = search.gridEvaluations;
record.crossingOrderTestUsed = search.crossingOrderTestUsed;
end
