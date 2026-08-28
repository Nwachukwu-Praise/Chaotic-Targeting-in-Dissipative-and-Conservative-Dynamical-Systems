function [search, bestApprox, stats] = ...
    run_shinbrot_paper_bisection( ...
    search, bestApprox, stats, sourceState, target, params, control, ...
    pLowFull, pHighFull, numHorizons)
%RUN_SHINBROT_PAPER_BISECTION Paper-faithful Shinbrot interval targeting.
%
% This is the deliberately separate implementation of the bisection method
% described by Shinbrot, Ott, Grebogi and Yorke, Phys. Rev. A 45, 4165
% (1992).  It uses endpoint branch checks, discontinuity isolation by
% midpoint range reduction, and a fixed count of twenty-four final bracket
% refinements.  Later project extensions remain under other method names.

if ~isfield(control, 'bisectionIterations')
    control.bisectionIterations = 24;
end
if control.bisectionIterations ~= 24
    error('shinbrotPaperBisection requires control.bisectionIterations = 24.');
end
if ~isfield(control, 'maxDiscontinuityIsolationDepth')
    control.maxDiscontinuityIsolationDepth = 24;
end
if ~isfield(control, 'maxDiscontinuityIntervals')
    control.maxDiscontinuityIntervals = 256;
end
if ~isfield(control, 'effectiveExactTargetError')
    control.effectiveExactTargetError = ...
        max(1e-12, 100 * eps(max(1, abs(target.x))));
end

search.crossingOrderTestUsed = true;
search.paperAlgorithmPath = true;
search.sourceRestartPolicy = ...
    'Every parameter evaluation restarts from the original full 3-D source state.';
search.parameterApplication = ...
    'Candidate p is held constant throughout each finite-horizon propagation.';
search.discontinuityTieBreaking = ...
    ['Project-specific stack and closer-endpoint ordering is used only ', ...
    'when the paper does not specify which discontinuity half to retain.'];
search.paperRefinementRequirement = control.bisectionIterations;
search.firstToleranceHitIteration = NaN;
search.endpointToleranceHit = false;
search.endpointToleranceHitParameter = NaN;
search.endpointToleranceHitError = NaN;
search.refinementsCompleted = 0;
search.finalParameterResolution = NaN;
search.finalBracketBisectionIterations = 0;
search.paperReplayVerified = false;
search.paperReplayTargetX = NaN;
search.paperReplayTargetError = Inf;
search.preReplaySelectedP = NaN;
search.preReplayTargetX = NaN;
search.preReplayTargetError = Inf;
search.horizonsExamined = zeros(0, 1);
search.evaluationSourceState = sourceState(:);
search.discontinuitiesDetected = 0;
search.discontinuityReductions = 0;
search.continuousIntervalsExamined = 0;
search.continuousTargetBracketFound = false;
search.discontinuityIsolationDepthReached = false;
search.crossingSignatureEvaluations = 0;
search.continuousIntervalsDiscarded = 0;
search.intervalsExaminedByHorizon = zeros(numHorizons, 1);
search.evaluatedParametersByHorizon = cell(numHorizons, 1);
search.evaluatedImagesByHorizon = cell(numHorizons, 1);
search.imageIntervalByHorizon = NaN(numHorizons, 2);
search.searchResolutionSettings = struct( ...
    'paperRefinements', control.bisectionIterations, ...
    'maxDiscontinuityIsolationDepth', ...
    control.maxDiscontinuityIsolationDepth, ...
    'maxDiscontinuityIntervals', control.maxDiscontinuityIntervals, ...
    'effectiveExactTargetError', control.effectiveExactTargetError);
search.completenessNote = ...
    ['Same-branch intervals that do not bracket the target are discarded. ', ...
    'Unsuccessful horizons therefore do not certify absence of a solution.'];

horizon = 0;
cacheP = zeros(1, 0);
cacheRecords = empty_record_array();
intervalsExamined = 0;

for horizonIndex = 1:numHorizons
    horizon = horizonIndex;
    search.horizonsExamined(end + 1, 1) = horizon; %#ok<AGROW>
    cacheP = zeros(1, 0);
    cacheRecords = empty_record_array();
    intervalsExamined = 0;

    left = evaluate_at(pLowFull);
    right = evaluate_at(pHighFull);

    record_endpoint_tolerance_hit(left);
    record_endpoint_tolerance_hit(right);

    if is_effectively_exact(left)
        [verified, replayRecord] = replay_record(left.p);
        record_horizon(horizon);
        if verified
            store_success(replayRecord, horizon, left.p, right.p, ...
                0, 0, abs(right.p - left.p), left);
            return;
        end
    end
    if is_effectively_exact(right)
        [verified, replayRecord] = replay_record(right.p);
        record_horizon(horizon);
        if verified
            store_success(replayRecord, horizon, left.p, right.p, ...
                0, 0, abs(right.p - left.p), right);
            return;
        end
    end

    [bracketFound, bracketLeft, bracketRight] = ...
        isolate_same_branch_target_bracket(left, right);
    if bracketFound
        search.continuousTargetBracketFound = true;
        [hitFound, selected, finalLeft, finalRight, ...
            finalIterations, firstHitIteration, finalResolution, ...
            preReplay] = refine_target_bracket(bracketLeft, bracketRight);
        record_horizon(horizon);
        if hitFound
            store_success(selected, horizon, finalLeft.p, finalRight.p, ...
                finalIterations, firstHitIteration, finalResolution, ...
                preReplay);
            return;
        end
    else
        record_horizon(horizon);
    end
end

% ========================= nested helpers ===============================

    function [bracketFound, bracketLeft, bracketRight] = ...
            isolate_same_branch_target_bracket(initialLeft, initialRight)
        bracketFound = false;
        bracketLeft = initialLeft;
        bracketRight = initialRight;
        stack = make_item(initialLeft, initialRight, 0);

        while ~isempty(stack)
            if intervalsExamined >= control.maxDiscontinuityIntervals
                search.discontinuityIsolationDepthReached = true;
                return;
            end

            [current, stack] = pop_best_item(stack);
            intervalsExamined = intervalsExamined + 1;
            intervalLeft = current.left;
            intervalRight = current.right;

            if ~intervalLeft.success || ~intervalRight.success
                continue;
            end

            if same_branch(intervalLeft, intervalRight)
                search.continuousIntervalsExamined = ...
                    search.continuousIntervalsExamined + 1;
                if contains_target(intervalLeft, intervalRight)
                    bracketFound = true;
                    bracketLeft = intervalLeft;
                    bracketRight = intervalRight;
                    return;
                end
                search.continuousIntervalsDiscarded = ...
                    search.continuousIntervalsDiscarded + 1;
                continue;
            end

            search.discontinuitiesDetected = ...
                search.discontinuitiesDetected + 1;
            if current.depth >= control.maxDiscontinuityIsolationDepth
                search.discontinuityIsolationDepthReached = true;
                continue;
            end

            [midpoint, degenerate] = midpoint_of(intervalLeft, intervalRight);
            if degenerate
                search.discontinuityIsolationDepthReached = true;
                continue;
            end
            stats.bisectionIterations = stats.bisectionIterations + 1;
            search.discontinuityReductions = ...
                search.discontinuityReductions + 1;

            stack(end + 1) = make_item( ...
                intervalLeft, midpoint, current.depth + 1); %#ok<AGROW>
            stack(end + 1) = make_item( ...
                midpoint, intervalRight, current.depth + 1); %#ok<AGROW>
        end
    end

    function [found, selected, intervalLeft, intervalRight, ...
            refinementsCompleted, firstHitIteration, finalResolution, ...
            preReplay] = refine_target_bracket(intervalLeft, intervalRight)
        found = false;
        selected = intervalLeft;
        firstHitIteration = NaN;
        refinementsCompleted = 0;
        exactRootReached = false;

        for iteration = 1:control.bisectionIterations
            [midpoint, degenerate] = midpoint_of(intervalLeft, intervalRight);
            if degenerate
                break;
            end
            stats.bisectionIterations = stats.bisectionIterations + 1;
            refinementsCompleted = iteration;

            midpointError = abs(midpoint.xFinal - target.x);
            if is_hit(midpoint) && isnan(firstHitIteration)
                firstHitIteration = iteration;
            end
            if midpoint.success && ...
                    midpointError <= control.effectiveExactTargetError
                intervalLeft = midpoint;
                intervalRight = midpoint;
                exactRootReached = true;
                break;
            end

            leftValid = same_branch(intervalLeft, midpoint) && ...
                contains_target(intervalLeft, midpoint);
            rightValid = same_branch(midpoint, intervalRight) && ...
                contains_target(midpoint, intervalRight);

            if leftValid && rightValid
                leftError = min(abs([intervalLeft.xFinal, ...
                    midpoint.xFinal] - target.x));
                rightError = min(abs([midpoint.xFinal, ...
                    intervalRight.xFinal] - target.x));
                if leftError <= rightError
                    intervalRight = midpoint;
                else
                    intervalLeft = midpoint;
                end
            elseif leftValid
                intervalRight = midpoint;
            elseif rightValid
                intervalLeft = midpoint;
            else
                [recovered, newLeft, newRight] = ...
                    isolate_same_branch_target_bracket(intervalLeft, ...
                    intervalRight);
                if ~recovered
                    break;
                end
                intervalLeft = newLeft;
                intervalRight = newRight;
            end
        end

        finalResolution = abs(intervalRight.p - intervalLeft.p);
        if exactRootReached
            preReplay = intervalLeft;
        else
            selectedP = 0.5 * (intervalLeft.p + intervalRight.p);
            preReplay = evaluate_at(selectedP);
        end

        [found, selected] = replay_record(preReplay.p);
        found = found && abs(selected.xFinal - target.x) <= target.tolerance;
    end

    function [verified, record] = replay_record(pValue)
        signatureResult = evaluate_shinbrot_crossing_signature( ...
            sourceState, pValue, horizon, params);
        record = record_from_signature(pValue, signatureResult);
        stats.parameterEvaluations = stats.parameterEvaluations + 1;
        stats.crossingPropagations = stats.crossingPropagations + ...
            record.numberOfPlaneCrossings;
        stats.crossingSignatureEvaluations = ...
            stats.crossingSignatureEvaluations + 1;
        search.crossingSignatureEvaluations = ...
            search.crossingSignatureEvaluations + 1;
        update_best(record);
        verified = record.success && ...
            abs(record.xFinal - target.x) <= target.tolerance;
    end

    function [midpoint, degenerate] = midpoint_of(recordLeft, recordRight)
        pMid = 0.5 * (recordLeft.p + recordRight.p);
        degenerate = (pMid == recordLeft.p) || (pMid == recordRight.p);
        if degenerate
            midpoint = recordLeft;
            return;
        end
        midpoint = evaluate_at(pMid);
    end

    function record = evaluate_at(pValue)
        cachedIndex = find(cacheP == pValue, 1);
        if ~isempty(cachedIndex)
            record = cacheRecords(cachedIndex);
            return;
        end

        signatureResult = evaluate_shinbrot_crossing_signature( ...
            sourceState, pValue, horizon, params);
        record = record_from_signature(pValue, signatureResult);

        cacheP(end + 1) = pValue;
        cacheRecords(end + 1) = record;

        stats.parameterEvaluations = stats.parameterEvaluations + 1;
        stats.crossingPropagations = stats.crossingPropagations + ...
            record.numberOfPlaneCrossings;
        stats.crossingSignatureEvaluations = ...
            stats.crossingSignatureEvaluations + 1;
        search.crossingSignatureEvaluations = ...
            search.crossingSignatureEvaluations + 1;
        update_best(record);
    end

    function record = record_from_signature(pValue, signatureResult)
        record.p = pValue;
        record.xFinal = signatureResult.xFinal;
        record.success = signatureResult.success;
        record.crossingSignature = signatureResult.crossingSignature;
        record.numberOfPlaneCrossings = ...
            signatureResult.numberOfPlaneCrossings;
        record.numberOfAcceptedCrossings = ...
            signatureResult.numberOfAcceptedCrossings;
    end

    function update_best(record)
        if record.success && isfinite(record.xFinal)
            approximationError = abs(record.xFinal - target.x);
            if approximationError < bestApprox.error
                bestApprox.p = record.p;
                bestApprox.horizon = horizon;
                bestApprox.x = record.xFinal;
                bestApprox.error = approximationError;
            end
        end
    end

    function tf = same_branch(first, second)
        tf = first.success && second.success && ...
            isequal(first.crossingSignature, second.crossingSignature);
    end

    function tf = contains_target(first, second)
        tf = first.success && second.success && ...
            min(first.xFinal, second.xFinal) <= target.x && ...
            target.x <= max(first.xFinal, second.xFinal);
    end

    function tf = is_hit(record)
        tf = record.success && ...
            abs(record.xFinal - target.x) <= target.tolerance;
    end

    function tf = is_effectively_exact(record)
        tf = record.success && ...
            abs(record.xFinal - target.x) <= ...
            control.effectiveExactTargetError;
    end

    function record_endpoint_tolerance_hit(record)
        if is_hit(record) && ~is_effectively_exact(record) && ...
                ~search.endpointToleranceHit
            search.endpointToleranceHit = true;
            search.endpointToleranceHitParameter = record.p;
            search.endpointToleranceHitError = ...
                abs(record.xFinal - target.x);
        end
    end

    function item = make_item(recordLeft, recordRight, depth)
        item.left = recordLeft;
        item.right = recordRight;
        item.depth = depth;
        item.priority = interval_priority(recordLeft, recordRight, depth);
    end

    function value = interval_priority(recordLeft, recordRight, depth)
        endpointErrors = [abs(recordLeft.xFinal - target.x), ...
            abs(recordRight.xFinal - target.x)];
        endpointErrors(~isfinite(endpointErrors)) = Inf;
        value = min(endpointErrors) + 1e-6 * depth;
    end

    function [item, stack] = pop_best_item(stack)
        [~, index] = min([stack.priority]);
        item = stack(index);
        stack(index) = [];
    end

    function record_horizon(horizonToRecord)
        search.intervalsExaminedByHorizon(horizonToRecord) = ...
            intervalsExamined;
        if isempty(cacheRecords)
            return;
        end
        pValues = [cacheRecords.p];
        xValues = [cacheRecords.xFinal];
        successFlags = logical([cacheRecords.success]);
        keep = successFlags & isfinite(xValues);
        search.evaluatedParametersByHorizon{horizonToRecord} = pValues(keep);
        search.evaluatedImagesByHorizon{horizonToRecord} = xValues(keep);
        if any(keep)
            search.imageIntervalByHorizon(horizonToRecord, :) = ...
                [min(xValues(keep)), max(xValues(keep))];
        end
    end

    function store_success(record, selectedHorizon, bracketLow, ...
            bracketHigh, finalIterations, firstHitIteration, ...
            finalResolution, preReplay)
        finalError = abs(record.xFinal - target.x);
        search.found = finalError <= target.tolerance;
        search.selectedP = record.p;
        search.horizon = selectedHorizon;
        search.targetX = record.xFinal;
        search.bestDistance = finalError;
        search.finalTargetError = finalError;
        search.selectionMethod = ...
            'Shinbrot paper bisection with crossing-order discontinuity reduction';
        search.bracketLow = bracketLow;
        search.bracketHigh = bracketHigh;
        search.finalBracketBisectionIterations = finalIterations;
        search.bisectionIterationsUsed = finalIterations;
        search.refinementsCompleted = finalIterations;
        search.firstToleranceHitIteration = firstHitIteration;
        search.finalParameterResolution = finalResolution;
        search.paperReplayVerified = search.found;
        search.paperReplayTargetX = record.xFinal;
        search.paperReplayTargetError = finalError;
        search.preReplaySelectedP = preReplay.p;
        search.preReplayTargetX = preReplay.xFinal;
        search.preReplayTargetError = abs(preReplay.xFinal - target.x);
        search.plotHorizon = selectedHorizon;
        search.failureCode = 'none';
        search.failureReason = '';
    end
end

function records = empty_record_array()
records = struct('p', {}, 'xFinal', {}, 'success', {}, ...
    'crossingSignature', {}, 'numberOfPlaneCrossings', {}, ...
    'numberOfAcceptedCrossings', {});
end
