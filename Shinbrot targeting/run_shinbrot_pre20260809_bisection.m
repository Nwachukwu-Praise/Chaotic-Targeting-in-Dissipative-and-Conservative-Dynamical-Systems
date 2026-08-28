function [search, bestApprox, stats] = ...
    run_shinbrot_pre20260809_bisection( ...
    search, bestApprox, stats, sourceState, target, params, control, ...
    pLowFull, pHighFull, numHorizons)
%RUN_SHINBROT_PRE20260809_BISECTION Reconstruction of the pre-9-August search.
%
% RECONSTRUCTION, NOT A RECOVERED FILE.  The original was overwritten on
% 9 August 2026 and no backup or version history survives.  This file is
% rebuilt from the behavioural description preserved in the header of
% run_shinbrot_discontinuity_aware_bisection.m:
%
%   "The previous version tested only whether X_t lay between the two
%    endpoint values X_n(p_L) and X_n(p_R), and discarded any same-signature
%    interval that failed that test without subdividing it. ... Subdivision
%    was driven entirely by the location of discontinuities, so the search
%    had no resolution parameter of its own."
%
% It is therefore the current search with the post-9-August machinery
% removed:
%
%   * no certified sweep, so no minCertifiedDepth and no declared
%     parameter resolution;
%   * no best-first Phase 2 traversal;
%   * no estimated-Lipschitz pruning;
%   * bracket refinement limited to control.bisectionIterations, not to
%     max(bisectionIterations, 60);
%   * a same-signature interval whose endpoint images do not bracket X_t is
%     DISCARDED rather than subdivided.
%
% VERIFICATION CRITERION.  This reconstruction is correct if and only if it
% reproduces the values reported in Table 2 of the report from the source at
% accepted crossing 150:
%
%       p = 0.086719,  final target error = 0.0069203,  horizon n = 3.
%
% If it does not, it is not the original and must not be used.  Check with
% verify_pre20260809_reconstruction.
%
% STATUS.  Retained for comparison only.  Subsection 3.1.5 of the report
% describes the post-9-August search, including equation (12), so reverting
% to this version operationally would falsify the method chapter.  Its
% intended use is to demonstrate the failure mode that motivated the change:
% a horizon at which a solution exists inside a continuous, non-monotone
% piece of X_n(p) is reported as unsuccessful here and found by the current
% search.
%
% The traversal order is depth-first with the left half taken first.  The
% original order is not recorded in the surviving description; because a
% successful horizon returns at the first hit, order can affect which of
% several admissible perturbations is returned.  This is the one point where
% the reconstruction may differ from the original, and it is the first thing
% to vary if the Table 2 check fails.

% ---- optional settings, matching the surviving pre-change names ---------
if ~isfield(control, 'maxDiscontinuityIsolationDepth')
    control.maxDiscontinuityIsolationDepth = 24;
end
if ~isfield(control, 'maxDiscontinuityIntervals')
    control.maxDiscontinuityIntervals = 256;
end

% ---- diagnostics that existed before the change -------------------------
search.crossingOrderTestUsed = true;
search.discontinuitiesDetected = 0;
search.discontinuityReductions = 0;
search.continuousIntervalsExamined = 0;
search.continuousTargetBracketFound = false;
search.discontinuityIsolationDepthReached = false;
search.crossingSignatureEvaluations = 0;
stats.crossingSignatureEvaluations = 0;

% ---- fields introduced after the change, held at neutral values so that
% ---- downstream code that reads them does not error ---------------------
search.certifiedParameterResolution = NaN;
search.verifiedCandidatesThisHorizon = 0;
search.continuousIntervalsSubdivided = 0;
search.continuousIntervalsPrunedByBound = 0;
search.continuousIntervalsAbandonedAtResolution = 0;
search.continuousIntervalsDiscarded = 0;
search.parameterResolutionFloorReached = false;
search.intervalBudgetExhausted = false;
search.intervalsExaminedByHorizon = zeros(numHorizons, 1);
search.lipschitzEstimateByHorizon = NaN(numHorizons, 1);
search.imageIntervalByHorizon = NaN(numHorizons, 2);
search.evaluatedParametersByHorizon = cell(numHorizons, 1);
search.evaluatedImagesByHorizon = cell(numHorizons, 1);
search.searchResolutionSettings = struct( ...
    'minCertifiedDepth', NaN, ...
    'maxContinuousSubdivisionDepth', NaN, ...
    'maxIntervalsExamined', NaN, ...
    'minParameterInterval', NaN, ...
    'lipschitzSafetyFactor', NaN, ...
    'maxBracketRefinements', control.bisectionIterations, ...
    'selectSmallestParameter', false);
search.completenessNote = ...
    ['Subdivision is driven entirely by the location of discontinuities. ', ...
    'A same-signature interval whose endpoint images do not bracket X_t ', ...
    'is discarded, so the search has no resolution parameter of its own ', ...
    'and an unsuccessful horizon carries no completeness claim.'];

% ---- per-horizon state, shared with the nested functions ----------------
horizon = 0;
cacheP = zeros(1, 0);
cacheRecords = empty_record_array();
intervalsExamined = 0;

for horizonIndex = 1:numHorizons
    horizon = horizonIndex;
    cacheP = zeros(1, 0);
    cacheRecords = empty_record_array();
    intervalsExamined = 0;

    left = evaluate_at(pLowFull);
    right = evaluate_at(pHighFull);

    if is_hit(left)
        store_success(left, horizon, left.p, right.p, 0);
        record_horizon(horizon);
        return;
    end
    if is_hit(right)
        store_success(right, horizon, left.p, right.p, 0);
        record_horizon(horizon);
        return;
    end

    [hitFound, hitRecord, bracketLow, bracketHigh, refineIterations] = ...
        explore_horizon(left, right);
    record_horizon(horizon);

    if hitFound
        store_success(hitRecord, horizon, ...
            bracketLow, bracketHigh, refineIterations);
        return;
    end
end

% ======================= nested helper functions ========================

    function [hitFound, hitRecord, bracketLow, bracketHigh, ...
            refineIterations] = explore_horizon(initialLeft, initialRight)
        % Depth-first traversal.  Subdivision happens only where the two
        % endpoints lie on different branches; a same-branch interval is
        % either refined, if its endpoint images bracket the target, or
        % discarded.  This is the behaviour the 9 August change replaced.
        hitFound = false;
        hitRecord = initialLeft;
        bracketLow = initialLeft.p;
        bracketHigh = initialRight.p;
        refineIterations = 0;

        stack = make_item(initialLeft, initialRight, 0);

        while ~isempty(stack)
            if intervalsExamined >= control.maxDiscontinuityIntervals
                search.discontinuityIsolationDepthReached = true;
                break;
            end

            current = stack(end);
            stack(end) = [];
            intervalsExamined = intervalsExamined + 1;

            intervalLeft = current.left;
            intervalRight = current.right;

            if is_hit(intervalLeft)
                hitFound = true;
                hitRecord = intervalLeft;
                bracketLow = intervalLeft.p;
                bracketHigh = intervalRight.p;
                return;
            end
            if is_hit(intervalRight)
                hitFound = true;
                hitRecord = intervalRight;
                bracketLow = intervalLeft.p;
                bracketHigh = intervalRight.p;
                return;
            end

            if ~intervalLeft.success && ~intervalRight.success
                continue;
            end

            if same_branch(intervalLeft, intervalRight)
                search.continuousIntervalsExamined = ...
                    search.continuousIntervalsExamined + 1;

                if ~contains_target(intervalLeft, intervalRight)
                    % The discarded case.  X_n(p) may still pass through
                    % X_t in the interior of this interval, which is
                    % precisely the failure mode the current search fixes.
                    search.continuousIntervalsDiscarded = ...
                        search.continuousIntervalsDiscarded + 1;
                    continue;
                end

                search.continuousTargetBracketFound = true;
                [found, selected, finalLeft, finalRight, iterations] = ...
                    refine_target_bracket(intervalLeft, intervalRight);
                if found
                    hitFound = true;
                    hitRecord = selected;
                    bracketLow = finalLeft.p;
                    bracketHigh = finalRight.p;
                    refineIterations = iterations;
                    return;
                end
                continue;
            end

            % Different ordered crossing signatures: a discontinuity of the
            % return map separates the endpoints, so the interval is reduced
            % rather than refined.  Footnote [11] of the paper.
            search.discontinuitiesDetected = search.discontinuitiesDetected + 1;

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
            search.discontinuityReductions = search.discontinuityReductions + 1;

            if is_hit(midpoint)
                hitFound = true;
                hitRecord = midpoint;
                bracketLow = intervalLeft.p;
                bracketHigh = intervalRight.p;
                return;
            end

            % Pushed right-half first so the left half is popped first.
            stack(end + 1) = make_item(midpoint, intervalRight, current.depth + 1); %#ok<AGROW>
            stack(end + 1) = make_item(intervalLeft, midpoint, current.depth + 1); %#ok<AGROW>
        end
    end

    function [found, selected, intervalLeft, intervalRight, ...
            iterationsUsed] = refine_target_bracket(intervalLeft, intervalRight)
        % Bisection inside one continuous bracket, limited to the declared
        % control.bisectionIterations.  The half whose endpoint images still
        % contain X_t is retained.
        found = false;
        selected = intervalLeft;
        iterationsUsed = 0;

        for iteration = 1:control.bisectionIterations
            iterationsUsed = iteration;
            [midpoint, degenerate] = midpoint_of(intervalLeft, intervalRight);
            if degenerate
                break;
            end
            stats.bisectionIterations = stats.bisectionIterations + 1;

            if is_hit(midpoint)
                found = true;
                selected = midpoint;
                return;
            end

            leftValid = same_branch(intervalLeft, midpoint) && ...
                contains_target(intervalLeft, midpoint);
            rightValid = same_branch(midpoint, intervalRight) && ...
                contains_target(midpoint, intervalRight);

            if leftValid
                intervalRight = midpoint;
            elseif rightValid
                intervalLeft = midpoint;
            else
                break;
            end
        end

        endpoints = [intervalLeft, intervalRight];
        errors = abs([endpoints.xFinal] - target.x);
        [~, bestIndex] = min(errors);
        selected = endpoints(bestIndex);
        found = is_hit(selected);
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

    function item = make_item(recordLeft, recordRight, depth)
        item.left = recordLeft;
        item.right = recordRight;
        item.depth = depth;
    end

    function record = evaluate_at(pValue)
        cachedIndex = find(cacheP == pValue, 1);
        if ~isempty(cachedIndex)
            record = cacheRecords(cachedIndex);
            return;
        end

        signatureResult = evaluate_shinbrot_crossing_signature( ...
            sourceState, pValue, horizon, params);

        record.p = pValue;
        record.xFinal = signatureResult.xFinal;
        record.success = signatureResult.success;
        record.crossingSignature = signatureResult.crossingSignature;
        record.numberOfPlaneCrossings = signatureResult.numberOfPlaneCrossings;
        record.numberOfAcceptedCrossings = ...
            signatureResult.numberOfAcceptedCrossings;

        cacheP(end + 1) = pValue;
        cacheRecords(end + 1) = record;

        stats.parameterEvaluations = stats.parameterEvaluations + 1;
        stats.crossingPropagations = stats.crossingPropagations + ...
            record.numberOfPlaneCrossings;
        stats.crossingSignatureEvaluations = ...
            stats.crossingSignatureEvaluations + 1;
        search.crossingSignatureEvaluations = ...
            search.crossingSignatureEvaluations + 1;

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

    function record_horizon(horizonToRecord)
        search.intervalsExaminedByHorizon(horizonToRecord) = intervalsExamined;
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

    function store_success(record, selectedHorizon, ...
            bracketLow, bracketHigh, finalIterations)
        finalError = abs(record.xFinal - target.x);
        search.found = finalError <= target.tolerance;
        search.selectedP = record.p;
        search.horizon = selectedHorizon;
        search.targetX = record.xFinal;
        search.bestDistance = finalError;
        search.finalTargetError = finalError;
        search.selectionMethod = ...
            'Shinbrot discontinuity-aware interval search, pre-2026-08-09';
        search.bracketLow = bracketLow;
        search.bracketHigh = bracketHigh;
        search.finalBracketBisectionIterations = finalIterations;
        search.bisectionIterationsUsed = finalIterations;
        search.plotHorizon = selectedHorizon;
    end
end

function records = empty_record_array()
records = struct('p', {}, 'xFinal', {}, 'success', {}, ...
    'crossingSignature', {}, 'numberOfPlaneCrossings', {}, ...
    'numberOfAcceptedCrossings', {});
end
