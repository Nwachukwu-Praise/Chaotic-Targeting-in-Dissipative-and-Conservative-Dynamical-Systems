function [search, bestApprox, stats] = ...
    run_shinbrot_coverage_bisection( ...
    search, bestApprox, stats, sourceState, target, params, control, ...
    pLowFull, pHighFull, numHorizons)
%RUN_SHINBROT_COVERAGE_BISECTION Shinbrot's criterion by Shinbrot's procedure.
%
% This is the 9 August search with everything that has no counterpart in the
% paper removed, and with the paper's own horizon test put back in front.
%
% WHAT IS KEPT FROM THE 9 AUGUST VERSION
%
%   * Ordered L/R crossing signatures to separate branches, following
%     Shinbrot et al. footnote [11].
%   * Subdivision of a same-signature interval whose endpoint images do not
%     bracket X_t.  Endpoint bracketing is only a sufficient test for the
%     paper's condition X_t in DeltaX_n, and the two agree only where
%     X_n(p) is monotone, so discarding such intervals can report failure at
%     a horizon where a solution exists.  This was the correct fix.
%   * The declared resolution of Subsection 3.1.5.  Every continuous feature
%     of X_n(p) wider than (pHighFull - pLowFull)/2^minCertifiedDepth is
%     sampled at least once, within any horizon that is actually searched.
%   * Proximity ordering of the remaining intervals below that width, and
%     the declared interval and depth budgets.
%
% WHAT IS REMOVED
%
%   * Estimated-Lipschitz pruning.  An empirical maximum secant slope is not
%     a certified bound on an unsampled interior slope, safety factor or
%     not, so it could exclude an interval containing a solution.  It also
%     appears nowhere in the paper and nowhere in Subsection 3.1.5.
%   * maxBracketRefinements, which silently raised the declared
%     control.bisectionIterations to at least 60.  Bracket refinement is
%     again limited to control.bisectionIterations.
%   * The conflation of exploratory subdivision with bracket bisection.
%     stats.bisectionIterations now counts only halvings performed inside a
%     target bracket.  Exploratory splits are reported separately in
%     search.exploratorySubdivisions.
%
% WHAT IS ADDED, AND WHY IT IS THE FAITHFUL PART
%
%   The paper's targeting condition is set coverage: the accessible image
%   set of Equation (5),
%
%       R_n(X_s) = { F^n(X_s; p_0 + dp) : |dp| <= dp_max },
%
%   grows with n until it encompasses X_t, after which some admissible p
%   must reach the target and only the estimate of that value remains to be
%   refined.  A horizon at which the image set does not reach X_t therefore
%   requires no search at all.
%
%   Each horizon now begins with a dyadic probe of the admissible interval,
%   from which the observed span of X_n is recorded as a lower bound on
%   DeltaX_n.  If X_t lies outside that span the horizon is rejected on the
%   paper's own criterion and the search advances, instead of exhausting the
%   interval budget to prove a negative.  This is where the 9 August version
%   spent most of its time: horizons with no solution were searched to
%   exhaustion.
%
%   The probe is a lower bound, not the true DeltaX_n, so the rejection is
%   conservative only to the extent that the probe resolves the image.  The
%   probe depth is declared in control.imageProbeDepth and reported in
%   search.imageProbeSamples, and every rejected horizon is flagged in
%   search.horizonRejectedByImageSpan so that a rejection is never confused
%   with an exhausted search.
%
% Field names and the calling signature are unchanged, so
% search_parameter_to_target, assert_verified_bisection_search and
% test_corrected_targeting_methods continue to work.  gridEvaluations is
% never incremented and no grid metadata is written.

% ---- optional settings, defaulted so existing callers are unaffected -----
if ~isfield(control, 'maxDiscontinuityIsolationDepth')
    control.maxDiscontinuityIsolationDepth = 24;
end
if ~isfield(control, 'maxDiscontinuityIntervals')
    control.maxDiscontinuityIntervals = 256;
end
if ~isfield(control, 'minCertifiedDepth')
    control.minCertifiedDepth = 8;
end
if ~isfield(control, 'maxContinuousSubdivisionDepth')
    control.maxContinuousSubdivisionDepth = 32;
end
if ~isfield(control, 'maxIntervalsExamined')
    control.maxIntervalsExamined = ...
        max(4096, 16 * control.maxDiscontinuityIntervals);
end
if ~isfield(control, 'minParameterInterval')
    control.minParameterInterval = ...
        abs(pHighFull - pLowFull) * 2^(-control.maxContinuousSubdivisionDepth);
end
if ~isfield(control, 'imageProbeDepth')
    % 2^5 + 1 = 33 dyadic samples.  Chosen to resolve the image span at a
    % cost far below a full sweep; it is a project-specific choice and the
    % paper prescribes none.
    control.imageProbeDepth = 5;
end
if ~isfield(control, 'imageSpanMargin')
    % A rejection requires X_t to lie outside the sampled span by more than
    % the target tolerance, so a marginal horizon is searched rather than
    % skipped.
    control.imageSpanMargin = target.tolerance;
end
if ~isfield(control, 'selectSmallestParameter')
    control.selectSmallestParameter = false;
end

% ---- diagnostics retained from previous versions ------------------------
search.crossingOrderTestUsed = true;
search.discontinuitiesDetected = 0;
search.discontinuityReductions = 0;
search.continuousIntervalsExamined = 0;
search.continuousTargetBracketFound = false;
search.discontinuityIsolationDepthReached = false;
search.crossingSignatureEvaluations = 0;
stats.crossingSignatureEvaluations = 0;

% ---- resolution and completeness diagnostics ----------------------------
search.certifiedParameterResolution = ...
    abs(pHighFull - pLowFull) * 2^(-control.minCertifiedDepth);
search.verifiedCandidatesThisHorizon = 0;
search.continuousIntervalsSubdivided = 0;
search.continuousIntervalsPrunedByBound = 0;   % retained at zero: no pruning
search.continuousIntervalsAbandonedAtResolution = 0;
search.parameterResolutionFloorReached = false;
search.intervalBudgetExhausted = false;
search.intervalsExaminedByHorizon = zeros(numHorizons, 1);
search.lipschitzEstimateByHorizon = NaN(numHorizons, 1);
search.imageIntervalByHorizon = NaN(numHorizons, 2);
search.evaluatedParametersByHorizon = cell(numHorizons, 1);
search.evaluatedImagesByHorizon = cell(numHorizons, 1);

% ---- separated iteration counters ---------------------------------------
search.bracketBisectionIterations = 0;
search.exploratorySubdivisions = 0;

% ---- coverage-gate diagnostics ------------------------------------------
search.imageProbeSamples = 2^control.imageProbeDepth + 1;
search.horizonRejectedByImageSpan = false(numHorizons, 1);
search.horizonSearched = false(numHorizons, 1);
search.imageSpanDistanceByHorizon = NaN(numHorizons, 1);

search.searchResolutionSettings = struct( ...
    'minCertifiedDepth', control.minCertifiedDepth, ...
    'maxContinuousSubdivisionDepth', control.maxContinuousSubdivisionDepth, ...
    'maxIntervalsExamined', control.maxIntervalsExamined, ...
    'minParameterInterval', control.minParameterInterval, ...
    'bracketRefinementLimit', control.bisectionIterations, ...
    'imageProbeDepth', control.imageProbeDepth, ...
    'imageSpanMargin', control.imageSpanMargin, ...
    'selectSmallestParameter', control.selectSmallestParameter);

search.completenessNote = sprintf( ...
    ['A horizon is searched only when the probed image span of X_n reaches ', ...
    'X_t, which is the paper''s condition X_t in DeltaX_n. Within a ', ...
    'searched horizon every continuous feature of X_n(p) wider than %.6g ', ...
    'in p is sampled at least once, and no interval is excluded by an ', ...
    'estimated bound. A horizon rejected on image span is flagged ', ...
    'separately from one exhausted at the declared resolution.'], ...
    search.certifiedParameterResolution);

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

    % ---- Equation (5): probe the accessible image set at this horizon ---
    probeCount = 2^control.imageProbeDepth + 1;
    probeP = linspace(pLowFull, pHighFull, probeCount);
    probeRecords = empty_record_array();
    directHit = [];
    for probeIndex = 1:probeCount
        record = evaluate_at(probeP(probeIndex));
        probeRecords(end + 1) = record; %#ok<AGROW>
        if is_hit(record)
            directHit = record;
            break;
        end
    end

    if ~isempty(directHit)
        record_horizon(horizon);
        search.horizonSearched(horizon) = true;
        store_success(directHit, horizon, pLowFull, pHighFull, 0);
        return;
    end

    [spanLow, spanHigh, anySuccess] = image_span(probeRecords);
    if anySuccess
        search.imageSpanDistanceByHorizon(horizon) = ...
            max([spanLow - target.x, target.x - spanHigh, 0]);
    end

    if ~anySuccess || ...
            target.x < spanLow - control.imageSpanMargin || ...
            target.x > spanHigh + control.imageSpanMargin
        % The image set does not reach the target at this horizon.  The
        % paper's condition fails, so no admissible perturbation exists
        % within the resolution of the probe and the horizon is skipped.
        record_horizon(horizon);
        search.horizonRejectedByImageSpan(horizon) = true;
        continue;
    end

    search.horizonSearched(horizon) = true;

    [hitFound, hitRecord, bracketLow, bracketHigh, refineIterations] = ...
        explore_horizon(probeP, probeRecords);
    record_horizon(horizon);

    if hitFound
        store_success(hitRecord, horizon, ...
            bracketLow, bracketHigh, refineIterations);
        return;
    end
end

% ======================= nested helper functions ========================

    function [hitFound, hitRecord, bracketLow, bracketHigh, ...
            refineIterations] = explore_horizon(probeP, probeRecords)
        % Traversal seeded from the probe subdivision.  Intervals shallower
        % than minCertifiedDepth are ordered breadth first, which delivers
        % the declared resolution; below that they are ordered by proximity
        % of their sampled images to X_t, as Subsection 3.1.5 states.  No
        % interval is excluded by an estimated bound.
        hitFound = false;
        hitRecord = probeRecords(1);
        bracketLow = probeP(1);
        bracketHigh = probeP(end);
        refineIterations = 0;
        collectMode = control.selectSmallestParameter;
        candidates = empty_candidate_array();

        queue = empty_item_array();
        for k = 1:(numel(probeRecords) - 1)
            queue(end + 1) = make_item(probeRecords(k), probeRecords(k + 1), ...
                control.imageProbeDepth); %#ok<AGROW>
        end

        while ~isempty(queue)
            if intervalsExamined >= control.maxIntervalsExamined
                search.intervalBudgetExhausted = true;
                break;
            end

            [~, bestIndex] = min([queue.priority]);
            current = queue(bestIndex);
            queue(bestIndex) = [];
            intervalsExamined = intervalsExamined + 1;

            intervalLeft = current.left;
            intervalRight = current.right;

            if is_hit(intervalLeft)
                candidates = add_candidate(candidates, intervalLeft, ...
                    intervalLeft.p, intervalRight.p, 0);
                if ~collectMode
                    hitFound = true;
                    hitRecord = intervalLeft;
                    return;
                end
            end
            if is_hit(intervalRight)
                candidates = add_candidate(candidates, intervalRight, ...
                    intervalLeft.p, intervalRight.p, 0);
                if ~collectMode
                    hitFound = true;
                    hitRecord = intervalRight;
                    return;
                end
            end

            if ~intervalLeft.success && ~intervalRight.success
                search.continuousIntervalsAbandonedAtResolution = ...
                    search.continuousIntervalsAbandonedAtResolution + 1;
                continue;
            end

            width = abs(intervalRight.p - intervalLeft.p);
            atFloor = width <= control.minParameterInterval;

            if same_branch(intervalLeft, intervalRight)
                search.continuousIntervalsExamined = ...
                    search.continuousIntervalsExamined + 1;

                if contains_target(intervalLeft, intervalRight)
                    search.continuousTargetBracketFound = true;
                    [found, selected, finalLeft, finalRight, iterations] = ...
                        refine_target_bracket(intervalLeft, intervalRight);
                    if found
                        candidates = add_candidate(candidates, selected, ...
                            finalLeft.p, finalRight.p, iterations);
                        if ~collectMode
                            hitFound = true;
                            hitRecord = selected;
                            bracketLow = finalLeft.p;
                            bracketHigh = finalRight.p;
                            refineIterations = iterations;
                            return;
                        end
                    end
                    % Inconclusive refinement, usually a fold splitting the
                    % bracket.  Fall through and subdivide.
                end

                if atFloor || ...
                        current.depth >= control.maxContinuousSubdivisionDepth
                    search.parameterResolutionFloorReached = true;
                    search.continuousIntervalsAbandonedAtResolution = ...
                        search.continuousIntervalsAbandonedAtResolution + 1;
                    continue;
                end

                [midpoint, degenerate] = midpoint_of(intervalLeft, intervalRight);
                if degenerate
                    search.parameterResolutionFloorReached = true;
                    search.continuousIntervalsAbandonedAtResolution = ...
                        search.continuousIntervalsAbandonedAtResolution + 1;
                    continue;
                end

                search.exploratorySubdivisions = ...
                    search.exploratorySubdivisions + 1;
                search.continuousIntervalsSubdivided = ...
                    search.continuousIntervalsSubdivided + 1;

                if is_hit(midpoint)
                    candidates = add_candidate(candidates, midpoint, ...
                        intervalLeft.p, intervalRight.p, 0);
                    if ~collectMode
                        hitFound = true;
                        hitRecord = midpoint;
                        return;
                    end
                end

                queue = push_halves(queue, intervalLeft, midpoint, ...
                    intervalRight, current.depth + 1);
                continue;
            end

            % Different ordered crossing signatures: footnote [11].
            search.discontinuitiesDetected = search.discontinuitiesDetected + 1;

            if atFloor || ...
                    current.depth >= control.maxDiscontinuityIsolationDepth
                search.discontinuityIsolationDepthReached = true;
                continue;
            end

            [midpoint, degenerate] = midpoint_of(intervalLeft, intervalRight);
            if degenerate
                search.discontinuityIsolationDepthReached = true;
                continue;
            end

            search.exploratorySubdivisions = search.exploratorySubdivisions + 1;
            search.discontinuityReductions = search.discontinuityReductions + 1;

            if is_hit(midpoint)
                candidates = add_candidate(candidates, midpoint, ...
                    intervalLeft.p, intervalRight.p, 0);
                if ~collectMode
                    hitFound = true;
                    hitRecord = midpoint;
                    return;
                end
            end

            queue = push_halves(queue, intervalLeft, midpoint, ...
                intervalRight, current.depth + 1);
        end

        if ~isempty(candidates)
            [~, order] = sortrows( ...
                [abs([candidates.p]).', [candidates.error].'], [1 2]);
            best = candidates(order(1));
            hitFound = true;
            hitRecord = best.record;
            bracketLow = best.bracketLow;
            bracketHigh = best.bracketHigh;
            refineIterations = best.iterations;
            search.verifiedCandidatesThisHorizon = numel(candidates);
        end
    end

    function [found, selected, intervalLeft, intervalRight, ...
            iterationsUsed] = refine_target_bracket(intervalLeft, intervalRight)
        % Bisection inside one continuous bracket, limited to the declared
        % control.bisectionIterations.  These are the only halvings counted
        % as bisection iterations.
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
            search.bracketBisectionIterations = ...
                search.bracketBisectionIterations + 1;

            if is_hit(midpoint)
                found = true;
                selected = midpoint;
                return;
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
                break;
            end
        end

        endpoints = [intervalLeft, intervalRight];
        errors = abs([endpoints.xFinal] - target.x);
        [~, bestIndex] = min(errors);
        selected = endpoints(bestIndex);
        found = is_hit(selected);
    end

    function [spanLow, spanHigh, anySuccess] = image_span(records)
        spanLow = NaN;
        spanHigh = NaN;
        anySuccess = false;
        if isempty(records)
            return;
        end
        values = [records.xFinal];
        flags = logical([records.success]);
        keep = flags & isfinite(values);
        if ~any(keep)
            return;
        end
        spanLow = min(values(keep));
        spanHigh = max(values(keep));
        anySuccess = true;
    end

    function candidates = empty_candidate_array()
        candidates = struct('record', {}, 'p', {}, 'error', {}, ...
            'bracketLow', {}, 'bracketHigh', {}, 'iterations', {});
    end

    function items = empty_item_array()
        items = struct('left', {}, 'right', {}, 'depth', {}, 'priority', {});
    end

    function candidates = add_candidate(candidates, record, ...
            bracketLow, bracketHigh, iterations)
        for index = 1:numel(candidates)
            if candidates(index).p == record.p
                return;
            end
        end
        item.record = record;
        item.p = record.p;
        item.error = abs(record.xFinal - target.x);
        item.bracketLow = bracketLow;
        item.bracketHigh = bracketHigh;
        item.iterations = iterations;
        candidates(end + 1) = item;
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

    function queue = push_halves(queue, recordLeft, midpoint, recordRight, depth)
        queue(end + 1) = make_item(recordLeft, midpoint, depth);
        queue(end + 1) = make_item(midpoint, recordRight, depth);
    end

    function item = make_item(recordLeft, recordRight, depth)
        item.left = recordLeft;
        item.right = recordRight;
        item.depth = depth;
        item.priority = interval_priority(recordLeft, recordRight, depth);
    end

    function value = interval_priority(recordLeft, recordRight, depth)
        if depth < control.minCertifiedDepth
            value = depth;
            return;
        end
        distance = interval_distance_to_target(recordLeft, recordRight);
        value = control.minCertifiedDepth + distance / (1 + distance);
    end

    function value = interval_distance_to_target(recordLeft, recordRight)
        value = Inf;
        if recordLeft.success && isfinite(recordLeft.xFinal)
            value = min(value, abs(recordLeft.xFinal - target.x));
        end
        if recordRight.success && isfinite(recordRight.xFinal)
            value = min(value, abs(recordRight.xFinal - target.x));
        end
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
            'Shinbrot interval search, image-coverage gated, no estimated pruning';
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
