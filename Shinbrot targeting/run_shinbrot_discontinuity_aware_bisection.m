function [search, bestApprox, stats] = ...
    run_shinbrot_discontinuity_aware_bisection( ...
    search, bestApprox, stats, sourceState, target, params, control, ...
    pLowFull, pHighFull, numHorizons)
%RUN_SHINBROT_DISCONTINUITY_AWARE_BISECTION Crossing-order interval search.
%
% This implementation is independent of the hybrid method's grid evaluation
% and bracket-refinement helpers. It separates trajectories with different
% ordered L/R crossing signatures, following Shinbrot et al. footnote [11],
% and refines continuous target brackets by bisection, following the
% refinement described in the body of the paper.
%
% WHAT CHANGED RELATIVE TO THE PREVIOUS VERSION
%
% The paper's targeting condition is X_t in DeltaX_n, where DeltaX_n is the
% image of the whole admissible dp interval after n iterations. The previous
% version tested only whether X_t lay between the two endpoint values
% X_n(p_L) and X_n(p_R), and discarded any same-signature interval that
% failed that test without subdividing it. Those two conditions agree only
% when X_n(p) is monotone on the interval. Because X_n(p) is generically
% non-monotone within a single continuous branch, and increasingly so as n
% grows, the previous version could report failure at a horizon where a
% solution exists in the interior of a continuous piece. Subdivision was
% driven entirely by the location of discontinuities, so the search had no
% resolution parameter of its own.
%
% This version subdivides every same-signature interval that does not
% bracket the target, rather than discarding it, and states its resolution
% explicitly:
%
%   Phase 1, certified sweep. Every interval shallower than
%   control.minCertifiedDepth is subdivided unconditionally, breadth first.
%   Consequently every continuous feature of X_n(p) wider than
%   (pHighFull - pLowFull) / 2^minCertifiedDepth in p is sampled at least
%   once. This is the search's declared resolution and is reported in
%   search.certifiedParameterResolution.
%
%   Phase 2, adaptive descent. Below that depth the remaining intervals are
%   examined best first, ordered by how close their sampled images come to
%   X_t, and are pruned only when an estimated Lipschitz bound on X_n(p)
%   rules the target out. The estimate is measured from same-branch pairs
%   already evaluated at this horizon and inflated by
%   control.lipschitzSafetyFactor. Pruning is therefore credible at the
%   declared resolution, not certified. Termination is reported through
%   search.parameterResolutionFloorReached, search.intervalBudgetExhausted
%   and the interval counters, so an inconclusive horizon is distinguishable
%   from an exhausted one.
%
% A failed refinement no longer abandons the horizon. The traversal
% continues, so a second bracket in the same continuous piece is still
% reachable.
%
% Every p value evaluated at each horizon is recorded, together with the
% observed span of X_n. search.imageIntervalByHorizon is a lower bound on
% DeltaX_n and can be compared against the delta*exp(lambda*n) growth that
% Equation (3) of the paper rests on.
%
% Backwards compatibility. The signature, the returned field names used by
% search_parameter_to_target, assert_verified_bisection_search and
% test_corrected_targeting_methods, and the meaning of every previously
% reported diagnostic are unchanged. gridEvaluations is still never
% incremented and no grid metadata is written. New settings are optional
% with defaults, and none of them is copied by
% build_bisection_noise_config, so saved noise-sweep configuration
% signatures remain valid.

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
if ~isfield(control, 'lipschitzSafetyFactor')
    control.lipschitzSafetyFactor = 4;
end
if ~isfield(control, 'minEvaluationsBeforeExclusion')
    control.minEvaluationsBeforeExclusion = 6;
end
if ~isfield(control, 'maxBracketRefinements')
    control.maxBracketRefinements = max(control.bisectionIterations, 60);
end
if ~isfield(control, 'selectSmallestParameter')
    % false reproduces the paper, which stops at the first perturbation that
    % reaches the target. true exhausts the horizon and then applies the
    % lexicographic rule (smallest |p|, then smallest target error), which is
    % what the hybrid grid method applies. Set it to true when comparing the
    % two methods on returned perturbation magnitude, otherwise the
    % comparison measures selection rules rather than search quality.
    control.selectSmallestParameter = false;
end

% ---- diagnostics retained from the previous version ---------------------
search.crossingOrderTestUsed = true;
search.discontinuitiesDetected = 0;
search.discontinuityReductions = 0;
search.continuousIntervalsExamined = 0;
search.continuousTargetBracketFound = false;
search.discontinuityIsolationDepthReached = false;
search.crossingSignatureEvaluations = 0;
stats.crossingSignatureEvaluations = 0;

% ---- diagnostics describing the resolution of the new search ------------
search.certifiedParameterResolution = ...
    abs(pHighFull - pLowFull) * 2^(-control.minCertifiedDepth);
search.verifiedCandidatesThisHorizon = 0;
search.continuousIntervalsSubdivided = 0;
search.continuousIntervalsPrunedByBound = 0;
search.continuousIntervalsAbandonedAtResolution = 0;
search.parameterResolutionFloorReached = false;
search.intervalBudgetExhausted = false;
search.intervalsExaminedByHorizon = zeros(numHorizons, 1);
search.lipschitzEstimateByHorizon = NaN(numHorizons, 1);
search.imageIntervalByHorizon = NaN(numHorizons, 2);
search.evaluatedParametersByHorizon = cell(numHorizons, 1);
search.evaluatedImagesByHorizon = cell(numHorizons, 1);
search.searchResolutionSettings = struct( ...
    'minCertifiedDepth', control.minCertifiedDepth, ...
    'maxContinuousSubdivisionDepth', control.maxContinuousSubdivisionDepth, ...
    'maxIntervalsExamined', control.maxIntervalsExamined, ...
    'minParameterInterval', control.minParameterInterval, ...
    'lipschitzSafetyFactor', control.lipschitzSafetyFactor, ...
    'maxBracketRefinements', control.maxBracketRefinements, ...
    'selectSmallestParameter', control.selectSmallestParameter);
search.completenessNote = sprintf( ...
    ['Every continuous feature of X_n(p) wider than %.6g in p is sampled ', ...
    'at least once. Below that width, intervals are examined best first ', ...
    'and pruned using an estimated Lipschitz bound inflated by %g; such ', ...
    'pruning is credible at the declared resolution, not certified.'], ...
    search.certifiedParameterResolution, control.lipschitzSafetyFactor);

% ---- per-horizon state, shared with the nested functions ----------------
horizon = 0;
cacheP = zeros(1, 0);
cacheRecords = empty_record_array();
lipschitzEstimate = 0;
evaluationsThisHorizon = 0;
intervalsExamined = 0;

for horizonIndex = 1:numHorizons
    horizon = horizonIndex;
    cacheP = zeros(1, 0);
    cacheRecords = empty_record_array();
    lipschitzEstimate = 0;
    evaluationsThisHorizon = 0;
    intervalsExamined = 0;

    left = evaluate_at(pLowFull);
    right = evaluate_at(pHighFull);
    update_lipschitz(left, right);

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
        % Best-first traversal of the admissible parameter interval. Phase 1
        % intervals carry a priority equal to their depth, which orders them
        % breadth first and guarantees the declared resolution. Phase 2
        % intervals carry a priority above every phase 1 depth, ordered by
        % proximity of the sampled images to the target.
        hitFound = false;
        hitRecord = initialLeft;
        bracketLow = initialLeft.p;
        bracketHigh = initialRight.p;
        refineIterations = 0;
        collectMode = control.selectSmallestParameter;
        candidates = empty_candidate_array();

        queue = make_item(initialLeft, initialRight, 0);

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
                % Neither endpoint produced a usable image. Nothing can be
                % concluded and nothing can be refined here.
                search.continuousIntervalsAbandonedAtResolution = ...
                    search.continuousIntervalsAbandonedAtResolution + 1;
                continue;
            end

            width = abs(intervalRight.p - intervalLeft.p);
            atFloor = width <= control.minParameterInterval;
            inCertifiedSweep = current.depth < control.minCertifiedDepth;

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
                    % Refinement was inconclusive, usually because a fold or
                    % a discontinuity splits the bracket. Fall through and
                    % subdivide instead of abandoning the horizon.
                elseif ~inCertifiedSweep && ...
                        target_excluded_by_bound(intervalLeft, intervalRight, width)
                    search.continuousIntervalsPrunedByBound = ...
                        search.continuousIntervalsPrunedByBound + 1;
                    continue;
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
                stats.bisectionIterations = stats.bisectionIterations + 1;
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

            % Different ordered crossing signatures. A discontinuity of the
            % return map separates the two endpoints, so the interval is
            % reduced rather than refined. This is the paper's footnote [11]
            % criterion.
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
            stats.bisectionIterations = stats.bisectionIterations + 1;
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

        % Queue exhausted. In collect mode, apply the lexicographic rule.
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
        % Bisection inside one continuous bracket, following the paper: the
        % half whose endpoint images still contain X_t is retained.
        found = false;
        selected = intervalLeft;
        iterationsUsed = 0;

        for iteration = 1:control.maxBracketRefinements
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

            if leftValid && rightValid
                % Both halves contain a solution by the intermediate value
                % theorem, so either choice is sound. Descend toward the
                % closer endpoint.
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
                % A fold or a discontinuity has split the bracket. Return so
                % the outer traversal subdivides this region instead.
                break;
            end
        end

        candidates = [intervalLeft, intervalRight];
        errors = abs([candidates.xFinal] - target.x);
        [~, bestIndex] = min(errors);
        selected = candidates(bestIndex);
        found = is_hit(selected);
    end

    function candidates = empty_candidate_array()
        candidates = struct('record', {}, 'p', {}, 'error', {}, ...
            'bracketLow', {}, 'bracketHigh', {}, 'iterations', {});
    end

    function candidates = add_candidate(candidates, record, ...
            bracketLow, bracketHigh, iterations)
        % Store one verified perturbation. Duplicates in p are skipped so the
        % lexicographic tie-break is not decided by how often a value was
        % revisited through the cache.
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

    function tf = target_excluded_by_bound(recordLeft, recordRight, width)
        % Exclude the interval only when an estimated Lipschitz bound on
        % X_n(p) puts the whole image farther from X_t than the tolerance.
        % For p in the left half, |X(p) - X_t| >= dLeft - L*width/2, and
        % symmetrically on the right, so min(dLeft, dRight) - L*width/2 is a
        % bound on the closest approach over the interval.
        % A zero estimate is admissible rather than disqualifying. Pruning is
        % unreachable above the certified depth, so by the time it is applied
        % every feature wider than the declared resolution has been sampled;
        % a flat sample there is evidence, not ignorance. Disqualifying a
        % zero estimate instead drives a constant X_n(p) to the interval
        % budget at every horizon.
        tf = false;
        if evaluationsThisHorizon < control.minEvaluationsBeforeExclusion
            return;
        end
        if ~recordLeft.success || ~recordRight.success || ...
                ~isfinite(lipschitzEstimate)
            return;
        end
        distanceLeft = abs(recordLeft.xFinal - target.x);
        distanceRight = abs(recordRight.xFinal - target.x);
        excursion = control.lipschitzSafetyFactor * ...
            lipschitzEstimate * width / 2;
        tf = (min(distanceLeft, distanceRight) - excursion) > target.tolerance;
    end

    function [midpoint, degenerate] = midpoint_of(recordLeft, recordRight)
        pMid = 0.5 * (recordLeft.p + recordRight.p);
        degenerate = (pMid == recordLeft.p) || (pMid == recordRight.p);
        if degenerate
            midpoint = recordLeft;
            return;
        end
        midpoint = evaluate_at(pMid);
        update_lipschitz(recordLeft, midpoint);
        update_lipschitz(midpoint, recordRight);
    end

    function update_lipschitz(recordLeft, recordRight)
        % Largest observed |dX/dp| across same-branch pairs at this horizon.
        if ~same_branch(recordLeft, recordRight)
            return;
        end
        width = abs(recordRight.p - recordLeft.p);
        if width <= 0
            return;
        end
        slope = abs(recordRight.xFinal - recordLeft.xFinal) / width;
        if isfinite(slope)
            lipschitzEstimate = max(lipschitzEstimate, slope);
        end
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
            % Breadth first through the certified sweep.
            value = depth;
            return;
        end
        distance = interval_distance_to_target(recordLeft, recordRight);
        % Map the distance into [0,1) so that every phase 2 interval sorts
        % after every phase 1 interval regardless of scale.
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
        % Cache exact dyadic parameter values within the current horizon.
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
        record.numberOfPlaneCrossings = ...
            signatureResult.numberOfPlaneCrossings;
        record.numberOfAcceptedCrossings = ...
            signatureResult.numberOfAcceptedCrossings;

        cacheP(end + 1) = pValue;
        cacheRecords(end + 1) = record;
        evaluationsThisHorizon = evaluationsThisHorizon + 1;

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
        search.lipschitzEstimateByHorizon(horizonToRecord) = lipschitzEstimate;
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
            'Shinbrot discontinuity-aware interval search, image-resolved';
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
