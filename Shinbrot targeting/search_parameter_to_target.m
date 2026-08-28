function search = search_parameter_to_target( ...
    sourceState, target, params, control, methodName)
%SEARCH_PARAMETER_TO_TARGET Find a directly verified bounded parameter.
%
% Supported methods:
%   shinbrotPaperBisection
%       Paper-faithful Shinbrot bisection with crossing-order checks.
%   shinbrotDiscontinuityAwareBisection
%       Later project extension with additional image-resolution machinery.
%   hybridGridBisection
%       Complete grid, all direct hits and adjacent brackets, followed by
%       local bisection inside every adjacent bracket.
%
% The selected p is held constant in the Lorenz equations from the original
% section source state through the requested accepted crossings. A search
% succeeds only when abs(X_n(p) - X_t) <= target.tolerance is evaluated
% directly.

if nargin < 5 || isempty(methodName)
    if isfield(control, 'searchMethod') && ~isempty(control.searchMethod)
        methodName = control.searchMethod;
    else
        error(['search_parameter_to_target requires an explicit ', ...
            'methodName or a nonempty control.searchMethod field.']);
    end
end

methodName = normalize_method_name(methodName);
pLowFull = -control.deltaP;
pHighFull = control.deltaP;
numHorizons = control.maxSearchCrossings;

search = initialize_search(numHorizons, methodName);
bestApprox = initialize_best_approximation();
stats = initialize_stats();
timerStart = tic;

switch methodName
    case 'shinbrotPaperBisection'
        [search, bestApprox, stats] = run_shinbrot_paper_bisection( ...
            search, bestApprox, stats, sourceState, target, params, ...
            control, pLowFull, pHighFull, numHorizons);
    case 'shinbrotDiscontinuityAwareBisection'
        [search, bestApprox, stats] = ...
            run_shinbrot_discontinuity_aware_bisection( ...
            search, bestApprox, stats, sourceState, target, params, ...
            control, pLowFull, pHighFull, numHorizons);
    case 'hybridGridBisection'
        pGrid = linspace(pLowFull, pHighFull, control.numParameterSamples);
        search = initialize_grid_storage(search, pGrid, numHorizons);
        [search, bestApprox, stats] = run_hybrid_grid_bisection( ...
            search, bestApprox, stats, sourceState, target, params, ...
            control, pGrid, numHorizons);
    otherwise
        error('Unknown search method: %s', methodName);
end

search.searchTime = toc(timerStart);
search = finalize_search(search, bestApprox, stats, target);
end

function methodName = normalize_method_name(methodName)
%NORMALIZE_METHOD_NAME Accept concise aliases for retained methods.

if isa(methodName, 'string')
    methodName = char(methodName);
end

switch lower(methodName)
    case {'shinbrot', 'shinbrotbisection', 'shinbrotpaper', ...
            'shinbrotpaperbisection', 'paper', 'paperbisection', ...
            'pureshinbrot', 'pureshinbrotbisection'}
        methodName = 'shinbrotPaperBisection';
    case {'shinbrotdiscontinuityawarebisection', ...
            'fullinterval', 'fullintervalbisection', ...
            'fullintervalbisectionbaseline'}
        methodName = 'shinbrotDiscontinuityAwareBisection';
    case {'hybrid', 'hybridgrid', 'hybridgridbisection'}
        methodName = 'hybridGridBisection';
    otherwise
        error('Unknown search method: %s', methodName);
end
end

function [search, bestApprox, stats] = run_hybrid_grid_bisection( ...
    search, bestApprox, stats, sourceState, target, params, control, ...
    pGrid, numHorizons)
%RUN_HYBRID_GRID_BISECTION Grid-first bracket discovery and refinement.

for horizon = 1:numHorizons
    [xValues, distanceValues, residualValues, stats, bestApprox] = ...
        evaluate_grid_horizon(sourceState, pGrid, horizon, params, ...
        target, stats, bestApprox);
    search.xByHorizon(horizon, :) = xValues;
    search.distanceByHorizon(horizon, :) = distanceValues;
    search.gridHorizonsEvaluated = horizon;

    features = find_hybrid_grid_features( ...
        pGrid, residualValues, target.tolerance);
    search.directGridHitsByHorizon(horizon) = ...
        numel(features.directHitIndices);
    search.adjacentBracketsByHorizon(horizon) = ...
        size(features.bracketIndexPairs, 1);

    candidates = empty_candidates();
    for index = features.directHitIndices(:).'
        candidates(end + 1) = make_candidate( ... %#ok<AGROW>
            pGrid(index), xValues(index), target, 'direct grid hit', ...
            pGrid(index), pGrid(index), 0);
    end

    for bracketNumber = 1:size(features.bracketIndexPairs, 1)
        indices = features.bracketIndexPairs(bracketNumber, :);
        lowIndex = indices(1);
        highIndex = indices(2);
        [found, pHit, xHit, iterationsUsed, stats, bestApprox] = ...
            refine_bracket(sourceState, target, params, horizon, ...
            pGrid(lowIndex), pGrid(highIndex), ...
            xValues(lowIndex), xValues(highIndex), control, ...
            stats, bestApprox);
        if found
            candidates(end + 1) = make_candidate( ... %#ok<AGROW>
                pHit, xHit, target, 'adjacent-grid bracket bisection', ...
                pGrid(lowIndex), pGrid(highIndex), iterationsUsed);
        end
    end

    search.verifiedCandidatesByHorizon(horizon) = numel(candidates);
    if isempty(candidates)
        continue;
    end

    [selected, uniqueCandidates] = ...
        select_verified_hybrid_candidate(candidates);
    search.verifiedCandidatesByHorizon(horizon) = numel(uniqueCandidates);
    search = record_hybrid_success(search, selected, horizon, target);
    return;
end
end

function candidates = empty_candidates()
%EMPTY_CANDIDATES Candidate schema shared with the selection helper.

candidates = struct('p', {}, 'x', {}, 'error', {}, 'origin', {}, ...
    'bracketLow', {}, 'bracketHigh', {}, 'iterations', {});
end

function candidate = make_candidate( ...
    pValue, xValue, target, origin, bracketLow, bracketHigh, iterations)
%MAKE_CANDIDATE Package one directly verified candidate.

candidate.p = pValue;
candidate.x = xValue;
candidate.error = abs(xValue - target.x);
candidate.origin = origin;
candidate.bracketLow = bracketLow;
candidate.bracketHigh = bracketHigh;
candidate.iterations = iterations;
end

function [xValues, distanceValues, residualValues, stats, bestApprox] = ...
    evaluate_grid_horizon(sourceState, pGrid, horizon, params, target, ...
    stats, bestApprox)
%EVALUATE_GRID_HORIZON Evaluate every grid point from the original source.

xValues = NaN(size(pGrid));
for i = 1:numel(pGrid)
    [xValues(i), ~, stats, bestApprox] = evaluate_parameter( ...
        sourceState, pGrid(i), horizon, params, target, stats, ...
        bestApprox, true);
end
residualValues = xValues - target.x;
distanceValues = abs(residualValues);
end

function [found, pHit, xHit, iterationsUsed, stats, bestApprox] = ...
    refine_bracket(sourceState, target, params, horizon, ...
    pLow, pHigh, xLow, xHigh, control, stats, bestApprox)
%REFINE_BRACKET Bisect one adjacent grid bracket and verify each trial.

found = false;
pHit = NaN;
xHit = NaN;
iterationsUsed = 0;
gLow = xLow - target.x;
gHigh = xHigh - target.x;
if ~isfinite(gLow) || ~isfinite(gHigh) || gLow * gHigh >= 0
    return;
end

for iteration = 1:control.bisectionIterations
    iterationsUsed = iteration;
    stats.bisectionIterations = stats.bisectionIterations + 1;
    pMid = 0.5 * (pLow + pHigh);
    [xMid, success, stats, bestApprox] = evaluate_parameter( ...
        sourceState, pMid, horizon, params, target, stats, ...
        bestApprox, false);
    if ~success || ~isfinite(xMid)
        return;
    end

    gMid = xMid - target.x;
    if abs(gMid) <= target.tolerance
        found = true;
        pHit = pMid;
        xHit = xMid;
        return;
    end

    if gLow * gMid < 0
        pHigh = pMid;
        gHigh = gMid;
    elseif gMid * gHigh < 0
        pLow = pMid;
        gLow = gMid;
    else
        % A sign-changing discontinuity need not contain a root.
        return;
    end
end
end

function [xFinal, success, stats, bestApprox] = evaluate_parameter( ...
    sourceState, pControl, horizon, params, target, stats, bestApprox, isGrid)
%EVALUATE_PARAMETER Evaluate X_n(p) and update benchmark counters.

stats.parameterEvaluations = stats.parameterEvaluations + 1;
if isGrid
    stats.gridEvaluations = stats.gridEvaluations + 1;
end

[xFinal, success, crossingsUsed] = x_after_n_crossings( ...
    sourceState, pControl, horizon, params);
stats.crossingPropagations = stats.crossingPropagations + crossingsUsed;
if success && isfinite(xFinal)
    bestApprox = update_best_approximation( ...
        bestApprox, pControl, horizon, xFinal, target);
end
end

function [xFinal, success, crossingsUsed] = x_after_n_crossings( ...
    sourceState, pControl, horizon, params)
%X_AFTER_N_CROSSINGS Hold p constant for the requested accepted crossings.

currentState = sourceState(:);
xFinal = NaN;
success = false;
crossingsUsed = 0;
for crossing = 1:horizon
    segment = next_valid_section_crossing(currentState, params, pControl);
    if ~segment.success
        return;
    end
    currentState = segment.eventState;
    crossingsUsed = crossingsUsed + 1;
end
xFinal = currentState(1);
success = true;
end

function search = initialize_search(numHorizons, methodName)
%INITIALIZE_SEARCH Create common outputs without assuming a parameter grid.

search.method = methodName;
search.selectionMethod = 'none';
search.found = false;
search.selectedP = NaN;
search.horizon = Inf;
search.targetX = NaN;
search.bestDistance = Inf;
search.finalTargetError = Inf;
search.failureCode = 'searchNotCompleted';
search.failureReason = '';

search.bracketLow = NaN;
search.bracketHigh = NaN;
search.finalBracketBisectionIterations = 0;
search.totalBisectionIterations = 0;
search.bisectionIterationsUsed = 0;
search.parameterEvaluations = 0;
search.crossingPropagations = 0;
search.gridEvaluations = 0;
search.searchTime = NaN;
search.controlledReplayCrossings = NaN;
search.replayHit = false;
search.replayFinalError = NaN;
search.replayTime = NaN;

search.bestApproxP = NaN;
search.bestApproxHorizon = Inf;
search.bestApproxTargetX = NaN;
search.bestApproxError = Inf;

search.pGrid = [];
search.xByHorizon = NaN(numHorizons, 0);
search.distanceByHorizon = NaN(numHorizons, 0);
search.plotHorizon = 1;
search.gridDistance = [];
search.gridMetadata = 'notApplicable';
search.gridHorizonsEvaluated = 0;
search.directGridHitsByHorizon = zeros(numHorizons, 1);
search.adjacentBracketsByHorizon = zeros(numHorizons, 1);
search.verifiedCandidatesByHorizon = zeros(numHorizons, 1);

search.crossingOrderTestUsed = false;
search.discontinuitiesDetected = 0;
search.discontinuityReductions = 0;
search.continuousIntervalsExamined = 0;
search.continuousTargetBracketFound = false;
search.discontinuityIsolationDepthReached = false;
search.crossingSignatureEvaluations = 0;
end

function search = initialize_grid_storage(search, pGrid, numHorizons)
%INITIALIZE_GRID_STORAGE Allocate grid diagnostics only for grid methods.

numP = numel(pGrid);
search.pGrid = pGrid;
search.xByHorizon = NaN(numHorizons, numP);
search.distanceByHorizon = NaN(numHorizons, numP);
search.gridDistance = NaN(1, numP);
search.gridMetadata = 'hybridGridBisection';
end

function stats = initialize_stats()
%INITIALIZE_STATS Counters common to both methods.

stats.parameterEvaluations = 0;
stats.crossingPropagations = 0;
stats.gridEvaluations = 0;
stats.bisectionIterations = 0;
stats.crossingSignatureEvaluations = 0;
end

function bestApprox = initialize_best_approximation()
%INITIALIZE_BEST_APPROXIMATION Retain closest finite image as a diagnostic.

bestApprox.p = NaN;
bestApprox.horizon = Inf;
bestApprox.x = NaN;
bestApprox.error = Inf;
end

function bestApprox = update_best_approximation( ...
    bestApprox, pValue, horizon, xValue, target)
%UPDATE_BEST_APPROXIMATION Retain the smallest finite target error.

errorValue = abs(xValue - target.x);
if errorValue < bestApprox.error
    bestApprox.p = pValue;
    bestApprox.horizon = horizon;
    bestApprox.x = xValue;
    bestApprox.error = errorValue;
end
end

function search = record_hybrid_success(search, candidate, horizon, target)
%RECORD_HYBRID_SUCCESS Store a verified candidate after full-horizon search.

search.found = candidate.error <= target.tolerance;
search.selectedP = candidate.p;
search.horizon = horizon;
search.targetX = candidate.x;
search.bestDistance = candidate.error;
search.finalTargetError = candidate.error;
search.selectionMethod = ['grid-first verified ', candidate.origin];
search.bracketLow = candidate.bracketLow;
search.bracketHigh = candidate.bracketHigh;
search.finalBracketBisectionIterations = candidate.iterations;
search.bisectionIterationsUsed = candidate.iterations;
search.plotHorizon = horizon;
search.gridDistance = search.distanceByHorizon(horizon, :);
search.failureCode = 'none';
search.failureReason = '';
end

function search = finalize_search(search, bestApprox, stats, target)
%FINALIZE_SEARCH Copy diagnostics and return structured numerical failure.

search.parameterEvaluations = stats.parameterEvaluations;
search.crossingPropagations = stats.crossingPropagations;
search.gridEvaluations = stats.gridEvaluations;
search.totalBisectionIterations = stats.bisectionIterations;
search.bestApproxP = bestApprox.p;
search.bestApproxHorizon = bestApprox.horizon;
search.bestApproxTargetX = bestApprox.x;
search.bestApproxError = bestApprox.error;

if search.found && search.finalTargetError <= target.tolerance
    search.failureCode = 'none';
    search.failureReason = '';
    return;
end

search.found = false;
search.selectedP = NaN;
search.horizon = Inf;
search.targetX = NaN;
search.finalTargetError = bestApprox.error;
search.bestDistance = bestApprox.error;
search.selectionMethod = 'no verified candidate';
search.failureCode = 'noVerifiedCandidate';
if strcmp(search.method, 'shinbrotPaperBisection')
    search.failureReason = ['No successful paper-bisection perturbation ', ...
        'was found within the selected targeting horizon, perturbation ', ...
        'interval, discontinuity-isolation limits and numerical tolerances.'];
else
    search.failureReason = ['No successful perturbation was found within ', ...
        'the selected targeting horizon, perturbation bound, grid resolution ', ...
        'and numerical tolerances.'];
end
if isfinite(bestApprox.horizon)
    search.plotHorizon = bestApprox.horizon;
    if bestApprox.horizon <= size(search.distanceByHorizon, 1)
        search.gridDistance = ...
            search.distanceByHorizon(bestApprox.horizon, :);
    end
end
end
