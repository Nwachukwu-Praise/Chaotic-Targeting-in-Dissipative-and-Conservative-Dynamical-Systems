function [connection, profile] = so_resolve_connection(zSource, targetComponents, cfg, profile, backwardCache)
%SO_RESOLVE_CONNECTION Resolved minimum-time forward-backward targeting.
%
% The search proceeds by increasing total time.  The returned time is the
% minimum over the resolved search domain, not an unconditional global
% optimum when earlier split holes exist.
if nargin < 5
    backwardCache = [];
end
sourceComponent = so_make_control_component(zSource, cfg);
connection = so_empty_connection();
if so_runtime_exceeded(cfg)
    connection = timeout_connection(connection);
    return;
end
maxF = cfg.maxForwardIterations;
maxB = cfg.maxBackwardIterations;
maxT = cfg.maxTotalTransferTime;

[direct, profile] = direct_target_containment(zSource, sourceComponent, targetComponents, cfg, profile);
if direct.success
    connection = direct;
    return;
end

[forwardFamily, profile] = so_build_curve_family(sourceComponent, 1, maxF, cfg, profile, 'forward');

if isempty(backwardCache)
    [backwardCache, profile] = so_build_backward_cache(targetComponents, cfg, profile);
else
    profile.backwardFamilyCacheReuses = profile.backwardFamilyCacheReuses + numel(targetComponents);
end

componentCount = numel(targetComponents);
splitStatus = strings(maxF + 1, maxB + 1, componentCount);
splitStatus(:, :, :) = "not_evaluated";
resolvablePairMask = false(maxF + 1, maxB + 1, componentCount);
failureRows = {};
best = so_empty_connection();
winningTotal = Inf;
earlierUnresolved = false;
winningDiagonalUnresolved = false;

for totalTime = 0:maxT
    if so_runtime_exceeded(cfg)
        connection = timeout_connection(connection);
        connection.splitStatus = splitStatus;
        connection.resolvablePairMask = resolvablePairMask;
        connection.resolutionFailures = make_failure_table(failureRows);
        connection.resolutionHoles = height(connection.resolutionFailures);
        return;
    end
    diagonalCandidates = {};
    diagonalUnresolved = false;
    for nForward = 0:totalTime
        nBackward = totalTime - nForward;
        if nForward > maxF || nBackward > maxB
            continue;
        end
        if so_runtime_exceeded(cfg)
            connection = timeout_connection(connection);
            connection.splitStatus = splitStatus;
            connection.resolvablePairMask = resolvablePairMask;
            connection.resolutionFailures = make_failure_table(failureRows);
            connection.resolutionHoles = height(connection.resolutionFailures);
            return;
        end
        fCurve = forwardFamily{nForward + 1};
        fResolved = strcmp(fCurve.resolutionStatus, "resolved");
        for c = 1:componentCount
            bCurve = backwardCache.families{c}{nBackward + 1};
            bResolved = strcmp(bCurve.resolutionStatus, "resolved");
            if fResolved && bResolved
                splitStatus(nForward + 1, nBackward + 1, c) = "resolved_no_crossing";
                resolvablePairMask(nForward + 1, nBackward + 1, c) = true;
                [candidate, profile] = test_resolved_split(sourceComponent, targetComponents{c}, ...
                    fCurve, bCurve, nForward, nBackward, cfg, profile);
                if candidate.success
                    splitStatus(nForward + 1, nBackward + 1, c) = "resolved_crossing";
                    diagonalCandidates{end + 1} = candidate; %#ok<AGROW>
                end
            else
                if ~fResolved && ~bResolved
                    status = "unresolved_both";
                elseif ~fResolved
                    status = "unresolved_forward_budget";
                else
                    status = "unresolved_backward_budget";
                end
                profile.unresolvedSplits = profile.unresolvedSplits + 1;
                splitStatus(nForward + 1, nBackward + 1, c) = status;
                diagonalUnresolved = true;
                failureRows(end + 1, :) = {nForward, nBackward, c, status, ...
                    string(fCurve.resolutionStatus), string(bCurve.resolutionStatus), ...
                    fCurve.pointCount, bCurve.pointCount, fCurve.maximumGap, bCurve.maximumGap, ...
                    fCurve.maximumMidpointDeviation, bCurve.maximumMidpointDeviation}; %#ok<AGROW>
            end
        end
    end

    if ~isempty(diagonalCandidates)
        for i = 1:numel(diagonalCandidates)
            cand = diagonalCandidates{i};
            if so_connection_better(cand, best)
                best = cand;
            end
        end
        winningTotal = totalTime;
        winningDiagonalUnresolved = diagonalUnresolved;
        break;
    end
    earlierUnresolved = earlierUnresolved || diagonalUnresolved;
end

if best.success
    connection = best;
    connection.tauResolved = winningTotal;
    connection.resolvablePairMask = resolvablePairMask;
    connection.splitStatus = splitStatus;
    connection.resolutionFailures = make_failure_table(failureRows);
    connection.resolvableForwardEnvelope = last_resolved_index(forwardFamily);
    connection.resolvableBackwardEnvelope = min_backward_envelope(backwardCache);
    connection.resolutionHoles = height(connection.resolutionFailures);
    connection.timeMinimumCertified = ~earlierUnresolved;
    connection.selectionCertified = ~earlierUnresolved && ~winningDiagonalUnresolved;
else
    connection.splitStatus = splitStatus;
    connection.resolvablePairMask = resolvablePairMask;
    connection.resolutionFailures = make_failure_table(failureRows);
    connection.resolvableForwardEnvelope = last_resolved_index(forwardFamily);
    connection.resolvableBackwardEnvelope = min_backward_envelope(backwardCache);
    connection.resolutionHoles = height(connection.resolutionFailures);
end
end

function [connection, profile] = direct_target_containment(zSource, sourceComponent, targetComponents, cfg, profile)
connection = so_empty_connection();
best = so_empty_connection();
for c = 1:numel(targetComponents)
    comp = targetComponents{c};
    controls = admissible_zero_time_controls(zSource, sourceComponent, comp, cfg);
    profile.directContainmentCandidates = profile.directContainmentCandidates + numel(controls);
    for q = 1:numel(controls)
        cand = make_direct_connection(zSource, controls(q), comp);
        if so_connection_better(cand, best)
            best = cand;
        end
    end
end
if best.success
    connection = best;
end
end

function controls = admissible_zero_time_controls(zSource, sourceComponent, comp, cfg)
analyticControls = [];
switch string(comp.type)
    case "circle"
        analyticControls = circle_zero_time_controls(zSource, comp, cfg);
    case "rectangle"
        analyticControls = rectangle_zero_time_controls(zSource, comp, cfg);
end

samples = linspace(-cfg.controlAmplitude, cfg.controlAmplitude, cfg.curve.initialControlSamples);
zSamples = sourceComponent.gamma0(samples);
inside = comp.contains(zSamples, cfg.containmentTolerance);
sampleControls = samples(inside);

controls = unique_tol([analyticControls(:); sampleControls(:)], 1e-13);
if isempty(controls)
    controls = zeros(0, 1);
end
end

function controls = circle_zero_time_controls(zSource, comp, cfg)
dx = so_wrap_diff_x(zSource(1) - comp.center(1));
rad2 = comp.radius^2 - dx^2;
if rad2 < -cfg.containmentTolerance
    controls = zeros(0, 1);
    return;
end
halfHeight = sqrt(max(0, rad2));
yLo = comp.center(2) - halfHeight - zSource(2);
yHi = comp.center(2) + halfHeight - zSource(2);
controls = minimum_abs_control_in_interval(yLo, yHi, cfg.controlAmplitude, cfg.containmentTolerance);
end

function controls = rectangle_zero_time_controls(zSource, comp, cfg)
rect = comp.rectangle;
x = so_wrap_x(zSource(1));
if rect.xMax >= rect.xMin
    xOK = x >= rect.xMin - cfg.containmentTolerance && x <= rect.xMax + cfg.containmentTolerance;
else
    xOK = x >= rect.xMin - cfg.containmentTolerance || x <= rect.xMax + cfg.containmentTolerance;
end
if ~xOK
    controls = zeros(0, 1);
    return;
end
yLo = rect.yMin - zSource(2);
yHi = rect.yMax - zSource(2);
controls = minimum_abs_control_in_interval(yLo, yHi, cfg.controlAmplitude, cfg.containmentTolerance);
end

function controls = minimum_abs_control_in_interval(yLo, yHi, controlAmplitude, tol)
lo = max(yLo, -controlAmplitude);
hi = min(yHi, controlAmplitude);
if lo > hi + tol
    controls = zeros(0, 1);
elseif lo <= 0 && hi >= 0
    controls = 0;
elseif abs(lo) < abs(hi)
    controls = lo;
else
    controls = hi;
end
end

function values = unique_tol(values, tol)
values = sort(values(:));
keep = true(size(values));
for i = 2:numel(values)
    keep(i) = abs(values(i) - values(find(keep(1:i-1), 1, 'last'))) > tol;
end
values = values(keep);
end

function connection = make_direct_connection(zSource, u, comp)
zPost = zSource + [0; u];
connection = so_empty_connection();
connection.success = true;
connection.nForward = 0;
connection.nBackward = 0;
connection.totalIterations = 0;
connection.tauResolved = 0;
connection.control = u;
connection.controlledInitialState = zPost;
connection.preControlState = zSource;
connection.postControlState = zPost;
connection.targetComponent = string(comp.id);
connection.targetPhasePointIndex = comp.phasePointIndex;
connection.targetChainID = string(comp.chainID);
connection.targetBoundaryParameter = NaN;
connection.targetBoundaryPoint = zPost;
connection.intersectionPoint = zPost;
connection.intersectionResidual = 0;
connection.forwardBackwardSplit = [0, 0];
connection.refinementBracket = [u, u, NaN, NaN];
connection.plannedPath = zPost;
connection.plannedSwitchIndex = 0;
connection.plannedSwitchState = zPost;
connection.plannedFinalState = zPost;
connection.finalState = zPost;
connection.targetContained = true;
connection.directContainment = true;
connection.zeroTimePretest = true;
connection.timeMinimumCertified = true;
connection.selectionCertified = true;
end

function [candidate, profile] = test_resolved_split(sourceComponent, targetComponent, fCurve, bCurve, ...
    nForward, nBackward, cfg, profile)
candidate = so_empty_connection();

if nBackward == 0
    inside = targetComponent.contains(fCurve.pointsLifted, cfg.containmentTolerance);
    if any(inside)
        idxCandidates = find(inside);
        [~, bestLocal] = min(abs(fCurve.parameters(idxCandidates)));
        idx = idxCandidates(bestLocal);
        u = fCurve.parameters(idx);
        refined = direct_forward_hit(sourceComponent, targetComponent, nForward, u, cfg);
        if refined.success
            candidate = connection_from_refined(refined, targetComponent, nForward, nBackward);
            return;
        end
    end
end

[hits, stats] = so_find_polyline_intersections_indexed(fCurve, bCurve, cfg);
profile.indexedCrossingQueries = profile.indexedCrossingQueries + stats.indexedQueries;
profile.exactSegmentTestsAfterSpatialFiltering = ...
    profile.exactSegmentTestsAfterSpatialFiltering + stats.exactTests;
profile.naiveComparisonsAvoided = profile.naiveComparisonsAvoided + stats.naiveAvoided;

for h = 1:numel(hits)
    profile.refinementCalls = profile.refinementCalls + 1;
    refined = so_refine_intersection(sourceComponent, targetComponent, nForward, nBackward, ...
        hits(h).sA, hits(h).sB, hits(h).shiftB, cfg, ...
        [hits(h).sALo, hits(h).sAHi], [hits(h).sBLo, hits(h).sBHi]);
    if refined.success
        cand = connection_from_refined(refined, targetComponent, nForward, nBackward);
        if so_connection_better(cand, candidate)
            candidate = cand;
        end
    end
end
end

function refined = direct_forward_hit(sourceComponent, targetComponent, nForward, u, cfg)
zControlled = sourceComponent.gamma0(u);
path = so_iterate(zControlled, nForward, cfg, 1, true);
finalState = path(:, end);
refined.success = targetComponent.contains(finalState, cfg.containmentTolerance);
refined.control = u;
refined.controlledInitialState = zControlled;
refined.intersectionPoint = finalState;
refined.intersectionResidual = 0;
refined.plannedPath = path;
refined.finalState = finalState;
refined.targetContained = refined.success;
refined.boundaryParameter = NaN;
refined.targetBoundaryPoint = finalState;
refined.bracket = [u, u, NaN, NaN];
refined.preControlState = sourceComponent.center;
refined.postControlState = zControlled;
end

function connection = connection_from_refined(refined, targetComponent, nForward, nBackward)
connection = so_empty_connection();
connection.success = true;
connection.nForward = nForward;
connection.nBackward = nBackward;
connection.totalIterations = nForward + nBackward;
connection.tauResolved = nForward + nBackward;
connection.control = refined.control;
connection.controlledInitialState = refined.controlledInitialState;
if isfield(refined, 'preControlState')
    connection.preControlState = refined.preControlState;
else
    connection.preControlState = [refined.controlledInitialState(1); refined.controlledInitialState(2) - refined.control];
end
connection.postControlState = refined.controlledInitialState;
connection.targetComponent = string(targetComponent.id);
connection.targetPhasePointIndex = targetComponent.phasePointIndex;
connection.targetChainID = string(targetComponent.chainID);
connection.targetBoundaryParameter = refined.boundaryParameter;
connection.targetBoundaryPoint = refined.targetBoundaryPoint;
connection.intersectionPoint = refined.intersectionPoint;
connection.intersectionResidual = refined.intersectionResidual;
connection.forwardBackwardSplit = [nForward, nBackward];
if isfield(refined, 'bracket')
    connection.refinementBracket = refined.bracket;
end
connection.plannedPath = refined.plannedPath;
connection.plannedSwitchIndex = nForward;
connection.plannedSwitchState = refined.plannedPath(:, nForward + 1);
connection.plannedFinalState = refined.finalState;
connection.finalState = refined.finalState;
connection.targetContained = refined.targetContained;
end

function tbl = make_failure_table(rows)
names = {'nForward','nBackward','componentIndex','status','forwardStatus','backwardStatus', ...
    'forwardPointCount','backwardPointCount','forwardMaxGap','backwardMaxGap', ...
    'forwardMaxMidpointDeviation','backwardMaxMidpointDeviation'};
if isempty(rows)
    tbl = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    tbl = cell2table(rows, 'VariableNames', names);
end
end

function n = last_resolved_index(family)
n = -1;
for i = 1:numel(family)
    if strcmp(family{i}.resolutionStatus, "resolved")
        n = i - 1;
    end
end
end

function n = min_backward_envelope(cache)
n = Inf;
for c = 1:numel(cache.families)
    n = min(n, last_resolved_index(cache.families{c}));
end
if isinf(n)
    n = -1;
end
end
