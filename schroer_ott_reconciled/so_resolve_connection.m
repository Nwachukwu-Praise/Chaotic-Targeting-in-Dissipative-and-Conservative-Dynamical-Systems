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
    diagonalCandidates = {};
    diagonalUnresolved = false;
    for nForward = 0:totalTime
        nBackward = totalTime - nForward;
        if nForward > maxF || nBackward > maxB
            continue;
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
for c = 1:numel(targetComponents)
    comp = targetComponents{c};
    if comp.contains(zSource, cfg.containmentTolerance)
        connection = make_direct_connection(zSource, 0, comp);
        return;
    end
    samples = linspace(-cfg.controlAmplitude, cfg.controlAmplitude, cfg.curve.initialControlSamples);
    zSamples = sourceComponent.gamma0(samples);
    inside = comp.contains(zSamples, cfg.containmentTolerance);
    if any(inside)
        u = samples(find(inside, 1, 'first'));
        connection = make_direct_connection(sourceComponent.gamma0(u), u, comp);
        return;
    end
end
end

function connection = make_direct_connection(z, u, comp)
connection = so_empty_connection();
connection.success = true;
connection.nForward = 0;
connection.nBackward = 0;
connection.totalIterations = 0;
connection.tauResolved = 0;
connection.control = u;
connection.controlledInitialState = z;
connection.targetComponent = string(comp.id);
connection.targetPhasePointIndex = comp.phasePointIndex;
connection.targetChainID = string(comp.chainID);
connection.targetBoundaryParameter = NaN;
connection.targetBoundaryPoint = z;
connection.intersectionPoint = z;
connection.intersectionResidual = 0;
connection.plannedPath = z;
connection.finalState = z;
connection.targetContained = true;
connection.directContainment = true;
connection.timeMinimumCertified = true;
connection.selectionCertified = true;
end

function [candidate, profile] = test_resolved_split(sourceComponent, targetComponent, fCurve, bCurve, ...
    nForward, nBackward, cfg, profile)
candidate = so_empty_connection();

if nBackward == 0
    inside = targetComponent.contains(fCurve.pointsLifted, cfg.containmentTolerance);
    if any(inside)
        idx = find(inside, 1, 'first');
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
    refined = so_refine_intersection(sourceComponent, targetComponent, nForward, nBackward, ...
        hits(h).sA, hits(h).sB, hits(h).shiftB, cfg);
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
connection.targetComponent = string(targetComponent.id);
connection.targetPhasePointIndex = targetComponent.phasePointIndex;
connection.targetChainID = string(targetComponent.chainID);
connection.targetBoundaryParameter = refined.boundaryParameter;
connection.targetBoundaryPoint = refined.targetBoundaryPoint;
connection.intersectionPoint = refined.intersectionPoint;
connection.intersectionResidual = refined.intersectionResidual;
connection.plannedPath = refined.plannedPath;
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
