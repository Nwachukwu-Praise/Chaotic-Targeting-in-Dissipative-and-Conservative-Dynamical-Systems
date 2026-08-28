function result = so_multistage_targeting(cfg, catalogue, route)
%SO_MULTISTAGE_TARGETING Schroer-Ott proxy targeting (reconciled build).
%
% The targeting logic is the fundamental-folder engine, unchanged: pass
% geometry is not an operational condition, and switches are chosen by
% minimizing J(j) = j + tau_res(w_j -> B_next) over resolved probes.
%
% Three additions, none of which alter the selected controls:
%   1. result.executionSegments records, per applied control, the state
%      before and after the y-kick and the replayed path.  The figures and
%      the noise validation consume this; nothing in the search reads it.
%   2. An optional wall-clock guard (cfg.runtime, see so_runtime_exceeded)
%      aborts a trial cleanly instead of stalling an ensemble.
%   3. A route with zero proxy chains is executed as a single direct
%      forward-backward stage onto the final target (horizontal case).
tStart = tic;
profile = so_new_performance_profile();
actualState = so_rectangle_center(cfg.sourceRectangle);
actualState(1) = so_wrap_x(actualState(1));
sourceState = actualState;

proxyCount = numel(route.targetComponents);
stagePlans = repmat(empty_stage(), 0, 1);
switchEvaluations = repmat(empty_switch_eval(), 0, 1);
executionSegments = repmat(empty_segment(), 0, 1);
executedControls = table();
executedSwitchIndices = [];
executedTrajectory = actualState;
skippedProxies = strings(0, 1);
resolutionFailures = table();
tieBreakLog = table();
failureDiagnostics = struct([]);

adoptedConnection = [];
adoptedTargetIndex = NaN;
totalExecutedIterations = 0;
controlRows = {};
stageIndex = 0;

for targetIndex = 1:proxyCount
    stageIndex = stageIndex + 1;
    stageTimer = tic;
    currentTarget = route.targetComponents{targetIndex};
    if ~isempty(adoptedConnection) && adoptedTargetIndex == targetIndex
        provisional = adoptedConnection;
    else
        [provisional, profile] = so_resolve_connection(actualState, currentTarget, cfg, profile);
    end
    stage = empty_stage();
    stage.stage = stageIndex;
    stage.targetIndex = targetIndex;
    stage.targetKind = "proxy";
    stage.targetID = string(route.chains(targetIndex).id);
    stage.rotationNumber = route.chains(targetIndex).omega;
    stage.provisionalConnection = provisional;
    stage.actualStateAtStart = actualState;
    if so_runtime_exceeded(cfg)
        stagePlans(end + 1) = stage; %#ok<AGROW>
        [result, failureDiagnostics] = abort_result('runtime_limit_exceeded', stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
            totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
            failureDiagnostics, tStart);
        return;
    end
    if ~provisional.success
        [result, failureDiagnostics] = abort_result('provisional_connection_failed', stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
            totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
            failureDiagnostics, tStart);
        return;
    end

    if targetIndex < proxyCount
        nextTarget = route.targetComponents{targetIndex + 1};
        nextTargetLabel = string(route.chains(targetIndex + 1).id);
    else
        nextTarget = route.finalTargetComponents;
        nextTargetLabel = "final-target";
    end

    [nextBackwardCache, profile] = so_build_backward_cache(nextTarget, cfg, profile);
    [bestProbe, probes, profile] = evaluate_switch_probes(actualState, provisional, ...
        nextTarget, nextBackwardCache, cfg, profile);
    sw = empty_switch_eval();
    sw.stage = stageIndex;
    sw.currentTargetID = stage.targetID;
    sw.nextTargetID = nextTargetLabel;
    sw.probes = probes;
    sw.selected = bestProbe;
    sw.pruned = sum([probes.pruned]);
    switchEvaluations(end + 1) = sw; %#ok<AGROW>

    if ~bestProbe.finite
        stagePlans(end + 1) = stage; %#ok<AGROW>
        % Distinguish "no switch works" from "we ran out of time looking",
        % otherwise an ensemble timeout is mis-reported as a real failure.
        if so_runtime_exceeded(cfg)
            category = 'runtime_limit_exceeded';
        else
            category = 'all_switch_probes_failed';
        end
        [result, failureDiagnostics] = abort_result(category, stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
            totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
            failureDiagnostics, tStart);
        return;
    end

    stage.switchProbeConnections = probes;
    stage.selectedSwitchIndex = bestProbe.j;
    stage.selectedNextConnection = bestProbe.connection;
    stage.runtimeSeconds = toc(stageTimer);

    if bestProbe.j == 0
        stage.skipped = true;
        skippedProxies(end + 1, 1) = stage.targetID; %#ok<AGROW>
        adoptedConnection = bestProbe.connection;
        adoptedTargetIndex = targetIndex + 1;
        stage.executedSegment = [];
        stage.actualStateAtEnd = actualState;
    else
        executed = execute_truncated_segment(actualState, provisional, bestProbe.j, cfg);
        if abs(provisional.control) > cfg.controlAmplitude + 10 * eps
            stagePlans(end + 1) = stage; %#ok<AGROW>
            [result, failureDiagnostics] = abort_result('control_bound_violation', stage, ...
                cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
                executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
                totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
                failureDiagnostics, tStart);
            return;
        end
        if executed.consistencyError > cfg.propagationConsistencyTolerance
            stagePlans(end + 1) = stage; %#ok<AGROW>
            [result, failureDiagnostics] = abort_result('propagation_consistency_failed', stage, ...
                cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
                executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
                totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
                failureDiagnostics, tStart);
            return;
        end
        executionSegments(end + 1) = make_segment(stageIndex, stage.targetID, ...
            actualState, provisional.control, bestProbe.j, executed, false); %#ok<AGROW>
        actualState = executed.finalState;
        totalExecutedIterations = totalExecutedIterations + bestProbe.j;
        executedTrajectory = [executedTrajectory, executed.path(:, 2:end)]; %#ok<AGROW>
        executedSwitchIndices(end + 1, 1) = bestProbe.j; %#ok<AGROW>
        controlRows(end + 1, :) = {stageIndex, string(stage.targetID), provisional.control, ...
            bestProbe.j, totalExecutedIterations, actualState(1), actualState(2), ...
            executed.consistencyError}; %#ok<AGROW>
        executedControls = make_control_table(controlRows);
        adoptedConnection = bestProbe.connection;
        adoptedTargetIndex = targetIndex + 1;
        stage.skipped = false;
        stage.executedSegment = executed;
        stage.actualStateAtEnd = actualState;
    end

    resolutionFailures = append_tables(resolutionFailures, provisional.resolutionFailures);
    resolutionFailures = append_tables(resolutionFailures, bestProbe.connection.resolutionFailures);
    stagePlans(end + 1) = stage; %#ok<AGROW>
end

stageIndex = stageIndex + 1;
stageTimer = tic;
if ~isempty(adoptedConnection) && adoptedTargetIndex == proxyCount + 1
    finalConnection = adoptedConnection;
else
    [finalConnection, profile] = so_resolve_connection(actualState, route.finalTargetComponents, cfg, profile);
end
finalStage = empty_stage();
finalStage.stage = stageIndex;
finalStage.targetIndex = proxyCount + 1;
finalStage.targetKind = "final";
finalStage.targetID = "final-target";
finalStage.rotationNumber = NaN;
finalStage.provisionalConnection = finalConnection;
finalStage.actualStateAtStart = actualState;
if ~finalConnection.success
    stagePlans(end + 1) = finalStage; %#ok<AGROW>
    [result, failureDiagnostics] = abort_result('final_connection_failed', finalStage, ...
        cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
        executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
        totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
        failureDiagnostics, tStart);
    return;
end

finalExecuted = execute_full_connection(actualState, finalConnection, cfg);
if finalExecuted.consistencyError > cfg.propagationConsistencyTolerance || ~finalExecuted.targetContained
    finalStage.executedSegment = finalExecuted;
    stagePlans(end + 1) = finalStage; %#ok<AGROW>
    [result, failureDiagnostics] = abort_result('final_independent_containment_failed', finalStage, ...
        cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
        executedSwitchIndices, executedTrajectory, executionSegments, sourceState, ...
        totalExecutedIterations, profile, resolutionFailures, tieBreakLog, ...
        failureDiagnostics, tStart);
    return;
end

executionSegments(end + 1) = make_segment(stageIndex, "final-target", actualState, ...
    finalConnection.control, finalConnection.totalIterations, finalExecuted, true);
totalExecutedIterations = totalExecutedIterations + finalConnection.totalIterations;
executedTrajectory = [executedTrajectory, finalExecuted.path(:, 2:end)];
executedSwitchIndices(end + 1, 1) = finalConnection.totalIterations;
controlRows(end + 1, :) = {stageIndex, "final-target", finalConnection.control, ...
    finalConnection.totalIterations, totalExecutedIterations, finalExecuted.finalState(1), ...
    finalExecuted.finalState(2), finalExecuted.consistencyError};
finalStage.executedSegment = finalExecuted;
finalStage.actualStateAtEnd = finalExecuted.finalState;
finalStage.runtimeSeconds = toc(stageTimer);
stagePlans(end + 1) = finalStage;
resolutionFailures = append_tables(resolutionFailures, finalConnection.resolutionFailures);

executedControls = make_control_table(controlRows);
result = completed_result(cfg, catalogue, route, stagePlans, switchEvaluations, ...
    executedControls, executedSwitchIndices, executedTrajectory, executionSegments, ...
    sourceState, totalExecutedIterations, finalExecuted.finalState, ...
    finalExecuted.targetContained, profile, resolutionFailures, tieBreakLog, ...
    skippedProxies, tStart);
end

function [bestProbe, probes, profile] = evaluate_switch_probes(actualState, provisional, ...
    nextTarget, nextBackwardCache, cfg, profile)
maxJ = min(provisional.totalIterations, cfg.maxTotalTransferTime);
probes = repmat(empty_probe(), 0, 1);
bestProbe = empty_probe();
Jbest = Inf;
for j = 0:maxJ
    if so_should_prune_probe(j, Jbest)
        pruned = empty_probe();
        pruned.j = j;
        pruned.pruned = true;
        probes(end + 1) = pruned; %#ok<AGROW>
        profile.jProbesPruned = profile.jProbesPruned + 1;
        continue;
    end
    if so_runtime_exceeded(cfg)
        break;
    end
    if j == 0
        wj = actualState;
    else
        wj = so_iterate(actualState + [0; provisional.control], j, cfg, 1, false);
    end
    [conn, profile] = so_resolve_connection(wj, nextTarget, cfg, profile, nextBackwardCache);
    probe = empty_probe();
    probe.j = j;
    probe.sourceState = wj;
    probe.connection = conn;
    probe.finite = conn.success;
    if conn.success
        probe.objectiveJ = j + conn.totalIterations;
    end
    probes(end + 1) = probe; %#ok<AGROW>
    if so_switch_candidate_better(probe, bestProbe)
        bestProbe = probe;
        Jbest = probe.objectiveJ;
    end
end
end

function executed = execute_truncated_segment(actualState, connection, j, cfg)
expected = so_iterate(actualState + [0; connection.control], j, cfg, 1, false);
path = so_iterate(actualState + [0; connection.control], j, cfg, 1, true);
executed.control = connection.control;
executed.iterations = j;
executed.path = path;
executed.finalState = path(:, end);
executed.expectedState = expected;
executed.consistencyError = norm(executed.finalState - expected);
executed.targetContained = false;
end

function executed = execute_full_connection(actualState, connection, cfg)
path = so_iterate(actualState + [0; connection.control], connection.totalIterations, cfg, 1, true);
executed.control = connection.control;
executed.iterations = connection.totalIterations;
executed.path = path;
executed.finalState = path(:, end);
executed.expectedState = connection.finalState;
executed.consistencyError = so_cylinder_distance(executed.finalState, connection.finalState);
executed.targetContained = route_final_contains(connection, executed.finalState, cfg);
end

function tf = route_final_contains(connection, state, cfg)
component = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
tf = connection.targetContained && component.contains(state, cfg.containmentTolerance);
end

function seg = make_segment(stageIndex, targetID, preState, control, iterations, executed, isFinal)
seg = empty_segment();
seg.stage = stageIndex;
seg.targetID = string(targetID);
seg.control = control;
seg.iterations = iterations;
seg.preControlState = preState;
seg.postControlState = preState + [0; control];
seg.path = executed.path;
seg.switchState = executed.path(:, end);
seg.finalState = executed.finalState;
seg.consistencyError = executed.consistencyError;
seg.isFinal = isFinal;
end

function seg = empty_segment()
seg.stage = NaN;
seg.targetID = "";
seg.control = NaN;
seg.iterations = NaN;
seg.preControlState = [NaN; NaN];
seg.postControlState = [NaN; NaN];
seg.path = zeros(2, 0);
seg.switchState = [NaN; NaN];
seg.finalState = [NaN; NaN];
seg.consistencyError = NaN;
seg.isFinal = false;
end

function tbl = make_control_table(rows)
names = {'stage','targetID','controlY','executedIterations','cumulativeIterations', ...
    'stateX','stateY','propagationConsistencyError'};
if isempty(rows)
    tbl = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    tbl = cell2table(rows, 'VariableNames', names);
end
end

function out = append_tables(a, b)
if isempty(b) || width(b) == 0 || height(b) == 0
    out = a;
elseif isempty(a) || width(a) == 0 || height(a) == 0
    out = b;
else
    out = [a; b];
end
end

function result = completed_result(cfg, catalogue, route, stagePlans, switchEvaluations, ...
    executedControls, executedSwitchIndices, executedTrajectory, executionSegments, ...
    sourceState, totalExecutedIterations, finalState, targetContained, profile, ...
    resolutionFailures, tieBreakLog, skippedProxies, tStart)
result.configuration = cfg;
result.orbitCatalogue = catalogue;
result.route = route;
result.stagePlans = stagePlans;
result.switchEvaluations = switchEvaluations;
result.executedControls = executedControls;
result.executedSwitchIndices = executedSwitchIndices;
result.executedTrajectory = executedTrajectory;
result.executionSegments = executionSegments;
result.sourceState = sourceState;
result.totalExecutedIterations = totalExecutedIterations;
result.numberOfControls = height(executedControls);
result.finalState = finalState;
result.targetContained = targetContained;
result.equationOneDiagnostics = so_equation_one_diagnostics(stagePlans, cfg);
result.tieBreakLog = tieBreakLog;
result.resolutionFailures = resolutionFailures;
result.failureDiagnostics = struct([]);
result.performanceProfile = profile;
result.skippedProxies = skippedProxies;
result.targetReached = targetContained;
result.runtimeSeconds = toc(tStart);
end

function [result, failureDiagnostics] = abort_result(category, stage, cfg, catalogue, route, ...
    stagePlans, switchEvaluations, executedControls, executedSwitchIndices, executedTrajectory, ...
    executionSegments, sourceState, totalExecutedIterations, profile, resolutionFailures, ...
    tieBreakLog, failureDiagnostics, tStart)
diag.stage = stage.stage;
diag.failureCategory = category;
diag.sourceState = stage.actualStateAtStart;
diag.targetID = stage.targetID;
diag.rotationNumber = stage.rotationNumber;
diag.iterationBounds = [cfg.maxForwardIterations, cfg.maxBackwardIterations, cfg.maxTotalTransferTime];
if isfield(stage, 'provisionalConnection')
    diag.connection = stage.provisionalConnection;
else
    diag.connection = so_empty_connection();
end
failureDiagnostics = [failureDiagnostics; diag];
result = completed_result(cfg, catalogue, route, stagePlans, switchEvaluations, ...
    executedControls, executedSwitchIndices, executedTrajectory, executionSegments, ...
    sourceState, totalExecutedIterations, executedTrajectory(:, end), false, profile, ...
    resolutionFailures, tieBreakLog, strings(0, 1), tStart);
result.failureDiagnostics = failureDiagnostics;
result.targetReached = false;
end

function stage = empty_stage()
stage.stage = NaN;
stage.targetIndex = NaN;
stage.targetKind = "";
stage.targetID = "";
stage.rotationNumber = NaN;
stage.actualStateAtStart = [NaN; NaN];
stage.actualStateAtEnd = [NaN; NaN];
stage.provisionalConnection = so_empty_connection();
stage.switchProbeConnections = repmat(empty_probe(), 0, 1);
stage.selectedSwitchIndex = NaN;
stage.selectedNextConnection = so_empty_connection();
stage.executedSegment = [];
stage.skipped = false;
stage.runtimeSeconds = NaN;
end

function sw = empty_switch_eval()
sw.stage = NaN;
sw.currentTargetID = "";
sw.nextTargetID = "";
sw.probes = repmat(empty_probe(), 0, 1);
sw.selected = empty_probe();
sw.pruned = 0;
end

function probe = empty_probe()
probe.j = NaN;
probe.sourceState = [NaN; NaN];
probe.connection = so_empty_connection();
probe.objectiveJ = Inf;
probe.finite = false;
probe.pruned = false;
end
