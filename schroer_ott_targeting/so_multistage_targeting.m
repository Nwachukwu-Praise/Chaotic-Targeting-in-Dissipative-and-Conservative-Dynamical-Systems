function result = so_multistage_targeting(cfg, catalogue, route)
%SO_MULTISTAGE_TARGETING Execute first-light Schroer-Ott proxy targeting.
%
% Pass geometry is not an operational condition here.  Switches are chosen
% by minimizing J(j)=j+tau_res(w_j -> B_next) over resolved probes.
tStart = tic;
profile = so_new_performance_profile();
actualState = so_rectangle_center(cfg.sourceRectangle);
actualState(1) = so_wrap_x(actualState(1));

proxyCount = numel(route.targetComponents);
stagePlans = repmat(empty_stage(), 0, 1);
switchEvaluations = repmat(empty_switch_eval(), 0, 1);
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
    if so_runtime_exceeded(cfg)
        stage = empty_stage();
        stage.stage = stageIndex + 1;
        stage.actualStateAtStart = actualState;
        [result, failureDiagnostics] = abort_result('runtime_limit_exceeded', stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
            resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
        result.failureDiagnostics = failureDiagnostics;
        return;
    end
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
    if ~provisional.success
        [result, failureDiagnostics] = abort_result('provisional_connection_failed', stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
            resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
        result.failureDiagnostics = failureDiagnostics;
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
    if so_runtime_exceeded(cfg)
        stagePlans(end + 1) = stage; %#ok<AGROW>
        [result, failureDiagnostics] = abort_result('runtime_limit_exceeded', stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
            resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
        result.failureDiagnostics = failureDiagnostics;
        return;
    end
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
        failureCategory = 'all_switch_probes_failed';
        if probes_include_runtime_timeout(probes)
            failureCategory = 'runtime_limit_exceeded';
        end
        [result, failureDiagnostics] = abort_result(failureCategory, stage, ...
            cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
            executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
            resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
        result.failureDiagnostics = failureDiagnostics;
        return;
    end

    stage.switchProbeConnections = probes;
    stage.selectedSwitchIndex = bestProbe.j;
    stage.selectedNextConnection = bestProbe.connection;

    if bestProbe.j == 0
        stage.skipped = true;
        skippedProxies(end + 1, 1) = stage.targetID; %#ok<AGROW>
        adoptedConnection = bestProbe.connection;
        adoptedTargetIndex = targetIndex + 1;
        stage.executedSegment = [];
        stage.actualStateAtEnd = actualState;
        stage.runtimeSeconds = toc(stageTimer);
        profile = record_stage_runtime(profile, stage.runtimeSeconds);
    else
        executed = execute_truncated_segment(actualState, provisional, bestProbe.j, cfg);
        profile = record_replay_runtime(profile, executed.replay.runtimeSeconds);
        if abs(provisional.control) > cfg.controlAmplitude + 10 * eps
            stagePlans(end + 1) = stage; %#ok<AGROW>
            [result, failureDiagnostics] = abort_result('control_bound_violation', stage, ...
                cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
                executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
                resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
            result.failureDiagnostics = failureDiagnostics;
            return;
        end
        if executed.consistencyError > cfg.propagationConsistencyTolerance
            stagePlans(end + 1) = stage; %#ok<AGROW>
            [result, failureDiagnostics] = abort_result('propagation_consistency_failed', stage, ...
                cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
                executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
                resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
            result.failureDiagnostics = failureDiagnostics;
            return;
        end
        actualState = executed.finalState;
        totalExecutedIterations = totalExecutedIterations + bestProbe.j;
        executedTrajectory = append_executed_path(executedTrajectory, executed); %#ok<AGROW>
        executedSwitchIndices(end + 1, 1) = bestProbe.j; %#ok<AGROW>
        controlRows(end + 1, :) = control_row(stageIndex, string(stage.targetID), ...
            executed, bestProbe.j, totalExecutedIterations); %#ok<AGROW>
        executedControls = make_control_table(controlRows);
        adoptedConnection = bestProbe.connection;
        adoptedTargetIndex = targetIndex + 1;
        stage.skipped = false;
        stage.executedSegment = executed;
        stage.actualStateAtEnd = actualState;
        stage.runtimeSeconds = toc(stageTimer);
        profile = record_stage_runtime(profile, stage.runtimeSeconds);
    end

    resolutionFailures = append_tables(resolutionFailures, provisional.resolutionFailures);
    resolutionFailures = append_tables(resolutionFailures, bestProbe.connection.resolutionFailures);
    stagePlans(end + 1) = stage; %#ok<AGROW>
end

stageIndex = stageIndex + 1;
stageTimer = tic;
if so_runtime_exceeded(cfg)
    finalStage = empty_stage();
    finalStage.stage = stageIndex;
    finalStage.targetKind = "final";
    finalStage.targetID = "final-target";
    finalStage.actualStateAtStart = actualState;
    stagePlans(end + 1) = finalStage; %#ok<AGROW>
    [result, failureDiagnostics] = abort_result('runtime_limit_exceeded', finalStage, ...
        cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
        executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
        resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
    result.failureDiagnostics = failureDiagnostics;
    return;
end
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
        executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
        resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
    result.failureDiagnostics = failureDiagnostics;
    return;
end

finalExecuted = execute_full_connection(actualState, finalConnection, cfg, route.finalTargetComponents);
profile = record_replay_runtime(profile, finalExecuted.replay.runtimeSeconds);
if finalExecuted.consistencyError > cfg.propagationConsistencyTolerance || ~finalExecuted.targetContained
    finalStage.executedSegment = finalExecuted;
    stagePlans(end + 1) = finalStage; %#ok<AGROW>
    [result, failureDiagnostics] = abort_result('final_independent_containment_failed', finalStage, ...
        cfg, catalogue, route, stagePlans, switchEvaluations, executedControls, ...
        executedSwitchIndices, executedTrajectory, totalExecutedIterations, profile, ...
        resolutionFailures, tieBreakLog, failureDiagnostics, tStart);
    result.failureDiagnostics = failureDiagnostics;
    return;
end

totalExecutedIterations = totalExecutedIterations + finalConnection.totalIterations;
executedTrajectory = append_executed_path(executedTrajectory, finalExecuted); %#ok<AGROW>
executedSwitchIndices(end + 1, 1) = finalConnection.totalIterations; %#ok<AGROW>
controlRows(end + 1, :) = control_row(stageIndex, "final-target", finalExecuted, ...
    finalConnection.totalIterations, totalExecutedIterations); %#ok<AGROW>
finalStage.executedSegment = finalExecuted;
finalStage.actualStateAtEnd = finalExecuted.finalState;
finalStage.runtimeSeconds = toc(stageTimer);
profile = record_stage_runtime(profile, finalStage.runtimeSeconds);
stagePlans(end + 1) = finalStage; %#ok<AGROW>
resolutionFailures = append_tables(resolutionFailures, finalConnection.resolutionFailures);

executedControls = make_control_table(controlRows);
result = completed_result(cfg, catalogue, route, stagePlans, switchEvaluations, ...
    executedControls, executedSwitchIndices, executedTrajectory, totalExecutedIterations, ...
    finalExecuted.finalState, finalExecuted.targetContained, profile, resolutionFailures, ...
    tieBreakLog, skippedProxies, tStart);
end

function [bestProbe, probes, profile] = evaluate_switch_probes(actualState, provisional, ...
    nextTarget, nextBackwardCache, cfg, profile)
maxJ = min(provisional.totalIterations, cfg.maxTotalTransferTime);
probes = repmat(empty_probe(), 0, 1);
bestProbe = empty_probe();
Jbest = Inf;
for j = 0:maxJ
    if so_runtime_exceeded(cfg)
        timeoutProbe = empty_probe();
        timeoutProbe.j = j;
        timeoutProbe.connection.diagnostics.runtimeLimitExceeded = true;
        timeoutProbe.connection.diagnostics.failureCategory = "runtime_limit_exceeded";
        probes(end + 1) = timeoutProbe; %#ok<AGROW>
        return;
    end
    if so_should_prune_probe(j, Jbest)
        pruned = empty_probe();
        pruned.j = j;
        pruned.pruned = true;
        probes(end + 1) = pruned; %#ok<AGROW>
        profile.jProbesPruned = profile.jProbesPruned + 1;
        continue;
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
replay = so_replay_stage(actualState, connection.control, connection.plannedPath, j, cfg, {});
if size(connection.plannedPath, 2) >= j + 1
    expected = connection.plannedPath(:, j + 1);
else
    expected = [NaN; NaN];
end
executed.control = connection.control;
executed.iterations = j;
executed.stageIterations = j;
executed.preControlState = replay.preControlState;
executed.postControlState = replay.postControlState;
executed.path = replay.replayPath;
executed.finalState = replay.finalState;
executed.switchState = replay.finalState;
executed.expectedState = expected;
executed.consistencyError = replay.maxPathwiseReplayError;
executed.maxPathwiseReplayError = replay.maxPathwiseReplayError;
executed.switchStateReplayError = replay.switchStateReplayError;
executed.finalStateReplayError = replay.finalStateReplayError;
executed.targetContained = false;
executed.replayPassed = replay.replayPassed;
executed.controlWithinBound = replay.controlWithinBound;
executed.replay = replay;
end

function executed = execute_full_connection(actualState, connection, cfg, targetComponents)
replay = so_replay_stage(actualState, connection.control, connection.plannedPath, ...
    connection.totalIterations, cfg, targetComponents);
executed.control = connection.control;
executed.iterations = connection.totalIterations;
executed.stageIterations = connection.totalIterations;
executed.preControlState = replay.preControlState;
executed.postControlState = replay.postControlState;
executed.path = replay.replayPath;
executed.finalState = replay.finalState;
executed.switchState = replay.finalState;
executed.expectedState = connection.finalState;
executed.consistencyError = replay.maxPathwiseReplayError;
executed.maxPathwiseReplayError = replay.maxPathwiseReplayError;
executed.switchStateReplayError = replay.switchStateReplayError;
executed.finalStateReplayError = replay.finalStateReplayError;
executed.targetContained = replay.targetContained;
executed.replayPassed = replay.replayPassed;
executed.controlWithinBound = replay.controlWithinBound;
executed.replay = replay;
end

function tf = route_final_contains(connection, state, cfg)
component = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
tf = connection.targetContained && component.contains(state, cfg.containmentTolerance);
end

function tbl = make_control_table(rows)
names = {'stage','targetID','preControlX','preControlY','postControlX','postControlY', ...
    'switchStateX','switchStateY','controlY','stageIterations','executedIterations', ...
    'cumulativeIterations','maxPathwiseReplayError','switchStateReplayError', ...
    'finalStateReplayError','targetContained','replayPassed','controlWithinBound', ...
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
    executedControls, executedSwitchIndices, executedTrajectory, totalExecutedIterations, ...
    finalState, targetContained, profile, resolutionFailures, tieBreakLog, skippedProxies, tStart)
result.configuration = cfg;
result.orbitCatalogue = catalogue;
result.route = route;
result.sourceProxyAudit = so_source_proxy_audit(cfg, route);
result.stagePlans = stagePlans;
result.switchEvaluations = switchEvaluations;
result.executedControls = executedControls;
result.executedSwitchIndices = executedSwitchIndices;
result.executedTrajectory = executedTrajectory;
result.executionSegments = collect_execution_segments(stagePlans);
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

function row = control_row(stageIndex, targetID, executed, stageIterations, cumulativeIterations)
row = {stageIndex, targetID, executed.preControlState(1), executed.preControlState(2), ...
    executed.postControlState(1), executed.postControlState(2), ...
    executed.switchState(1), executed.switchState(2), executed.control, ...
    stageIterations, stageIterations, cumulativeIterations, ...
    executed.maxPathwiseReplayError, executed.switchStateReplayError, ...
    executed.finalStateReplayError, executed.targetContained, executed.replayPassed, ...
    executed.controlWithinBound, executed.switchState(1), executed.switchState(2), ...
    executed.consistencyError};
end

function trajectory = append_executed_path(trajectory, executed)
trajectory = [trajectory, [NaN; NaN], executed.preControlState, ...
    executed.postControlState, executed.path(:, 2:end)];
end

function profile = record_stage_runtime(profile, seconds)
profile.stageRuntimeSeconds(end + 1, 1) = seconds;
profile.runtimeByStage = profile.stageRuntimeSeconds;
end

function profile = record_replay_runtime(profile, seconds)
profile.replayRuntimeSeconds(end + 1, 1) = seconds;
end

function segments = collect_execution_segments(stagePlans)
segments = struct([]);
for i = 1:numel(stagePlans)
    if isfield(stagePlans(i), 'executedSegment') && ~isempty(stagePlans(i).executedSegment)
        if isempty(segments)
            segments = stagePlans(i).executedSegment;
        else
            segments(end + 1) = stagePlans(i).executedSegment; %#ok<AGROW>
        end
    end
end
end

function tf = probes_include_runtime_timeout(probes)
tf = false;
for i = 1:numel(probes)
    if isfield(probes(i).connection, 'diagnostics') && ...
            isfield(probes(i).connection.diagnostics, 'runtimeLimitExceeded') && ...
            probes(i).connection.diagnostics.runtimeLimitExceeded
        tf = true;
        return;
    end
end
end

function [result, failureDiagnostics] = abort_result(category, stage, cfg, catalogue, route, ...
    stagePlans, switchEvaluations, executedControls, executedSwitchIndices, executedTrajectory, ...
    totalExecutedIterations, profile, resolutionFailures, tieBreakLog, failureDiagnostics, tStart)
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
if isfield(diag.connection, 'diagnostics') && isfield(diag.connection.diagnostics, 'runtimeLimitExceeded') && ...
        diag.connection.diagnostics.runtimeLimitExceeded
    category = 'runtime_limit_exceeded';
    diag.failureCategory = category;
end
failureDiagnostics = [failureDiagnostics; diag];
result = completed_result(cfg, catalogue, route, stagePlans, switchEvaluations, ...
    executedControls, executedSwitchIndices, executedTrajectory, totalExecutedIterations, ...
    executedTrajectory(:, end), false, profile, resolutionFailures, tieBreakLog, strings(0, 1), tStart);
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
