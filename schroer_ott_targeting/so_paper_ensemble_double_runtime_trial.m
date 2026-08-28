function [record, result] = so_paper_ensemble_double_runtime_trial(sourceIndex, sourceState, ...
    baseCfg, paper, catalogue, runtimeBudgetSeconds, seed, configurationFingerprint)
%SO_PAPER_ENSEMBLE_DOUBLE_RUNTIME_TRIAL Run one preserved paper source.
record = so_empty_paper_ensemble_record();
record.sourceIndex = sourceIndex;
record.sourceState = sourceState(:);
record.sourceX = sourceState(1);
record.sourceY = sourceState(2);
record.runtimeBudgetSeconds = runtimeBudgetSeconds;
record.seed = seed;
record.configurationFingerprint = string(configurationFingerprint);

cfg = baseCfg;
cfg.runtimeLimitSecondsPerTrial = runtimeBudgetSeconds;
cfg.sourceRectangle = point_rectangle(sourceState, sprintf('paper-source-%02d', sourceIndex));

trialTimer = tic;
try
    cfg.runtime.startTime = trialTimer;
    cfg.runtime.limitSeconds = runtimeBudgetSeconds;
    route = so_construct_omega_route(catalogue, cfg, paper.expectedRotationNumbers, ...
        'estimated Figure-3 seven-resonance route');
    result = so_multistage_targeting(cfg, catalogue, route);
    record.runtimeSeconds = toc(trialTimer);
    record = populate_record_from_result(record, result, paper);
catch err
    result = err;
    record.runtimeSeconds = toc(trialTimer);
    record.terminationReason = "numerical_failure";
    record.legacyFailureCategory = string(err.identifier);
    record.errorIdentifier = string(err.identifier);
    record.errorMessage = string(err.message);
    record.finalTargetDistance = NaN;
end
end

function rect = point_rectangle(z, id)
rect = struct('xMin', z(1), 'xMax', z(1), 'yMin', z(2), 'yMax', z(2), 'id', id);
end

function record = populate_record_from_result(record, result, paper)
record.found = isfield(result, 'targetReached') && result.targetReached;
[record.terminationReason, record.legacyFailureCategory] = classify_result(result);

if isfield(result, 'route') && isfield(result.route, 'rotationNumbers')
    record.selectedRoute = result.route.rotationNumbers;
end
if isfield(result, 'totalExecutedIterations')
    record.iterations = result.totalExecutedIterations;
end
if isfield(result, 'numberOfControls')
    record.controls = result.numberOfControls;
end
if isfield(result, 'finalState')
    record.finalState = result.finalState(:);
    record.finalX = record.finalState(1);
    record.finalY = record.finalState(2);
    record.finalTargetDistance = rectangle_distance(record.finalState, paper.target);
end
if isfield(result, 'targetContained')
    record.targetContained = result.targetContained;
end
if isfield(result, 'performanceProfile') && isfield(result.performanceProfile, 'unresolvedSplits')
    record.unresolvedSplits = result.performanceProfile.unresolvedSplits;
end

if isfield(result, 'executedControls') && ~isempty(result.executedControls) && height(result.executedControls) > 0
    controls = result.executedControls.controlY(:).';
    record.controlMagnitudes = controls;
    record.maxAbsControl = max(abs(controls));
    record.totalAbsControl = sum(abs(controls));
    record.switchPoints = [result.executedControls.stateX(:).'; result.executedControls.stateY(:).'];
    record.maxPathwiseReplayError = max(result.executedControls.maxPathwiseReplayError);
else
    record.controlMagnitudes = zeros(1, 0);
    record.switchPoints = zeros(2, 0);
end

if isfield(result, 'stagePlans') && ~isempty(result.stagePlans)
    residuals = [];
    certified = true;
    for i = 1:numel(result.stagePlans)
        conn = result.stagePlans(i).provisionalConnection;
        if isfield(conn, 'intersectionResidual') && isfinite(conn.intersectionResidual)
            residuals(end + 1) = conn.intersectionResidual; %#ok<AGROW>
        end
        if isfield(conn, 'timeMinimumCertified')
            certified = certified && conn.timeMinimumCertified;
        end
    end
    if ~isempty(residuals)
        record.maxIntersectionResidual = max(residuals);
    end
    record.minimumCertified = certified;
end
end

function [reason, legacy] = classify_result(result)
if isfield(result, 'targetReached') && result.targetReached
    reason = "success";
    legacy = "";
    return;
end

legacy = "not_reached";
if isfield(result, 'failureDiagnostics') && ~isempty(result.failureDiagnostics)
    diagnostics = result.failureDiagnostics;
    if isfield(diagnostics, 'failureCategory')
        legacy = string(diagnostics(end).failureCategory);
    end
end

switch legacy
    case "runtime_limit_exceeded"
        reason = "time_budget_exhausted";
    case {"all_switch_probes_failed", "switch_probe_budget_exhausted"}
        reason = "candidate_budget_exhausted";
    case {"provisional_connection_failed", "final_connection_failed", "not_reached"}
        reason = "no_eligible_route";
    case {"propagation_consistency_failed", "final_independent_containment_failed", ...
            "control_bound_violation"}
        reason = "invalid_replay";
    otherwise
        reason = "numerical_failure";
end
end

function d = rectangle_distance(z, rect)
if any(~isfinite(z))
    d = NaN;
    return;
end
x = so_wrap_x(z(1));
if rect.xMax >= rect.xMin
    xInside = x >= rect.xMin && x <= rect.xMax;
else
    xInside = x >= rect.xMin || x <= rect.xMax;
end
if xInside
    dx = 0;
else
    dx = min(abs(so_wrap_diff_x(x - rect.xMin)), abs(so_wrap_diff_x(x - rect.xMax)));
end
if z(2) < rect.yMin
    dy = rect.yMin - z(2);
elseif z(2) > rect.yMax
    dy = z(2) - rect.yMax;
else
    dy = 0;
end
d = hypot(dx, dy);
end
